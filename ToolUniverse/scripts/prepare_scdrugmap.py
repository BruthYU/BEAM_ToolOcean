"""Build a compact read-only scDrugMap index beside existing extracted matrices (Python 3.9+)."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile


GEO_SHA256 = "6a228029df713006fb465f02cc3e1e17540b0ddaf4b69374b116aeb135cdf73f"
BCC_SHA256 = "8fad3be2c9777f9b56f7dcade4a54d1d71d6157e8a9f86711e6e75b509bd5235"
GEO_URLS = {
    "GSE120575": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz",
    "GSE123813": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123813/suppl/GSE123813_bcc_all_metadata.txt.gz",
}
REGIMENS = {
    "GSE120575_anti_PD1": "anti-PD1",
    "GSE120575_anti_CTLA4": "anti-CTLA4",
    "GSE120575_anti_CTLA4_PD1": "anti-CTLA4+PD1",
}
UNKNOWN = (None, None, None, None, None)
DATASET_FIELDS = ("dataset_id", "collection", "accession", "subgroup", "source", "tissue", "cancer",
                  "treatment_type", "regimen", "platform", "n_cells", "n_genes")
SCHEMA = """
PRAGMA user_version=1;
CREATE TABLE datasets(
 dataset_id TEXT PRIMARY KEY, collection TEXT, accession TEXT, subgroup TEXT, source TEXT, tissue TEXT,
 cancer TEXT, treatment_type TEXT, regimen TEXT, platform TEXT, n_cells INTEGER, n_genes INTEGER,
 n_patients INTEGER, n_samples INTEGER, metadata_json TEXT NOT NULL, matrix_path TEXT NOT NULL,
 matrix_size INTEGER NOT NULL, matrix_mtime_ns INTEGER NOT NULL, orientation TEXT NOT NULL
);
CREATE TABLE genes(dataset_id TEXT NOT NULL, gene_id TEXT NOT NULL, ordinal INTEGER NOT NULL,
 row_offset INTEGER, PRIMARY KEY(dataset_id,gene_id)) WITHOUT ROWID;
CREATE TABLE cells(dataset_id TEXT NOT NULL,cell_id TEXT NOT NULL,ordinal INTEGER NOT NULL,
 label TEXT NOT NULL,patient_id TEXT,sample_id TEXT,timepoint TEXT,row_offset INTEGER,
 response_original TEXT,sample_regimen TEXT,PRIMARY KEY(dataset_id,cell_id)) WITHOUT ROWID;
"""


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stable_barcode(barcode):
    return barcode.replace("-", ".")


def verified_geo(path, is_bcc):
    """Use only the exact published mapping version for which the cell/therapy joins were audited."""
    if path is None:
        return {}
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    require(digest == (BCC_SHA256 if is_bcc else GEO_SHA256), "Unrecognized GEO mapping SHA256")
    result = {}
    with gzip.open(path, "rt", encoding="utf-8" if is_bcc else "latin1", newline="") as handle:
        if is_bcc:
            for row in csv.DictReader(handle, delimiter="\t"):
                cell = row["cell.id"]
                require(cell not in result, "GSE123813 barcode collision")
                result[cell] = (row["patient"], row["treatment"])
        else:
            for row in csv.reader(handle, delimiter="\t"):
                if not row or not re.fullmatch(r"Sample \d+", row[0]):
                    continue
                cell = stable_barcode(row[1])
                require(cell not in result, "GSE120575 barcode collision")
                sample, response, regimen = row[4:7]
                require(re.fullmatch(r"(Pre|Post)_P\d+(_\d+)?", sample), "Unknown sample format")
                result[cell] = (sample.split("_")[1], sample, sample.split("_")[0].lower(), response, regimen)
    require(result, "Empty clinical mapping")
    return result


def clinical(record, cell, label, maps):
    acc = record["accession"]
    if not maps.get(acc):
        return UNKNOWN
    if acc == "GSE120575":
        item = maps[acc][stable_barcode(cell)]
        require(REGIMENS.get(record["subgroup"]) == item[4], "GEO therapy and archive slice differ")
        require({"Responder":"sensitive","Non-responder":"resistant"}[item[3]] == label, "GEO and archive labels differ")
        return item
    patient, timepoint = maps[acc][cell]
    row = record["biopsy_library_mapping"][cell.split("_",1)[0]]
    require((patient,timepoint) == (row["patient_id"],row["timepoint"]), "BCC biopsy identity mismatch")
    require({"Yes":"sensitive","No":"resistant"}[row["response_original"]] == label, "BCC patient label differs")
    return patient, row["sample_id"], timepoint, row["response_original"], row["regimen"]


def within(root, name):
    path = (root / name).resolve()
    require(root in path.parents and path.is_file(), f"Missing or invalid extracted matrix member: {name}")
    return path


def parse_csv_line(raw, delimiter=","):
    return next(csv.reader([raw.decode("utf-8-sig")],delimiter=delimiter))


def external_cell_labels(annotation, header):
    with annotation.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.reader(handle)
        require(next(rows) == ["", "condition"], "Unexpected annotation header")
        labels = {}
        for cell, label in rows:
            key = stable_barcode(cell)
            require(cell and key not in labels and label, "Empty or colliding annotation barcode or label")
            labels[key] = (cell, label)
    require(len(header) == len(labels), "Matrix and annotation cell counts differ")
    # Match by exact normalized identity; preserve matrix positions and public annotation IDs.
    for cell in header:
        key = stable_barcode(cell)
        require(key in labels, "Matrix barcode cannot be joined uniquely to the annotation")
        yield labels.pop(key)
    require(not labels, "Unmatched cell labels")


def primary_label(raw):
    prefix = raw.split(b"\t",2)[:2]
    require(len(prefix)==2, "Missing primary cell identifier or condition")
    return parse_csv_line(b"\t".join(prefix), "\t")


def gene_and_separator(raw):
    # Author CSVs have quoted and unquoted identifiers. A quoted comma inside a gene
    # must remain part of the identifier. Numeric columns are unquoted delimiter-separated.
    sep = raw.find(b'",') + 1 if raw.startswith(b'"') else raw.find(b",")
    require(sep > 0, "Missing CSV gene identifier or delimiter")
    gene = parse_csv_line(raw[:sep])[0]
    return gene, sep


def batch_store(db, table, records):
    size = 4 if table == "genes" else 10
    db.executemany(f"INSERT INTO {table} VALUES ({','.join('?' for _ in range(size))})", records)
    records.clear()


def index_one(db, root, rec, maps):
    dataset_id, kind = rec["dataset_id"], rec["collection"]
    require(kind in ("primary","external") and dataset_id.startswith(kind + ":"), "Unknown collection or ID")
    primary = kind == "primary"
    path = within(root,rec["matrix"])
    stat = path.stat()
    digest = hashlib.sha256()
    labels = collections.Counter()
    patients, samples = set(), set()
    n_cells = 0
    n_genes = 0
    mapped = 0
    cells, genes = [], []
    def add_gene(name, pos):
        nonlocal n_genes
        require(name != "", "Blank gene identifier")
        genes.append((dataset_id,name,n_genes,pos))
        n_genes += 1
        if len(genes) == 4096:
            batch_store(db,"genes",genes)
    def add_cell(name, label, pos):
        nonlocal n_cells, mapped
        require(bool(name and label), "Blank cell identifier or archive label")
        patient,sample,tp,response,rx = clinical(rec,name,label,maps)
        if patient is not None:
            patients.add(patient)
            samples.add(sample)
            mapped += 1
        labels[label] += 1
        cells.append((dataset_id,name,n_cells,label,patient,sample,tp,pos,response,rx))
        n_cells += 1
        if len(cells) == 4096:
            batch_store(db,"cells",cells)
    with path.open("rb", buffering=4*1024**2) as handle:
        raw_header = handle.readline()
        digest.update(raw_header)
        offset = len(raw_header)
        hdr = parse_csv_line(raw_header, "\t" if primary else ",")
        if primary:
            require(hdr[:2] == ["Cell_barcode","Condition"], "Unexpected primary header")
            for name in hdr[2:]:
                add_gene(name,None)
        else:
            require(len(hdr)>1, "Missing external expression header")
            for name,label in external_cell_labels(within(root,rec["annotation"]),hdr[1:]):
                add_cell(name,label,None)
        for raw in handle:
            digest.update(raw)
            if raw not in (b"\n",b"\r\n"):
                if primary:
                    require(raw.count(b"\t") == len(hdr)-1, "Irregular primary expression row")
                    name,label=primary_label(raw)
                    add_cell(name,label,offset)
                else:
                    name, sep=gene_and_separator(raw)
                    require(raw.count(b",",sep) == len(hdr)-1, "Irregular external expression row")
                    add_gene(name,offset)
            offset += len(raw)
    batch_store(db,"genes",genes)
    batch_store(db,"cells",cells)
    require((n_cells,n_genes)==(rec["n_cells"],rec["n_genes"]),f"{dataset_id} matrix dimensions differ")
    require(dict(labels)==rec["label_counts"], f"{dataset_id} deposited labels differ")
    require(path.stat().st_mtime_ns==stat.st_mtime_ns and offset==stat.st_size, "Source matrix changed during indexing")
    require(mapped in (0,n_cells), "Clinical linkage must be complete for each supported slice")
    accession=rec["accession"]
    metadata = {k:rec.get(k) for k in ("source_study","pmid","archive","source_record","clinical_response_source")}
    label_unit = {"GSE120575":"sampled_lesion_response", "GSE123813":"patient_clinical_response"}.get(accession,"archive_cell_condition")
    metadata.update(matrix_member=rec["matrix"],annotation_member=rec.get("annotation"),matrix_sha256=digest.hexdigest(),
                    mapped_cells=mapped,label_unit=label_unit,expression_scale="as_deposited",
                    label_counts=dict(labels),mapping_source=GEO_URLS.get(accession) if mapped else None,
                    mapping_sha256=(GEO_SHA256 if accession=="GSE120575" else BCC_SHA256) if mapped else None)
    if accession == "GSE120575":
        metadata["timepoint_definition"] = "pre=baseline; post=on-treatment (not necessarily after completion)"
    if accession == "GSE243649":
        metadata["cross_slice_overlap"] = "Identical archived expression matrices in external:Data9 and external:Data10"
    values = [rec[k] for k in DATASET_FIELDS]+[len(patients) if mapped else None,len(samples) if mapped else None,
        json.dumps(metadata,ensure_ascii=False),str(path),stat.st_size,stat.st_mtime_ns,"cell_by_gene" if primary else "gene_by_cell"]
    db.execute("INSERT INTO datasets VALUES ("+",".join("?" for _ in values)+")",values)
    return {"dataset_id":dataset_id,"n_cells":n_cells,"n_genes":n_genes,"mapped_cells":mapped}


def prepare(source, catalog, output, geo=None, bcc=None):
    output=Path(output).resolve()
    require(not output.exists(),"Choose a new output directory; an existing registry is never replaced")
    root=(Path(source)/"extracted").resolve()
    catalog=Path(catalog)
    records=json.loads(catalog.read_text(encoding="utf8"))
    require(isinstance(records,list) and records and len({r["dataset_id"] for r in records})==len(records),"Invalid catalog or nonunique IDs")
    output.parent.mkdir(parents=True,exist_ok=True)
    staging=Path(tempfile.mkdtemp(prefix="."+output.name+".building-",dir=output.parent))
    db=None
    try:
        maps={"GSE120575":verified_geo(geo,False),"GSE123813":verified_geo(bcc,True)}
        db=sqlite3.connect(staging/"registry.sqlite")
        db.executescript("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF;"+SCHEMA)
        totals=[]
        for rec in records:
            result=index_one(db,root,rec,maps)
            totals.append(result)
            print(f"{result['dataset_id']}: {result['n_cells']} cells, {result['n_genes']} genes",flush=True)
        db.executescript("CREATE UNIQUE INDEX cell_page ON cells(dataset_id,ordinal);"
                         "CREATE UNIQUE INDEX gene_position ON genes(dataset_id,ordinal);"
                         "CREATE INDEX cell_sample ON cells(dataset_id,sample_id,patient_id,timepoint);")
        db.commit()
        require(db.execute("PRAGMA quick_check").fetchall()==[("ok",)], "SQLite index integrity check failed")
        db.close()
        db=None
        summary={"schema_version":1,"catalog_sha256":hashlib.sha256(catalog.read_bytes()).hexdigest(),
                 "dataset_count":len(totals),"cell_records":sum(r["n_cells"] for r in totals),
                 "mapped_cell_records":sum(r["mapped_cells"] for r in totals),"datasets":totals}
        (staging/"manifest.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
        # A concurrent successful builder leaves a nonempty directory, so rename cannot replace it.
        os.rename(staging,output)
        return summary
    finally:
        if db is not None:
            db.close()
        if staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,required=True,help="scDrugMap collection root with extracted/<archive member>")
    default_catalog=Path(__file__).resolve().parent/"resources"/"scdrugmap_collection_catalog.json"
    parser.add_argument("--catalog",type=Path,default=default_catalog)
    parser.add_argument("--output",type=Path,required=True,help="New directory to publish the finished registry")
    parser.add_argument("--geo",type=Path,help="Published GSE120575_patient_ID_single_cells.txt.gz")
    parser.add_argument("--bcc",type=Path,help="Published GSE123813_bcc_all_metadata.txt.gz")
    args=parser.parse_args()
    prepare(args.source,args.catalog,args.output,args.geo,args.bcc)
