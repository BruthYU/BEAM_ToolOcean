"""Read-only access to the indexed scDrugMap collection. Indices are made by prepare_scdrugmap.py."""
import csv
import json
import math
import os
import sqlite3
from contextlib import closing
from pathlib import Path

from .base_tool import BaseTool
from .tool_registry import register_tool


DATASET_FILTERS = ("collection", "accession", "subgroup", "source", "tissue", "cancer",
                   "treatment_type", "regimen", "platform")
CASE_FILTERS = ("label", "patient_id", "sample_id", "timepoint", "response_original", "sample_regimen")
CATALOG_FIELDS = ("dataset_id", *DATASET_FILTERS, "n_cells", "n_genes", "n_samples", "n_patients")


class QueryError(ValueError):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


def require(condition, code, message):
    if not condition:
        raise QueryError(code, message)


def filter_sql(filters, allowed):
    require(isinstance(filters, dict), "invalid_filter", "filters must map field names to lists of exact values.")
    clauses, bindings = [], []
    for name, values in filters.items():
        require(name in allowed, "invalid_filter", f"Unsupported filter {name}. Allowed: {', '.join(allowed)}.")
        require(isinstance(values, list) and 1 <= len(values) <= 50 and
                all(isinstance(x, str) and 0 < len(x) <= 512 for x in values),
                "invalid_filter", f"{name} needs 1–50 nonempty strings.")
        clauses.append(name + " IN (" + ",".join("?" for _ in values) + ")")
        bindings.extend(values)
    return (" AND " + " AND ".join(clauses)) if clauses else "", bindings


def page(db, selection, params, order, args):
    limit, offset = args.get("limit", 50), args.get("offset", 0)
    require(type(limit) is int and 1 <= limit <= 200 and type(offset) is int and 0 <= offset <= 10**9,
            "invalid_pagination", "limit must be 1–200; offset must be a nonnegative integer up to 1 billion.")
    total = db.execute("SELECT COUNT(*) FROM (" + selection + ")", params).fetchone()[0]
    rows = [dict(r) for r in db.execute(selection + " ORDER BY " + order + " LIMIT ? OFFSET ?",
                                         [*params, limit, offset])]
    return {"rows":rows, "total":total, "limit":limit, "offset":offset,
            "next_offset":offset+len(rows) if offset+len(rows)<total else None}


def requested_ids(values, key, maximum):
    require(isinstance(values, list) and 1 <= len(values) <= maximum and
            all(isinstance(x, str) and 0 < len(x) <= 512 for x in values) and len(set(values))==len(values),
            "invalid_ids", f"{key} requires 1–{maximum} distinct, nonempty exact IDs.")
    return values


def locations(db, table, key, dataset_id, values):
    sql = f"SELECT {key},ordinal,row_offset FROM {table} WHERE dataset_id=? AND {key} IN (" + ",".join("?" for _ in values) + ")"
    by_id = {r[key]:r for r in db.execute(sql, [dataset_id, *values])}
    missing = [val for val in values if val not in by_id]
    require(not missing, "unknown_gene" if key=="gene_id" else "unknown_cell",
            f"Unknown {key} in {dataset_id}: {missing[:10]}.")
    return [by_id[x] for x in values]


def scalar(value):
    if value in ("", "NA", "NaN", "nan", "null"):
        return None
    if value.isdigit():
        return int(value)
    try:
        number = float(value)
    except ValueError as exc:
        raise QueryError("invalid_expression", "Source contains a non-numeric expression value.") from exc
    return number if math.isfinite(number) else None


def expression(db, group, args):
    did = group["dataset_id"]
    genes = requested_ids(args.get("genes"), "genes", 64)
    cells = requested_ids(args.get("cell_ids"), "cell_ids", 128)
    gene_records = locations(db, "genes", "gene_id", did, genes)
    cell_records = locations(db, "cells", "cell_id", did, cells)
    primary = group["orientation"]=="cell_by_gene"
    require(primary or group["orientation"]=="gene_by_cell", "invalid_registry","Unknown matrix orientation.")
    record_ids, other, key, shift, delimiter = ((cell_records, gene_records, "cell_id", 2, "\t") if primary
                                                 else (gene_records, cell_records, "gene_id", 1, ","))
    values = []
    with open(group["matrix_path"], "rb") as stream:
        stat = os.fstat(stream.fileno())
        require((stat.st_size, stat.st_mtime_ns)==(group["matrix_size"],group["matrix_mtime_ns"]),
                "source_changed","Source size or modification time changed. Rebuild the index before reading.")
        for record in record_ids:
            stream.seek(record["row_offset"])
            line = stream.readline().decode("utf-8-sig")
            fields = next(csv.reader([line], delimiter=delimiter))
            width = (group["n_genes"] if primary else group["n_cells"]) + shift
            require(len(fields)==width and fields[0]==record[key], "source_changed",
                    "Indexed row identifier or width differs from the source. Rebuild the index.")
            values.append([scalar(fields[r["ordinal"]+shift]) for r in other])
    if primary:
        values = [list(row) for row in zip(*values)]
    meta = json.loads(group["metadata_json"])
    return {"dataset_id":did, "genes":genes, "cell_ids":cells, "values":values, "axis_order":["gene","cell"],
            "expression_scale":meta["expression_scale"],"missing_values":sum(v is None for row in values for v in row),
            "matrix_sha256":meta.get("matrix_sha256")}


def query(db, op, args):
    if op=="list":
        cond, vals = filter_sql(args.get("filters",{}), DATASET_FILTERS)
        sql = "SELECT "+",".join(CATALOG_FIELDS)+" FROM datasets WHERE 1=1"+cond
        return page(db, sql, vals, "collection DESC,length(dataset_id),dataset_id", args)
    did = args.get("dataset_id")
    require(isinstance(did,str) and 0<len(did)<=128, "invalid_dataset_id",
            "Use a collection-qualified dataset_id from scDrugMap_list_datasets.")
    group=db.execute("SELECT * FROM datasets WHERE dataset_id=?", (did,)).fetchone()
    require(group is not None,"unknown_dataset","Unknown dataset_id. Use scDrugMap_list_datasets; a GEO accession may have multiple groups.")
    if op=="info":
        out={k:group[k] for k in CATALOG_FIELDS}
        meta=json.loads(group["metadata_json"])
        out.update({k:meta[k] for k in ("label_unit","mapped_cells","expression_scale","label_counts","matrix_member","archive",
                                         "mapping_source","source_study","matrix_sha256","annotation_member","mapping_sha256",
                                         "source_record","pmid","clinical_response_source","timepoint_definition","cross_slice_overlap")
                    if k in meta})
        out.update(orientation=group["orientation"],sample_mapping_available=group["n_samples"] is not None and group["n_samples"]>0)
        return out
    if op=="genes":
        prefix=args.get("prefix","")
        require(isinstance(prefix,str) and len(prefix)<=512, "invalid_prefix","prefix must be a string of at most 512 characters.")
        return page(db,"SELECT gene_id,ordinal FROM genes WHERE dataset_id=? AND substr(gene_id,1,length(?))=?",
                    [did,prefix,prefix],"ordinal",args)
    if op=="expression":
        return expression(db,group,args)
    cond, vals=filter_sql(args.get("filters",{}), CASE_FILTERS)
    if op=="cells":
        return page(db,"SELECT cell_id,label,patient_id,sample_id,timepoint,response_original,sample_regimen FROM cells WHERE dataset_id=?"+cond,
                    [did,*vals],"ordinal",args)
    if op=="samples":
        require(group["n_samples"] is not None and group["n_samples"]>0,"sample_mapping_unavailable",
                "No verified cell-to-sample mapping for this group. Cell rows remain available through scDrugMap_query_cells.")
        # A filter on a clinical label should never change the cell denominator of a biopsy.
        grouped="""SELECT sample_id, patient_id, timepoint, count(*) AS n_cells,
                   CASE WHEN count(DISTINCT label)=1 THEN min(label) END AS label,
                   CASE WHEN count(DISTINCT response_original)=1 THEN min(response_original) END AS response_original,
                   CASE WHEN count(DISTINCT sample_regimen)=1 THEN min(sample_regimen) END AS sample_regimen
                   FROM cells WHERE dataset_id=? AND sample_id IS NOT NULL
                   GROUP BY patient_id,sample_id,timepoint"""
        return page(db,"SELECT * FROM ("+grouped+") WHERE 1=1"+cond,[did,*vals],"patient_id,sample_id,timepoint",args)
    raise QueryError("unknown_operation","Unknown registered operation.")


@register_tool("ScDrugMapQueryTool")
class ScDrugMapQueryTool(BaseTool):
    def run(self, arguments):
        try:
            require(isinstance(arguments,dict),"invalid_arguments","Arguments must be a JSON object.")
            op=self.tool_config["fields"]["operation"]
            allowed=self.tool_config["parameter"]["properties"]
            require(all(k in allowed for k in arguments),"invalid_arguments","Unexpected argument for this tool.")
            require(arguments.get("operation",op)==op,"invalid_arguments","Operation does not match the registered tool.")
            deployment=os.environ.get("TOOLUNIVERSE_SCDRUGMAP_DATA")
            default=Path(__file__).resolve().parents[4]/"dataset"/"scDrugMap"/"query_data"
            db_file=(Path(deployment) if deployment else default)/"registry.sqlite"
            require(db_file.is_file(),"snapshot_unavailable","Set TOOLUNIVERSE_SCDRUGMAP_DATA to the directory containing registry.sqlite.")
            with closing(sqlite3.connect(db_file.resolve().as_uri()+"?mode=ro",uri=True)) as db:
                db.row_factory=sqlite3.Row
                require(db.execute("PRAGMA user_version").fetchone()[0]==1, "unsupported_snapshot",
                        "This reader requires index format version 1.")
                return {"status":"success","data":query(db,op,arguments)}
        except QueryError as exc:
            return {"status":"error","code":exc.code,"error":str(exc)}
        except (OSError, sqlite3.DatabaseError, UnicodeError, csv.Error, json.JSONDecodeError) as exc:
            return {"status":"error","code":"snapshot_read_error","error":f"Could not read the registered snapshot ({type(exc).__name__})."}
