"""Read-only queries over a versioned, locally deployed IMvigor210 snapshot."""
import json
import os
from pathlib import Path
import sqlite3
from .base_tool import BaseTool
from .tool_registry import register_tool

@register_tool("IMvigor210QueryTool")
class IMvigor210QueryTool(BaseTool):
    def run(self, arguments):
        try:
            return self._run(arguments)
        except Exception as exc:
            return {"status": "error", "error": str(exc)}

    def _run(self, args):
        import numpy as np
        directory = os.environ.get("TOOLUNIVERSE_IMVIGOR210_DATA")
        if not directory:
            raise ValueError("Set TOOLUNIVERSE_IMVIGOR210_DATA to the prepared query_data directory")
        root = Path(directory).resolve()
        manifest = json.loads((root / "manifest.json").read_text())
        operation = self.tool_config["fields"]["operation"]
        if operation == "info":
            return {"status": "success", "data": manifest}
        con = sqlite3.connect((root / "imvigor210.sqlite").as_uri() + "?mode=ro", uri=True)
        con.row_factory = sqlite3.Row
        try:
            fields = manifest["fields"]
            def quote(field):
                if field not in fields:
                    raise ValueError("Unknown field: " + str(field))
                return '"' + field.replace('"', '""') + '"'
            clauses, values = [], []
            filters = args.get("filters", {})
            if not isinstance(filters, dict):
                raise ValueError("filters must be an object of field: [accepted values]")
            for field, accepted in filters.items():
                column = quote(field)
                if not isinstance(accepted, list) or not accepted or len(accepted) > 500:
                    raise ValueError("Filter values must be nonempty lists of at most 500 values")
                nonnull = [x for x in accepted if x is not None]
                if any(not isinstance(x, (str, int, float)) or isinstance(x, bool) for x in nonnull):
                    raise ValueError("Filter values must be strings, numbers or null")
                alternatives = []
                if nonnull:
                    alternatives.append(column + " IN (" + ",".join("?" for _ in nonnull) + ")")
                    values.extend(nonnull)
                if None in accepted:
                    alternatives.append(column + " IS NULL")
                clauses.append("(" + " OR ".join(alternatives) + ")")
            where = " WHERE " + " AND ".join(clauses) if clauses else ""
            selected = [dict(x) for x in con.execute('SELECT * FROM samples' + where + ' ORDER BY sample_id', values)]
            common = {"dataset_version": manifest["version"], "matched_samples": len(selected),
                      "matched_patients": len({x['patient_id'] for x in selected if x['patient_id'] is not None})}
            if operation == "samples":
                requested = args.get("fields", ["sample_id", "patient_id", "response", "immune_phenotype"])
                if not requested or len(requested) != len(set(requested)):
                    raise ValueError("fields must be nonempty and unique")
                for field in requested: quote(field)
                limit, offset = args.get("limit", 50), args.get("offset", 0)
                if type(limit) is not int or not 1 <= limit <= 500 or type(offset) is not int or offset < 0:
                    raise ValueError("limit must be 1..500 and offset a nonnegative integer")
                page = selected[offset:offset+limit]
                data = dict(common, rows=[{k:r[k] for k in requested} for r in page], returned_samples=len(page),
                            next_offset=offset+limit if offset+limit < len(selected) else None)
            elif operation == "summary":
                group = args.get("group_by", "response"); quote(group)
                groups = {}
                for row in selected:
                    value = row[group]
                    key = json.dumps(value)
                    entry = groups.setdefault(key, {"value": value, "samples": 0, "patient_ids": set()})
                    entry['samples'] += 1
                    if row['patient_id'] is not None: entry['patient_ids'].add(row['patient_id'])
                data = dict(common, group_by=group, groups=[{"value":r['value'], "samples":r['samples'], "patients":len(r['patient_ids'])} for r in groups.values()],
                            missing_values={f:sum(r[f] is None for r in selected) for f in fields},
                            note="Patient counts across groups may overlap. No statistical association or survival estimate is computed.")
            elif operation == "expression":
                if args.get("expression_type", "raw_counts") != "raw_counts":
                    raise ValueError("Only raw_counts is available; counts are not TPM or FPKM")
                genes, ids = args.get("genes"), args.get("sample_ids")
                if not isinstance(genes, list) or not 1 <= len(genes) <= 100 or any(not isinstance(g,str) or not g for g in genes):
                    raise ValueError("genes must contain 1..100 nonempty gene IDs or exact symbols")
                if not isinstance(ids, list) or not 1 <= len(ids) <= 500 or any(not isinstance(i,str) for i in ids) or len(ids)!=len(set(ids)):
                    raise ValueError("sample_ids must contain 1..500 unique sample IDs")
                all_samples = [x[0] for x in con.execute('SELECT sample_id FROM samples ORDER BY rowid')]
                index = {s:i for i,s in enumerate(all_samples)}
                unknown = [s for s in ids if s not in index]
                if unknown: raise ValueError("Unknown sample IDs: " + ', '.join(unknown))
                allowed = {s['sample_id'] for s in selected}
                if any(s not in allowed for s in ids): raise ValueError("Some sample_ids do not satisfy filters")
                try:
                    matrix = np.load(root/'counts.npy', mmap_mode='r', allow_pickle=False)
                except OSError:
                    matrix = np.load(root/'counts.npy', allow_pickle=False)
                if matrix.shape != (manifest['genes'], manifest['samples']): raise ValueError("Snapshot matrix dimensions do not match manifest")
                result, missing = [], []
                for gene in dict.fromkeys(genes):
                    found = con.execute('SELECT * FROM genes WHERE gene_id=?', (gene,)).fetchall()
                    if not found: found = con.execute('SELECT * FROM genes WHERE symbol=?', (gene,)).fetchall()
                    if not found: missing.append(gene); continue
                    for row in found:
                        result.append({"query":gene,"gene_id":row['gene_id'],"symbol":row['symbol'],
                                       "values":[int(matrix[row['gene_index'],index[s]]) for s in ids]})
                data = dict(common, expression_type='raw_counts', units='un-normalized gene read counts', sample_ids=ids,
                            genes=result, missing_genes=missing, note="Ambiguous symbols return all matching rows; values follow sample_ids order.")
            else:
                raise ValueError("Unknown operation")
            return {"status":"success", "data":data}
        finally:
            con.close()
