"""Run on a compute worker: export published model artifacts and small official fixtures.

Usage: python prepare_model_bundles.py --upstream DIR --output DIR
Inputs are downloaded by the development machine. No network access is needed.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import shutil


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest(root, model, source, files, **metadata):
    doc = dict(model=model, source=source, files={key: {'path': str(Path(path).relative_to(root)), 'sha256': sha(path)} for key, path in files.items()}, **metadata)
    (root/'manifest.json').write_text(json.dumps(doc, indent=2)+'\n')


def prepare_rna(upstream, output):
    import numpy as np
    import pandas as pd
    import rdata
    import rdata.conversion._conversion as conversion
    import anndata as ad
    original = conversion.convert_array
    # Preserve R named coefficient vectors that default rdata drops the names of.
    def named_array(obj, attrs=None):
        values = original(obj, attrs)
        if attrs and 'names' in attrs and 'dim' not in attrs:
            return pd.Series(np.asarray(values), index=attrs['names'])
        return values
    conversion.convert_array = named_array
    source = upstream/'RNAAgeCalc-master'
    raw = rdata.read_rda(source/'R/sysdata.rda')
    root = output/'RNAAgeCalc'
    root.mkdir(parents=True, exist_ok=True)
    coeff = raw['genelist_all']['all_tissues']['GTExAge']
    if not isinstance(coeff, pd.Series):
        raise TypeError(f'Coefficient names lost: {type(coeff)}')
    path = root/'coefficients.csv'
    pd.DataFrame({'gene': coeff.index[1:], 'coefficient': coeff.values[1:]}).to_csv(path, index=False)
    shutil.copy2(source/'R/sysdata.rda', root/'sysdata.rda')
    shutil.copy2(source/'R/predict_age.R', root/'predict_age.R')
    manifest(root, 'RNAAgeCalc', 'https://github.com/reese3928/RNAAgeCalc',
             {'coefficients': path, 'original_coefficients': root/'sysdata.rda', 'upstream_predictor':root/'predict_age.R'},
             tissue='all_tissues', signature='GTExAge', stype='all', intercept=float(coeff.iloc[0]),
             license='GPL-2', source_archive_sha256=sha(upstream/'RNAAgeCalc.tar.gz'))
    fpkm = rdata.read_rda(source/'data/fpkmExample.RData')['fpkm']
    if not isinstance(fpkm, pd.DataFrame):
        fpkm = fpkm.to_pandas()
    fixture = ad.AnnData(fpkm.T.astype(float))
    fixture.write_h5ad(output/'rnaagecalc_example.h5ad')
    # Independent direct vector application corresponding to the R apply() expression.
    matrix = fpkm.loc[coeff.index[1:]].to_numpy()
    reference = [sum(np.log2(matrix[:, i]+1)*coeff.values[1:])+coeff.iloc[0] for i in range(matrix.shape[1])]
    pd.DataFrame({'sample_id':fpkm.columns.astype(str), 'predicted_age':reference}).to_csv(output/'rnaagecalc_reference.csv', index=False)
    conversion.convert_array = original


def prepare_sencid(upstream, output):
    import pandas as pd
    import anndata as ad
    source = upstream/'SenCID-main/SenCID'
    root = output/'SenCID'
    root.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(source/'resource/seneset.txt', header=None, sep='\t')[0].astype(str).tolist()
    (root/'features.json').write_text(json.dumps(features))
    tree = ast.parse((source/'Pred.py').read_text())
    assignment = next(node for node in ast.walk(tree) if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id=='feature_dict' for t in node.targets))
    (root/'sid_features.json').write_text(json.dumps(ast.literal_eval(assignment.value)))
    shutil.copy2(upstream/'SenCID-main/LICENSE',root/'LICENSE')
    files = {'features':root/'features.json','sid_features':root/'sid_features.json','license':root/'LICENSE'}
    for path in (source/'model').glob('*.pkl'):
        dest = root/path.name
        shutil.copy2(path, dest)
        files[path.stem] = dest
    manifest(root, 'SenCID', 'https://github.com/JackieHanLab/SenCID', files,
             license='MIT', denoising=False, source_archive_sha256=sha(upstream/'SenCID.tar.gz'))
    fixture = pd.read_csv(source/'demo/demo/origin_matrix_GSE94980.txt', sep='\t', index_col=0).iloc[:, :24]
    ad.AnnData(fixture.T.astype(float)).write_h5ad(output/'sencid_example.h5ad')


def prepare_scage(upstream, output):
    import pandas as pd
    import anndata as ad
    root = output/'scAgeClock'
    root.mkdir(parents=True, exist_ok=True)
    features = pd.read_csv(upstream/'h5ad_var.tsv', sep='\t')['h5ad_var'].astype(str).tolist()
    (root/'features.json').write_text(json.dumps(features))
    # Retain the exact official inference classes, excluding training-only relative imports.
    source = (upstream/'GatedMultiheadAttention.py').read_text()
    tree = ast.parse(source)
    selected = [ast.get_source_segment(source, node) for node in tree.body if isinstance(node, ast.ClassDef) and node.name in ('GatedMultiheadAttentionFCNet','FeatureGate','MultiheadAttention')]
    (root/'gma_inference.py').write_text('import torch\nfrom torch import nn\nimport numpy as np\n\n'+'\n\n'.join(selected)+'\n')
    shutil.copy2(upstream/'GatedMultiheadAttention.py',root/'upstream_architecture.py')
    shutil.copy2(upstream/'scAgeClock_GMA_model_state_dict.pth', root/'weights.pth')
    mapping = {}
    for name in features[:4]:
        frame = pd.read_csv(upstream/f'{name}_numeric_index.tsv',sep='\t')
        mapping[name] = dict(zip(frame['categorical_value'].astype(str), frame['numeric_index'].astype(int)))
    (root/'categories.json').write_text(json.dumps(mapping))
    shutil.copy2(upstream/'scage_LICENSE',root/'LICENSE')
    files = {'license':root/'LICENSE','features':root/'features.json','architecture_source':root/'gma_inference.py','original_architecture':root/'upstream_architecture.py','weights':root/'weights.pth','categories':root/'categories.json'}
    manifest(root,'scAgeClock','https://github.com/gangcai/scageclock',files,
             license='BSD-3-Clause',architecture={'cat_cardinalities':[21,664,52,3],'num_numeric_features':len(features)-4,'projection_dim':512,'prediction_hidden_sizes':[256,128],'num_heads':8})
    fixture = ad.read_h5ad(upstream/'scage_example.h5ad')[:32].copy()
    if fixture.n_vars != len(features):
        raise ValueError('Official fixture width differs from official feature metadata.')
    fixture.var_names = features
    fixture.write_h5ad(output/'scageclock_example.h5ad')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--upstream',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--models',nargs='+',default=['rna','sencid','scage'])
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    for model in args.models:
        globals()['prepare_'+model](args.upstream,args.output)
        print(f'Prepared {model}',flush=True)
