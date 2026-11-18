"""Execute real published models on their author-provided fixtures on a worker."""
import argparse, ast, importlib.util, json, sys, types
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from tooluniverse.aging_models_tool import RNAAgeCalcTool, ScAgeClockTool, SenCIDTool


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def run_models(args,engine):
    from jsonschema import validate
    configs=engine.all_tool_dict
    def invoke(name, request):
        result=engine.run_one_function({'name':name,'arguments':{'request':request}})
        assert result.get('operation')==request['operation'],result
        validate(result,configs[name]['return_schema'])
        return result
    results={}
    root=args.bundles;out=args.output;out.mkdir(parents=True,exist_ok=True)
    specifications=[('rna','RNAAgeCalc',RNAAgeCalcTool,'rnaagecalc_example.h5ad',dict(operation='predict_age',expression_type='FPKM',gene_id_type='SYMBOL',tissue='all_tissues',data_kind='bulk')),
                    ('sencid','SenCID',SenCIDTool,'sencid_example.h5ad',dict(operation='score_senescence',expression_type='counts',gene_id_type='SYMBOL',min_gene_coverage=0.8)),
                    ('scage','scAgeClock',ScAgeClockTool,'scageclock_example.h5ad',dict(operation='predict_age',expression_type='scageclock_formatted',device='cuda',batch_size=8))]
    for key,name,cls,fixture,extra in specifications:
        if key not in args.models:continue
        request=dict(input_path=str(root/fixture),model_dir=str(root/name),output_dir=str(out),species='human',**extra)
        inspection=invoke('AgingIO',dict(operation='inspect_dataset',input_path=request['input_path']))
        result=invoke(name,request)
        result['input_inspection']=inspection
        artifact=result.get('predictions',result.get('scores'))
        prediction=pd.read_csv(artifact['path'],dtype=str,keep_default_na=False)
        data=ad.read_h5ad(root/fixture)
        id_key='sample_id' if key=='rna' else 'cell_id'
        assert prediction[id_key].tolist()==data.obs_names.astype(str).tolist()
        if key=='rna':
            reference=pd.read_csv(root/'rnaagecalc_reference.csv')
            np.testing.assert_allclose(prediction.predicted_age.astype(float),reference.predicted_age,rtol=1e-12,atol=1e-10)
            result['parity']='exact coefficient formula, independent scalar evaluation; not R-runtime parity'
        if key=='sencid':
            source=args.upstream/'SenCID-main/SenCID'
            # Execute the original no-DCA preprocessing after removing only the unused
            # eager dca import. No classifier/preprocessing outputs are mocked.
            tree=ast.parse((source/'DataPro.py').read_text())
            tree.body=[n for n in tree.body if not (isinstance(n,ast.ImportFrom) and n.module=='dca.api')]
            module=types.ModuleType('sencid_upstream_datapro');module.__file__=str(source/'DataPro.py')
            sys.path.insert(0,str(source))
            exec(compile(tree,str(source/'DataPro.py'),'exec'),module.__dict__)
            scaled,_=module.ScaleData(data.copy(),denoising=False)
            upstream_pred=load_module('sencid_upstream_pred',source/'Pred.py')
            for sid in range(1,7):
                reference=upstream_pred.Pred(scaled,sid)
                np.testing.assert_allclose(prediction[f'SID{sid}_score'].astype(float),reference.SID_Score,rtol=1e-10,atol=1e-10)
            result['parity']='official ScaleData(denoising=False) and Pred: six SID score vectors match'
        if key=='scage':
            import torch
            manifest=json.loads((root/name/'manifest.json').read_text())
            module=load_module('gma_reference',root/name/'gma_inference.py')
            model=module.GatedMultiheadAttentionFCNet(**manifest['architecture'])
            model.load_state_dict(torch.load(root/name/'weights.pth',map_location='cpu',weights_only=True));model.eval()
            with torch.inference_mode():
                reference=model(torch.tensor(data.X.toarray(),dtype=torch.float32)).flatten().numpy()
            np.testing.assert_allclose(prediction.predicted_age.astype(float),reference,rtol=1e-4,atol=1e-3)
            result['parity']='official GMA full-batch CPU vs wrapped batched CUDA'
            # Separate synthetic engineering fixture: this only tests formatting and
            # serial execution. It is not biological validation or study evidence.
            from scipy.sparse import csr_matrix
            features=json.loads((root/name/'features.json').read_text())
            rng=np.random.default_rng(23)
            raw=ad.AnnData(csr_matrix(rng.poisson(3,(6,len(features)-4)).astype(float)),
                          obs=pd.DataFrame(index=[f'engineering_{i}' for i in range(6)]),
                          var=pd.DataFrame(index=features[4:]))
            cats=data.X[:6,:4].toarray()
            for i,col in enumerate(features[:4]): raw.obs[col+'_index']=cats[:,i].astype(int)
            raw_path=out/'scage_engineering_raw.h5ad';raw.write_h5ad(raw_path)
            before=raw.X.copy()
            preparation=dict(operation='prepare_input',input_path=str(raw_path),model_dir=str(root/name),
                             output_dir=str(out),species='human',expression_type='counts',
                             category_columns={c:c+'_index' for c in features[:4]},categories_are_indices=True)
            prepared=invoke(name,preparation)
            formatted=ad.read_h5ad(prepared['output']['path'])
            reference_formatter=load_module('scage_original_formatter',args.upstream/'scage_formatting.py')
            reference_formatted=reference_formatter.format_anndata_multiple(raw.copy(),model_genes=features[4:],
                                        cat_cols=[c+'_index' for c in features[:4]],normalize=True)
            np.testing.assert_allclose(formatted.X.toarray(),reference_formatted.X.toarray(),rtol=1e-6,atol=1e-6)
            assert (ad.read_h5ad(raw_path).X!=before).nnz==0
            assert formatted.obs_names.tolist()==raw.obs_names.tolist()
            engineering_request=dict(request,input_path=prepared['output']['path'])
            engineering_prediction=invoke(name,engineering_request)
            validate(engineering_prediction,configs[name]['return_schema'])
            result['engineering_preparation']={'fixture':'synthetic_poisson_counts_with_official_category_indices',
                'biological_validation':False,'formatting_parity':'official format_anndata_multiple',
                'preparation':prepared,'prediction':engineering_prediction}

        results[name]=result
        (out/'model_smoke.json').write_text(json.dumps(results,indent=2))
        print(f'{name}: PASS ({len(prediction)} observations)',flush=True)

def main(args):
    from tooluniverse import ToolUniverse
    engine=ToolUniverse(workspace=str(args.output/'workspace'),load_workspace=False)
    try:
        engine.load_tools(categories=['Aging'],clear_existing=True)
        run_models(args,engine)
    finally:
        engine.close()

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--bundles',type=Path,required=True);parser.add_argument('--upstream',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--models',nargs='+',default=['rna','sencid','scage']);main(parser.parse_args())
