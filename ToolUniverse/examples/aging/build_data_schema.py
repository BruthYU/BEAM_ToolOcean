"""Generate the data-tool schema fragment used to assemble aging_tools.json."""
import json
from pathlib import Path
S={'type':'string','minLength':1}
A={'type':'array','items':S,'minItems':1,'uniqueItems':True}
I={'type':'integer','minimum':0}
B={'type':'boolean'}
def enum(*x):return {'type':'string','enum':list(x)}
def branch(op,props,required):
    return {'type':'object','properties':{'operation':{'type':'string','enum':[op]},**props},'required':['operation',*required],'additionalProperties':False}
artifact={'type':'object','properties':{'path':S,'shape':{'type':'array','items':{'type':'integer','minimum':0}}},'required':['path','shape']}
outputs={
'inspect_dataset':{'shape':{'type':'array'},'obs_columns':{'type':'array'},'var_columns':{'type':'array'},'layers':{'type':'array'},'integer_valued':B,'recorded_metadata':{'type':'object'},'zero_total_observations':I},
'import_expression':{'expression':artifact},'attach_metadata':{'expression':artifact,'extra_metadata_rows':I},
'export_results':{'table':artifact},'filter_qc':{'expression':artifact,'audit':artifact,'retained':I,'excluded':I},
'align_gene_ids':{'expression':artifact,'audit':artifact},
'normalize_expression':{'expression':artifact,'layer':S,'target_sum':{'type':'number','exclusiveMinimum':0},'method':S},
'prepare_model_input':{'expression':artifact,'coverage':{'type':['object','null']},'preprocessing':S},
'aggregate_counts':{'expression':artifact,'table':artifact,'metadata':artifact},'summarize_cell_counts':{'table':artifact}}
def config(name,ops):
    return {'name':name,'type':name+'Tool','description':{'AgingIO':'Inspect, import, annotate and export local expression data with explicit biological metadata.','AgingPreprocess':'Prepare model-specific aging inputs without silently changing gene IDs or expression scales.','AgingPseudobulk':'Aggregate raw single-cell counts within donor, sample and cell type; preserve biological replicate identities.'}[name],
    'parameter':{'type':'object','properties':{'request':{'oneOf':ops}},'required':['request'],'additionalProperties':False},
    'return_schema':{'oneOf':[{'type':'object','properties':{'operation':{'const':op['properties']['operation']['enum'][0]},'provenance':{'type':'object'},**outputs[op['properties']['operation']['enum'][0]]},'required':['operation','provenance',*outputs[op['properties']['operation']['enum'][0]]]} for op in ops]}}
base={'input_path':S,'output_dir':S}; layer={'layer':S}
ops=[branch('inspect_dataset',{'input_path':S,**layer},['input_path']),
branch('import_expression',{**base,'orientation':enum('observations_by_genes','genes_by_observations'),'species':S,'tissue':S,'expression_scale':enum('counts','fpkm','tpm','cpm','log_normalized'),'gene_id_type':S,'data_kind':enum('bulk','single_cell','pseudobulk'),'metadata_path':S,'id_key':S},[*base,'orientation','species','tissue','expression_scale','gene_id_type','data_kind']),
branch('attach_metadata',{**base,'metadata_path':S,'id_key':S},[*base,'metadata_path','id_key']),
branch('export_results',{**base,**layer,'obs_columns':A,'genes':A},list(base))]
configs=[config('AgingIO',ops)]
ops=[branch('filter_qc',{**base,**layer,'min_counts':I,'min_genes':I,'max_counts':I},list(base)),
branch('align_gene_ids',{**base,'mapping_path':S,'target_gene_id_type':S,'strip_ensembl_version':B},list(base)),
branch('normalize_expression',{**base,**layer,'target_sum':{'type':'number','exclusiveMinimum':0},'output_layer':S},list(base)),
branch('prepare_model_input',{**base,**layer,'model':enum('RNAAgeCalc','scAgeClock','SenCID'),'species':enum('human'),'expression_scale':enum('counts','fpkm','tpm','cpm'),'data_kind':enum('bulk','single_cell','pseudobulk'),'feature_path':S,'min_gene_coverage':{'type':'number','minimum':0,'maximum':1},'categorical_keys':A,'gene_length_key':S},[*base,'model','species','expression_scale','data_kind'])]
configs.append(config('AgingPreprocess',ops))
props={**base,**layer,'donor_key':S,'sample_key':S,'cell_type_key':S,'metadata_columns':A,'min_cells':{'type':'integer','minimum':1}}
configs.append(config('AgingPseudobulk',[branch(op,props,list(base)) for op in ('aggregate_counts','summarize_cell_counts')]))
Path(__file__).with_name('schema_data.json').write_text(json.dumps(configs,indent=2)+'\n')
