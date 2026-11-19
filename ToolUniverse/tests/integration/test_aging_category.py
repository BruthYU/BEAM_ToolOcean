"""Aging SDK/MCP and a real donor-level serial analysis on small count data.

The fixture is synthetic with an explicit known signal. It verifies software
execution and statistical direction, not biological model validity.
"""
import asyncio
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from jsonschema import Draft202012Validator
from tooluniverse import ToolUniverse

pytestmark = pytest.mark.integration
NAMES = {'AgingIO','AgingPreprocess','AgingPseudobulk','RNAAgeCalc','scAgeClock','SenCID',
         'AgingJointSummary','AgingAssociation','AgingComposition',
         'AgingDifferentialExpression','AgingPathway','AgingValidation'}


@pytest.fixture
def universe(tmp_path):
    tu=ToolUniverse(workspace=str(tmp_path/'workspace'),load_workspace=False)
    try:
        tu.load_tools(categories=['Aging'],clear_existing=True)
        yield tu
    finally:
        tu.close()


def call(tu,name,operation,**request):
    args={'request':{'operation':operation,**request}}
    result=tu.run_one_function({'name':name,'arguments':args})
    assert result.get('operation')==operation,result
    Draft202012Validator(tu.all_tool_dict[name]['return_schema']).validate(result)
    return result


def test_aging_discovery_and_schemas(universe,tmp_path):
    assert set(universe.all_tool_dict)==NAMES
    for config in universe.all_tool_dict.values():
        Draft202012Validator.check_schema(config['parameter'])
        Draft202012Validator.check_schema(config['return_schema'])
    error=universe.run_one_function({'name':'AgingDifferentialExpression','arguments':{'request':{
        'operation':'compare_groups','input_path':'does_not_exist','output_dir':str(tmp_path/'invalid')}}})
    assert 'error' in error
    assert not (tmp_path/'invalid').exists()


def test_aging_serial_pseudobulk_de_pathway(universe,tmp_path):
    import anndata as ad
    rng=np.random.default_rng(19)
    genes=[f'GENE{i}' for i in range(100)]
    rows=[];expression=[]
    # Six independent donors per arm, two cell types, ten cells each.
    for donor in range(12):
        baseline=rng.uniform(8,12,len(genes))
        for cell_type in ['T','B']:
            mu=baseline.copy()
            if donor>=6 and cell_type=='T':mu[:10]*=5
            for cell in range(10):
                rows.append({'cell_id':f'd{donor}_{cell_type}_{cell}','donor_id':f'd{donor:02d}',
                             'sample_id':f's{donor:02d}','cell_type':cell_type,
                             'group':'treated' if donor>=6 else 'control'})
                expression.append(rng.negative_binomial(20,20/(20+mu)))
    metadata=pd.DataFrame(rows)
    pd.DataFrame(expression,index=metadata.cell_id,columns=genes).to_csv(tmp_path/'counts.csv',index_label='cell_id')
    metadata.to_csv(tmp_path/'cells.csv',index=False)
    out=str(tmp_path/'outputs')
    imported=call(universe,'AgingIO','import_expression',input_path=str(tmp_path/'counts.csv'),
                  metadata_path=str(tmp_path/'cells.csv'),id_key='cell_id',orientation='observations_by_genes',
                  species='human',tissue='blood',expression_scale='counts',gene_id_type='symbol',data_kind='single_cell',output_dir=out)
    inspected=call(universe,'AgingIO','inspect_dataset',input_path=imported['expression']['path'])
    assert inspected['shape']==[240,100]
    qc=call(universe,'AgingPreprocess','filter_qc',input_path=imported['expression']['path'],min_counts=1,min_genes=10,output_dir=out)
    pb=call(universe,'AgingPseudobulk','aggregate_counts',input_path=qc['expression']['path'],
            donor_key='donor_id',sample_key='sample_id',cell_type_key='cell_type',metadata_columns=['group'],min_cells=5,output_dir=out)
    bulk=ad.read_h5ad(pb['expression']['path'])
    assert bulk.shape==(24,100)
    original=np.asarray(expression)
    for _,row in bulk.obs.reset_index().iterrows():
        selector=(metadata.donor_id==row.donor_id)&(metadata.cell_type==row.cell_type)
        np.testing.assert_allclose(bulk[row.pseudobulk_id].X.toarray().ravel(),original[selector].sum(axis=0))
    # The complete multi-cell-type matrix must not be treated as 24 independent donors.
    invalid=universe.run_one_function({'name':'AgingDifferentialExpression','arguments':{'request':{
        'operation':'compare_groups','input_path':pb['expression']['path'],'expression_scale':'counts',
        'compare_groups':['control','treated'],'output_dir':out}}})
    assert 'error' in invalid
    de=call(universe,'AgingDifferentialExpression','compare_groups',input_path=pb['expression']['path'],
            expression_scale='counts',cell_type='T',compare_groups=['control','treated'],output_dir=out)
    assert de['n_donors']==12
    table=pd.read_csv(de['results']['path'])
    assert table.loc[table.gene_id.isin(genes[:10]),'log2_fold_change'].median()>1
    sets={'known_induced':genes[:10],'background_program':genes[40:60]}
    gsea=call(universe,'AgingPathway','run_preranked_gsea',input_path=de['results']['path'],
              gene_sets=sets,gene_set_source='known-signal synthetic fixture',n_permutations=50,output_dir=out)
    enrichment=pd.read_csv(gsea['results']['path'])
    assert enrichment.loc[enrichment.gene_set=='known_induced','es'].item()>0
    ora=call(universe,'AgingPathway','run_overrepresentation',selected_genes=genes[:10],background_genes=genes,
             gene_sets=sets,gene_set_source='known-signal synthetic fixture',output_dir=out)
    assert Path(ora['results']['path']).exists()


def test_aging_mcp_discovery(universe):
    from fastmcp import Client
    from tooluniverse.smcp import SMCP
    async def check():
        server=SMCP(tooluniverse_config=universe,tool_categories=['Aging'],search_enabled=False,auto_expose_tools=True)
        try:
            async with Client(server) as client:
                tools={tool.name:tool for tool in await client.list_tools()}
                assert set(tools)>=NAMES
                for name in NAMES:
                    assert tools[name].inputSchema==universe.all_tool_dict[name]['parameter']
        finally:
            await server.close()
    asyncio.run(check())
