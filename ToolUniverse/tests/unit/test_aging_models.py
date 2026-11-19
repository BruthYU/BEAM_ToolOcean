"""Contract failures and numerical parity are tested separately from real model smoke."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import anndata as ad

from tooluniverse.aging_models_tool import RNAAgeCalcTool, ScAgeClockTool, SenCIDTool, _bundle
from tooluniverse.exceptions import ToolValidationError


def bundle(tmp_path):
    root=tmp_path/'rna';root.mkdir()
    path=root/'coefficients.csv'
    pd.DataFrame({'gene':['A','B'],'coefficient':[2.,-1.]}).to_csv(path,index=False)
    doc={'model':'RNAAgeCalc','source':'test fixture, not a published model','tissue':'all_tissues','intercept':5.,'files':{'coefficients':{'path':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}}}
    (root/'manifest.json').write_text(json.dumps(doc))
    return root


def request(tmp_path, values, genes=('A','B')):
    path=tmp_path/'data.h5ad'
    data=ad.AnnData(np.asarray(values),obs=pd.DataFrame(index=['001','NA']),var=pd.DataFrame(index=list(genes)))
    data.write_h5ad(path)
    return {'input_path':str(path),'output_dir':str(tmp_path),'model_dir':str(bundle(tmp_path)),'species':'human','expression_type':'FPKM','gene_id_type':'SYMBOL','tissue':'all_tissues','data_kind':'bulk'}


def test_rna_exact_formula_and_id_preservation(tmp_path):
    r=request(tmp_path,[[3.,7.],[7.,3.]])
    result=object.__new__(RNAAgeCalcTool).predict_age(r)
    values=pd.read_csv(result['predictions']['path'],dtype={'sample_id':str},keep_default_na=False)
    assert values.sample_id.tolist()==['001','NA']
    np.testing.assert_allclose(values.predicted_age,[6.,9.])


def test_rna_incomplete_gene_coverage_fails(tmp_path):
    r=request(tmp_path,[[3.],[7.]],genes=['A'])
    with pytest.raises(ToolValidationError,match='coverage'):
        object.__new__(RNAAgeCalcTool).predict_age(r)


def test_model_bundle_checksum_rejects_corruption(tmp_path):
    root=bundle(tmp_path)
    (root/'coefficients.csv').write_text('corrupted')
    with pytest.raises(ToolValidationError,match='checksum'):
        _bundle(root,'RNAAgeCalc')


def test_rna_rejects_tissue_mismatch(tmp_path):
    r=request(tmp_path,[[3.,7.],[7.,3.]])
    r['tissue']='brain'
    with pytest.raises(ToolValidationError,match='tissue'):
        object.__new__(RNAAgeCalcTool).predict_age(r)


def test_sencid_rejects_normalized_expression(tmp_path):
    with pytest.raises(ToolValidationError,match='raw counts'):
        object.__new__(SenCIDTool).score_senescence({'species':'human','expression_type':'lognorm','gene_id_type':'SYMBOL'})


def test_models_reject_cross_species(tmp_path):
    for tool,operation in [(RNAAgeCalcTool,'predict_age'),(ScAgeClockTool,'predict_age'),(SenCIDTool,'score_senescence')]:
        with pytest.raises(ToolValidationError,match='human only'):
            getattr(object.__new__(tool),operation)({'species':'mouse'})
