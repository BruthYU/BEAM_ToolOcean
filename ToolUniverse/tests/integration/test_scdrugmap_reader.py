"""ToolUniverse integration tests for the registered scDrugMap reader.

For the additional release-data tests, set TEST_SCDRUGMAP_DATA to the prepared
scDrugMap index directory. These tests read only local files and use no network.
"""
import asyncio
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from jsonschema import Draft202012Validator
from tooluniverse import ToolUniverse

NAMES = {
    'list_datasets', 'get_dataset_info', 'query_samples', 'query_cells',
    'list_genes', 'get_expression',
}
ANTI_PD1 = 'external:Data20'


def build_fixture(root):
    """Two matrix orientations, shared IDs, a clinical unknown and a real-valued readout."""
    source = root / 'collection'
    primary = source / 'extracted/lxndt_filter'
    external = source / 'extracted/ext_data'
    primary.mkdir(parents=True)
    external.mkdir()
    a = primary / 'replicate_countsMatrix.txt'
    a.write_text('Cell_barcode\tCondition\tA_1\tAX1\tX\n'
                 'B\tresistant\t2\t7\t0\nA\tsensitive\t9\t0\t4\n')
    b = external / 'replicate.rawcount.csv'
    b.write_text(',C.1,A\nX,0,2.5\nA_1,0.125,4\nAX1,NA,9\n')
    annotation = external / 'replicate.annotation.csv'
    annotation.write_text(',condition\nA,sensitive\nC-1,persistent\n')
    descriptors = []
    for kind, matrix, counts in [
        ('primary', a, {'resistant': 1, 'sensitive': 1}),
        ('external', b, {'persistent': 1, 'sensitive': 1}),
    ]:
        descriptor = dict(dataset_id=kind + ':Data1', collection=kind, accession='GSE_TOY',
                          subgroup='replicate', source='patients', tissue='Tumor tissue',
                          cancer='Melanoma', regimen='anti-PD1', treatment_type='Immunotherapy',
                          platform='Smart-seq2', n_cells=2, n_genes=3, label_counts=counts,
                          matrix=str(matrix.relative_to(source/'extracted')))
        if kind == 'external':
            descriptor['annotation'] = str(annotation.relative_to(source/'extracted'))
        descriptors.append(descriptor)
    catalog = root / 'catalog.json'
    catalog.write_text(json.dumps(descriptors))
    script = Path(__file__).resolve().parents[2] / 'scripts/prepare_scdrugmap.py'
    if not script.is_file():
        script = Path(__file__).with_name('prepare_scdrugmap.py')
    index = root / 'query_data'
    subprocess.run([sys.executable, str(script), '--source', str(source), '--catalog', str(catalog),
                    '--output', str(index)], check=True, capture_output=True)
    return index


class ReaderAssertions:
    def open_universe(self, index):
        env = patch.dict(os.environ, TOOLUNIVERSE_SCDRUGMAP_DATA=str(index))
        env.start()
        self.addCleanup(env.stop)
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.tu = ToolUniverse(workspace=directory.name, load_workspace=False)
        self.addCleanup(self.tu.close)
        self.tu.load_tools(categories=['scdrugmap'], clear_existing=True)
        loaded = set(self.tu.all_tool_dict)
        self.assertEqual({name for name in loaded if name.startswith('scDrugMap_')},
                         {'scDrugMap_' + name for name in NAMES})

    def call(self, tool, **kwargs):
        name = 'scDrugMap_' + tool
        spec = self.tu.all_tool_dict[name]
        Draft202012Validator(spec['parameter']).validate(kwargs)
        response = self.tu.run({'name': name, 'arguments': kwargs}, use_cache=False, verbose=False)
        Draft202012Validator(spec['return_schema']).validate(response)
        return response

    def data(self, tool, **kwargs):
        result = self.call(tool, **kwargs)
        self.assertEqual(result.get('status'), 'success', result)
        return result['data']


class TestLocalMatrices(ReaderAssertions, unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.index = build_fixture(Path(temp.name))
        self.open_universe(self.index)

    def test_matrix_orientation_namespace_and_order(self):
        primary = self.data('get_expression', dataset_id='primary:Data1',
                            genes=['X', 'A_1'], cell_ids=['A', 'B'])
        external = self.data('get_expression', dataset_id='external:Data1',
                             genes=['X', 'A_1'], cell_ids=['A', 'C-1'])
        self.assertEqual(primary['values'], [[4, 0], [9, 2]])
        self.assertEqual(external['values'], [[2.5, 0], [4, 0.125]])
        self.assertEqual(external['cell_ids'], ['A', 'C-1'])
        self.assertEqual(external['genes'], ['X', 'A_1'])
        missing = self.data('get_expression', dataset_id='external:Data1', genes=['AX1'], cell_ids=['C-1'])
        self.assertEqual((missing['values'], missing['missing_values']), ([[None]], 1))

    def test_pagination_and_filters(self):
        ids = set()
        for offset in (0, 1):
            data = self.data('list_datasets', limit=1, offset=offset)
            self.assertEqual(data['total'], 2)
            self.assertEqual(data['next_offset'], 1 if offset == 0 else None)
            ids.add(data['rows'][0]['dataset_id'])
        self.assertEqual(ids, {'primary:Data1', 'external:Data1'})
        page = self.data('query_cells', dataset_id='external:Data1', filters={'label':['persistent']})
        self.assertEqual(page['total'], 1)
        self.assertEqual(page['rows'][0]['cell_id'], 'C-1')
        self.assertIsNone(page['rows'][0]['patient_id'])
        genes = self.data('list_genes', dataset_id='external:Data1', prefix='A_')
        self.assertEqual(genes['total'], 1)
        self.assertEqual(genes['rows'][0]['gene_id'], 'A_1')

    def test_unknown_sample_mapping_and_invalid_identifiers(self):
        error = self.call('query_samples', dataset_id='primary:Data1')
        self.assertEqual(error.get('code'), 'sample_mapping_unavailable', error)
        for values in [
                {'dataset_id':'primary:Data1', 'genes':['MISSING'], 'cell_ids':['A']},
                {'dataset_id':'external:Data1', 'genes':['X'], 'cell_ids':['MISSING']},
                {'dataset_id':'Data1', 'genes':['X'], 'cell_ids':['A']}]:
            with self.subTest(values=values):
                result = self.call('get_expression', **values)
                self.assertEqual(result.get('status'), 'error', result)

    def test_changed_source_rejects_indexed_lookup(self):
        matrix = self.index.parent / 'collection/extracted/ext_data/replicate.rawcount.csv'
        with matrix.open('ab') as handle:
            handle.write(b'outdated,1,2\n')
        error = self.call('get_expression', dataset_id='external:Data1', genes=['X'], cell_ids=['A'])
        self.assertEqual(error.get('code'), 'source_changed', error)

    def test_non_coding_mcp_client_can_discover_and_read_expression(self):
        from fastmcp import Client
        from tooluniverse.smcp import SMCP

        async def request():
            server = SMCP(tooluniverse_config=self.tu, tool_categories=['scdrugmap'],
                          search_enabled=False, auto_expose_tools=True)
            try:
                async with Client(server) as client:
                    catalog = {tool.name: tool for tool in await client.list_tools()}
                    for short in NAMES:
                        name = 'scDrugMap_' + short
                        discovered = catalog[name].inputSchema
                        expected = self.tu.all_tool_dict[name]['parameter']
                        self.assertEqual(set(discovered['properties']), set(expected['properties']))
                        self.assertEqual(set(discovered.get('required', [])), set(expected['required']))
                    result = await client.call_tool('scDrugMap_get_expression', {
                        'dataset_id': 'external:Data1', 'genes': ['X'], 'cell_ids': ['A']})
                    self.assertFalse(result.is_error)
                    reply = json.loads(result.content[0].text)
                    self.assertEqual(reply['data']['values'], [[2.5]])
            finally:
                await server.close()

        asyncio.run(request())


@unittest.skipUnless(os.environ.get('TEST_SCDRUGMAP_DATA'),
                     'Set TEST_SCDRUGMAP_DATA to run independently verified release-data checks')
class TestPublishedCollection(ReaderAssertions, unittest.TestCase):
    def setUp(self):
        self.open_universe(Path(os.environ['TEST_SCDRUGMAP_DATA']))

    def test_dataset_identity_and_collection_dimensions(self):
        page = self.data('list_datasets', limit=2, offset=0)
        self.assertEqual(page['total'], 60)
        self.assertEqual(page['next_offset'], 2)
        next_page = self.data('list_datasets', limit=2, offset=page['next_offset'])
        self.assertTrue({r['dataset_id'] for r in page['rows']}.isdisjoint(
                        r['dataset_id'] for r in next_page['rows']))
        for dataset_id, cells in [('primary:Data1', 1529), ('external:Data1', 4754),
                                  (ANTI_PD1, 11652)]:
            with self.subTest(dataset_id=dataset_id):
                info = self.data('get_dataset_info', dataset_id=dataset_id)
                self.assertEqual((info['dataset_id'], info['n_cells']), (dataset_id, cells))

    def test_source_labels_and_missing_patient_mapping(self):
        cells = self.data('query_cells', dataset_id='external:Data6',
                          filters={'label': ['persistent']}, limit=2)
        self.assertEqual(cells['total'], 673)
        self.assertEqual({r['label'] for r in cells['rows']}, {'persistent'})
        missing = self.call('query_samples', dataset_id='primary:Data1')
        self.assertEqual(missing['code'], 'sample_mapping_unavailable', missing)

    def test_baseline_patients_do_not_count_cells_as_patients(self):
        biopsies = self.data('query_samples', dataset_id=ANTI_PD1,
                             filters={'timepoint': ['pre']}, limit=50)
        rows = biopsies['rows']
        self.assertEqual(biopsies['total'], 12)
        self.assertEqual(len({r['patient_id'] for r in rows}), 12)
        self.assertEqual(sum(r['n_cells'] for r in rows), 3795)
        self.assertEqual(sum(r['label'] == 'sensitive' for r in rows), 4)
        cells = self.data('query_cells', dataset_id=ANTI_PD1,
                          filters={'timepoint': ['pre']}, limit=1)
        self.assertEqual(cells['total'], 3795)

    def test_stored_decimals_and_requested_axis_order(self):
        genes = ['NKG7', 'TCF7']
        cells = ['A10_P4_M15', 'A10_P1_M15']
        result = self.data('get_expression', dataset_id=ANTI_PD1, genes=genes, cell_ids=cells)
        self.assertEqual(result['genes'], genes)
        self.assertEqual(result['cell_ids'], cells)
        self.assertEqual(result['values'], [[12.46, 10.38], [0, 2.67]])
        self.assertTrue(result['expression_scale'])
        failed = self.call('get_expression', dataset_id=ANTI_PD1, genes=['MISSING_GENE_QA'], cell_ids=cells)
        self.assertEqual(failed.get('status'), 'error', failed)

    def test_bcc_technical_libraries_and_multiple_biopsies(self):
        info = self.data('get_dataset_info', dataset_id='external:Data17')
        self.assertEqual(info['label_unit'], 'patient_clinical_response')
        self.assertEqual((info['n_cells'], info['n_samples'], info['n_patients']), (39134, 18, 8))
        pre = self.data('query_samples', dataset_id='external:Data17', filters={'timepoint':['pre']})
        self.assertEqual(pre['total'], 9)
        self.assertEqual(sum(row['n_cells'] for row in pre['rows']), 14019)
        su001 = [row for row in pre['rows'] if row['patient_id']=='su001']
        self.assertEqual({r['sample_id'] for r in su001}, {'su001_pre1', 'su001_pre2'})
        self.assertEqual(sum(r['n_cells'] for r in su001), 1149)
        self.assertEqual({r['response_original'] for r in su001}, {'Yes'})


if __name__ == '__main__':
    unittest.main()
