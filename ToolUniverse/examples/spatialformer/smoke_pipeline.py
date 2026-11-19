"""Tiny serial SDK smoke with real checkpoints and explicitly synthetic counts.

Run on the GPU worker. This checks interfaces, not biological accuracy.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from jsonschema import Draft202012Validator
from tooluniverse import ToolUniverse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((args.runtime/'single/manifest.json').read_text())
    vocabulary = json.loads((args.runtime/'single'/manifest['files']['tokenizer']['path']).read_text())
    genes = [g for g, token in vocabulary.items() if token >= 27 and not g.startswith('<')][:48]
    rng = np.random.default_rng(41)
    events = []
    tu = ToolUniverse(load_workspace=False)
    tu.load_tools(categories=['spatial_transcriptomics'], clear_existing=True)

    def call(name, operation, **kwargs):
        result = tu.run_one_function({'name': name, 'arguments': {'request': {
            'operation': operation, 'output_dir': str(args.output / 'results'), **kwargs}}})
        if 'error' in result:
            raise RuntimeError(result)
        Draft202012Validator(tu.all_tool_dict[name]['return_schema']).validate(result)
        events.append({'tool': name, 'operation': operation, 'result': result})
        (args.output/'events.json').write_text(json.dumps(events, indent=2))
        print('PASS', name, operation, flush=True)
        return result

    try:
        samples = []
        for i in range(4):
            raw = pd.DataFrame(rng.poisson(4, (4, len(genes))) + 1, columns=genes)
            raw.insert(0, 'cell_id', ['001', '002', '003', '004'])
            raw['x'] = [0, 5, 100, 105]
            raw['y'] = 0
            path = args.output/f'sample{i}.csv'
            raw.to_csv(path, index=False)
            imported = call('SpatialSingleCellIO', 'import_cells', input_path=str(path),
                            platform='Xenium', sample_id=f's{i}', gene_columns=genes,
                            coordinate_columns=['x', 'y'], coordinate_unit='micrometer', coordinate_frame='global')
            annotations = pd.DataFrame({'cell_id': raw.cell_id, 'cell_type': ['T', 'B', 'T', 'B']})
            ann = args.output/f'annotations{i}.csv'
            annotations.iloc[::-1].to_csv(ann, index=False)
            labeled = call('SpatialSingleCellIO', 'attach_annotations', input_path=imported['output_path'],
                           annotation_path=str(ann), label_column='cell_type', input_id_key='source_cell_id')
            samples.append(ad.read_h5ad(labeled['output_path']))
        cells = ad.concat(samples, merge='same', uns_merge='same')
        cells_path = args.output/'cells.h5ad'
        cells.write_h5ad(cells_path)
        reference = call('SpatialSingleCellIO', 'prepare_reference', input_path=str(cells_path), cell_type_key='cell_type')
        # The prepared reference is accepted by the existing cell2location reference stage.
        call('cell2location', 'fit_reference', reference_path=reference['reference_path'],
             cell_type_key='cell_type', count_layer='counts', max_epochs=1)
        common = dict(input_path=str(cells_path), count_layer='counts', assay='Xenium',
                      tissue='Lung', condition='Disease', device='cuda', batch_size=2,
                      num_workers=0, max_genes_per_cell=48)
        embedded = call('SpatialFormer', 'extract_embeddings', model_bundle_path=str(args.runtime/'single'), **common)
        z = ad.read_h5ad(embedded['output_path']).obsm['X_SpaF']
        assert z.shape[0] == len(cells) and np.isfinite(z).all()
        benchmark = call('SpatialPairValidation', 'build_pair_benchmark', input_path=str(cells_path),
                         sample_key='sample_id', frame_key='sample_id', coordinate_unit='micrometer',
                         radius=10, pairs_per_class=2, random_seed=0)
        predicted = call('SpatialFormer', 'predict_pairs', model_bundle_path=str(args.runtime/'pair'),
                         pairs_path=benchmark['pairs_path'], **common)
        call('SpatialPairValidation', 'evaluate_pairs', prediction_path=predicted['predictions_path'],
             benchmark_path=benchmark['benchmark_path'], pretraining_overlap='unknown')
        summary = call('SpatialPairSummary', 'summarize_type_pairs', prediction_path=predicted['predictions_path'],
                       input_path=str(cells_path), cell_type_key='cell_type', sample_key='sample_id',
                       frame_key='sample_id', benchmark_path=benchmark['benchmark_path'],
                       sampling_description='Synthetic smoke, two requested pairs per class and sample')
        meta = args.output/'patients.csv'
        pd.DataFrame({'sample_id':[f's{i}' for i in range(4)],'patient_id':[f'p{i}' for i in range(4)],
                      'group':['A','A','B','B']}).to_csv(meta,index=False)
        call('SpatialPairSummary','compare_groups',summary_path=summary['summary_path'],metadata_path=str(meta),
             compare_groups=['A','B'],n_permutations=19,random_seed=0)
        (args.output/'SUCCESS.json').write_text(json.dumps({'status':'passed','n_operations':len(events),
            'data':'synthetic counts, not a biological benchmark','embedding_shape':list(z.shape)}))
    finally:
        tu.close()


if __name__ == '__main__':
    main()
