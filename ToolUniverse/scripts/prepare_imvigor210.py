"""Download and prepare the official IMvigor210 1.0.1 snapshot (one-time setup)."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import tarfile
import tempfile
import urllib.request

URL = 'http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/IMvigor210CoreBiologies_1.0.1.tar.gz'
PACKAGE_SHA256 = '890b80da68a837db6edbc1d348c53de4c7c8c7e9eddf3000bf26df5ef2b4668a'
SOURCE_SHA256 = '290b48cb854a3ef0f510118f0e34693fdf7d6b88ec544fc6f6fe7dbac658e268'


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def pairs(obj):
    result = {}
    while obj is not None and obj.value:
        tag = obj.tag
        if tag.info.type.name == 'REF':
            tag = tag.referenced_object
        result[tag.value.value.decode()] = obj.value[0]
        obj = obj.value[1]
    return result


def prepare(source, output):
    import numpy as np
    import rdata
    parsed = rdata.parser.parse_file(source)
    attrs = pairs(parsed.object.value[0].attributes)
    converter = rdata.conversion.SimpleConverter()
    samples = converter.convert(pairs(attrs['phenoData'].attributes)['data'])
    genes = converter.convert(pairs(attrs['featureData'].attributes)['data'])
    counts = None
    for bucket in attrs['assayData'].value.hash_table.value:
        if bucket.info.type.name == 'LIST':
            entry = pairs(bucket)
            if 'counts' in entry:
                counts = converter.convert(entry['counts'])
    if counts is None:
        raise ValueError('Official object has no counts assay')
    matrix = np.asarray(counts.values)
    if (list(counts.coords['dim_0'].values.astype(str)) != list(genes.index.astype(str))
            or list(counts.coords['dim_1'].values.astype(str)) != list(samples.index.astype(str))):
        raise ValueError('Gene/sample ordering does not match annotations')
    if matrix.shape != (31286, 348) or not np.issubdtype(matrix.dtype, np.integer) or np.any(matrix < 0):
        raise ValueError('Unexpected raw-count matrix')
    aliases = {'ANONPT_ID': 'patient_id', 'Best Confirmed Overall Response': 'response',
               'Immune phenotype': 'immune_phenotype', 'os': 'os_months', 'censOS': 'os_event',
               'IC Level': 'pd_l1_ic', 'TC Level': 'pd_l1_tc',
               'FMOne mutation burden per MB': 'tmb_per_mb'}
    samples = samples.rename(columns=aliases).copy()
    samples.index.name = 'sample_id'
    if not samples.index.is_unique:
        raise ValueError('Duplicate sample IDs')
    genes = genes[['entrez_id', 'symbol']].copy()
    genes['entrez_id'] = genes['entrez_id'].astype(str)
    genes['gene_index'] = np.arange(len(genes))
    with sqlite3.connect(output / 'imvigor210.sqlite') as connection:
        samples.reset_index().to_sql('samples', connection, index=False)
        genes.rename(columns={'entrez_id': 'gene_id'}).to_sql('genes', connection, index=False)
    np.save(output / 'counts.npy', matrix, allow_pickle=False)
    manifest = {'dataset': 'IMvigor210', 'version': 'IMvigor210CoreBiologies_1.0.1',
                'samples': len(samples), 'patients': int(samples.patient_id.nunique()), 'genes': len(genes),
                'expression_types': ['raw_counts'], 'fields': list(samples.reset_index().columns),
                'field_mapping': aliases, 'source': URL, 'source_sha256': SOURCE_SHA256,
                'package_sha256': PACKAGE_SHA256,
                'units': {'raw_counts': 'un-normalized gene read counts', 'os_months': 'months',
                          'os_event': '0=alive,1=dead', 'tmb_per_mb': 'mutations/Mb'}}
    manifest['files_sha256'] = {name: sha256(output / name) for name in ['imvigor210.sqlite', 'counts.npy']}
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf8')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path, help='New query_data directory; never overwritten')
    parser.add_argument('--archive', type=Path, help='Existing official tar.gz; skip network download')
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error('Output already exists; choose a new directory')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.imvigor210-', dir=output.parent) as temporary:
        temporary = Path(temporary)
        archive = args.archive.resolve() if args.archive else temporary / 'official.tar.gz'
        if not args.archive:
            print('Downloading official package:', URL, flush=True)
            with urllib.request.urlopen(URL, timeout=120) as response, archive.open('wb') as target:
                shutil.copyfileobj(response, target)
        if sha256(archive) != PACKAGE_SHA256:
            raise ValueError('Official package SHA256 mismatch; refusing to prepare unknown data')
        source = temporary / 'cds.RData'
        # Extract only the expected member, never arbitrary archive paths.
        with tarfile.open(archive) as package:
            member = package.getmember('IMvigor210CoreBiologies/data/cds.RData')
            if not member.isfile():
                raise ValueError('cds.RData is not a regular file')
            with package.extractfile(member) as stream, source.open('wb') as target:
                shutil.copyfileobj(stream, target)
        if sha256(source) != SOURCE_SHA256:
            raise ValueError('cds.RData SHA256 mismatch')
        staging = temporary / 'query_data'
        staging.mkdir()
        manifest = prepare(source, staging)
        staging.rename(output)
    print(json.dumps({'output': str(output), 'samples': manifest['samples'], 'patients': manifest['patients'],
                      'genes': manifest['genes']}))
    print('Configure TOOLUNIVERSE_IMVIGOR210_DATA=' + str(output))


if __name__ == '__main__':
    main()
