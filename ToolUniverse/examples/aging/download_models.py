"""Download the reviewed upstream artifacts on the networked development host.

Checksums lock the exact reviewed snapshot even though source URLs use branch names.
Usage: python download_models.py --output /shared/cluster-runtime/aging/upstream
Run prepare_model_bundles.py on a compute worker afterwards.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    records=json.loads((Path(__file__).parent/'model_sources.json').read_text())
    for entry in records:
        path=args.output/entry['filename']
        if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256']:
            continue
        temporary=path.with_suffix(path.suffix+'.download')
        with urllib.request.urlopen(entry['url'],timeout=900) as response,temporary.open('wb') as stream:
            while True:
                chunk=response.read(1024*1024)
                if not chunk:break
                stream.write(chunk)
        if hashlib.sha256(temporary.read_bytes()).hexdigest()!=entry['sha256']:
            temporary.unlink();raise RuntimeError(f'Upstream changed: checksum mismatch for {path.name}')
        temporary.replace(path)
    for filename in ['RNAAgeCalc.tar.gz','SenCID.tar.gz']:
        with tarfile.open(args.output/filename) as archive:
            archive.extractall(args.output,filter='data')
    print('All model assets verified. Prepare bundles on the compute worker.')

if __name__=='__main__':main()
