"""Assemble reviewed schema fragments into the category registry."""
import json
from pathlib import Path
root=Path(__file__).resolve().parents[2]
configs=[]
for part in ('data','models','statistics','mechanism'):
    configs.extend(json.loads((Path(__file__).parent/f'schema_{part}.json').read_text()))
for config in configs:
    config['mcp_schema_mode'] = 'passthrough'
    config['parameter']['properties']['request'].setdefault('type', 'object')
assert len(configs)==12
assert len({x['name'] for x in configs})==12
(root/'src/tooluniverse/data/aging_tools.json').write_text(json.dumps(configs,indent=2)+'\n')
