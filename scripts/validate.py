#!/usr/bin/env python3
"""app.yml against the Pennsieve app manifest schema, and CITATION.cff with cffconvert.

Needs jsonschema >= 4 and pyyaml; cffconvert (which pins an old jsonschema) is looked for on PATH
or as CFFCONVERT=/path/to/cffconvert. The schema is read from SCHEMA (default: the sibling
workflow-single-cell repo).
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parent.parent
schema_path = Path(os.environ.get("SCHEMA", ROOT.parent / "workflow-single-cell/docs/app-manifest.v1.schema.json"))
schema = json.loads(schema_path.read_text())
errors = list(jsonschema.validators.validator_for(schema)(schema).iter_errors(yaml.safe_load((ROOT / "app.yml").read_text())))
for e in errors:
    print(f"app.yml: {e.message} at {list(e.path)}")
print(f"app.yml: {'valid' if not errors else f'{len(errors)} error(s)'}")
cff = os.environ.get("CFFCONVERT") or shutil.which("cffconvert")
rc = 0
if cff:
    rc = subprocess.run([cff, "--validate", "-i", str(ROOT / "CITATION.cff")]).returncode
else:
    print("cffconvert not found; CITATION.cff not validated")
sys.exit(1 if errors or rc else 0)
