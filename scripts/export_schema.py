"""
export_schema.py — Export DesignPackage JSON Schema

Writes the FactorForge DesignPackage JSON Schema to a file for
documentation or downstream validation use.

Usage:
    python scripts/export_schema.py
    python scripts/export_schema.py --output schema.json
"""

import json
from pathlib import Path

from factorforge.schemas import DesignPackage


output = Path("src/factorforge/schemas/design_package.schema.json")
output.write_text(json.dumps(DesignPackage.model_json_schema(), indent=2))
print(f"Schema exported to {output}")
