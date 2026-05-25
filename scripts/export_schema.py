"""Export the DesignPackage JSON Schema to a file."""

import json
from pathlib import Path

from factorforge.schemas import DesignPackage


output = Path("src/factorforge/schemas/design_package.schema.json")
output.write_text(json.dumps(DesignPackage.model_json_schema(), indent=2))
print(f"Schema exported to {output}")
