import re

with open(r'api/optimize.py', 'r', encoding='utf-8') as f:
    text = f.read()

# The local imports are at line ~553
# We want to remove 'import tempfile', 'import os', 'import json' inside do_POST
lines = text.split('\n')
out = []
for line in lines:
    if line.strip() in ['import tempfile', 'import os', 'import json']:
        # Only skip if it is heavily indented (inside do_POST)
        if line.startswith('                    import'):
            continue
    out.append(line)

with open(r'api/optimize.py', 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
