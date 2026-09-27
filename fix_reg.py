import re

with open(r'src/factorforge/rules/registry.py', 'r', encoding='utf-8') as f:
    text = f.read()

text = text.replace('f"{r.rule_id}:{r.version}:{f\\"{self.resolved_enforcement(r).value}:{self.resolved_authorized_action(r).value}\\"}:{r.authority.authority_type.value}"', 'f"{r.rule_id}:{r.version}:{self.resolved_enforcement(r).value}:{self.resolved_authorized_action(r).value}:{r.authority.authority_type.value}"')
text = text.replace('f"{r.rule_id}:{r.version}:{f"{self.resolved_enforcement(r).value}:{self.resolved_authorized_action(r).value}"}:{r.authority.authority_type.value}"', 'f"{r.rule_id}:{r.version}:{self.resolved_enforcement(r).value}:{self.resolved_authorized_action(r).value}:{r.authority.authority_type.value}"')

with open(r'src/factorforge/rules/registry.py', 'w', encoding='utf-8') as f:
    f.write(text)
