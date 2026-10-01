"""Public repository must not schedule an absent private learning runner."""

from pathlib import Path


def test_public_learning_workflow_is_manual_only():
    workflow = (Path(__file__).resolve().parents[1] / '.github' / 'workflows'
                / 'agentos_continuous_learning.yml').read_text(encoding='utf-8')
    assert 'workflow_dispatch:' in workflow
    assert 'schedule:' not in workflow
    assert 'AgentOS runner is not installed in this public repository' in workflow
