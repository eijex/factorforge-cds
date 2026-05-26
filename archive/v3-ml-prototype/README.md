# v3 ML Prototype Archive

This directory preserves the historical ML-oriented prototype work, including
early training, evaluation, decoder, Kaggle notebook, and model-card materials.

The archived prototype is not installed, registered, or used by the FactorForge
v3.1.x production workflow. Current public releases use deterministic,
constraint-aware design paths and the `profile` engine namespace.

Large generated experiment outputs, especially long CSV loss logs and decoded
comparison tables under the former `experiments/results/alpha_run*` paths, are
intentionally excluded from this archive to keep the public repository compact.
The retained code and reports are enough to document the research direction
without presenting the prototype as production functionality.

The historical Kaggle training notebook is intentionally excluded from the
public archive because it contains external workspace and checkpoint handling
that is not part of the supported release.
