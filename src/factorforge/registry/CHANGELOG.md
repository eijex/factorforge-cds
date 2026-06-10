# Registry Changelog

## [patch] 2026-06-07 - 091 canonical label normalization

- job: 091 / PR: #49 / merge commit: `88aa900`
- change type: patch / canonical label normalization
- behavior change: no
- before: `["BsaI", "BsmBI", "BbsI"]` -> after: `["BsaI", "BpiI", "BsmBI"]`
- rationale: align registry label to established FactorForge production convention
- claim boundary: no biological scanning target change, no recognition-site behavior change, no benchmark scoring change
