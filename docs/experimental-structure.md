# Position-conditioned DP research module v0.1.0

`factorforge.engines.structural_pause.StructuralPauseDPOptimizer` is an opt-in
research interface. It is available through an explicit public browser testing route; it is not a
profile default. Supply explicit one-based residue positions, a structure reference,
target relative codon adaptiveness and nonnegative penalty weight. The target
is uncalibrated; no universal 0.35–0.50 biological interval is established.

The additive objective is log-adaptiveness minus a position-local absolute
distance penalty. Output CAI is recomputed independently from selected codons;
`score` is the penalized objective. The existing position/GC/automaton state
coordinates are computational dimensions, not a protein-folding simulator.
Infeasible GC bands are rejected in this experimental module, while the
unchanged baseline DP API can report a nearest-band fallback as infeasible.

No confirmed PDB/AlphaFold domain mapping is inferred automatically. A reviewed explicit
mapping is required. Missing annotation does not support a folding claim.
Synthetic exhaustive comparisons verify encoded objective and constraints;
they do not show expression, soluble-yield or aggregation improvement.

Sent designs remain historical. Regeneration requires a new identity and
parameter/structure manifest, independent translation/GC/motif/flank checks,
an old/new comparison and human biological review before synthesis replacement.
The experimental YieldPredict 14D contract requires independent evaluation and
retraining; default 13D checkpoints must not be silently repurposed.

## Direct 3DP browser testing (web build; Jobs 389–397)

The FactorForge form supports explicit 1DP/3DP selection. 3DP uploads a PDB/mmCIF file, selects chain and one-based file model order, and requires reviewed one-based protein positions. The selected observed chain must exactly match the complete protein; partial or homologous structures require the retained local reviewed mapping workflow. Direct requests use `/api/optimize` with `design_mode=paired_3dp` and `structure_context={content,format,chain,model,positions,positions_reviewed,target_adaptiveness,penalty_weight}`.

Limits: 3–300 amino acids, structure content 250 KB, HTTP JSON body 1 MB. Both candidates use the configured codon reference, requested GC and canonical motif exclusions plus custom additions. Invalid inputs and infeasible candidates fail without a baseline fallback. Existing sequence-only settings and API requests remain supported. 3DP uses its explicit DP settings and does not apply Kozak, TpA, template or seed transformations. SOP checks are reported for both candidates; regeneration is manual and never silently changes mode. Watermark stays disabled.

This is stateless computational testing, with no ResearchOS approval or ValidationHub evidence recording. Those existing integrations remain available for later operational use. The structure file is not retained in the response or browser history; normal design results may be retained in local browser history. Submit only public-safe inputs. Structural conformity metrics are codon-weight proxies, not folding/expression measurements. Web deployment revisions are tracked separately from a package release.

## Observation integration (Job 397)

Both browser modes attach `scientific_review` to the primary output and each
paired candidate: sense-CDS motif scan (`SEQUENCE_SCAN`) and initiation-window
MFE (`MODEL_PREDICTION`). The scan reports matched patterns, not confirmed RNA
processing. A zero-hit scan is not proof of no biological risk.

MFE uses local ViennaRNA with per-call model details (37 Celsius by default),
a 48-nt CDS start window, and at most four explicitly supplied upstream bases.
No upstream sequence is invented. Coordinates are zero-based relative to CDS
start, end-exclusive. Missing dependencies, invalid results or non-finite values
are unavailable/error, never PASS. A required unavailable policy observation
returns HTTP 422; the default advisory does not change candidate generation.
The retained -3 kcal/mol advisory is not a calibrated plant efficacy cutoff.

POST `/api/optimize` with `action=recommend_structure`, sequence and a matching
structure context suggests positions using declared AlphaFold pLDDT or relative
experimental B-factor observations. It does not assign DSSP secondary structure,
validate physical pause/linker sites or choose a winner. Unknown sources, missing
CA observations and flat unsupported regions do not create midpoint suggestions.
Returned `reviewed=false` requires explicit user acceptance and separate review.
Changing sequence or structure invalidates review. Results, HTML reports and local
history retain observations; uploaded structure files are not retained.

Scope guard: a returned sequence must translate exactly to the input protein before it is treated as sense CDS for these observations. Expanded constructs without explicit transcript/CDS mapping report NOT_AVAILABLE rather than folding promoter/vector flanks as an initiation region.

## Local NB-Fold reviewed CLI pipeline (Jobs 398–399; Phase A/B/C)

`factorforge nbfold` provides a two-step local reviewed CLI workflow (`factorforge.nbfold`) for pre-computed AlphaFold 2/3 or experimental PDB/mmCIF structure files without running remote folding servers or protein language models:

1. `factorforge nbfold annotate`: validates a single-record 20-AA protein FASTA (`<= 16 KB`, `3–300 aa`) and local `.pdb`/`.cif`/`.mmcif` file (`<= 250 KB`) with explicit chain and 1-based file model order, checks 100% exact residue sequence match without gaps/disorder/modified residues, extracts `CA`-atom confidence (`CA_atom_confidence` for AF3 vs `residue_plddt` for AF2) plus optional verified confidence/PAE/InterPro/UniProt sidecars (`<= 2 MB` each, total input `<= 5 MB`), runs Phase C structural context observers, and atomically writes `annotation.json` with `reviewed=false`.
2. `factorforge nbfold optimize`: verifies an explicit user-prepared `review.json` (`review_type=LOCAL_USER_ACKNOWLEDGEMENT`, `positions_reviewed=true`, matching protein/structure/sidecar/annotation/configuration SHA-256 digests and 1-based unique positions), checks that any active structural annotation conflicts are explicitly acknowledged (`conflicts_acknowledged=true` plus non-empty `conflict_resolution_rationale`), runs paired `baseline_1dp` (`DPV2Optimizer`) and `structure_3dp` (`StructuralPauseDPOptimizer`) under identical resolved host reference, GC bounds, canonical Type IIS motifs, and flanks, evaluates independent hard QC, sense-CDS plant RNA pattern scans (`SEQUENCE_SCAN`), and 5' initiation MFE (`MODEL_PREDICTION`), and atomically writes `paired_sequences.fasta`, `observations.json`, `run_manifest.json`, and `NB_FOLD_REPORT.md` (or sanitized `error_manifest.json` on failure).

Phase C structural context annotations (`factorforge.nbfold.annotations`):
- `domain_boundary_assignment` (`domain_caller`): segments domains and linker/tail regions from verified PAE matrices ($\ge 20\text{ aa}$, 20-residue sliding window) and/or verified InterPro JSON sidecars (`--interpro-sidecar`), reporting unannotated terminal intervals and gaps between selected domain intervals; these are heuristic context labels, not established disorder or physical linkers and computing `intra_domain_mean_pae` vs `inter_domain_mean_pae` in `domain_relative_confidence`. Boundary divergences (`> 10 aa`) between PAE and InterPro are flagged as `PAE_INTERPRO_BOUNDARY_CONFLICT`.
- `transmembrane_topology` (`tm_signal_caller`): maps 1-based Signal peptide, Transmembrane helix, Topological domain, Disulfide bond, and Glycosylation features from a verified local UniProt JSON sidecar (`--uniprot-sidecar`).
- `backbone_secondary_structure`: observes local `H` (helix), `E` (strand), and `C` (coil/loop) geometry from backbone dihedral angles ($\phi, \psi$), $i \to i+4$/$i \to i+3$ $\text{C}=\text{O} \cdots \text{N}$ oxygen-nitrogen distances, when full `N, CA, C, O` backbone coordinates are present; returns `NOT_AVAILABLE` on `CA`-only structures without inventing `H` or `E`. External binary DSSP (`secondary_structure_dssp`) remains separately reported as `NOT_AVAILABLE`.
- When candidate or reviewed positions overlap a Signal peptide, Transmembrane segment, Disulfide/Glycosylation site, rigid domain interior core, or rigid `H`/`E` secondary-structure core (or when PAE/InterPro or SP/TM sidecar conflicts exist), `optimize` rejects unacknowledged reviews with `REVIEW_CONFLICTS_UNACKNOWLEDGED` (exit code `3`).

Boundaries: `LOCAL_USER_ACKNOWLEDGEMENT` is not a cryptographic signature or server approval. Phase C structural annotations provide auxiliary review context and do not alter `auto_annotator` candidate recommendation rules or the `1DP`/`3DP` optimization objectives (`DPV2Optimizer`, `StructuralPauseDPOptimizer`). `codon_weight_target_conformity` (`structure_fit_*`) measures proximity of selected codon relative adaptiveness weights to `target_adaptiveness` at reviewed positions, not physical folding accuracy or soluble expression yield. 14D `YieldPredict` remains `NOT_AVAILABLE` (`value=null`), and `ESMFold`/`ESM-2` capabilities remain disabled (`available=false`, failing closed with `UNSUPPORTED_BACKEND`).

Changing provider/confidence interpretation or export metadata invalidates the previous annotation review. Positions differing from recommendations, including a `NO_CANDIDATES` result, require matching `user_curated_positions` in a rehashed annotation or a non-empty `modification_rationale` in the explicit review manifest.

The current timeout is checked after computational stages; it rejects overdue runs but does not interrupt an optimizer while it is running. Staged files are replaced individually via per-file atomic replacement (`os.replace`, 파일별 원자적 교체); export is not a transaction across all four files. A filesystem failure can leave partial artifacts without a completed manifest. Treat a run as complete only with `run_manifest.json` and matching artifact hashes.

Phase C intake requires at least one valid sequence/hash/length binding field; malformed supplied fields are rejected. Length-only binding checks coordinate compatibility and does not establish sequence identity. Supply the exact sequence or its SHA-256 for identity verification. Unknown UniProt source versions remain null. Missing or degenerate dihedrals are `U` (unassigned), with `PARTIAL`/`NOT_AVAILABLE` status when appropriate. The local H/E/C observer is an unvalidated heuristic, separate from DSSP; synthetic geometry tests do not establish agreement on real structures. PAE segmentation and domain-core warnings remain experimental review heuristics.

## Phase C specification alignment & synthetic reconstruction checks (Job 400)

- **Baseline preservation (`v1.2.1`):** For `v1.2.1`, errors found in the Job 399 review were fixed and related tests passed (이번 리뷰에서 발견한 오류를 수정하고 관련 테스트를 통과했다); `v1.2.1` is preserved as the operational baseline rather than claiming all possible runtime bugs are closed or calling the thresholds a universal biological standard.
- **Explicit PAE parameter semantics (`PAE_PARAMETER_PROFILES`):**
  - `min_protein_aa` (`20 aa`): minimum total protein length to attempt PAE segmentation.
  - `min_subdomain_aa`: initial 20-residue sliding-window seed length (`8 aa` in `v1_2_1_baseline` vs `10 aa` in `prd_v1_0_0_candidate`).
  - `min_domain_aa`: post-merge global domain minimum length (`8 aa` in `v1_2_1_baseline` vs `40 aa` in `prd_v1_0_0_candidate`).
  - `merge_cutoff_angstrom`: cross-subdomain merge threshold (`12.0 Å` in `v1_2_1_baseline` vs `15.0 Å` in `prd_v1_0_0_candidate`).
- **InterPro `entry_type` separation:** `domain_caller` preserves all normalized InterPro entries in `interpro_annotations`, filters `interpro_domains` and `consensus_domains` strictly to `entry_type == "domain"` (`used_for_domain_consensus=True`, plus untyped legacy entries for backward compatibility), and records `family`, `homologous_superfamily`, `repeat`, and `site` entries in `interpro_non_domain_annotations` (`used_for_domain_consensus=False`) so full-length Family/Superfamily spans or short Repeat/Site motifs never collapse multi-domain boundaries or erase inter-domain linkers. Overlapping `Domain` entries record `merged_domain_ids` and `merged_accessions`.
- **Synthetic reconstruction checks:** current Job400 fixtures generate backbone coordinates from hand-entered H/E/C labels and compare them with those same labels. PAE and InterPro intervals are constructed fixtures. They verify code behavior, not independent DSSP agreement, real public-coordinate performance, or biological superiority of the default PAE profile. Genuine public-coordinate/provider/DSSP evaluation remains open.

Web structure handling requires an exact-match local upload and explicit human-reviewed positions. Remote inference, private reference cache lookup, artificial coordinate/confidence fallback, and automatic approval are disabled. Results without bound structure data hide the 3D viewer. Uploaded sources remain provider-unverified unless separately established; panels annotate synonymous-codon choices on the same coordinates.
