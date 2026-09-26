# FactorForge SOP Quick Guide

## 1. What is a FactorForge SOP?
A FactorForge Standard Operating Procedure (SOP) file is a YAML document that tells the system two things:
1. **Design Settings**: How to generate the DNA sequence (e.g., optimization method, restriction enzymes to avoid).
2. **Laboratory Review Policy**: How to classify findings in the generated sequence (e.g., whether a cryptic splice site is a critical failure or an acceptable warning).

## 2. Finding vs. Policy
- **Finding**: A computational observation (e.g., "BsaI site detected at nt 100"). The FactorForge backend scanner generates these autonomously.
- **Policy**: Your laboratory's rule for that finding (e.g., "BsaI sites are `HARD_FAIL`"). The SOP defines the policy.

## 3. Policy Enforcement Levels
Each rule in the SOP can be assigned one of three levels:
- **`HARD_FAIL`**: The sequence is rejected immediately. Use this for strict technical constraints (e.g., synthesis limits, essential cloning restriction sites).
- **`WARNING`**: The sequence is flagged for human review. It is not blocked, but the scientist must explicitly approve the finding before proceeding.
- **`IGNORE`**: The finding is suppressed and not reported in the review dashboard.

## 4. Design Settings vs. Review Policy
- **Design Settings** (`workflow:`) control the active generation of the sequence. For example, setting `design_method: feasibility_best` tells the optimizer to balance biological constraints mathematically.
- **Review Policy** (`rules:`) controls how the final sequence is judged.

## 5. Modifying the SOP
To update your laboratory policy:
1. Go to the FactorForge **Sequence Policy** panel.
2. Click **Download Example SOP** to get a baseline YAML file.
3. Open the file in any text editor.
4. Modify the values under the `LABORATORY SETTINGS - EDIT BELOW` section (e.g., changing `rna.cryptic_splice.v1: WARNING` to `IGNORE`).
5. Save the file and upload it back to FactorForge.

## 6. Upload, Apply, and Export
- **Upload**: Select your `.yaml` file to load it into the browser.
- **Apply SOP**: Validates the SOP and sets it as the active policy for your browser session. (Note: Invalid SOPs are rejected to ensure safety).
- **Export YAML**: Downloads your currently active, validated SOP.

## 7. Important Boundaries
**The Conservative Template is Not an Approved SOP**
The default templates provided by FactorForge (e.g., "Conservative Plant Expression") are in-silico starting points. They do not constitute an approved, validated biological SOP for any specific laboratory. You must review and approve these settings internally.

**Privacy and Data Security**
Do not include sensitive information in your SOP file. Keep raw DNA sequences, protein sequences, internal construct IDs, patient data, and collaborator confidential information out of the YAML file.
