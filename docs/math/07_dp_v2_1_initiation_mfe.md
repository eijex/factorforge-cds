# Math Specification: DP v2.1 (Initiation-Aware Piecewise Objective Exact DP)

**Module:** `factorforge.engines.dp_v2_1`  
**Specification Version:** 2.1.0-RC1  
**Status:** Canonical Engineering Baseline  
**Theoretical Classification:** Exact Discrete Optimization over Finite-State & Bounded-Sum Constrained Search Spaces with Position-Dependent Piecewise Log-Fitness Objectives

---

## 1. Problem Definition & Scope

Given:
* A target amino acid sequence $P = (a_1, a_2, \dots, a_N)$ of length $N$,
* An upstream transcript context $U$ (e.g., $5'$ UTR, Kozak consensus, or construct flank),
* A specified stop codon $c_{\text{stop}} \in \{\text{TAA}, \text{TAG}, \text{TGA}\}$ (or stop design space),
* A downstream transcript context $R$ (e.g., construct flank or terminator),
* Host-specific codon fitness weight vector $w(c) \in (0, 1]$,
* Target global GC content bounds $[\text{GC}_{\min}, \text{GC}_{\max}] \in [0, 1]$,
* A forbidden recognition motif set $\mathcal{F} = \bigcup_{m \in M} \{m, \text{RC}(m)\}$ (e.g. Type IIS BsaI, BsmBI, BpiI),
* An initiation window length $K \in [10, 20]$ (default $K = 15$ codons, $45\text{ nt}$),
* Reference codon-preference ramp profile $r_i \in [0, 1]$ for $1 \le i \le K$,

Design an optimal Coding DNA Sequence (CDS) $X = (c_1, c_2, \dots, c_N)$ that:

1. **Translational Identity (Hard Invariant):** Translates exactly to $P$ ($\forall i, \text{translate}(c_i) = a_i$) with valid open reading frame and zero internal stop codons.
2. **Full-Construct Motif Exclusion (Hard Invariant):** Zero occurrences of any motif in $\mathcal{F}$ across intra-codon and inter-codon junctions, construct flanks, and stop junctions ($\text{count}(\mathcal{F}, U + X + c_{\text{stop}} + R) \equiv 0$).
3. **Hard Global GC Band Compliance (Hard Invariant / UNSAT):** Cumulative GC count $H(X) = \sum_{i=1}^N g(c_i)$ strictly satisfies $[\lceil 3N \cdot \text{GC}_{\min} \rceil, \lfloor 3N \cdot \text{GC}_{\max} \rfloor]$. If no reachable state satisfies the band, the engine raises `UnsatisfiableDesignError` (unless exploratory fallback is explicitly requested).
4. **Piecewise Objective Optimization (Exact Optimization):** Maximizes the composite discrete objective $J_{\text{v2.1}}(X)$ balancing $5'$ open-topology structural relaxation, codon-profile harmonization, and downstream elongation CAI maximization with full-path deterministic LexRank tie-breaking.
5. **Thermodynamic Post-Validation (Independent Proof):** True $5'$ transcript minimum free energy ($\Delta G_{5'\text{-window}}$) evaluated independently via ViennaRNA/RNAfold over the contiguous sequence $U_{\text{tail}(15\text{nt})} + c_1 + \dots + c_K$.

---

## 2. Mathematical Formalism & State Space

### 2.1. Definitions
* **Amino Acid Sequence:** $P = (a_1, \dots, a_N), \quad a_i \in \Sigma_{\text{AA}}$
* **Synonymous Codon Mapping:** $\Gamma(a_i) \subset \Sigma_{\text{NT}}^3$ where $\Sigma_{\text{NT}} = \{\text{A}, \text{C}, \text{G}, \text{T}\}$
* **Host Codon Adaptiveness:** $w(c) \in (0, 1]$, relative adaptiveness of codon $c$ in the target host.
* **Host-Relative Codon Preference Proxy:** $p(c) = w(c) \in (0, 1]$.
* **Codon GC Contribution:** $g(c) \in \{0, 1, 2, 3\}$, count of G and C nucleotides in codon $c$.
* **Codon AU Contribution:** $\text{AU}(c) = 3 - g(c) \in \{0, 1, 2, 3\}$.
* **Aho-Corasick Automaton:** DFA $\mathcal{A} = (Q, \Sigma_{\text{NT}}, \delta, q_0, T)$ tracking motif overlaps with $\mathcal{F}$.

### 2.2. Piecewise Position-Dependent Log-Fitness Function $\ell_i(c)$

For position $i \in \{1, \dots, N\}$ and candidate codon $c \in \Gamma(a_i)$:

$$\ell_i(c) = \begin{cases} 
\alpha \cdot H(i, c) + \beta \cdot \text{Pen}_{\text{open}}(c), & \text{if } 1 \le i \le K \quad \text{(Stage 1: Initiation Ramp)} \\
\ln w(c), & \text{if } K < i \le N \quad \text{(Stage 2: Elongation Body)}
\end{cases}$$

Where:
1. **Position-Dependent Codon-Profile Harmonization:**
   $$H(i, c) = - \left| p(c) - r_i \right|$$
   where $r_i = r_{\text{start}} + (r_{\text{end}} - r_{\text{start}}) \cdot \frac{i - 1}{\max(1, K - 1)}$ with default $r_{\text{start}} = 0.50$, $r_{\text{end}} = 0.90$.
2. **Open-Topology Linear Structural Proxy:**
   $$\text{Pen}_{\text{open}}(c) = - \gamma \cdot g(c) + \delta \cdot \text{AU}(c)$$
   with default $\gamma = 0.5$, $\delta = 0.5$.
3. **Hyperparameters:** $\alpha \ge 0, \beta \ge 0$ (defaults: $\alpha = 1.0, \beta = 1.0$).

### 2.3. Sufficient State Property

Let the DP state at position $i \in \{0, \dots, N\}$ be:
$$(i, h, q)$$
* $i \in \{0, \dots, N\}$: Processed codon count.
* $h \in \{0, \dots, 3N\}$: Cumulative GC nucleotide count.
* $q \in Q \setminus T$: Current safe state in the Aho-Corasick automaton.

> **Theorem (Sufficient State for Piecewise DP):**  
> Under position-dependent additive scoring $\ell_i(c)$ and additive GC accumulation, $(i, h, q)$ contains all necessary and sufficient information for future transitions. Therefore, the Bellman optimality principle holds unconditionally over $\prod_{i=1}^N \Gamma(a_i)$.

---

## 3. Recurrence Relation & Full-Path LexRank Tie-Breaking

### 3.1. Deterministic Multi-Tier Ranking

$$\text{maximize } \Phi_{\text{v2.1}}(X) = \left( J_{\text{v2.1}}(X), \; -\left| \sum_{i=1}^N g(c_i) - 3N \cdot \text{GC}_{\text{mid}} \right|, \; -\text{LexRank}(X) \right)$$

Where:
$$J_{\text{v2.1}}(X) = \sum_{i=1}^N \ell_i(c_i)$$
$$\text{GC}_{\text{mid}} = \frac{\text{GC}_{\min} + \text{GC}_{\max}}{2}$$
$$\text{LexRank}(X) = \text{Canonical lexicographical ordering over full synonymous codon paths } (c_1, \dots, c_N)$$

### 3.2. Base Case & Construct Flank Initialization
Given upstream context $U$:
1. Verify $\text{Safe}(q_0, U) == 1$. If $0$, raise `ValueError("Upstream context contains forbidden motif")`.
2. $q_{\text{start}} = \delta^*(q_0, U)$.
3. Initialize DP table:
   $$D_0(h, q) = \begin{cases} 0.0, & \text{if } h = 0 \text{ and } q = q_{\text{start}} \\ -\infty, & \text{otherwise} \end{cases}$$

### 3.3. Forward Recurrence with Canonical Path-Rank Tracking
For $i = 1 \dots N$, for each active $(h, q)$ at step $i-1$ with predecessor path rank $R_{i-1}(h, q)$, and for each $c \in \Gamma(a_i)$:
1. If $\text{Safe}(q, c) == 1$:
   $$q' = \delta^*(q, c), \quad h' = h + g(c)$$
   $$\text{cand\_score} = D_{i-1}(h, q) + \ell_i(c)$$
   $$\text{cand\_key} = (R_{i-1}(h, q), \text{Index}(c))$$
2. Relax into $D_i(h', q')$:
   - If $\text{cand\_score} > D_i(h', q') + \epsilon$: accept candidate.
   - If $|\text{cand\_score} - D_i(h', q')| \le \epsilon$: accept candidate if $\text{cand\_key} < \text{existing\_key}$ (True Full-Path LexRank).
3. Establish canonical layer path rank $R_i(h', q')$ by sorting active states by winning $(R_{i-1}, \text{Index}(c))$.

### 3.4. Termination & Full Construct Tail Verification
Let $H_{\text{valid}} = [\lceil 3N \cdot \text{GC}_{\min} \rceil, \lfloor 3N \cdot \text{GC}_{\max} \rfloor]$ and $S_{\text{tail}} = c_{\text{stop}} + R$:

$$J^* = \max_{\substack{h \in H_{\text{valid}} \\ q \in Q \setminus T \\ \text{Safe}(q, S_{\text{tail}}) = 1}} D_N(h, q)$$

If no state satisfies $h \in H_{\text{valid}}$ with $\text{Safe}(q, S_{\text{tail}}) == 1$:
* Raise `UnsatisfiableDesignError` (P0 Hard Invariant Enforcement).
* Otherwise, backtrack using deterministic $\text{Parent}$ pointers to emit the unique, global optimum sequence $X^*$.

---

## 4. Complexity & Theoretical Bounds

* **Time Complexity:** $O(N \cdot \Delta H \cdot |Q| \cdot |\Gamma|_{\max})$ where $|\Gamma|_{\max} \le 6$.
* **Backtracking Storage:** $O(N \cdot \Delta H \cdot |Q|)$ stored as sparse layer predecessors.
* **Empirical Latency:** $< 25\text{ ms}$ for $N=235\text{ AA}$ (Humira LC), $< 50\text{ ms}$ for $N=472\text{ AA}$ (Humira HC) on a single CPU core.

---

## 5. Thermodynamic Validation Protocol (Independent Verification)

1. **Folding Window:** $\text{Window} = U_{\text{tail}(15\text{nt})} + c_1 + \dots + c_K$ (Total length: $45 \sim 60\text{ nt}$).
2. **Thermodynamic Model:** ViennaRNA RNAfold (`RNA.fold` at $37^\circ\text{C}$, Turner energy parameters).
3. **Evaluation Outcomes:**
   - $\Delta G \ge -12.0\text{ kcal/mol} \implies \text{PASS}$ (relaxed, low secondary structure).
   - $\Delta G < -12.0\text{ kcal/mol} \implies \text{WARNING}$ (elevated initiation hairpin risk).
   - ViennaRNA unavailable $\implies \text{INDETERMINATE}$ with `enforcement = WARNING`.
4. **Provenance Package:** SHA-256 digest of $(U, c_1 \dots c_K, c_{\text{stop}}, R)$ recorded in `Metrics.context_digest`.
