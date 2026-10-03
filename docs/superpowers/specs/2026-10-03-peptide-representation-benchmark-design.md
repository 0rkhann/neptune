# Peptide representation benchmark — design spec

**Status:** design approved in conversation; not yet an implementation plan.
**Date:** 2026-10-03
**Evidence base:** `docs/superpowers/research/` — 12 cited research files, ~7,100 lines.

---

## 1. What this paper is

**The field measures validity and treats it as search quality. They are
different properties, and the evidence says they come apart.**

Two facts are only compatible if that is true:

- Tokenizer choice moves **validity** by ~28 points in chemistry — Atom-in-SMILES
  81.26 against NPBPE30k 53.77 [^1].
- Representation choice moves **optimization outcomes** barely at all —
  NovoMolGen across >30,000 experiments reports "only modest performance
  differences overall" [^2]; PMO states "SELFIES does not seem to offer an
  immediate benefit in optimization performance compared to SMILES except in
  GA" [^3]; and two pilots in this work (§4) found no confirmed effect.

Evolutionary computation has had the distinction for twenty years — **locality**
is neighbourhood preservation under a single move, not validity [^4] — and the
molecular generative literature cites none of it. Four independent searches
returned zero citations of Rothlauf, locality, redundancy or causality in
cheminformatics generative work [research: `2026-10-02-representation-locality…`].

Rothlauf & Oetzel already showed the failure mode: grammatical evolution, a
grammar-based genotype-to-phenotype map, "has problems with locality as many
neighboring genotypes do not correspond to neighboring phenotypes", and
"leads to lower performance for mutation-based search approaches" [^5].
**SELFIES is structurally the same kind of object and has never been analysed
this way.**

### Secondary claim, independently supported

**The apparatus could not detect a search difference even if one existed.**
See §3.

---

## 2. Scope and domain

**Primary domain:** peptidomimetics and cyclic peptides.
**AMP** supplies documented foreclosure evidence and one tier-2 oracle.
**Catalysis** appears as a one-paragraph motivating example only — the authors'
own objective requires exactly one phosphazene or TMG superbase, which FASTA
cannot express, so a FASTA arm scores zero on their real task by construction.
No catalysis experiments: ORD holds 1,430 organocatalysis reactions of 2.4M, and
no public superbase dataset exists [research: `catalysis-oracles-and-cost.md`].

---

## 3. Established findings (no further experiments needed)

### 3.1 The oracle layer is representation-locked

- **102 of 322** HELMCore peptide monomers have no canonical analog; the
  remaining 220 collapse onto 20 letters at a mean of **11.0 per letter**,
  worst case F at 30. Computed from `HELMCoreLibrary.json` (MIT).
- **RDKit's `MolFromHELM` accepts 46 of 649 monomers** and returns `None` for
  the rest, silently. All 649 probed [research: `assembly-stack-and-licences.md`].
- **pyPept** replaces a non-canonical residue with its closest natural analogue
  at Tanimoto ≥ 0.5, "otherwise, the non-natural amino acid is replaced by an
  alanine" [^6]. **No generative paper states whether its metrics inherit this.**
- **99.8%** of CycPeptMPDB peptides carry ≥1 nonstandard monomer [^7].
- The **AMP Challenge** (submissions Sept 2026, wet lab 2027) bans
  "noncanonical amino acids, stapled peptides, peptidomimetics, and chemically
  modified variants", FASTA only.

### 3.2 The field cannot compare itself

Survey of 8 peptide generative-design papers
[research: `peptide-evaluation-conventions.md`]:

- **0 of 8** state an oracle-call or compute budget.
- **5 mutually incompatible definitions of "valid"** in simultaneous use,
  including one (PepThink-R1) that folds oracle applicability into validity.
- **3 of 7** generative papers compare against baselines in a different
  representation without controlling for it. Only HELM-GPT names the confound.
- **2 of 17** methods accept a user-supplied objective without editing core code
  [research: `baseline-availability-audit.md`].
- **No peptide generative paper has run a controlled representation ablation.**
- Novelty is string identity in 6 of 7 papers; PepThink-R1 does set arithmetic
  on **uncanonicalised** SMILES.
- MOSES metrics are **undefined** on peptides, not merely untuned: `Scaf`
  returns NaN (`min_rings=2`), `Filters` returns 0 (≥8-atom ring rule), `Frag`
  degenerates under BRICS.

### 3.3 The metric determines the winner

PMO's full per-call results (figshare DOI 10.6084/m9.figshare.20123453, 2,875
runs) were downloaded and the top-10 curves recomputed from scratch, **validated
to max absolute error 0.0016 against all 23 published REINVENT values in Table 2**
[research: `pmo-saturation-and-metrics.md`]:

- **AUC and calls-to-threshold disagree on the winner on 3 of 4 saturated
  oracles.** drd2: SynNet best AUC, REINVENT fastest. albuterol: GP BO best AUC,
  REINVENT fastest. isomers_c7h8n2o2: SMILES-GA best AUC, STONED fastest.
- `qed`: 20 of 25 methods within **0.0093** of each other while needing
  **300 to 9,400 calls** — a 31× speed spread, invisible at the endpoint.
- For REINVENT on `qed`, **94.1% of the AUC integral accrues after the curve
  flattened**.
- **Correction to a common assumption:** PMO is *not* generally saturated.
  Only 2 of 23 oracles reach 99% of final by 3,000 calls; 16 of 23 are still
  improving past 7,000. The problem is per-oracle, not systemic.
- **GuacaMol has no budget argument at all** in `GoalDirectedGenerator`, and
  measured ~**200× spread** in actual oracle calls between rows of the same
  table (SMILES LSTM HC ~20k against Graph MCTS ~4M per task).

COCO states the problem verbatim: "if the fixed budget is too large, the
algorithm might solve the function before the budget is exceeded… The remedy is
to define a final target value and measure the runtime" [^8].

**The metric this situation needs already exists**: **Oracle Burden** — calls to
produce N molecules above threshold T, censoring as "Failed" [^9]. PMO logs no
calls-to-threshold metric at all.

### 3.4 The field's one representation metric rests on n=1

SELFIES' validity-under-mutation, resolved across three document versions
against the PDFs [research: `selfies-validity-numbers.md`]:

| Document | Scope | SMILES @1 mut | SELFIES |
|---|---|---|---|
| arXiv v1 (2019) | **QM9, population** | 18.1% | **88.1%** |
| arXiv v2 (2020) | MDMA, n=1 | 9.9% | 100% |
| **Published** [^10] | MDMA, n=1 | **26.6%** | 100% |

**Cite the published set: 26.6 / 9.0 / 3.7.** Three caveats travel with it:
**one molecule (MDMA)**; the operator is a **bitflip restricted to tokens
already present**; **no sample size is stated anywhere**. The only
population-scale measurement was in the 2019 preprint and was **dropped before
publication**.

---

## 4. Experiments run, and what they returned

Full data in `pilot_length_matched/`. Pre-registrations in that directory were
written before the corresponding runs existed.

### 4.1 The original signal was a configuration artifact

The pre-existing `l_UNCONSTRAINED_results` logP comparison showed monomer-level
tokenization beating character-level decisively (best logP 33.77 against 23.65).
**The two arms differed in four ways, not one:**

| | monomer | submonomer |
|---|---|---|
| `max_sequence_length` | 128 tokens ≈ 128 residues | 128 tokens ≈ 8 residues |
| vocabulary | 1190 | 34 |
| `diversity_filter` | `MonomerFingerprint` | `SubclassFingerprint` |
| replay `scaffold_type` | `monomer` | `subclass` |

Also: the committed config lists a `peptide_length` component absent from the
output, while the output carries `synthesizability_factor` absent from the
config. **Those results were not produced by the configs in the repository.**

### 4.2 Pilot A — matched arms, n=3 (`PREREGISTRATION.md`)

Configs byte-identical apart from the prior. Pre-registered primary: top-10 AUC
over the full budget.

**Result: D = 0.0171, R = 0.0229 → not detectable.** The direction also
**reversed** — submonomer slightly ahead.

A **post-hoc** trajectory analysis then found complete separation on
calls-to-threshold at 0.9, 0.95 and 0.99.

### 4.3 Pilot B — confirmatory, n=5 fresh seeds (`PREREGISTRATION_CONFIRM.md`)

Mann-Whitney U one-tailed on seeds 45–49 only; seeds 42–44 excluded as
hypothesis-generating.

**Result: U = 5.0, p = 0.0754 at all three thresholds → NOT CONFIRMED at
α = 0.05.**

Mechanism of the Pilot A false positive is visible: at a 1,000-call budget
**D barely moved (+0.0528 → +0.0566) while R quadrupled (0.0311 → 0.1427)**.
Three seeds could not estimate the spread.

### 4.4 What the data does show, descriptively

```
calls to top-10 mean >= 0.99
Pilot A  submonomer 1454, 1522, 1560          monomer 1962, 1963, 2327
Pilot B  submonomer 1078, 1446, 1451, 1694    monomer 1763, 1770, 2095, 2276, 2278
         submonomer 2838  <- the exception
```

**7 of 8 submonomer seeds are faster than every one of the 8 monomer seeds**,
consistently across three thresholds and two independent pilots.

And the arms differ in **reliability**, not only speed:

```
spread of calls to 0.95:  monomer 1329-1525 (196)   submonomer 763-1965 (1202)
```

**Both observations are post-hoc and neither has passed a clean test.**
Resolving them requires n=10 fresh seeds, pre-registered on speed *and* spread.

### 4.5 Confounds not addressed by any experiment so far

- **Embedding capacity**: vocabulary 1190 against 34, tied embeddings.
- **One oracle**, logP, which is degenerate — unbounded, decomposable, measured
  at r = 0.870 against heavy-atom count in these runs; PMO excludes it [^3].
- **Representation was never varied.** Only tokenization, within HELM.
- Wall clock: ~11 min (monomer) against ~56 min (submonomer) per run,
  identical oracle budget, sequential runs.

---

## 5. Benchmark design

### 5.1 Three factors, partially crossed

| Representation | Tokenizations |
|---|---|
| FASTA | residue |
| HELM | monomer, sub-monomer |
| CHUCKLES | residue-template, character |
| SMILES | atom-regex, character |

**Seven cells.** Three within-representation tokenization contrasts; one
four-way character-level row isolating representation at matched granularity.
Oracle nature is the third factor, spanning measured hackability (logP r = 0.870
against heavy atoms; permeability r = +0.09) and Ehrlich's ruggedness knob `q`.

**Graph is excluded as an arm.** It requires a different architecture —
Graph-Mamba must impose a degree ordering "because SSMs are recurrent models and
require ordered input" — which reintroduces the confound. Graph GA appears in
tier-2 as a labelled *uncontrolled* baseline.

### 5.2 Protocol

**Identity.** Assemble to an RDKit molecule via pinned `helmkit` (MIT, one
dependency, verified on macrocycles, disulfides, branches, staples, custom
libraries; raises `ValueError` on unknown monomers, never `None`), then key on
canonical isomeric SMILES with stereo preserved, one documented protonation
rule, **no tautomer normalisation** — RDKit's `Canonicalize()` erases every
α-carbon stereocentre by documented design (`removeSp3Stereo=true`).
**RegistrationHash is excluded**: its hashes changed in 2023.03, 2024.03,
2024.09, 2025.09 and 2026.03.6.

**Budget.** Two fixed budgets — unique scorable molecules (PMO-compatible) and
total generations. Run ends when either is exhausted; the consumed ratio is
reported. **Length matched in molecules, not tokens.**

**Validity.** One test, arm-independent: parses in its notation and assembles to
an RDKit-sanitisable molecule.

**Metrics.** **Oracle Burden** [^9] as the headline, not endpoint AUC.
Wilson or Jeffreys intervals for success rate (Brown, Cai & DasGupta: Wald
"should not be used"; Wilson or Jeffreys for n ≤ 40) [^11]. IQM with stratified
bootstrap CIs for continuous metrics, pooled over tasks to clear Agarwal's N=10
threshold [^12]. **Multiple Comparison Matrix** for ordering — pool-invariant,
unlike critical-difference diagrams, which are "open to both inadvertent and
intentional manipulation" [^13]. **IQM is excluded for success rate**: on
pooled 0/1 data it saturates at q ≤ 0.25 and q ≥ 0.75.
**Seed spread reported alongside every central estimate.**

**Objectives are bounded with a peak, not a plateau.** A saturating sigmoid is
not a bound: in these runs 208 molecules exceeded logP 20 at reward exactly
1.0000, so growing from logP 20 to 43.74 earned +0.0000 and cost nothing.
Use target windows. **No length cap of any kind** — structural caps are the
one mechanism argued against on principle, and the damage provably scales with
tokenization granularity [^14][^15].

**Disabled for all arms:** SMILES-only augmentation, the GraphGA mutator that
emits SMILES into non-SMILES runs, and beam enumeration (hardcodes `^`).

**Per-arm hyperparameter tuning is mandatory**, disclosed PMO-style. An untuned
arm is an invalid comparison, not a weak one [^16].

### 5.3 Oracles

**Tier-1 — constructed, known optimum.** Chemical Ehrlich functions [^17]:
structural SMARTS predicates at controlled topological distances on the
assembled molecule, paired on whether the motif requires a non-canonical
monomer. **No decomposition**, so an invented building block is judged on
whether it bears the feature. Optimum planted by construction and exhibited as a
witness molecule. Generalises `oracles/peptidomimetic/topological_potential_field.py`.

**Holdout by oracle score, not similarity.** Every corpus molecule scoring above
a stated threshold is removed from pretraining. Exact, metric-free. Similarity
is used only for the post-hoc leakage audit, with **count-based ECFP plus
topological torsion** — state of the art across 126 peptide datasets — against a
random-pair baseline measured in-corpus. **MAP4C is excluded**: independent
evidence is thin for MAP4 and absent for MAP4C; `pip install map4` 404s;
`mapchiral`'s LICENSE is 0 bytes; and v0.0.7 silently ignores its documented
`seed` argument.

**Tier-2 — realistic, externally validated only.**
Cyclic-peptide permeability as **logPapp regression under a scaffold split**,
inside bounded windows. Measured hack-resistant by adversarial probe against
PepINVENT's shipped XGBoost oracle: grafting alkyl to logP 34.7 *lowered* the
score; a greedy hill-climb converged to logP −3.55. The classifier is **not**
used — 51.7% of CycPeptMPDB already scores >0.9 under it.
AMP selectivity using BATTLE-AMP's **37,565 pre-generated composition-preserving
shuffled negatives** (MIT) [^18], paired with hemolysis via **Macrel 1.6.1**
(MIT, 1,000 peptides in 0.68 s, both probabilities on one row).
**QM and docking are excluded** — docking fails validity (Vina-family pose
accuracy collapses above ~10 rotatable bonds; peptides have 15–25; measured ~25%
top-pose success on 2–5-residue peptides [^19]), determinism (AutoDock-GPU
maintainers state fixed seeds do not guarantee identical results), and objective
validity.

### 5.4 Capacity and architecture

Three total-parameter-matched sizes over ~10×. Kaplan's non-embedding convention
is **not available** — it breaks when vocabulary is the variable under study, and
no scaling law covers V = 50 [^20][^21]. Tied-versus-untied as a crossed factor
at the reference size; no matched-parameter tying ablation exists in modern work.
Mamba for the full matrix, **transformer repeat at the reference size** on a
`q`-spanning subset — chemistry puts architecture span at ≈1.7 points against a
tokenizer span of ≈28 [^1], but genomics documents an SSM-versus-attention
granularity reversal.

### 5.5 Artifact

Standalone pip-installable repo, CPU-only, no `mamba-ssm`. NEPTUNE is one
entrant via an adapter. `helmkit` pinned with a conformance test suite.
Clean-room CHUCKLES implementation — PepFoundry is AGPL-3.0. **STONED excluded:
no licence file at all.** Ship HELMCoreLibrary, DBAASP, DRAMP general/clinical
(excluding the patent subset); fetch CycPeptMPDB by script — no licence
statement exists on its site.

---

## 6. Open hypotheses, in priority order

**H1.** Single-edit **validity** is statistically independent of single-edit
**molecular displacement**. If true, the field has been measuring the wrong
property. CPU, days.

**H2.** Locality differs *less* across representations than validity does.
Same computation as H1. Directly explains the nulls in §4.

**H3.** Is SMILES redundancy synonymous or non-synonymous under token edits?
Rothlauf's theory says that distinction decides whether redundancy is harmless
or destroys guided search [^4].

**H4.** n=10 fresh seeds per arm, pre-registered on **both** speed and spread,
to settle §4.4. ~10 hours on existing hardware.

**H5.** Does a speed difference appear on a **non-degenerate** oracle
(permeability, hack-resistant) when it did not on logP (hackable)? This is the
representation × oracle interaction.

**H6.** Capacity — parameter-matched retrained priors. Only if H1–H3 find a
mechanism.

**Note on H1/H2 for this system:** because NEPTUNE assembles from a library of
valid building blocks, **validity is pinned at 1.0 by construction** — observed
in all 12 runs. That makes it the controlled setting the field lacks: validity
constant, only locality free to vary.

---

## 7. Known limitations, stated rather than mitigated

- Capacity confounded throughout §4.
- One oracle, degenerate, in all experiments to date.
- Representation never varied; only tokenization, within HELM.
- §4.4's two observations are post-hoc.
- Hemolysis has **low predictability under 60%-identity homology splits**
  (QMAP, Sci Rep 2026) — acceptable for a fixed benchmark target, not for a
  biological claim.
- xTB's own HOMO correlates r = 0.769 with ωB97X-D/def2-SVP across 665k
  drug-sized molecules (QMugs, 2022).
- The research corpus was produced by AI agents; `docs/superpowers/research/README.md`
  lists eight claims that were wrong on first pass and corrected.

---

## 8. Timeline

NeurIPS **Evaluations & Datasets** track (renamed from Datasets & Benchmarks).
2026 cycle closed May 2026. Target 2027, expect early May. **Code and data must
be final at submission** — the track treats them as part of the paper, not
supplementary material.

---

## Citations

[^1]: Wang, Sultan, Volkamer & Klakow, arXiv:2602.13958 (2026). 3 architectures × 8 tokenizers × 2 splits.
[^2]: NovoMolGen, arXiv:2508.13408 (2025).
[^3]: Gao, Fu, Sun & Coley, "Sample Efficiency Matters", NeurIPS 2022 Datasets and Benchmarks, arXiv:2206.12411.
[^4]: Rothlauf, *Representations for Genetic and Evolutionary Algorithms*, Springer, https://link.springer.com/book/10.1007/3-540-32444-5
[^5]: Rothlauf & Oetzel, "On the Locality of Grammatical Evolution", https://link.springer.com/chapter/10.1007/11729976_29
[^6]: pyPept, J Cheminform 15:79 (2023).
[^7]: HELM-BERT, JCIM (2026).
[^8]: Hansen et al., COCO, arXiv:1605.03560 (2016).
[^9]: Guo & Schwaller, "Beam Enumeration", ICLR 2024, arXiv:2309.13957 — Oracle Burden.
[^10]: Krenn, Häse, Nigam, Friederich & Aspuru-Guzik, *Mach. Learn.: Sci. Technol.* 1 045024 (2020), DOI 10.1088/2632-2153/aba947.
[^11]: Brown, Cai & DasGupta, *Statistical Science* 16:101–133 (2001).
[^12]: Agarwal et al., "Deep RL at the Edge of the Statistical Precipice", NeurIPS 2021, arXiv:2108.13264.
[^13]: Ismail-Fawaz et al. (2023), Multiple Comparison Matrix.
[^14]: Grammar-Aligned Decoding, NeurIPS 2024.
[^15]: DOMINO, ICML 2024 — constrained-decoding damage scales with tokenization granularity.
[^16]: Tripp & Hernández-Lobato, arXiv:2310.09267 (2023). The follow-up is the **AI-for-Science workshop**, not ICML main track.
[^17]: Stanton et al., Ehrlich functions, arXiv:2407.00236 (2024); code `prescient-design/holo-bench`, PyPI `pytorch-holo`.
[^18]: BATTLE-AMP (2026), `szczurek-lab/battleamp-snakemake`, MIT.
[^19]: Rentzsch & Renard (2015), flexible-ligand Vina on 2–5-residue peptides.
[^20]: Kaplan et al., arXiv:2001.08361 (2020).
[^21]: Tao et al., "Scaling Laws with Vocabulary", arXiv:2407.13623 (2024).

**Further sources** — GuacaMol (Brown et al., JCIM 59:1096–1108, 2019);
MOSES (Polykovskiy et al., Front Pharmacol 11:565644, 2020);
Renz, Van Rompaey, Wegner, Hochreiter & Klambauer, Drug Discov Today Technol
(2019/2020), doi 10.1016/j.ddtec.2020.09.003; Langevin, Vuilleumier & Bianciotto,
J Cheminform 14:20 (2022); Ciepliński et al., JCIM 63:3238–3247 (2023),
arXiv:2006.16955 — **cite the range −0.01 to −0.79, not a single value**;
Guo & Schwaller, Augmented Memory, JACS Au 4:2160–2172 (2024);
Blaschke et al., J Cheminform 12:68 (2020); Demšar, JMLR 7:1–30 (2006);
Jelassi et al., ICML 2024, arXiv:2402.01032; MambaByte, arXiv:2401.13660 (2024);
Skinnider, Nat Mach Intell (2024); MAP4, Capecchi, Probst & Reymond,
J Cheminform 12:43 (2020); MAP4C, J Cheminform 16:53 (2024);
Venkatraman et al., Pharmaceuticals 17:992 (2024); Boldini et al.,
J Cheminform 16:35 (2024); Sidorczuk et al., Brief Bioinform (2022);
CycPeptMPDB, Li et al., JCIM 63:2240–2250 (2023); helmkit, J Cheminform (2026),
PMC13449311; PepINVENT, arXiv:2409.14040; PepEVOLVE, arXiv:2511.16912 (2025);
PepFoundry, JCIM 66:1264–1273 (2026); CHUCKLES, Siani, Weininger & Blaney,
JCICS 34:588–593 (1994); BetterBench, NeurIPS 2024 D&B.

**Verify before citing.** `docs/superpowers/research/README.md` records eight
claims corrected mid-investigation, including venue attributions and one numeric
correlation that does not exist in the form it was first reported.
