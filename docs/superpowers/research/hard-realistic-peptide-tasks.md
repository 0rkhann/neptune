# Hard, realistic oracle tasks for a peptide generative-design benchmark

Research dossier. Every claim carries a URL and year. Numbers marked **[MEASURED HERE]** were
computed in this session against the real artefacts (CycPeptMPDB v1.2 CSV; PepINVENT's shipped
XGBoost permeability model), not taken from an abstract. Scripts and data are in the same
scratchpad directory (`an1.py`, `an2.py`, `oracle.py`, `probe.py`, `probe2.py`, `cyc_all.csv`,
`PepInvent/`).

Reproduction notes: `cyc_all.csv` = `http://cycpeptmpdb.com/static//download/peptides/CycPeptMPDB_Peptide_All.csv`
(8,466 rows, 7,991 unique SMILES, v1.2, retrieved 2026-10-02). The PepINVENT oracle pickle needs
`xgboost==1.7.5` + `scikit-learn==1.2.2` + `numpy<2` on Python 3.10, plus the attribute shim in
`oracle.py` (the pickle predates several ctor params of the modern `XGBClassifier`).

---

## SECTION 1 — Is cyclic-peptide permeability hackable?

### 1.0 Bottom line

**No, not by the size/lipophilicity inflation route.** Measured directly against PepINVENT's
shipped oracle, grafting an alkyl tail onto cyclosporin A drives Crippen logP from 3.3 to 34.7 —
squarely into the pilot run's logP 33–43 regime — and the predicted permeability score goes
**down**, 0.816 → 0.772, then flat. A greedy hill-climb against the same oracle converged to a
peptide at Crippen logP **−3.55**. The degenerate move is actively anti-correlated with reward.

The real failure mode of this oracle is a different one: **saturation**. 51.7% of CycPeptMPDB
already scores >0.9 and 15.9% scores >0.99, so the objective is nearly solved by the training
distribution before an optimiser starts. See §1.6.

### 1.1 Why it is not hackable: the labels themselves are not a function of size or grease

Correlations of **experimental** log Papp against RDKit descriptors, over the whole of
CycPeptMPDB v1.2. **[MEASURED HERE]** (`an1.py`). Censored values at the −10 detection floor
dropped (n = 8,193):

| descriptor | Pearson r | Spearman ρ | range in DB |
|---|---|---|---|
| Crippen MolLogP | **+0.092** | +0.126 | −4.0 … 7.3 |
| MolWt | −0.034 | −0.051 | 342 … 1778 |
| HeavyAtomCount | −0.029 | −0.051 | 24 … 127 |
| NumHDonors | −0.152 | −0.187 | 0 … 13 |
| NumHAcceptors | −0.088 | −0.111 | 4 … 19 |
| TPSA | −0.094 | −0.153 | 73 … 436 |
| NumRotatableBonds | −0.005 | +0.036 | 3 … 59 |

For contrast, within this same set `r(MolLogP, HeavyAtomCount) = +0.240`. So logP and size are
coupled to each other but **neither is coupled to the measured property**. This is the structural
reason the task resists the degenerate strategy: the label is not a monotone function of the thing
an optimiser can trivially inflate.

The pooled numbers understate the real signal, because the database is a union of congeneric
series. Within-source correlations, sources with n ≥ 100 **[MEASURED HERE]** (`an2.py`):

| source | n | r(logP) | r(HAC) | r(HBD) | r(TPSA) |
|---|---|---|---|---|---|
| 2020_Townsend | 2881 | +0.243 | −0.128 | −0.421 | −0.486 |
| 2021_Kelly | 1519 | +0.449 | +0.233 | −0.523 | −0.517 |
| 2013_CHUGAI | 878 | **−0.397** | −0.477 | −0.326 | −0.410 |
| 2016_Furukawa | 680 | **−0.199** | −0.278 | +0.096 | −0.044 |
| 2023_Ohta | 584 | +0.491 | +0.151 | −0.343 | −0.192 |
| 2018_CHUGAI | 374 | −0.062 | −0.249 | −0.285 | −0.311 |
| 2024_Faris | 207 | +0.227 | +0.167 | — | — |
| 2022_Bhardwaj | 136 | +0.332 | −0.171 | −0.033 | −0.266 |
| **n-weighted mean** | | **+0.173** | **−0.085** | **−0.360** | **−0.403** |

Read this carefully — it is the single most useful table here:

- **logP flips sign across series** (−0.40 to +0.49). Lipophilicity is not a direction, it is a
  window whose optimum depends on scaffold. An optimiser cannot ride it.
- **HBD and TPSA are the only consistent drivers**, and both are *negative*, i.e. the oracle
  rewards *removing* polar groups. That direction is **bounded below** — HBD cannot go below 0,
  TPSA cannot go below the backbone minimum. Unlike logP, it cannot be inflated without limit.
- **Heavy-atom count is near-zero to negative.** Growing the molecule does not help.

Ratio-form descriptors behave better than either numerator or denominator alone **[MEASURED HERE]**:
`r(TPSA/HeavyAtomCount, perm) = −0.217` vs `r(TPSA) = −0.094`; `r(MolLogP/HAC) = +0.156` vs
`r(MolLogP) = +0.092`. Relevant to §5.

### 1.2 How much of the oracle can a size+grease-only model reproduce?

5-fold random-split HistGradientBoosting on CycPeptMPDB, n = 8,193 **[MEASURED HERE]** (`an2.py`):

| feature set | Pearson r | R² | MAE |
|---|---|---|---|
| Crippen logP alone | 0.389 | 0.152 | 0.551 |
| logP + MW | 0.525 | 0.276 | 0.506 |
| logP + MW + HBD + TPSA | 0.577 | 0.333 | 0.483 |
| 7 simple 2D descriptors | 0.603 | 0.363 | 0.470 |
| **published CycPeptMP fusion model** | **0.883** | — | **0.355** |

A pure bulk-property model tops out at r = 0.60 / MAE 0.47. CycPeptMP reaches 0.883 / 0.355
(Li, Yanagisawa & Akiyama, *Brief Bioinform* 25(5):bbae417, 2024,
https://doi.org/10.1093/bib/bbae417). So roughly **two-thirds of the explainable variance lives
outside bulk size and lipophilicity** — in stereochemistry, N-methylation pattern and conformation.

### 1.3 Direct adversarial probe of PepINVENT's shipped oracle

PepINVENT (Geylan et al., *Chem Sci*, 2025, https://doi.org/10.1039/D4SC07642G;
repo https://github.com/MolecularAI/PepInvent) ships `data/models/predictive_model.pckl` +
`feature_scalar.pckl`. Reading `pepinvent/scoring_function/scoring_components/predictive_model.py`:
the featurisation is `AllChem.GetMorganFingerprint(radius=4, useChirality=True, useCounts=True)`,
folded modulo 2048, MinMax-scaled, then `XGBClassifier.predict_proba()[:,1]`. Loaded model:
**450 trees, max_depth 18, 2048 features**. Trained on CycPeptMPDB PAMPA; paper reports balanced
accuracy 0.78, MCC 0.59.

**Oracle output vs descriptors across all 7,991 unique CycPeptMPDB peptides [MEASURED HERE]**
(`probe.py`):

| | Pearson r | Spearman ρ |
|---|---|---|
| P(permeable) vs Crippen MolLogP | **+0.165** | +0.215 |
| P(permeable) vs MolWt | +0.203 | +0.233 |
| P(permeable) vs HeavyAtomCount | +0.204 | +0.234 |
| P(permeable) vs NumHDonors | −0.097 | −0.055 |
| P(permeable) vs TPSA | +0.104 | +0.077 |
| P(permeable) vs NumRotatableBonds | +0.179 | +0.282 |
| *(sanity)* P(permeable) vs measured log Papp | +0.632 | — |

Compare the pilot's `r(logP, heavy atoms) = 0.87`. Here the oracle–logP coupling is **0.165**.

**Adversarial growth test [MEASURED HERE].** Greasy tail grafted onto cyclosporin A, lengthening:

| added carbons | Crippen logP | heavy atoms | MW | **P(permeable)** |
|---|---|---|---|---|
| 0 (CsA itself) | 3.27 | 85 | 1203 | **0.8162** |
| 2 | 4.30 | 87 | 1233 | 0.8024 |
| 5 | 5.47 | 90 | 1275 | 0.7828 |
| 10 | 7.42 | 95 | 1345 | 0.7720 |
| 20 | 11.32 | 105 | 1485 | 0.7720 |
| 40 | 19.12 | 125 | 1766 | 0.7720 |
| 80 | **34.72** | 165 | 2327 | **0.7720** |

Monotone decrease, then a hard plateau. At the exact logP the pilot run reached (33–43), the
oracle score is *below* the unmodified parent.

Growing a homo-oligomer of the single most permeability-favourable residue type (N-Me-Leu)
**[MEASURED HERE]**:

| residues | Crippen logP | heavy atoms | HBD | TPSA | **P(permeable)** |
|---|---|---|---|---|---|
| 3 | 2.79 | 29 | 1 | 81 | **0.7776** |
| 5 | 4.54 | 47 | 1 | 122 | 0.7200 |
| 8 | 7.16 | 74 | 1 | 183 | 0.6106 |
| 15 | 13.27 | 137 | 1 | 325 | 0.5885 |
| 30 | 26.37 | 272 | 1 | 630 | 0.5885 |
| 50 | **43.83** | 452 | 1 | 1036 | **0.5885** |

Stacking the "best" residue is punished. Pure lipophilic blobs score at or below chance:
C40 alkane (logP 15.85) → 0.5005; C100 alkane (logP 39.26) → 0.5005; polystyrene 10-mer
(logP 21.81) → 0.4770; squalene-like (logP 10.60) → 0.3332. Median over the real database is 0.911.

**What the oracle's top-scored molecules actually look like [MEASURED HERE]:** median Crippen logP
of the top-50 is **1.56**, *below* the database median of 2.84; median heavy-atom count 77 vs 62;
the 15 highest-scoring entries sit at logP 0.7–2.3, MW 998–1197, HBD 4–5, TPSA 239–288 Å². The
oracle's optimum is a mid-polarity, mid-size macrocycle — which is also where the real oral bRo5
drugs sit (§2).

**Simulated optimiser [MEASURED HERE]** (`probe2.py`). Greedy hill-climb over residue substitution,
chain growth and chain shrinkage, from a random 8-mer, 13-residue vocabulary. Converged in 7 steps
at score 0.939 to a **10-residue peptide, Crippen logP −3.55, MW 947, HBD 7, TPSA 299 Å²**,
composition `{Sar:2, meAla:2, Thr, Gly, meVal, Phe, Ser, Pro}`. The optimiser went *away* from
lipophilicity and did not inflate length. Random linear peptides show the same: mean oracle score
falls monotonically with length, 0.673 at 6 residues → 0.355 at 30 residues.

### 1.4 Published feature-attribution and ablation work

- **CycPeptMP ablation** (Li et al., *Brief Bioinform* 2024, https://doi.org/10.1093/bib/bbae417):
  all three feature levels (atom / monomer / peptide) contribute. Stated plainly in the paper:
  SVM "could partially predict permeability by using lipophilicity descriptors, such as LogP, which
  are largely dependent on molecular weight" — i.e. the authors explicitly identify the logP/size
  route as the *weak* baseline (MAE 0.488), not the strong one. They also note cyclic-peptide
  permeation "negatively correlated with molecule size".
- **Systematic 13-model benchmark** (*J Cheminform* 17, 2025, https://doi.org/10.1186/s13321-025-01083-4):
  tried logP and TPSA as **auxiliary prediction tasks**. The auxiliary tasks were learned well
  (R² > 0.9) but gave "limited or no benefit" to permeability prediction. Their stated reason is
  directly on point: logP and TPSA in CycPeptMPDB are "calculated using fragment-based additive
  methods, which sum contributions from predefined chemical groups" and "do not account for 3D
  structure or the burial of polar surfaces". **The additively-decomposable descriptors carry
  almost none of the permeability signal.** Also: scaffold split yields substantially *lower*
  generalisability than random split; DMPNN is the best representation.
- **MultiCycPermea Grad-CAM** (Wang et al., *BMC Biol* 23:63, 2025,
  https://doi.org/10.1186/s12915-025-02166-2): on permeability-cliff pairs the model attends to
  **number of aromatic rings** (penalising), **phenol / polar functional groups** (penalising) and
  **all N-methylation sites** (rewarding), plus side-chain alkylation. In-distribution MSE 0.16 vs
  Multi_CycGT 0.29; filtered ~80% of low-permeability peptides in a downstream design study.
- **CPMP** (Frontiers in Bioinformatics 2025, https://doi.org/10.3389/fbinf.2025.1566174;
  https://github.com/panda1103/CPMP): MAT-based, R² 0.67 PAMPA / 0.75 Caco-2 / 0.62 RRCK / 0.73 MDCK.
  Stratified error analysis: best for MW 800–900 (R² 0.71) and **>1100 Da (R² 0.76)**, worse at
  ≤700 Da. Highest R² at **high logP > 4.0 (0.798)**, collapses for highly polar peptides
  (TPSA 350–400 Å²: R² 0.148). Accuracy *improves* for heavily modified peptides
  (modified-AA ratio > 0.6: R² 0.772–0.890).
- **PeptideCLM** (Feller & Wilke, *JCIM* 65:571, 2025; https://github.com/AaronFeller/PeptideCLM;
  weights https://huggingface.co/aaronfeller): SMILES transformer pretrained with ncAAs,
  D-α-carbons, N-methylation, PEGylation, all five cyclisation types. Permeability ROC-AUC 0.820
  on k-means holdout of CycPeptMPDB, beating ChemBERTa-2 (0.770). MSE 0.551 regression,
  ROC-AUC 0.781 classification at a −5.5 threshold.
- Other models for completeness: Multi_CycGT (*J Med Chem* 67:1888, 2024) acc 0.820 / ROC-AUC 0.865;
  PharmPapp MAE 0.317 R² 0.672 (KNIME-only); CyclePermea MAE 0.334; MuCoCP R² 0.503;
  MSF-CPMP acc 0.9062 / AUROC 0.9546 (https://github.com/wanglabhku/MSF-CPMP);
  CYCLOPS web tool acc 0.824, regression MAE 0.477 / R² 0.44 (*Digital Discovery* 2025,
  https://doi.org/10.1039/D4DD00375F).

### 1.5 Med-chem structure–permeability relationships (Lokey, Jacobson, Kessler, Kihlberg)

The mechanistic literature says the same thing as the numbers: permeability is a **conformational**
property, not a compositional one.

- **Permeability cliffs from stereochemistry alone.** Hewitt, Leung, Pye, …, Jacobson & Lokey,
  *JACS* 137, 2015 (https://pubmed.ncbi.nlm.nih.gov/25517352/): among cyclic hexapeptide
  diastereomers **6.6 and 6.7, differing by a single stereocentre at Leu5, Caco-2 permeability
  differs 60-fold.** Identical molecular formula, identical MW, identical HBD/HBA count, identical
  Crippen logP. No additive per-atom function can distinguish them. Their conclusion: 3D SASA and
  ΔG_desolv correlate highly with Caco-2 permeability, and "steric occlusion of polar groups from
  solvent can outweigh intramolecular hydrogen bonding".
- **N-methylation.** White, Renzelman, Rand, …, Jacobson & Lokey, *Nat Chem Biol* 7:810, 2011
  (https://doi.org/10.1038/nchembio.664): on-resin N-methylation of a 755 Da cyclic hexapeptide with
  **three N-methyl groups gave 28% oral bioavailability in rat**. Critically, selectivity and degree
  of N-methylation are "dependent on backbone stereochemistry" — conformation dictates where
  methylation even happens. Ovadia/Kessler had earlier shown that **both the number and, more
  importantly, the relative positions** of N-methyl groups drive permeability. Lokey's own library
  data shows the counter-examples: "some compounds with 4 N-Me groups showed very low permeability,
  while some non-N-methylated compounds appeared to be highly permeable". **[MEASURED HERE]** across
  CycPeptMPDB the raw count of N-methylated monomers correlates with permeability at only
  r = **+0.094** — confirming the count is not the driver, the pattern is.
- **Intramolecular H-bonding / chameleonicity.** Rezai, Bock, …, Lokey & Jacobson, *JACS*
  128:14073, 2006 (https://doi.org/10.1021/ja063076p) — conformational flexibility + internal
  H-bonding predicts relative permeability. Rossi Sebastiano et al., *J Med Chem* 61:4189, 2018
  (https://doi.org/10.1021/acs.jmedchem.8b00347) — dynamically exposed polarity governs both
  permeability and solubility. Danelius et al., *Chem Eur J* 26:5231, 2020. Chamelogk
  (Caron & Ermondi, *J Med Chem* 2023): median chameleonicity 0.19 for Ro5 vs 0.77 for bRo5.
- **The lipophilicity window, stated as a window.** Over/Ahlbach et al., "Charting islands of
  permeability", 2015 (https://europepmc.org/article/MED/25974856): 62 cyclic hexapeptides,
  Caco-2 + PAMPA. Model: "peptides with very high permeability have high lipophilicity **and** few
  solvent hydrogen bond interactions, whereas peptides with very low permeability have low
  lipophilicity **or** many solvent interactions." Two conditions, conjunctive — not a single axis.
  Hewitt et al. 2015 likewise report that side-chain substitution preserves permeability only
  "within a reasonable lipophilicity window", and that the highly lipophilic Nal substitution gave
  erratic PAMPA values "most likely due to solubility or aggregation issues" — the solubility wall
  (§2.4) bites before the lipophilicity benefit saturates.

### 1.6 The real risk with this oracle class: saturation, not inflation

**[MEASURED HERE]** fraction of CycPeptMPDB that PepINVENT's own oracle already scores above
threshold:

| threshold | peptides | fraction |
|---|---|---|
| > 0.90 | 4,135 / 7,991 | **51.7%** |
| > 0.95 | 3,307 | 41.4% |
| > 0.99 | 1,270 | **15.9%** |
| > 0.999 | 19 | 0.2% |

A classifier probability is a bounded, saturating objective. Half the training distribution is
already at >0.9, so a "best-of-dataset" baseline nearly solves it — exactly the failure GuacaMol
diagnoses for its trivial benchmarks (§5.3). Three mitigations, all cheap: (a) use a **regression**
oracle on log Papp, not a classifier probability — the *J Cheminform* 2025 benchmark found
"regression generally outperforms classification" on ROC-AUC anyway; (b) score against a
**scaffold-split** model, which that benchmark shows is substantially harder; (c) require the
score **jointly** with a solubility and a novelty constraint so the saturated region is not
by itself a win.

Second caveat, stated plainly: the plateaus in §1.3 (0.7720, 0.5885, 0.5005 repeated exactly) are
the signature of a **tree ensemble extrapolating**. Off-distribution, the oracle is not wrong in a
hackable direction, it is simply piecewise-constant and uninformative. That makes it hack-*resistant*
but it also means a benchmark built on it is only meaningful inside a chemically sane region; pair
it with a validity/alert filter as PepINVENT itself does (§6).

---

## SECTION 2 — The real bRo5 property space

### 2.1 Doak / Kihlberg: the measured outer limits

Doak, Over, Giordanetto & Kihlberg, "Oral Druggable Space beyond the Rule of 5",
*Chem Biol* 21(9):1115–1142, 2014 (https://doi.org/10.1016/j.chembiol.2014.08.013). 226 orally
administered drugs and clinical candidates with MW > 500.

| set | criteria | coverage |
|---|---|---|
| Lipinski Ro5 (1997) | MW ≤ 500, 0 ≤ cLogP ≤ 5, HBA ≤ 10, HBD ≤ 5 | reference |
| **extended Ro5** | MW ≤ 700, 0 ≤ cLogP ≤ 7.5, HBD ≤ 5, PSA ≤ 200 Å², NRotB ≤ 20 | 62% (141/226) |
| **outer limits of oral bRo5** | **MW ≤ 1000, −2 ≤ cLogP ≤ 10, HBD ≤ 6, HBA ≤ 15, PSA ≤ 250 Å², NRotB ≤ 20** | **93% (211/226)** |

The asymmetry is the whole point: relative to Ro5, the MW ceiling doubles (500 → 1000) and the
cLogP ceiling doubles (5 → 10), PSA goes up 1.8× (140 → 250), but **HBD rises only 1.2× (5 → 6)**.
HBD is the hard constraint. It is also the one an optimiser cannot inflate (floor at 0) — which is
why §1.1's HBD-dominated signal is good news for benchmark design.

Also from the same paper, and directly relevant to §5: **PSA scaled by MW is invariant across the
size range.** 95% of compounds < 700 Da have PSA/MW between 0.10 and 0.35; 92% of those > 700 Da
have PSA/MW between 0.15 and 0.30. A **ratio** is the size-invariant quantity, the raw PSA is not.
CNS penetration is tighter still: PSA ≤ 90 Å², MW ≤ 450, HBD ≤ 3.

bRo5 compounds score QED 0.16 (95% CI 0.14–0.17) vs 0.31 (0.29–0.32) for extended-Ro5 — both far
below the 0.49 QED of molecules medicinal chemists labelled "unattractive". **QED is actively
misleading in this chemical space and must not be used as a guard here.**

Follow-ups: Poongavanam, Doak & Kihlberg, *Curr Opin Chem Biol* 44:23–29, 2018
(https://pubmed.ncbi.nlm.nih.gov/29803972/); Caron, Kihlberg, Goetz, Ratkova, Poongavanam & Ermondi,
*ACS Med Chem Lett* 2021 (https://pmc.ncbi.nlm.nih.gov/articles/PMC7812602/).

### 2.2 Published multi-parameter scores for bRo5

**AB-MPS** — DeGoey, Chen, Cox & Wendt, *J Med Chem* 61(7):2636, 2018
(https://pubmed.ncbi.nlm.nih.gov/28926247/), derived from 1,116 AbbVie bRo5 compounds with PK data:

```
AB-MPS = |cLogD − 3| + NAR + NRotB
```

NAR = number of aromatic rings, NRotB = number of rotatable bonds, cLogD at pH 7.4.
**AB-MPS ≤ 14** separates acceptable (F > 27.3%) from low/moderate oral bioavailability;
**≤ 15** separates orals from parenterals.

This formulation is worth copying for three reasons:
1. The lipophilicity term is `|cLogD − 3|` — an **absolute value around an optimum**, so inflating
   logD *increases* the penalty. Non-monotone by construction.
2. NAR and NRotB both **grow with molecule size** and both carry a positive sign in a quantity
   being minimised. Size inflation is penalised twice.
3. It is the opposite of a "sum of favourable per-atom contributions".

Benchmarked on macrocyclic drugs by García Jiménez, Poongavanam & Kihlberg,
*J Med Chem* 2023 (https://exa.ai/library/publication/qlnyw1xrctw): AB-MPS ≤ 15 classified oral
macrocycles at 79%/61% sensitivity and parenterals at 71%/79% specificity across train/test sets.

**ETR = EPSA / TPSA** — "polarity reduction" descriptor, *J Med Chem* 2024
(https://pubmed.ncbi.nlm.nih.gov/38498697/), validated on ~1,000 compounds with human absorption
data plus ~10,000 AbbVie tool compounds (~1,000 PROTACs, ~7,000 Ro5, ~2,000 bRo5). A **ratio of a
measured 3D polarity to a computed 2D polarity** — the quantity is explicitly a measure of how much
polarity the molecule hides. Cannot be written as a sum over atoms (§5).

### 2.3 Cyclic-peptide-specific limits

- CycPeptMPDB v1.2 spans MW 342–1778 and TPSA 73–702 Å². The permeability threshold used
  throughout the field is log Papp ≥ **−6.0** (1×10⁻⁶ cm/s) for "permeable"; PepINVENT's classifier
  and Multi_CycGT use −6, PeptideCLM uses −5.5 (CycPeptMPDB, *JCIM* 63:2240, 2023,
  https://doi.org/10.1021/acs.jcim.2c01573).
- **Over 99.6% of CycPeptMPDB peptides contain non-natural amino acids** (CycPeptMP 2024). A
  benchmark on canonical residues only is not this chemical space.
- **[MEASURED HERE]** the top-5% most permeable peptides in CycPeptMPDB have median MW 1225
  (vs 885 overall), median heavy atoms 88 (vs 64), median HBD 5 (vs 4), median TPSA 285 Å²
  (vs 186), median Crippen logP 2.51 (vs 2.84). Note that this cut is partly source-confounded, but
  the direction is still informative: **the best real permeable cyclic peptides are large and
  relatively polar, and slightly *less* lipophilic than the database median.**

### 2.4 Permeability and solubility are documented as antagonistic, with numbers

- **Whitty, Zhong, Viarengo-Baker, Beglov, Hall & Vajda, *Drug Discov Today* 21(5):712–717, 2016**
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC5821503/). The cleanest quantitative statement in the
  literature, from 20–22 approved oral macrocyclic drugs:
  - good aqueous solubility requires **TPSA ≥ 0.2 × MW**
  - good passive permeability requires nonpolar-conformer **3D PSA ≤ 140 Å²**
  - "**one or other of these limits is almost invariably violated for compounds with MW > 600 Da**"
  - a compound at MW 1100 with TPSA 220 Å² must bury ≈80 Å², about one third of its polar surface,
    to permeate
  - **none of the approved oral macrocyclic drugs has TPSA ≤ 140 Å²**

  This is a genuine, published, quantified **frontier**: for MW > 600 the two objectives cannot both
  be satisfied by a static conformer, only by chameleonic behaviour. A benchmark task that scores
  solubility and permeability jointly at MW > 600 is therefore hard *for a stated physical reason*,
  not by arbitrary weighting.

- **Townsend / Lokey et al., "drug-like properties in macrocycles above MW 1000", 2020**
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC7719619/): matched cyclic decapeptide pairs. At
  ALogP 0.4 solubility is 1900 mg/mL; at ALogP 3.7 the rigid isomer falls to **7 mg/mL** — a ~270×
  solubility loss over 3.3 ALogP units, exactly where permeability peaks. The chameleonic isomer
  retains 210 mg/mL at the same ALogP. LogD–LogPe correlation R² = 0.93 in the soluble regime vs
  0.16 in the insoluble regime.
- **Naylor et al., *J Med Chem* 61:11169–11182, 2018** — "Lipophilic Permeability Efficiency
  Reconciles the **Opposing Roles of Lipophilicity** in Membrane Permeability and Aqueous Solubility".
- Cyclosporin A aqueous solubility is **23 µM at 25 °C** (Caron/Ermondi, *J Med Chem* 2023) — the
  canonical permeable macrocycle is barely soluble.
- PepINVENT's own paper frames this as the core design tension and reports that in its RL runs
  "the solubility component was learned first, in the first 100 steps, followed by improving the
  permeability in the soluble peptide space … **permeability was harder to optimize**".

---

## SECTION 3 — Solubility and aggregation oracles for peptides with non-canonical residues

**Headline: exactly one established tool natively handles ncAAs (CamSol-PTM, web-only, no code),
and the SMILES-based small-molecule solubility models are demonstrably broken on peptide-sized
molecules.** For a tight loop you will almost certainly have to train or fine-tune your own.

### 3.1 CamSol family (Sormanni / Vendruscolo, Cambridge)

- **CamSol-PTM** — Oeller et al., *Nat Commun* 14:7475, **2023**,
  https://doi.org/10.1038/s41467-023-42940-w (PMC10656490). First sequence-based method handling
  non-natural residues. Extends CamSol's hydrophobicity / charge / α-helix / β-sheet scales to
  modified amino acids; **hydrophobicity of a new modified AA is derived from RDKit Crippen**, so
  new residues are added from SMILES. Authors' own stated restriction: "currently restricted to
  modifications that are of similar size to canonical amino acids" — not lipids, not glycans.
  They scored all 10,000 modified AAs from Amarasinghe et al. and 40,000 Nrf2 variants, so
  throughput exists in principle.
- **Input**: sequence (one-letter + modified-AA codes). Not HELM. **No published support for
  backbone N-methylation, D-amino acids (chirality is invisible to the scales), or cyclisation.**
  That is a serious gap for cyclic-peptide work specifically.
- **Code: none released.** Web server only, https://www-cohsoftware.ch.cam.ac.uk/index.php/camsolptm,
  behind a registration/login wall (verified 2026).
- **Licence**: academic use via the registered server; commercial licence sold through Cambridge
  Enterprise, https://licensing.enterprise.cam.ac.uk/product/camsol-method (1-year single-site,
  price on application).
- **Accuracy**: Pearson r = 0.78 (PYY variants), 0.81 (18A), 0.58 (GLP-1); r = 0.6 for extended PYY
  with 4 unseen modified AAs; **average r = 0.72 over 37 synthesised variants**. Relative
  (intrinsic) solubility only, not absolute. Validation peptides are 18–36 residues
  (MW ≈ 2000–4300) — at or above the top of the MW 500–2000 window; **no reported error on
  MW 500–2000 macrocycles**.
- **VERDICT: not loop-usable.** No weights, no documented public API, registration-gated,
  academic-use ToS. Viable only as an offline ranking of a few thousand designs, or as a label
  source to distil into a local surrogate.
- Note: PepINVENT nevertheless uses CamSol-PTM as its solubility scoring component (§6), which means
  **the published PepINVENT MPO is not reproducible end-to-end without that licence.**

### 3.2 AGGRESCAN family

- **AGGRESCAN** (Conchillo-Solé et al., *BMC Bioinformatics* 8:65, 2007, http://bioinf.uab.es/aggrescan/)
  — FASTA only, 20 canonical residues, no code download.
- **A3D / A3D 2.0** (*NAR* 43:W306, 2015; *NAR* gkz321, 2019; https://biocomp.chem.uw.edu.pl/A3D2/)
  — input is a **PDB structure**, not sequence. Free, no login, documented RESTful API.
- **A3D standalone**: https://bitbucket.org/lcbio/aggrescan3d — **MIT licence**, Python 2.7,
  pip/conda/Docker (*Bioinformatics* btz143, 2019). Dynamic mode and mutation modelling need
  FoldX/Modeller (separately licensed).
- **Aggrescan4D** (*NAR* 52:W170–W175, **2024**, https://biocomp.chem.uw.edu.pl/a4d/) — adds
  pH 4.0–9.0, mmCIF/UniProt input, ~500k precalculated jobs. Free, no login.
- **Non-canonical support: none.** The scale is an experimentally derived aggregation-propensity
  scale for the 20 natural residues; HETATM residues have no scale value. (No tool named
  "AGGRESCAN-BIO" exists; the lineage is AGGRESCAN → A3D → A3D 2.0 → A3D-MODB → A4D.)
- **VERDICT: no** for ncAA peptides. MIT standalone is computationally loop-friendly and you could
  in principle add your own residue scale values, but it needs a 3D structure per call.

### 3.3 TANGO

Fernandez-Escamilla et al., *Nat Biotechnol* 22:1302, 2004. Sequence only, 20 canonical residues,
assumes a fully denatured chain at fixed 1 mM concentration. **Distributed as a compiled binary
only — the licensing page states they do not provide or sell source code**
(https://tango.crg.es/licensing-and-services). Free academic licence agreement required; commercial
licence via CRG. Validated on 179 + 71 literature peptides (2004); no ncAA capability.
**VERDICT: no** (no ncAAs). Fast enough to loop if canonical-only were acceptable.

### 3.4 WALTZ / WALTZ-DB 2.0

*NAR* 2020 (PMC6943037), http://waltzdb.switchlab.org/. **1,416 hexapeptides, 512 amyloid-forming**,
with TEM/ThT/FTIR annotations plus WALTZ/TANGO/PASTA scores and CORDAX steric-zipper models.
Downloadable CSV/Excel/JSON. **All entries are canonical L-hexapeptides** (N-acetylated /
C-amidated termini only); no D-residues, no ncAAs. WALTZ itself is a PSSM served from the SwitchLab
server, not redistributable.
**VERDICT: no as an oracle; yes as training data for a hexapeptide aggregation classifier.**

### 3.5 Protein-scale solubility predictors — all FASTA-only, all canonical-only

| tool | year | code / weights | accuracy |
|---|---|---|---|
| Protein-Sol (*Bioinformatics* 33:3098) | 2017 | web only | ACC 0.516, MCC 0.03 on SoluProt independent set |
| SOLart (*Bioinformatics* 36:1445) | 2020 | web server, structure-based | regression, R² typically 0.4–0.6 for this class |
| SoluProt (*Bioinformatics* 37:23) | 2021 | web + standalone | **ACC 58.5%**, AUC 0.62, MCC 0.17 |
| **NetSolP** (*Bioinformatics* 38:941) | 2022 | **BSD-3-Clause**, https://github.com/TviNet/NetSolP-1.0, 5.63 GB weights | ESM1b-based, SOTA on strictly partitioned sets |
| DeepSoluE (*BMC Biol* 21:12) | 2023 | web only | ACC 0.595, MCC 0.19 |
| GraphSol / **GATSol** (*BMC Bioinformatics* 2024) | 2021/2024 | GATSol public, needs AlphaFold structures | **R² 0.517** (eSOL), 0.424 (yeast_108) |

Every one is trained on *E. coli* heterologous-expression solubility of **whole proteins**, uses a
20-letter alphabet, and has **no reported figure on MW 500–2000**. NetSolP is the only one with a
permissive licence plus downloadable weights, so it is the only realistic fine-tuning backbone.
**VERDICT: no as-is; NetSolP = with caveats as a backbone.**

### 3.6 SMILES-based aqueous solubility models — and the quantified reason they fail here

- **AqSolDB** (Sorkun et al., *Sci Data* 6:143, **2019**, https://doi.org/10.1038/s41597-019-0151-1).
  MW distribution of the full 9,982-compound set, computed from the data:
  mean **266.7**, median 228.7, p90 422.4, p95 552.7, p99 959.2, max 5299.5.
  **MW 500–2000: 613 compounds = 6.1%. MW ≥ 1000: 85 compounds = 0.85%. 93.8% is below MW 500.**
  Anything trained on AqSolDB is a sub-MW-500 model.
- **ESOL** (Delaney, *J Chem Inf Comput Sci* 44:1000, **2004**): 2,874 compounds, 4-descriptor linear
  regression, AAE 0.75 log overall / 0.83 on a 528-compound blind set. Its own "Large" subset has
  **mean MW 341 ± 71** — ESOL has never seen a peptide.
- **SolTranNet** (Francoeur & Koes, *JCIM* 61:2530, **2021**): https://github.com/gnina/SolTranNet,
  **Apache-2.0**, `pip install soltrannet`, weights ship with the package, CLI + Python API.
  Scaffold-split CV RMSE **1.459** on AqSolDB; **1.711** on a withheld SC2+FreeSolv set; 94.8%
  sensitivity as an insolubility filter. Trained purely on AqSolDB, so it inherits the MW
  distribution above. **No reported accuracy on peptides.**
- **Solubility Challenges**: SC1 best models RMSE 0.7–1.1 log despite ~0.05 log experimental error.
  SC2 (2019, 132 drug-like molecules): mean of 37 submissions RMSE **1.14** (tight) / **1.62**
  (loose); best post-hoc model RMSE 0.86, R² 0.54 (Conn et al., *JCIM* 2023,
  https://doi.org/10.1021/acs.jcim.2c01189). Authors conclude the limiting factor is "the deficiency
  of QSPR methods", not experimental uncertainty.
- **The key citation for extrapolation failure** — **Avdeef & Kansy, "Flexible-Acceptor General
  Solubility Equation for beyond Rule of 5 Drugs", *Mol Pharm*, 2020**,
  https://doi.org/10.1021/acs.molpharmaceut.0c00689. On a 32-molecule bRo5 test set including
  cyclosporine A, **gramicidin A, leuprolide, nafarelin, oxytocin** and vancomycin:
  **traditional GSE RMSE = 3.0 log units**; their flexible-acceptor GSE = **1.10**; RF regression = 1.07.
  Corroborated by Ermondi, Poongavanam, Vallaro, Kihlberg & Caron, *ADMET & DMPK* 8(3):207–214, 2020,
  https://doi.org/10.5599/admet.834: GSE and the Abraham solvation equation "failed to predict the
  solubility of the larger compounds in bRo5 space"; bRo5 solubility modelling "is in its infancy".
- **VERDICT**: SolTranNet/Chemprop-class models are **technically** loop-ready (permissive licence,
  local weights, ms per molecule) but **scientifically out of domain by roughly 2 log units** on
  peptide-sized molecules. Use as a crude insolubility tripwire only, never as the objective.

### 3.7 Peptide-specific / ncAA-aware models — the realistic options

- **PeptideCLM** (Feller & Wilke, *JCIM* 65:571, **2025**; https://github.com/AaronFeller/PeptideCLM;
  weights https://huggingface.co/aaronfeller; pretraining data Zenodo 10.5281/zenodo.14194470).
  SMILES-in RoFormer/BERT, 44M params, 586-token SMILES-pair vocabulary, pretrained on 23M molecules
  including 10M generated peptides with ncAAs (SwissSidechain), **10% D-α-carbons, N-methylation on
  20% of backbone amines, PEGylation, all five cyclisation types**. Fine-tuned in the paper on
  **permeability, not solubility**. **VERDICT: yes, with caveats** — best available encoder for this
  exact input class; you supply solubility labels and the head.
- **PepLand** (arXiv 2311.04419, 2023; *Brief Bioinform* 2025; https://github.com/zhangruochi/PepLand).
  Multi-view heterogeneous GNN (atom + fragment + junction), two-stage pretraining
  (canonical → non-canonical). **Solubility is a reported downstream task**, alongside permeability,
  binding, synthesizability. **VERDICT: yes, with caveats** — cheap forward passes; check the repo
  licence.
- **PeptiVerse** (bioRxiv 2025.12.31.697180, **Jan 2026**, PMC12773018, CC-BY-4.0; HF org
  https://huggingface.co/ChatterjeeLab/PeptiVerse). Unified therapeutic-peptide property platform
  that **accepts either amino-acid sequence or chemically modified peptide SMILES**, covering
  hemolysis, **solubility**, non-fouling, toxicity, permeability, half-life, binding. Non-canonical
  peptides split by Morgan-fingerprint Tanimoto; canonical by MMseqs2. Explicitly positioned as
  "guidance oracles within generative modeling workflows". **VERDICT: the most promising drop-in**
  — verify per-task metrics and solubility-dataset provenance before trusting it. Preprint.
- **SolPepBench / PepSol2000** (https://github.com/Ascaris-Equi/PepSol2000, MIT code / CC-BY-4.0
  data, currently an anonymous review artefact). Solvent-conditioned peptide solubility: 1,337 unique
  sequences, 4,405 peptide–solvent records, 7 solvent conditions, binary labels, official
  pair-stratified / sequence-disjoint / solvent-held-out splits. Useful **training data**; baselines
  deliberately weak. Linear peptides.
- **Peptide self-assembly (aggregation)**: Martini CG-MD **aggregation propensity (AP)** score —
  dipeptides (Frederix et al., *JPCL* 2:2380, 2011), all 8,000 tripeptides + hydrophobicity-corrected
  **APH** (Frederix et al., *Nat Chem* 7:30, **2015**); active-learning extension to tetra/penta/hexa
  (Van Teijlingen & Tuttle, *JCTC* 17:3221, 2021, ~50× speedup, filters on logP < 0 for solubility);
  **AI-expert** MCTS + RF + MD (*Nat Chem* 14, **2022**, https://doi.org/10.1038/s41557-022-01055-3).
  Accepts any residue parameterisable in Martini, which includes D-residues and many ncAAs — but each
  ground-truth call is a CG-MD run (minutes–hours). **VERDICT: with caveats** — the RF/SVM surrogate
  is the oracle; the MD is the label generator.

### 3.8 Does CycPeptMPDB carry solubility?

**No.** v1.2 = 7,991 cyclic peptides, 385 monomers, 56 sources. Stored per peptide: **membrane
permeability only** (PAMPA 7,298 / Caco-2 1,332 / RRCK 186 / MDCK 64), HELM, SMILES, 3D conformers,
plus *computed* RDKit descriptors and measured EPSA where reported. MW 342–1778 — exactly the right
window, the wrong property. No cyclic-peptide database carrying experimental aqueous solubility was
found; SolPepBench is the closest and is linear-peptide, binary, solvent-conditioned.

---

## SECTION 4 — Are AMP activity predictors hackable?

**Yes, demonstrably — far more so than cyclic-peptide permeability.** The strongest single result:
**BATTLE-AMP (bioRxiv 2026)** shows **9 of 21 benchmarked model variants have false-positive rate
> 0.80 on shuffled decoys with identical amino-acid composition**, and **no model exceeds MCC 0.17**
on that task. Caveat: BATTLE-AMP and GenPept-2025 are preprints, not yet peer-reviewed.

### 4.1 Documented degenerate / adversarial generation

- **BATTLE-AMP**, Szymczak et al., bioRxiv 2026,
  https://www.biorxiv.org/content/10.1101/2026.06.19.733349v1 (doi 10.64898/2026.06.19.733349) —
  the central oracle-hacking evidence, detailed in §4.2 and §4.5.
- **Beierle et al.**, bioRxiv 2025 / 2026, https://www.biorxiv.org/content/10.1101/2025.10.29.685317
  (PubMed 42106831). They must **explicitly filter generated sequences for "no recurrence of
  individual amino acids more than seven times"** — i.e. homopolymer runs — plus length 3–36. The
  fraction removed differs by model. They also report that "evaluating AMP predictors on data biased
  toward specific peptide properties revealed strong, model-specific preferences" and warn in-silico
  prediction "should not be used as a sole indicator for candidate selection".
- **AMP-Designer**, *Sci Adv* 2024, https://www.science.org/doi/10.1126/sciadv.ads8932 — the
  conditional-token baseline "produces invalid sequences, such as '[mask]'"; full fine-tuning
  overfits and "generates sequences that exhibit repetition either within themselves or to the real
  AMPs".
- **ProDCARL**, arXiv 2602.00157 (2026), https://arxiv.org/html/2602.00157v1 — explicit reward-hacking
  discussion: entropy regularisation + early stopping used "to preserve diversity and reduce reward
  hacking". Mean pAMP only 0.081 → 0.178; joint hit rate 6.3% (pAMP > 0.7, pTox < 0.3); no wet-lab
  validation.
- **ApexAmphion / deep-RL antibiotic platform**, arXiv 2509.18153 (2025),
  https://arxiv.org/pdf/2509.18153 — PPO against an ApexMIC oracle produces measurable compositional
  drift: **lysine frequency +0.088**, Gly and Leu each > +0.03, Glu −0.030, D/Q/S/T each > −0.02 vs
  natural AMPs. Cationic-residue inflation driven by the reward, caught only because they measured it.
- **Krishnan et al., IEEE SPW 2026** (doi 10.1109/SPW72489.2026.00013) and *Attack-as-Design*,
  https://github.com/armandyam/attack-as-design — adversarial search against a surrogate *is* the
  design loop. They constrain net charge, hydrophobicity, aromatic/charged fraction, sequence
  complexity and **longest homopolymer run** to natural ranges precisely because unconstrained search
  escapes them (on GFP: 100% of designs stay in realistic biophysical range vs 5–19% for baselines).
- **Genotypic Triggers**, arXiv 2608.06779 (2026), https://arxiv.org/html/2608.06779 — backdoored AMP
  generators raise predicted immunogenicity risk for HLA-allele carriers by **743% on average** while
  retaining high predicted potency and low predicted toxicity, **passing conventional safety screens**.
- **cdGAN**, *Brief Bioinform* 2025, https://pubmed.ncbi.nlm.nih.gov/41137855/ — FBGAN's
  classifier-feedback loop "may introduce bias and constrain the diversity and quality of the
  generated peptides".

### 4.2 Feature over-reliance, with numbers

- **BATTLE-AMP 2026**: a **hand-written two-feature filter** (net charge ≥ +3 AND hydrophobic moment
  above library median) achieves **Precision@100 = 0.12, LR+ = 7.5** on the leakage-free SLAY library
  (n ≈ 438,000, ~1.8% actives) — **outperforming trained binary classifiers**, which fall to LR+ ≈ 1.0
  there (AMPlify: 5 actives per 100 top-ranked; MBC-Attention: 4). The authors name the problem
  "compositional shortcuts".
- **Mooney et al., *PLoS ONE* 2012**,
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0045012 — AntiBP2 and CAMP cannot
  separate bioactive peptides from their own scrambles: **MCC −0.01 to 0.06**. PeptideRanker MCC 0.11
  (short) / 0.34 (long).
- **Porto et al., "Sense the Moment", *BBA Gen Subj* 2020**,
  https://doi.org/10.1016/j.bbagen.2020.129707 — "any arbitrary sequence with a similar composition to
  an AMP would be predicted as an AMP". A **single scalar** (geometric hydrophobic moment, no ML)
  gives > 80% accuracy and > 90% sensitivity, exceeding deep-learning methods on their set; only
  24/388 (6.19%) shuffles increased GHM.
- **Lobanov et al. 2023** — amino-acid-composition-only method performs "at the level of the best
  methods" (~80% accuracy/specificity) vs BERT-MLP/BERT-LA.
  Server http://bioproteom.protres.ru/antimicrob/
- **Lin et al., *iScience* 2023**, https://www.cell.com/iscience/fulltext/S2589-0042(23)02327-1 —
  forward feature selection: for the antibacterial class "all the features were charge-based except
  the PAAC-related amino acids"; discriminative residues C and E; AUC 0.849–0.913.
- **AmPEP**, *Sci Rep* 2018, https://www.nature.com/articles/s41598-018-19752-w — reducing 105
  descriptors to 23 gives comparable performance on all metrics except some precision loss.

### 4.3 Benchmark bias

**Sidorczuk et al., *Brief Bioinform* 23(5):bbac343, 2022**,
https://academic.oup.com/bib/article/23/5/bbac343/6672903 (PMC9487607). 660 models = 12 architectures
× 1 positive set × 11 negative-sampling methods × 5 runs. Findings:

- Models score much better when train and benchmark negatives come from the **same sampling method**.
  Mean AUC uplift: SVM-LZ +2.3%, AmpGram +3.6%, CS-AMPPred +4.4%, AMAP/MACREL +7.5%,
  ampir/AmPEPpy ~+9.5%, **all remaining architectures > +10%**. AmPEP and iAMP-2L were exceptions and
  were simply bad (mean AUC 0.65).
- Mean AUC correlates negatively with train/benchmark amino-acid-composition difference
  (**Spearman ρ = −0.53**, p < 2.2e-16) and with median-length difference (**ρ = −0.44**).
- Only **8 of 26** surveyed models met the minimal "bronze" reproducibility standard (~70%
  non-reproducible). Server: http://BioGenies.info/AMPBenchmark

Related:
- **BATTLE-AMP 2026**: of **48 published methods surveyed, fewer than 25% were reproducible**.
  HomologySplit (CD-HIT at 80/60/40% identity): median MCC drops 0.37 → 0.32; AMPredictor 0.54 → 0.42;
  HydrAMP_AMP 0.38 → 0.28. On 1–10 aa peptides most models fall below MCC 0.30. Thresholds used:
  active MIC ≤ 32 µg/mL, inactive MIC ≥ 128 µg/mL.
- **GenPept-Curated-2025**, bioRxiv 2026,
  https://www.biorxiv.org/content/10.64898/2026.04.25.720793v1 — 11,000 sequences (5,500/5,500),
  IPG-based deduplication; documents length-dependent class imbalance (AMP fraction **14.2%** in
  10–50 aa bins vs **77.1%** in 101–150 aa) as a distribution shift benchmarks ignore.
- **AMPlify**, *BMC Genomics* 2022, https://doi.org/10.1186/s12864-022-08310-4 — authors deliberately
  length-match negatives "so that the model did not learn to distinguish classes based on sequence
  lengths", and concede "all the publicly available AMP prediction tools face difficulty in
  differentiating between AMPs and non-AMPs that are highly similar in their sequences".
- **Gabere & Noble**, *Bioinformatics* 2017, https://pmc.ncbi.nlm.nih.gov/articles/PMC5860510/ —
  independent cross-tool evaluation; AntiBP outperformed its successor AntiBP2.

### 4.4 Wet-lab validation rates

| study | generated / screened | synthesised | active | threshold |
|---|---|---|---|---|
| Das et al., *Nat Biomed Eng* 2021 (CLaSS), https://www.nature.com/articles/s41551-021-00689-x | ~90k sampled, filtered | 20 | **2 potent**; HydrAMP's recount 3/21 = **14%** | MIC ≤ 128 µg/mL |
| Szymczak et al., *Nat Commun* 2023 (HydrAMP), https://www.nature.com/articles/s41467-023-36994-z | 900/prototype | 7 (no preselection) | 2 = **29%** | active vs *E. coli* |
| HydrAMP + classifier consensus + MD preselection | 180 → top 24 | 24 | 23 = **96%** @128; **9/24 = 38%** @ ≤32 µg/mL | — |
| Wan et al., *Nat Biomed Eng* 2024 (APEX de-extinction), https://pubmed.ncbi.nlm.nih.gov/38862735/ | 10,311,899 mined; 37,176 flagged | 69 | 41 = **59%** (prior scoring fn: **24%**) | MIC ≤ 128 µmol/L |
| Torres/Maus et al., *Nat Mach Intell* 2026 (ApexGO), https://www.nature.com/articles/s42256-026-01237-5 | 10 templates + intermediates | 100 | 86 = **86%**; 68% improved on template | MIC ≤ 64 µmol/L |
| Santos-Júnior et al. (AMPSphere), *Cell* 2024 | 863,498 c_AMPs | 50 | 27 = **54%** vs pathogens; 32 = 72% incl. commensals; **none vs MRSA** | complete growth inhibition |
| Ma et al., *Nat Biotechnol* 40:921, 2022, https://doi.org/10.1038/s41587-022-01226-0 | 2,349 candidates | 216 | 181 = **>83%** | — |
| Huang et al., *Nat Biomed Eng* 2023 (SMEP) | 10¹¹ 6–9-mers | 55 | 54 | — |
| PepGAN, ChemRxiv 2020, https://doi.org/10.26434/chemrxiv.12116136 | — | 6 | 1 highly active (MIC 3.1 µg/mL) | — |

**The failures matter more than the headlines.** Das et al.'s Supplementary Tables 7–8 list measured
MICs for all 20 AI-designed candidates: **most are > 1000 µg/mL** against both *S. aureus* and
*E. coli*, with only three showing anything. HydrAMP: all four generated Temporin-A "positive"
analogues were inactive (MIC ≥ 512 µg/mL). The high hit rates above are **after** aggressive
multi-classifier consensus, MD and expert handpicking; HydrAMP's unfiltered round-1 rate was 29%.

### 4.5 False-positive characterisation

- **BATTLE-AMP 2026, Fig. 4** — three decoy families.
  *SyntheticRandom* (uniform) is rejected by most models.
  *SyntheticRealistic* (matching global AMP residue frequencies) already breaks binary classifiers:
  **AMPlify FPR 0.79, AMPpred-MFA 1.00, SenseXAMPclf 0.84**.
  *SyntheticShuffled* (per-sequence permutations of real actives, **composition exactly preserved**):
  **9 of 21 variants exceed FPR 0.80**, including the top performers **MBC-Attention 0.87 and
  AMPredictor 0.89**; **no model achieves MCC > 0.17**. Most classifiers assign near-identical P_AMP
  to actives and their shuffles. ActivityCliff set (46 pairs, median 86% pairwise identity): ML
  regressors ρ = −0.13 to −0.25; MD descriptors near chance. Both ML and physics fail.
- **Bello-Madruga & Torrent Burgas**, *Comput Struct Biotechnol J* 23:972–981, 2024,
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10884422/ — 14 IDR-derived peptides with AMP-compatible
  composition, each predicted active by ≥ 3 of 5 tools (CAMPR4, AI4AMP, AmpGram, HydrAMP,
  sAMPpred-GAT): **12/14 showed no activity up to 100 µM**. Even the best predictor, HydrAMP,
  "failed to accurately predict 50% of the IDR-derived peptides". Cause: failure to fold, not
  composition.
- **Hemolysis is the real constraint, and it is driven by the same features the classifiers reward.**
  Chen et al. 2005 / Jiang et al., *Chem Biol Drug Des* 2011 — peptide D1(V13): MIC geomean 2.9 µM,
  MHC 5.2 µM → **therapeutic index 1.8**. A single V13K "specificity determinant" cut hemolysis
  > 32-fold and raised TI **90-fold** (Gram-neg) — with charge *increased*. "The vast majority of
  native AMPs are very hemolytic."
  Mourtada et al., *Nat Biotechnol* 2019 (StAMPs), https://pmc.ncbi.nlm.nih.gov/articles/PMC7437984/ —
  a "lyticity index" from i,i+3 / i,i+4 hydrophobic networks correlates with % hemolysis at
  **R² = 0.995**; LI < 600 captured all StAMPs with < 10% hemolysis at 25 µg/mL.
  Edwards et al., *ACS Infect Dis* 2016, https://pubs.acs.org/doi/full/10.1021/acsinfecdis.6b00045 —
  "Antimicrobial activity increased with amphipathicity, but unfortunately so did toxicity"; TIs span
  3 (thanatin) to 335 (tachyplesin-1).
  Greco et al., *Sci Rep* 2020, https://www.nature.com/articles/s41598-020-69995-9 — TIs 50–500 across
  24 peptides, with a caution that erythrocyte species choice inflates reported TIs.
- **Sequence order matters at fixed composition**: Wu et al., *Antibiotics* 14:1077, 2025,
  https://www.mdpi.com/2079-6382/14/11/1077 — 12 permutations of KKWRKWLKWLAKK with identical MW,
  net charge and hydrophobicity differ substantially in therapeutic index.

### 4.6 Non-canonical amino acids in AMP prediction

Almost nothing mainstream handles them.
- **DBAASP** (https://dbaasp.org/) is the only AMP database recording per-residue NCAA identity and
  position, but its own predictors require the 20 canonical residues. Das et al. 2021 explicitly
  stripped B/J/O/U/X/Z before training CLaSS.
- **AmpHGT, 2025**, https://pmc.ncbi.nlm.nih.gov/articles/PMC12217533/ — heterogeneous-graph model
  handling **> 600 NCAAs** plus the 20 CAAs via SMILES fragment views. Its dataset is itself a
  cautionary tale: NCAA-AMPs contained **486 unique NCAAs vs 84 in non-AMPs**, so the authors had to
  restrict training to the 26-NCAA intersection (340/340) **to stop the model cheating on NCAA
  identity alone**; CD-HIT could not be applied because of the expanded alphabet.
- **MODAN, 2023**, https://github.com/ycu-iil/MODAN — multi-objective Bayesian optimisation over
  Morgan/MACCS fingerprints, covering α,α-disubstituted NPAAs (Aib, Ac5c, Ac6c), Dab, Orn,
  side-chain stapling; 123,390 candidates, 7 synthesised in round 1, all low-hemolysis;
  starting models r = 0.70–0.88 on 82 peptides.
- **ADAPT/QLAPD-D, 2026** — first curated benchmark of paired L/D-substituted AMPs; before it "no
  publicly accessible benchmark links D-substituted AMP sequences with corresponding antimicrobial
  activity assays", and AlphaFold2 "yields unreliable geometries for D-chiral backbones". 42.5% of
  paired peptides improved on D-substitution.
- Standard PLM predictors (ESM2, ProtTrans) are structurally incapable here — UniRef, 20-letter
  alphabet.

### 4.7 Verdict for benchmark purposes

An AMP classifier is a **bad** "hard realistic oracle" on its own: it is gameable by composition
alone, its false-positive rate on composition-matched shuffles exceeds 0.8 for the best models, and
its wet-lab unfiltered hit rate is ~14–29%. If used at all, it must be paired with
(a) composition-matched shuffled decoys as hard negatives, (b) a hemolysis/therapeutic-index
counter-objective, and (c) a homology-aware split. BATTLE-AMP's SyntheticShuffled protocol is the
ready-made hack-resistance test.

---

## SECTION 5 — Non-decomposable objectives

The degenerate strategy the pilot hit — stack more atoms, each contributing positively — works
exactly when the objective is **`f(molecule) = Σ_atoms g(atom)`** or close to it, and is unbounded
above. Crippen logP is literally this: a sum of 68 atom-type contributions (Wildman & Crippen,
*JCIM* 39:868, 1999). So is TPSA (Ertl, sum of fragment contributions), molar refractivity, heavy
atom count, and most fragment-based counts.

The resistant objectives are those where that factorisation fails.

### 5.1 PMO's explicit statement about decomposability

Gao, Fu, Sun & Coley, "Sample Efficiency Matters: A Benchmark for Practical Molecular Optimization",
*NeurIPS* 35:21342–21357, 2022; arXiv 2206.12411; https://github.com/wenhao-gao/mol_opt.

The verbatim claim, in their oracle-landscape analysis:

> "Isomer-type oracles are summations of atomic contribution, while all other MPOs are mainly based
> on similarity measured by fingerprints"

and the consequence they measure: **string-based GAs (SMILES GA, STONED) reach superior relative
performance specifically on the isomer-type tasks** (`isomer_c7h8n2o2`, `isomer_c9h10n2o2pf2cl`,
`sitagliptin_mpo`, `zaleplon_mpo`), while similarity-based MPOs cluster separately and ML bioactivity
oracles (DRD2, GSK3β, JNK3) fall into the similarity cluster. **Decomposability of the oracle changes
which algorithm wins** — that is PMO's empirical evidence that the property matters.

PMO also notes that several oracles "have been described as 'trivial'", asserting this "is only true
when the number of oracle queries is not controlled", and finds QED, DRD2 and `osimertinib_mpo`
solvable within hundreds of calls.

### 5.2 Taxonomy of resistant objective forms, with concrete instances

**(a) Ratios and normalised quantities.** A ratio `A/B` where both numerator and denominator are
extensive cannot be increased by scaling the molecule, because both grow together.
- `PSA / MW`, invariant at 0.10–0.35 across the whole oral drug range (Doak et al. 2014, §2.1).
- **ETR = EPSA / TPSA** (*J Med Chem* 2024) — "polarity reduction", explicitly a ratio.
- Whitty's `TPSA ≥ 0.2 × MW` solubility criterion (2016) is a ratio threshold, not an absolute one.
- **[MEASURED HERE]** `TPSA/HeavyAtomCount` correlates with cyclic-peptide permeability at −0.217 vs
  −0.094 for raw TPSA; `MolLogP/HAC` at +0.156 vs +0.092 for raw logP. Normalising **improves** the
  signal while removing the inflation route.
- Ligand efficiency and LipE/LLE are the classical med-chem instances of the same idea.

**(b) Objectives with an interior optimum (absolute value / window / Gaussian).**
- **AB-MPS = |cLogD − 3| + NAR + NRotB** (DeGoey et al. 2018). Penalises deviation from logD 3 in
  *either* direction, and adds positive penalties for two size-correlated counts.
- GuacaMol's modifier vocabulary: `Gaussian(target, σ)`, `MinGaussian`, `MaxGaussian` — e.g.
  `logP → Gaussian(8, 2)` or `MinGaussian(5.0, 2)` in CNS MPO (Brown et al., *JCIM* 59:1096, 2019).
- REINVENT/PepINVENT `double_sigmoid` and `step` transforms (§6) — PepINVENT's own config uses
  `lipophilicity → step(low=1, high=4)`, i.e. a **window**, with weight 0.
- Any desirability function peaking at an interior value kills a monotone exploit by construction.

**(c) Properties depending on 3D conformation / intramolecular interactions, not composition.**
This is the strongest class, and the cyclic-peptide literature is the best-documented example.
- **Lokey's permeability cliffs**: compounds 6.6 and 6.7 differ **only in one stereocentre** —
  same formula, same MW, same HBD, same Crippen logP, **60-fold Caco-2 difference**
  (Hewitt et al., *JACS* 2015). No per-atom sum can represent this; a `useChirality=True` fingerprint
  can at least see it, which is exactly what PepINVENT uses.
- **3D PSA in a nonpolar environment** vs 2D TPSA: Whitty et al. 2016 require 3D PSA ≤ 140 Å² while
  TPSA ≥ 0.2×MW — the two are different functions of the same atoms, and their gap *is* chameleonicity.
- **EPSA** (chromatographically measured exposed polar surface area) and **Chamelogk**
  (Caron & Ermondi, *J Med Chem* 2023) are experimental, conformation-dependent, non-additive.
- The *J Cheminform* 2025 benchmark's finding that additive logP/TPSA auxiliary tasks give no benefit
  to permeability prediction (§1.4) is direct evidence that the permeability signal is in the
  non-additive part.
- Intramolecular H-bond count, ΔG_desolv, 3D SASA — all require a conformer ensemble.

**(d) Objectives with internal consistency constraints / conjunctions.**
- Whitty's frontier: solubility needs TPSA ≥ 0.2·MW **and** permeability needs 3D PSA ≤ 140 Å², and
  for MW > 600 "one or other of these limits is almost invariably violated". A **conjunction of two
  antagonistic constraints**, with a published MW at which it becomes infeasible without chameleonicity.
- Over/Ahlbach 2015: high permeability needs high lipophilicity **and** few solvent H-bond
  interactions; low permeability follows from low lipophilicity **or** many solvent interactions.
- Isomer tasks (exact molecular formula) are constraints, but PMO shows they are *also* additive and
  therefore easy for string GAs — a constraint alone is not sufficient.

**(e) Reference-dependent objectives (similarity, rediscovery, scaffold-hopping).**
PMO's "all other MPOs are mainly based on similarity measured by fingerprints". Tanimoto similarity
to a fixed reference is bounded at 1 and **decreases** when you add atoms the reference lacks — so
growth is self-penalising. GuacaMol's `deco_hop`, `scaffold_hop`, `valsartan_smarts` add substructure
presence/absence requirements.

**(f) Order/permutation-dependent objectives.** At fixed composition, permuting a sequence changes
the property: Wu et al. 2025 show 12 permutations of one 13-mer with identical MW, net charge and
hydrophobicity differ substantially in therapeutic index; BATTLE-AMP's SyntheticShuffled decoys use
exactly this to expose classifiers (§4.5). Any objective that distinguishes a sequence from its
permutations is, by definition, not a sum over residues.

### 5.3 GuacaMol's own treatment of trivial objectives

Brown, Fiscato, Segler & Vaucher, *JCIM* 59(3):1096–1108, 2019;
https://doi.org/10.1021/acs.jcim.8b00839; arXiv 1811.09621; https://github.com/BenevolentAI/guacamol.

Appendix 7.5/8.5, "Results for Trivial Goal-Directed Benchmarks", is the relevant passage. They
report that for `logP (target −1.0)`, `logP (target 8.0)`, `TPSA (target 150)`, `CNS MPO`, `QED`,
`C7H8N2O2` and `pioglitazone MPO`, **"best of dataset" alone scores 1.000 on four of seven**, and
all baselines reach ≈1.000. Their conclusion:

> "it is evident that such trivial benchmarks are not suitable for the assessment of generative
> models … The ChEMBL database alone achieves excellent scores already"

and the design criterion they state: the difficulty is "to find tasks that are hard enough, that is,
objectives that are not already satisfied by molecules in the training dataset." **This is exactly
the saturation problem measured in §1.6** — and it is a *different* failure from size inflation.
Note also that GuacaMol's logP tasks use `Gaussian(target, σ)` modifiers, i.e. they target a value
rather than maximise — the unbounded "penalized logP" task used elsewhere in the literature is the
one that produces absurd long-chain molecules.

GuacaMol states three further hack-resistance rationales explicitly:

- **Why geometric mean.** "With a geometric mean, all requirements … must be met, at least partly",
  whereas an arithmetic mean lets a molecule score well with one contribution at zero.
- **Why similarity-based rather than ML oracles.** Tanimoto's advantage is that "deficiencies in the
  ML models … cannot be as easily exploited by the generative models." Directly relevant: a learned
  permeability oracle is the exploitable kind, which is why §1.3's probe mattered.
- **Why no realism penalty term.** "A seemingly obvious solution is to add a contribution penalizing
  unrealistic molecules to the scoring function. Unfortunately, it is not trivial to encode what
  medicinal chemists would consider as molecules that are acceptable."

Instead GuacaMol applies **post-hoc quality filters** (Walters' `rd_filters` with
SureChEMBL/Glaxo/PAINS/in-house rule sets, https://github.com/PatWalters/rd_filters), explicitly
framed as "high precision, low recall". Top-100 pass rates: **SMILES LSTM 77%** (equal to the
best-of-ChEMBL baseline), **Graph GA 40%, SMILES GA 36%, Graph MCTS 22%**. The ranking inverts the
optimisation leaderboard — **the strongest optimisers are the worst hackers**, quantified.

Related critiques:
- **Renz, Van Rompaey, Wegner, Hochreiter & Klambauer, "On Failure Modes in Molecule Generation and
  Optimization", *Drug Discov Today Technol* 32–33:55–63, 2019**,
  https://doi.org/10.1016/j.ddtec.2020.09.003 (https://pubmed.ncbi.nlm.nih.gov/33386095/) — two
  results. (a) The **AddCarbon** model: adding a single carbon to a random training molecule scores
  perfectly on novelty, validity and uniqueness and highly on KL divergence. (b) **The control-model
  divergence test** — train an *independent* model on the same data; during goal-directed generation
  the optimisation score rises while the control score **diverges**, proving the generator is
  exploiting features specific to the particular scoring model rather than the underlying property.
  **This is the ready-made, general hack-detection protocol** and it applies directly to a
  permeability oracle: hold out a second CycPeptMPDB model (different architecture or split) and
  watch for divergence. Follow-up: Langevin et al., "Explaining and avoiding failure modes in
  goal-directed generation of small molecules", *J Cheminform* 14:20, 2022,
  https://doi.org/10.1186/s13321-022-00601-y.
- **Gao & Coley, "The Synthesizability of Molecules Proposed by Generative Models", *JCIM*
  60(12):5714, 2020**, https://doi.org/10.1021/acs.jcim.0c00174; arXiv 2002.07007 — the
  SA-score-as-guard paper. Compares post-hoc filtering, training-set biasing, and **multiplying the
  normalised objective by a normalised SA Score or SCScore**. Figure 3a shows SMILES GA on
  *Osimertinib MPO* emitting boron/thiol/polyketone garbage — a worked picture of an optimiser gaming
  a GuacaMol MPO.
- **Conformal-prediction PepINVENT** (van Weesep, …, Geylan, AstraZeneca, arXiv 2605.05770) — the
  clearest statement of this problem *in the peptide domain*: RL "can suggest molecules that lie
  outside the predictor's domain of applicability … potentially steering designs into high reward but
  also high uncertainty chemical spaces." Their fix pays reward only inside a conformal applicability
  domain. See §6.6.

---

## SECTION 6 — Published multi-parameter peptide design objectives

### 6.1 PepINVENT's actual shipped MPO configuration

From `PepInvent/data/experiment_configurations/config_crbp_peptide.json` in
https://github.com/MolecularAI/PepInvent (read directly, this session). Target: cyclic Rev-binding
peptide. **Aggregation: `geometric_mean`.** 1000 RL steps, batch 32.

| component | weight | transform |
|---|---|---|
| `camsol_solubility` | 1 | `sigmoid(low=0, high=0.5, k=0.5)` |
| `predictive_model` (PAMPA permeability) | **5** | `no_transformation` (raw P(permeable) ∈ [0,1]) |
| `maximum_ring_size` | 1 | `double_sigmoid(low=12, high=50, coef_div=10, coef_si=20, coef_se=20)` |
| `lipophilicity` | **0** | `step(low=1, high=4)` |
| `custom_alerts` | 1 | `no_transformation`, 12 SMARTS patterns |

Diversity filter: `IdenticalMurckoScaffold`, score_threshold 0.4, bucket_size 25,
similarity_threshold 0.4, **penalty 0.5**.

Five things are directly relevant to hack-resistance, and this is the most useful single artefact
in the dossier:

1. **Lipophilicity is weight 0 and `step(1, 4)`.** It is *tracked, not optimised*, and when used at
   all it is a **window** [1, 4], not a maximand. AstraZeneca did not treat logP as a reward. The
   paper confirms: "Lipophilicity was also reported as the Wildman–Crippen log P value during the
   runs **to visualize the shifts**".
2. **Geometric mean, not weighted sum.** `geometric_mean.py` computes `total = Π sᵢ^wᵢ` then
   `^(1/Σwᵢ)`. A geometric mean is zero if any component is zero, so no component can be traded away.
3. **Gates sit outside the mean.** `substructure_match` is **excluded** from the mean and multiplied
   in afterwards: `final = gmean · s_substructure`. `custom_alerts` returns `1 − any(SMARTS match)`,
   a 0/1 killer term — **a single alert zeroes the whole geometric mean.**
4. **Ring size is a `double_sigmoid` window [12, 50]**, not "maximise ring size". The topology
   configs show the full idiom: `sigmoid(12,50,k=0.4)` to maximise ring size, `double_sigmoid(12,30)`
   to pin head-to-tail macrocycles, `reverse_sigmoid(12,50,k=0.4)` to force linear peptides.
5. **12 custom-alert SMARTS** penalise implausible substructures: `[#8][#8]` (peroxide), `[#6;+]`
   (carbocation), `C#C`, `C(=[O,S])[O,S]`, N–N, N–S, plus six aminal/acetal/thioacetal patterns.

**Correction worth flagging:** `camsol_solubility` appears in the shipped config but is **not** in
the scoring-component factory registry (`scoring_component_factory.py`, which lists
`maximum_ring_size`, `molecular_weight`, `substructure_match`, `predictive_model`, `lipophilicity`,
`custom_alerts`). Combined with CamSol-PTM being licence-gated (§3.1), **the published PepINVENT MPO
cannot be run end-to-end from the public repo.** The permeability half can.

Permeability oracle: XGBoost v1.7.5 on CycPeptMPDB PAMPA, 90/10 split stratified on labels **and
data sources**, balanced accuracy 0.78, MCC 0.59; score = `predict_proba()[:,1]`, untransformed.

Transform formulas as implemented (`transformations/score_transformations.py`):
- `sigmoid`: `1/(1 + 10^(10k(x−(low+high)/2)/(low−high)))`
- `reverse_sigmoid`: `1/(1 + 10^(k(x−(high+low)/2)·10/(high−low)))`
- `double_sigmoid`: `A/B − C` with `A=10^(cse·x/cdiv)`, `B=A+10^(cse·low/cdiv)`,
  `C=10^(csi·x/cdiv)/(10^(csi·x/cdiv)+10^(csi·high/cdiv))`
- `step`/`left_step`/`right_step`: hard 0/1; plus `custom_interpolation` (scipy `interp1d`).

Reported outcome: solubility learned in the first ~100 steps, permeability after; "permeability was
harder to optimize"; a peak around step 200 then a drop attributed to the diversity-filter penalty;
backbone heterocycle incorporation emerged as a preferred design.

### 6.2 REINVENT 4 scoring machinery (PepINVENT's parent)

Loeffler et al., *J Cheminform* 16:20, 2024; https://github.com/MolecularAI/REINVENT4
(`configs/SCORING.md`, `reinvent/scoring/`).

**Transforms**: `sigmoid`, `reverse_sigmoid`, `double_sigmoid`, `right_step`, `left_step`, `step`,
`value_mapping`, `exponential_decay`. v4 rewrote the sigmoids as a numerically stable base-10
logistic (`x = v − (high+low)/2; k = 10k/(high−low)`); `high==low` degenerates to `hard_sigmoid`.
`double_sigmoid` splits at `x_center` with `k_left=coef_si/coef_div`, `k_right=coef_se/coef_div`
(defaults coef_div 100, coef_si/se 150).

**Aggregation** (`scoring/aggregators/means.py`): weights are **normalised** (`w/Σw`);
`arithmetic_mean = nansum(s·w/Σw)`; `geometric_mean = nanprod(s^(w/Σw))` with `s` clamped to ≥ 1e-8.
NaNs are masked out of the weight normalisation per molecule.

**Three component classes** (`scoring/config.py`): *scorers* enter the mean; *filters* produce a
validity mask applied **before** scoring; *penalties* are multiplied onto the aggregate
(`total_scores * penalties`). `custom_alerts` ships as a filter with **no weight** — "as a filter it
will be applied globally". Example `configs/scoring.toml`: geometric_mean over QED (w 0.25),
MolecularWeight (w 0.25, `double_sigmoid(low 200, high 500, coef_div 500, si/se 20)`),
TanimotoDistance (w 0.1).

**No permeability or peptide oracle ships with REINVENT4** — only SAScore, chemprop, qptuna,
dockstream, synthsense and RDKit physchem.

### 6.3 QED and CNS MPO — the template desirability formulations

**QED** — Bickerton, Paolini, Besnard, Muresan & Hopkins, *Nat Chem* 4:90–98, 2012,
https://doi.org/10.1038/nchem.1243. Reference implementation: `rdkit/Chem/QED.py`. Eight properties:
MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, **ALERTS** (count of 116 structural-alert SMARTS). Each mapped
by an **asymmetric double sigmoid**:

```
d(x) = [A + B/(1+exp(−(x−C+D/2)/E)) · (1 − 1/(1+exp(−(x−C−D/2)/F)))] / DMAX
```

(MW: A 2.817, B 392.575, C 290.749, D 2.420, E 49.223, F 65.372, DMAX 104.981.)
Aggregation = **weighted geometric mean**: `QED = exp(Σ wᵢ ln dᵢ / Σ wᵢ)`. Default `WEIGHT_MEAN`:
MW 0.66, ALOGP 0.46, HBA 0.05, HBD 0.61, PSA 0.06, ROTB 0.65, AROM 0.48, **ALERTS 0.95** — the
anti-garbage term carries the **highest weight of all eight**.

**But QED must not be used in bRo5 space**: Doak et al. 2014 measured mean QED 0.16 for oral bRo5
drugs vs 0.49 for molecules medicinal chemists rejected as unattractive (§2.1). QED is a
small-molecule prior that actively penalises the molecules this benchmark is about.

**CNS MPO** — Wager et al., *ACS Chem Neurosci* 1:435–449, 2010, https://doi.org/10.1021/cn100008c;
reassessment https://doi.org/10.1021/acschemneuro.6b00029 (2016). Six properties: ClogP, ClogD7.4,
MW, TPSA, HBD, most-basic pKa. **Piecewise-linear** desirability T0 ∈ [0,1]; monotonic decreasing for
ClogP, ClogD, MW, pKa, HBD; a **hump** for TPSA. **Unweighted sum → 0–6**, threshold ≥ 4 (74% of 119
marketed CNS drugs). Breakpoints: ClogP 1→0 over [3,5]; ClogD over [2,4]; MW over [360,500]; HBD over
[0.5,3.5]; pKa over [8,10]; TPSA = 1 on [40,90], 0 at ≤ 20 and ≥ 120.

*Caveat:* GuacaMol's `CNS_MPO_ScoringFunction` is **not** Wager's score — it uses 5 Gaussian terms
(TPSA twice as min- and max-Gaussian, plus HBD, logP, MW; no ClogD, no pKa) and returns
`0.2·(o1+…+o5)`.

### 6.4 bRo5 / cyclic-peptide-specific schemes

- **AB-MPS = |cLogD − 3| + NAR + NRB**, lower better, **≤ 14** flags acceptable oral PK (≤ 15 for
  oral vs parenteral) — DeGoey, Chen, Cox & Wendt, *J Med Chem* 61:2636–2651, 2018,
  https://doi.org/10.1021/acs.jmedchem.7b00717. An **unweighted sum of raw penalties, no
  normalisation**. Best-validated bRo5 MPO and intrinsically anti-inflation (§2.2, §5.2b).
- **Doak/Kihlberg bRo5 windows** — a conjunction of thresholds, not a scalar: MW ≤ 1000,
  −2 ≤ cLogP ≤ 10, HBD ≤ 6, HBA ≤ 15, TPSA ≤ 250 Å², NRotB ≤ 20 (*Chem Biol* 2014).
- **ETR = EPSA/TPSA**, *J Med Chem* 2024 — ratio-based polarity-reduction descriptor validated on
  ~1,000 human-absorption compounds plus ~10,000 AbbVie compounds.
- **Whitty's two-constraint frontier**: TPSA ≥ 0.2·MW (solubility) AND 3D PSA ≤ 140 Å² (permeability),
  *Drug Discov Today* 2016 — a ready-made conjunctive objective with a published infeasibility
  region (MW > 600).
- **LPE** (lipophilic permeability efficiency): `LPE = logD_dec/w − 1.06·ALogP + 5.47` (Naylor et al.
  2018).
- **HBD bi-property models** for oral vs parenteral macrocycles — García Jiménez, Poongavanam &
  Kihlberg, *J Med Chem* 2023: 83–94% sensitivity for orals, 67–79% specificity for parenterals,
  beating AB-MPS's 79%/61% and 71%/79%.
- **No published Wager-style desirability MPO exists specific to macrocycles.** The bRo5 field uses
  raw sums (AB-MPS) or hard windows. That is itself a gap and an opportunity.

### 6.5 Generative cyclic-peptide objectives, as actually coded

- **PepEVOLVE** (arXiv 2511.16912, 2025, https://arxiv.org/abs/2511.16912), Eq. 15 — explicit
  **weighted geometric mean**: `Score = (S_perm³ · S_ring · S_lipo · S_SMARTS)^(1/6)`. Permeability
  exponent 3, others 1; lipophilicity **targeted at ≈ −4.0** to match PepINVENT. GRPO-style
  group-relative advantage `(R − R̄)/σ_R`.
- **HELM-GPT** (https://github.com/charlesxu90/helm-gpt, 2024): `score_type ∈ {sum, product, weight}`,
  default `weight = (scores·w/Σw).sum()`. Permeability: `sigmoid(−8, −4, k=1)` on predicted logPapp —
  **a window on the log scale, saturating at −4**. KRAS Kd/IC50: `rsigmoid(0, 2, k=1)` on log10 nM.
  Dual objective = unweighted sum of two sigmoids (range 0–2). REINVENT loss, σ=60.
- **RFpeptides** (*Nat Chem Biol* 2025, https://www.nature.com/articles/s41589-025-01929-w): **no
  scalar reward at all** — a hard filter cascade. AfCycDesign normalised iPAE < 0.3 (MCL1/MDM2),
  < 0.13 (GABARAP), < 0.28 (RbtA); Cα RMSD < 1.5 Å; Rosetta ddG < −50/−30/−40 kcal/mol; contact
  molecular surface > 300 Å²; **SAP < 35** (an aggregation guard).
- **AfCycDesign** (*Nat Commun* 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC12095755/): weighted-sum
  hallucination loss `(1 − pLDDT) + PAE/31 + con/2`; selection pLDDT > 0.9, Rosetta Pnear > 0.9.
- **PepThink-R1** (arXiv 2508.14765, 2025): `R = dup_fac · (0.8·prop_smooth + 0.2·sim_fac)`, with
  `prop_smooth = (1/3)Σσ((xᵢ−tᵢ)/kᵢ)` over LogD / MRT_rat / SIF t½, `sim_fac = σ(α(s − 0.6))`
  Tanimoto, and `dup_fac = (1/max(1,n+1))^γ` — an explicit **multiplicative anti-duplication
  penalty**.
- **MODAN** (https://github.com/ycu-iil/MODAN, 2023) — multi-objective Bayesian optimisation for AMPs
  with α,α-disubstituted non-proteinogenic amino acids; objectives are activity **plus hemolysis**,
  over Morgan/MACCS fingerprints. Notable for using the toxicity counter-objective (§4.5).
- **CycPeptMP** (https://github.com/akiyamalab/cycpeptmp) — single-objective logPapp regressor used
  downstream as an oracle; no MPO code in the repo.

### 6.6 GuacaMol / PMO oracle definitions, and who discusses hack-resistance

**GuacaMol modifiers** (`score_modifier.py`): `Linear`, `Squared(1−c(t−x)²)`, `Absolute(1−|t−x|)`,
`Gaussian`, `MinGaussian`/`MaxGaussian` (half-Gaussians, 1.0 on the desired side), `ClippedScore`,
`SmoothClippedScore` (logistic, k=4/(upper−lower)), `ThresholdedLinear`, `Chained`. Exact MPOs
(`standard_benchmarks.py`), all geometric mean unless noted:

- **Osimertinib MPO** = gmean[FCFP4 Tanimoto·Clipped(upper 0.8); ECFP6 Tanimoto·MinGaussian(0.85,0.1)
  ("but_not_too_similar"); TPSA·MaxGaussian(100,10); logP·MinGaussian(1,1)]
- **Fexofenadine MPO** = gmean[AP Tanimoto·Clipped(0.8); TPSA·MaxGaussian(90,10); logP·MinGaussian(4,1)]
- **Perindopril MPO** = gmean[ECFP4 Tanimoto; #aromatic rings·Gaussian(2, 0.5)]
- **Amlodipine MPO** = gmean[ECFP4 Tanimoto; #rings·Gaussian(3, 0.5)]
- **Sitagliptin MPO** = gmean[ECFP4·**Gaussian(μ=0, σ=0.1)** (rewards *dissimilarity*);
  logP·Gaussian(target,0.2); TPSA·Gaussian(target,5); Isomer C16H15F6N5O]
- **Zaleplon MPO** = gmean[ECFP4 Tanimoto; Isomer C19H17N3O2]
- **Pioglitazone MPO** = gmean[ECFP4·Gaussian(0,0.1); MW·Gaussian(356.4,10); NRotB·Gaussian(2,0.5)]
- **Ranolazine MPO** = *arithmetic* mean[AP Tanimoto·Clipped(0.7); logP·MaxGaussian(7,1);
  F-count·Gaussian(1,1); #aromatic rings·MinGaussian(1,1)]
- **Physchem MPO** (arithmetic) = Bertz·MaxGaussian(1500,200), MW·MinGaussian(400,40),
  aromatic rings·MinGaussian(3,1), F·Gaussian(6,1)
- **Median molecules 1/2** = gmean of two Tanimotos (camphor/menthol ECFP4; tadalafil/sildenafil ECFP6)
- **Valsartan SMARTS** = gmean[SMARTS match; logP, TPSA, Bertz each Gaussian at sitagliptin's values]
- **Isomer** scorer = gmean of per-element AtomCounter·Gaussian(n,1) plus total-atom·Gaussian(N,2)

Benchmark score = weighted average over top-1/top-10/top-100.

Note the pattern across all of them: **every physicochemical term is Gaussian, MinGaussian,
MaxGaussian or Clipped — i.e. bounded above.** Not one is a raw maximised logP. GuacaMol's authors
designed out the exact failure the pilot run hit.

**PMO** reimplements 19 GuacaMol oracles in PyTDC (`tdc/chem_utils/oracle/oracle.py`) with identical
modifiers and gmean aggregation, plus QED, DRD2, GSK3β, JNK3. Metric: **AUC Top-10**, trapezoidal
area under (mean of top-10) vs oracle calls, normalised by `max_oracle_calls = 10,000`, logged every
100 calls. **SA score and internal diversity are logged alongside every run** — the degeneracy
diagnostics are built into the harness.

**Who explicitly discusses hack-resistance:**

| source | what it says |
|---|---|
| GuacaMol 2019 | trivial objectives; gmean rationale; Tanimoto chosen because ML oracles are exploitable; rd_filters pass rates (§5.3) |
| PMO 2022 | decomposability of isomer oracles; "trivial only when queries are uncontrolled" — the 10k budget *is* the device; atom-by-atom builders waste budget on unstable molecules |
| Renz et al. 2019 | AddCarbon; **control-model divergence test** (§5.3) |
| Gao & Coley 2020 | SA/SCScore as a multiplicative guard; pictures of GA-generated garbage |
| **Conformal PepINVENT**, arXiv 2605.05770 | RL "steering designs into high reward but also high uncertainty chemical spaces"; replaces raw XGBoost probability with Mondrian inductive conformal P-values, non-conformity `α = 0.5 − (P̂(yᵢ|x) − max_{y≠yᵢ}P̂(y|x))/2`; five reward variants including a **"harsh"** binary (1 iff P1>0.2 and P0<0.2) and a **"soft"** partial-credit version; the soft one converged in ~150 steps with the largest fraction of confidently-permeable designs. **The only applicability-domain reward native to peptides.** |
| ProDCARL 2026 | names "reward hacking"; entropy regularisation + early stopping |
| Krishnan / Attack-as-Design 2026 | constrains homopolymer run length, charge, hydrophobicity to natural ranges |
| BATTLE-AMP 2026 | composition-preserving shuffled decoys as hard negatives |
| AmpHGT 2025 | restricted NCAA vocabulary to the train/test intersection to stop identity-cheating |
| PepINVENT / PepEVOLVE / PepThink-R1 | not in those words, but the configs encode it: gmean, logP at weight 0 as a window, ring-size window, alert gates outside the mean, diversity/duplication penalties |

### 6.7 Design patterns that recur across every hack-resistant objective found

1. **Weighted geometric mean with normalised weights** is the dominant anti-degeneracy aggregator —
   QED, GuacaMol, REINVENT4, PepINVENT, PepEVOLVE. Arithmetic mean only where partial failure should
   be survivable (GuacaMol's Ranolazine, Physchem, Cobimetinib).
2. **Gates multiply, they do not average.** PepINVENT's `substructure_match`, REINVENT4's penalty
   class, QED's ALERTS at the highest weight of all eight terms.
3. **Windowed transforms beat monotone ones** wherever the optimiser could inflate without bound:
   PepINVENT's `double_sigmoid(12,30)` on ring size, CNS MPO's TPSA hump, AB-MPS's `|cLogD − 3|`,
   every GuacaMol physchem term.
4. **Two published hack-*detection* devices**: an oracle-call budget with AUC-over-budget scoring
   (PMO), and an independent control model watched for divergence (Renz et al.). The newest and only
   peptide-native one is the **conformal applicability-domain reward**.

---

## GAPS

- **Model-side, not label-side, generalisation is unmeasured.** §1.1–§1.2 are label correlations and
  §1.3 probes one oracle (PepINVENT's XGBoost). CycPeptMP, CPMP, MultiCycPermea and PeptideCLM were
  **not** probed adversarially here; CycPeptMP needs commercial MOE for feature generation, which
  blocks independent replication of the best-performing published model.
- **No solubility oracle exists that is simultaneously ncAA-aware, cyclisation-aware and
  weight-downloadable.** CamSol-PTM is the only ncAA-native tool and it is web-only, licence-gated,
  blind to chirality/N-methylation/cyclisation, and validated on 18–36mers rather than MW 500–2000
  macrocycles. PepINVENT's published MPO is therefore not externally reproducible end-to-end.
- **Permeability oracles are saturated, and nobody has published the saturation statistic.** 51.7% of
  CycPeptMPDB exceeds P = 0.9 under PepINVENT's own classifier. The literature discusses accuracy,
  never headroom. A regression-on-logPapp oracle under a scaffold split is the obvious fix but its
  difficulty has not been characterised as a generative benchmark.
- **No published cyclic-peptide benchmark reports what its optimiser converged to.** Chemical-property
  distributions of generated molecules (logP, MW, HBD) are reported as plots in PepINVENT and nowhere
  as the degeneracy diagnostic the pilot run needed; there is no accepted hack-resistance metric for
  peptide design analogous to BATTLE-AMP's SyntheticShuffled or Renz's control-model divergence test.
  Relatedly, **no Wager-style desirability MPO exists for macrocycles** — the bRo5 field still uses
  raw sums (AB-MPS) and hard windows, so there is no calibrated multi-property score to benchmark
  against.
- **The AMP-side evidence leans on three 2026 preprints** (BATTLE-AMP, GenPept-2025, ProDCARL,
  Genotypic Triggers) that are not peer-reviewed. The peer-reviewed core (Sidorczuk 2022,
  Mooney 2012, Porto 2020, Bello-Madruga 2024) supports the same conclusion but with weaker numbers.
- **The conformation-dependent quantities that make permeability hard (3D PSA, EPSA, ΔG_desolv,
  chameleonicity) have no fast differentiable oracle.** EPSA and Chamelogk are chromatographic
  measurements; MD-based AP scores cost minutes to hours per call. Using them as benchmark objectives
  requires training a surrogate first, which reintroduces the shortcut risk the surrogate was meant
  to avoid.
