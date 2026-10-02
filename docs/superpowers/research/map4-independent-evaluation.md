# MAP4 / MAP4C under independent scrutiny

Adversarial evidence review. Compiled 2 October 2026. Facts and measured numbers only.
Every claim carries a URL and a year. Evidence for and against is reported in both directions.

Primary sources under evaluation (authors' own work, listed for reference, not evaluated here):

- Capecchi, Probst, Reymond, *J. Cheminform.* **12**:43 (2020) — MAP4.
  https://doi.org/10.1186/s13321-020-00445-4 · full text https://pmc.ncbi.nlm.nih.gov/articles/PMC7291580/
- Orsi, Reymond, *J. Cheminform.* **16**:53 (2024) — MAP4C.
  https://doi.org/10.1186/s13321-024-00849-6 · full text https://pmc.ncbi.nlm.nih.gov/articles/PMC11090803/

---

## SECTION 1 — INDEPENDENT EVALUATIONS

Four third-party evaluations that actually include MAP4 were found. Two are peer-reviewed papers,
one is a blog post by a recognised cheminformatics source, one is an application paper. No
independent evaluation of **MAP4C** was found at all.

### 1.1 Venkatraman, Gaiser, Demekas, Roy, Xiong, Wheeler (2024) — largest independent VS benchmark that includes MAP4

*"Do Molecular Fingerprints Identify Diverse Active Drugs in Large-Scale Virtual Screening? (No)"*,
*Pharmaceuticals* **17**(8):992 (2024). NTNU Trondheim / University of Arizona / NIH-NIAID /
University of Montana — no Reymond-group authors.
https://doi.org/10.3390/ph17080992 · full text https://europepmc.org/articles/PMC11356940

Task: ligand-based virtual screening, 32 fingerprints, four benchmarks (DEKOIS, DUD-E, MUV,
LIT-PCBA), plus the MMV St. Jude malaria set (2,507 confirmed actives, 303,303 **experimentally
confirmed** inactives). MAP4 used at 2,048 bits. Metrics: ROC AUC and a new "decoy retention
factor" DRF at p = 0.1 (**lower DRF is better**).

Table 1 of that paper (selected rows, exact values):

| Fingerprint | AUC DEKOIS | AUC DUD-E | AUC MUV | AUC LIT-PCBA | DRF DEKOIS | DRF DUD-E | DRF MUV | DRF LIT-PCBA |
|---|---|---|---|---|---|---|---|---|
| **MAP4 (2048)** | **0.81** | **0.83** | 0.56 | 0.54 | **0.14** | **0.07** | 0.91 | 0.85 |
| MHFP (2048) | 0.81 | 0.81 | 0.54 | 0.53 | 0.17 | 0.10 | 0.97 | 0.94 |
| ECFP4 (1024) | 0.76 | 0.80 | 0.54 | 0.51 | 0.19 | 0.09 | 0.99 | 1.00 |
| ECFP6 (1024) | 0.75 | 0.78 | 0.54 | 0.52 | 0.20 | 0.10 | 0.99 | 0.98 |
| TT | 0.80 | 0.80 | **0.61** | 0.56 | 0.15 | 0.10 | **0.72** | 0.85 |
| AT2D (4096) | 0.78 | 0.79 | 0.58 | 0.55 | 0.20 | 0.12 | 0.86 | 0.83 |
| AP2D (4096) | 0.64 | 0.66 | 0.49 | 0.51 | 0.55 | 0.39 | 1.24 | 1.00 |
| AVALON (1024) | 0.72 | 0.73 | 0.60 | 0.55 | 0.30 | 0.18 | 0.85 | 0.97 |
| RDK7 (1024) | 0.62 | 0.63 | 0.58 | **0.59** | 0.74 | 0.70 | 0.98 | 0.85 |
| MACCS (166) | 0.71 | 0.75 | 0.55 | 0.54 | 0.33 | 0.18 | 0.99 | 0.93 |

**For MAP4:** it is the single best fingerprint of the 32 on DUD-E (AUC 0.83, DRF 0.07) and tied
best on DEKOIS (AUC 0.81, best DRF 0.14). This is the strongest piece of independent evidence in
MAP4's favour and it is from a paper actively hostile to fingerprints in general.

**Against MAP4:** on the two harder, less biased benchmarks it is mid-pack — MUV AUC 0.56 (TT 0.61,
AVALON 0.60, AT2D/RDK5/RDK7 0.58) and LIT-PCBA AUC 0.54 (RDK7 0.59, RDK6 0.58, RDK5 0.56). The
authors themselves state that DEKOIS/DUD-E performance "is offset by concerns of benchmark bias …
such as artificial enrichment … analogue bias … and false negative bias", so MAP4's lead sits on
exactly the benchmarks the paper distrusts.

**Specifically about MAP4 thresholds** (St. Jude malaria set, Section 2.4): "MAP4 shows an apparent
relative abundance of actives, but note enrichment is still only ∼10-fold; also note that <1.5% of
actives show Tanimoto similarity >0.1 to another active, raising concerns about how to establish
meaningful MAP4 thresholds."

The paper's global conclusion applies to MAP4 as to all 32: the best DRF(0.1) anywhere in the study
was 0.09 (~11-fold enrichment), and fingerprint similarity "do[es] not correlate with compound
potency" (their Section 2.5, Kendall τ heatmap, Figure 5).

### 1.2 Boldini, Ballabio, Consonni, Todeschini, Grisoni, Sieber (2024) — natural-products benchmark

*"Effectiveness of molecular fingerprints for exploring the chemical space of natural products"*,
*J. Cheminform.* **16**:35 (2024). TU Munich / Milano-Bicocca / TU Eindhoven — no Reymond authors.
https://doi.org/10.1186/s13321-024-00830-3 · full text https://pmc.ncbi.nlm.nih.gov/articles/PMC10964529/
Code: https://github.com/dahvida/NP_Fingerprints

20 fingerprints (MAP4 at 1024, radius 2, with the modified Jaccard for categorical MinHash vectors),
129,869 curated COCONUT natural products for similarity analysis plus 12 CMNPD bioactivity
classification datasets (Random Forest and DNN, scaffold split, 5 replicates).

*Similarity distribution* (their Table 4, COCONUT): MAP4 has the **narrowest and lowest** distribution
of all 20 — min 0.000, 25th pct 0.002, median **0.011**, 75th pct 0.026, max **0.067**. For comparison
MHFP median 0.052, ECFP 0.108, AP 0.184, MACCS 0.410. The authors attribute this partly to the
modified Jaccard (two positions match only if they hold the same integer) and partly to "MinHashing
paths rather than circular fragments". Practical consequence: across 50 × ~50 million natural-product
pairs, **no MAP4 pair exceeded 0.067**.

*Correlation:* MAP4 and MHFP correlate at Pearson R = 0.85 on COCONUT. MAP4's correlation with ECFP
is 0.67 and with FCFP 0.77 (MHFP: 0.77 and 0.88).

*QSAR:* the best performers named are RAD2D (best MCC, 0.506), LSTAR (best ROC-AUC, 0.900) and MHFP
(best PR-AUC, 0.669). MAP4 is **not** among the winners on any of the three headline metrics. The
only positive note for MAP4 is qualitative: "MAP4, MHFP and LSTAR tend to have less false positives".
Overall finding: "no encoding significantly outperforms all others across all QSAR datasets".

### 1.3 Walters (2020) — Practical Cheminformatics blog, regression benchmark

Pat Walters, *"Benchmarking 'One Molecular Fingerprint to Rule Them All'"* (21 March 2020).
http://practicalcheminformatics.blogspot.com/2020/03/benchmarking-one-molecular-fingerprint.html

24 ChEMBL IC50 datasets (203–5,207 molecules each, from a 2019 Bender-group paper), XGBoost
regression, 10 folds, R² on test, MAP4 vs Morgan. Result: Cohen's *d* was **negative in all 24
comparisons**, and **< −0.8 (large effect) in 20 of 24**. Conclusion as written: "I'll probably stick
with the Morgan fingerprints for my regression models." Caveat the post itself states: it used
Iwatobipen's patched MAP4 because of a bug in the original folded-fingerprint code.

This is the only head-to-head *regression* comparison found and it is clearly negative for MAP4.

### 1.4 López-López, Robles, Plisson, Medina-Franco (2023) — application, not a comparison

*"Mapping the structure–activity landscape of non-canonical peptides with MAP4 fingerprinting"*,
*Digital Discovery* **2**, 1494–1505 (2023). UNAM / CINVESTAV / Universidad Veracruzana.
https://doi.org/10.1039/D3DD00098B

Uses MAP4 + the MinHashed distance to build SAS maps and SALI values for 223 anti-MRSA peptides
(25,185 pairs). This is an independent *use* of MAP4, and it reports a positive correlation between
MAP4-based similarity and sequence-based similarity. It is **not** a benchmark: no other fingerprint
is run as a control, so it supplies no comparative evidence either way.

### 1.5 Benchmarks that OMIT MAP4

These are the major fingerprint/representation benchmarks published **after** MAP4 (2020) that could
have included it and did not. Checked by full-text search for "MAP4"/"MinHash":

- **MoleculeACE / activity cliffs** — van Tilborg, Alenicheva, Grisoni, *JCIM* **62**(23):5938–5951
  (2022). 24 ML approaches, 30 targets. Descriptors evaluated: ECFP, MACCS, WHIM, 11 physicochemical
  properties. **MAP4 absent.** https://doi.org/10.1021/acs.jcim.2c01073 ·
  https://github.com/molML/MoleculeACE
- **Molecular Fingerprints Are Strong Models for Peptide Function Prediction** — Adamczyk, Ludynia,
  Czech, arXiv:2501.17901 (2025, v3 checked). The largest peptide-fingerprint study found: 6
  benchmarks, 126 datasets, SOTA on LRGB Peptides-func/struct. Fingerprints evaluated: count ECFP,
  count Topological Torsion, count RDKit. **MAP4 absent from v1, v2 and v3** — notable because these
  authors maintain scikit-fingerprints, which ships a MAP implementation.
  https://arxiv.org/abs/2501.17901
- **Systematic benchmarking of 13 AI methods for cyclic peptide membrane permeability** — *J.
  Cheminform.* **17**:129 (2025), ~6,000 CycPeptMPDB peptides. Fingerprint arm = 210 RDKit
  descriptors. **MAP4 absent** (zero occurrences of "MAP4" or "MHFP" in the full text).
  https://doi.org/10.1186/s13321-025-01083-4
- **Riniker & Landrum** benchmarking platform, *J. Cheminform.* **5**:26 (2013) and **O'Boyle &
  Sayle**'s literature-based similarity benchmark, *J. Cheminform.* **8**:36 (2016) both predate MAP4,
  so their omission is chronological, not a judgement. MAP4's own evaluation is an extension of the
  Riniker/Landrum platform. https://doi.org/10.1186/1758-2946-5-26 ·
  https://doi.org/10.1186/s13321-016-0148-0

**Reimplementation (counts as adoption, not as evaluation):** scikit-fingerprints provides
`MAPFingerprint` and lists MAP4 among "new ones" it added (Adamczyk & Ludynia, arXiv:2407.13291,
2024; library docs https://scikit-fingerprints.readthedocs.io/latest/modules/fingerprints.html).
Shipping an implementation is not an endorsement of performance, and the same group's peptide paper
did not use it.

**Summary of Section 1.** Independent evidence on MAP4 amounts to two peer-reviewed benchmarks and
one blog post. It is split: strongly positive on DEKOIS/DUD-E retrieval (Venkatraman 2024),
mid-pack-to-negative on MUV/LIT-PCBA retrieval (same paper), not a winner on natural-product QSAR
(Boldini 2024), and clearly worse than Morgan for regression (Walters 2020). MAP4C has **zero**
independent evaluations.

---

## SECTION 2 — IS "RESOLUTION" A VALID ARGUMENT?

### 2.1 What the flagship statistic actually is

MAP4 paper, Table 4 (https://pmc.ncbi.nlm.nih.gov/articles/PMC7291580/): over 96,456 stereochemistry-
stripped HMDB metabolites, the number with an indistinguishable nearest neighbour (Jaccard distance
exactly 0):

| Fingerprint | zero-distance NNs | % |
|---|---|---|
| MAP4-1024 | 0 | 0.0% |
| **AP (RDKit, unhashed)** | **1,677** | **1.7%** |
| TT (unhashed) | 68,623 | 71.1% |
| MHFP6-1024 | ≈69,900 | 72.5% (percentage quoted in the paper's text) |
| ECFP4-1024 | 70,329 | 72.9% |

Two things follow directly from the authors' own table:

1. The dramatic "over 70%" contrast is against **substructure** fingerprints. Against the obvious
   non-MinHashed alternative — the plain RDKit atom-pair fingerprint — the gap is **1.7 percentage
   points**, not 70. AP achieves 98.3% resolution with no MinHashing at all.
2. Part of MAP4's zero-collision result is an artefact of the comparison rule, not of information
   content. For categorical MinHash vectors the modified Jaccard counts a match only when two
   positions hold the *same integer*, so exact ties are mechanically rarer than for binary vectors
   (Boldini et al. 2024, Methods, "Similarity metrics", https://pmc.ncbi.nlm.nih.gov/articles/PMC10964529/).

### 2.2 Resolution is not the accepted validation criterion; neighbourhood behaviour is

Patterson, Cramer, Ferguson, Clark, Weinberger, *"Neighborhood behavior: a useful concept for
validation of 'molecular diversity' descriptors"*, *J. Med. Chem.* **39**(16):3049–3059 (1996).
https://doi.org/10.1021/jm960290n · PDF https://gwern.net/doc/science/chemistry/1996-patterson.pdf

The validation test defined there is explicitly *not* about discrimination. It plots **differences in
descriptor value against differences in biological activity** over 20 literature datasets and asks
whether the points concentrate in a trapezoid ("neighbourhood enhancement"). Measured results: 2D
fingerprints of side chains and the two topomeric fields were significant in 43/51 = 85% of cases
with mean enhancements 1.59, 1.40, 1.47; whole-molecule 2D fingerprints, **atom pairs** and
autocorrelation vectors were significant in only 20/60 = 33% of cases with mean enhancements 1.21,
**1.15**, 1.13; connectivity indices, log P, MR and strain energy showed no significant neighbourhood
behaviour at all. A random number per molecule showed none in any of 20 datasets. The operating
threshold they adopt is an enhancement > 1.1.

The decisive point for this question: a random number per molecule has **perfect resolution**
(essentially zero collisions) and **zero** neighbourhood behaviour. Patterson's framework scores it
as worthless. Resolution therefore cannot be sufficient. *(This is reasoning from the cited result,
not a quotation.)* The same holds for the degenerate baseline of hashing a canonical SMILES string —
0% collisions, no usable similarity structure.

### 2.3 Direct measurements: more resolution buys very little

Greg Landrum (RDKit author), *"Colliding bits III, expanded"*, RDKit blog, 23 January 2023.
https://greglandrum.github.io/rdkit-blog/posts/2023-01-23-colliding-bits-iii-expanded.html
4 million PubChem compounds; Spearman R between folded-fingerprint Tanimoto and collision-free
Tanimoto:

| nBits | mfp2 | mfp3 | rdk5 | hashap | hashtt |
|---|---|---|---|---|---|
| 16384 | 0.999 | 0.999 | 1.000 | 0.997 | 0.999 |
| 1024 | 0.990 | 0.987 | 0.981 | 0.874 | 0.983 |
| 512 | 0.979 | 0.972 | 0.929 | 0.702 | 0.967 |
| 128 | 0.897 | 0.839 | 0.528 | 0.343 | 0.828 |

At 1,024 bits the mean similarity error from collisions is 0.014 (mfp2), with 90th-percentile |d| of
0.036. Stated conclusion: "Down to a fingerprint size of 1024 bits all the fingerprints except
hashap have an R of >0.90 and comparatively low deviations."

Greg Landrum, *"Colliding bits II, revisited"*, 25 December 2022.
https://greglandrum.github.io/rdkit-blog/posts/2022-12-25-colliding-bits-ii-revisited.html
1K vs 16K fingerprints across 37 ChEMBL target sets × 5 fingerprints × 5 learners. Median AUC change
from going 16× larger: between −0.087 and +0.031; for XGBoost the difference is "same" for mfp2 and
hashtt. rdk5 and hashap get *worse* with more bits under RF/BRF/XGB, which the post attributes to
overfitting. TL;DR as written: "a small, but real improvement".

*"Hash Collisions in Molecular Fingerprints: Effects on Property Prediction and Bayesian
Optimization"*, arXiv:2511.17078 (2025): exact (collision-free) fingerprints give "a small yet
consistent improvement in predictive accuracy" on five DOCKSTRING benchmarks.
https://arxiv.org/abs/2511.17078

Counter-evidence (resolution does sometimes help): O'Boyle & Sayle (2016) report that going from
1,024 to 16,384 bits improves VS on their literature-based similarity benchmark
(https://doi.org/10.1186/s13321-016-0148-0). But Venkatraman et al. (2024) re-tested exactly this on
MUV and LIT-PCBA and "found that longer fingerprints yield little to no gain in efficacy"
(their Section 2.1, Table S1).

### 2.4 Direct evidence that similarity *values* can be meaningless even at high resolution

Venkatraman et al. (2024) measured the thing resolution does not measure
(https://europepmc.org/articles/PMC11356940):

- active–active and active–decoy Tanimoto distributions "substantially overlap" for all 32
  fingerprints (their Figure 1, DEKOIS);
- where early enrichment does occur (GBA, OPRK1, PPARG in LIT-PCBA), the high-scoring actives are
  "almost entirely bioisosteres" sharing the query scaffold — for OPRK1, all 5 actives with Tc > 0.5
  are built on the query scaffold (their Figure 3);
- fingerprint similarity to a potent active does not correlate with the potency of other actives
  (their Section 2.5 and Figure 5).

Related methodological critiques of fingerprint evaluation itself: Sieg, Flachsenberg & Rarey,
*"In Need of Bias Control"*, *JCIM* (2019) https://doi.org/10.1021/acs.jcim.8b00712; Scior et al.,
*"Recognizing Pitfalls in Virtual Screening: A Critical Review"*, *JCIM* (2012)
https://doi.org/10.1021/ci200528d; *"Data Leakage and Redundancy in the LIT-PCBA Benchmark"*,
arXiv:2507.21404 (2025) https://arxiv.org/abs/2507.21404.

### 2.5 Verdict on Section 2

Resolution and meaningfulness are **different properties and are measured by different experiments**.
No source was found in which a low collision rate is proposed as evidence that a similarity measure
is good. The accepted validation criterion since Patterson et al. (1996) is the relationship between
descriptor distance and property distance. Empirically, large reductions in collisions (1K → 16K
bits, or exact fingerprints) change similarity values by ~0.01–0.04 and model AUC by a few percent at
most, in both directions. MAP4's zero-collision statistic is real and verifiable, but it is a
statement about resolution only, and against the fingerprint it most resembles (RDKit AP) the
improvement is 1.7 percentage points.

---

## SECTION 3 — MINHASH ESTIMATION ERROR

### 3.1 The canonical analysis

Broder, A. Z., *"On the resemblance and containment of documents"*, Proc. Compression and Complexity
of Sequences 1997, IEEE, pp. 21–29 (1998). doi:10.1109/SEQUEN.1997.666900. The unbiasedness theorem
and the sketch construction are in Broder, Glassman, Manasse, Zweig, *"Syntactic Clustering of the
Web"* (1997) http://bitsavers.trailing-edge.com/pdf/dec/tech_reports/SRC-TN-1997-015.pdf, and the
min-wise independence requirement in Broder, Charikar, Frieze, Mitzenmacher, STOC '98
https://eecs.harvard.edu/~michaelm/CS222/minwise.pdf.

The estimator is the mean of *k* i.i.d. Bernoulli(J) variables, so it is unbiased with expected error
O(1/√k) (https://en.wikipedia.org/wiki/MinHash, which notes "400 hashes would be required to estimate
J(A, B) with an expected error less than or equal to .05"). The same bound is restated for genomics
in Mash: "the error bound of the Jaccard estimate ε = O(1/√s) relies only on the sketch size and is
independent of genome size … the relative error can grow quite large for very small Jaccard values"
(Ondov et al., *Genome Biol.* **17**:132, 2016, https://doi.org/10.1186/s13059-016-0997-x).

### 3.2 Numbers at MAP4's permutation counts

SE(Ĵ) = √(J(1−J)/k). MAP4C (`mapchiral`) defaults to **k = 2048**; the Reymond `map4` reference
implementation defaults to **dimensions = 1024** (verified in source,
https://github.com/reymond-group/map4/blob/master/map4/map4.py, `MAP4Calculator.__init__`); the MAP4C
paper states "all fingerprints were used with 2,048-bits"; the PyPI `map4` README example uses
`dimensions=2048`.

| true J | SE at k=2048 | 95% CI half-width | relative 95% CI | SE at k=1024 |
|---|---|---|---|---|
| 0.90 | 0.0066 | ±0.013 | ±1.4% | 0.0094 |
| 0.70 | 0.0101 | ±0.020 | ±2.8% | 0.0143 |
| 0.50 | 0.0110 | ±0.022 | ±4.3% | 0.0156 |
| 0.30 | 0.0101 | ±0.020 | ±6.6% | 0.0143 |
| 0.10 | 0.0066 | ±0.013 | ±13% | 0.0094 |
| 0.05 | 0.0048 | ±0.009 | ±19% | 0.0068 |
| 0.011 | 0.0023 | ±0.0045 | **±41%** | 0.0033 |

The last row matters: 0.011 is the **median** MAP4 pairwise similarity across ~50 M natural-product
pairs (Boldini et al. 2024, Table 4). In the regime where MAP4 actually operates on large molecules,
the 95% confidence interval on a single similarity value is roughly ±40% of the value itself.

### 3.3 Measured confirmation (run for this review)

Script: `bench.py` in this scratchpad. Python 3.12, RDKit 2026.03.5, `mapchiral` 0.0.7, Intel Core
i7-14700HX, single-threaded. 12 molecule pairs drawn from 20 drug-like molecules + 10 peptides;
each pair's MinHash Jaccard re-estimated under **25 independently seeded permutation families**
(seeding the underlying `rdkit.Chem.rdMHFPFingerprint.MHFPEncoder(k, seed)` directly).

| Ĵ (mean over seeds) | observed SD, k=2048 | predicted √(J(1−J)/2048) | observed SD, k=1024 | predicted at k=1024 |
|---|---|---|---|---|
| 0.1484 | 0.0076 | 0.0079 | 0.0092 | 0.0111 |
| 0.1297 | 0.0088 | 0.0074 | 0.0119 | 0.0105 |
| 0.0188 | 0.0025 | 0.0030 | 0.0032 | 0.0043 |
| 0.0167 | 0.0023 | 0.0028 | 0.0035 | 0.0040 |
| 0.0131 | 0.0020 | 0.0025 | 0.0038 | 0.0036 |
| 0.0085 | 0.0021 | 0.0020 | 0.0024 | 0.0029 |
| 0.0040 | 0.0017 | 0.0014 | 0.0021 | 0.0020 |
| 0.0018 | 0.0008 | 0.0009 | 0.0014 | 0.0014 |

Observed spread tracks the binomial prediction within a factor of ~1.3 across three orders of
magnitude of J. The error is real, it is the size theory says, and it is invisible in any single run
because the permutations are fixed by a hard-coded seed.

### 3.4 Do the MAP4 papers discuss this?

**No.** Full-text searches of both papers (Europe PMC XML, PMC7291580 and PMC11090803) return **zero**
occurrences of "standard error", "variance", or "threshold". "Approximate" appears only in reference
to LSH-based approximate nearest-neighbour search, never to the Jaccard estimate itself. Both papers
cite Broder's MinHash lineage for the *construction* and refer the reader to the MHFP6 paper for
details, but neither states an error bar on a reported similarity. The MAP4 paper does discuss
collisions — but as an argument that MinHashing beats modulo folding, not as an error analysis.

### 3.5 What this implies for a hard similarity threshold

- A threshold applied to a single MAP4/MAP4C Jaccard value is applied to a random variable with
  SD ≈ 0.011 at J = 0.5 and k = 2048. Two different permutation seeds will put borderline pairs on
  opposite sides of any cut.
- Comparing two candidates: SE of the *difference* is √2 × SE ≈ 0.0156 at J ≈ 0.5, so a single
  pairwise ranking is not resolvable at 95% confidence unless the true gap exceeds ≈ 0.044. Ranking
  N candidates multiplies the number of such comparisons, so ties become pervasive.
- On large molecules the problem is worse in relative terms, not better: at the natural-product
  median J = 0.011 the relative 95% interval is ±41%.
- A fixed seed makes this *reproducible*, not *accurate*. Reproducibility across runs is not evidence
  that the value is close to the true Jaccard.
- Independent corroboration that no principled MAP4 threshold exists: Venkatraman et al. (2024),
  "<1.5% of actives show Tanimoto similarity >0.1 to another active, raising concerns about how to
  establish meaningful MAP4 thresholds."
- Mitigation that the literature supports: increase k (error falls only as 1/√k, so halving the error
  costs 4× the permutations and 4× the storage), or use a non-sketched fingerprint where the Jaccard
  is exact. The RDKit atom-pair and Morgan fingerprints compute exact Tanimoto on their (folded)
  feature sets — their error is bias from collisions, which is systematic and measurable
  (Section 2.3), not sampling noise that changes with a seed.

---

## SECTION 4 — PUBLISHED CRITICISM AND FAILURE CASES

### 4.1 Documented in the authors' own MAP4C paper: rotationally symmetric macrocycles

Orsi & Reymond (2024), Table 1 (https://pmc.ncbi.nlm.nih.gov/articles/PMC11090803/); all fingerprints
at 2,048 bits. Numbers are distinct fingerprint values obtained out of the total possible
stereoisomers:

| Query | N stereocentres / symmetry | Total stereoisomers | MAP6C | **MAP4C** | MAP2C | APC | ECFP6C | ECFP4C |
|---|---|---|---|---|---|---|---|---|
| Gramicidin S (11) | 10 / C2 | 528 | 528 | **504** | 334 | 25 | 448 | 243 |
| **Valinomycin (12)** | 12 / C3 | **1,376** | 1,250 | **714** | 416 | 112 | 616 | 27 |
| Nonactin (13) | 16 / C4 | 16,456 | 16,425 | **16,176** | 10,045 | 13,189 | 6,474 | 675 |
| NP213 (14) | 7 / C7 | 20 | 7 | **13** | 17 | 13 | 5 | 3 |
| Quinaldopeptin (9) | – / C2 | 136 | 136 | 136 | – | – | – | – |
| Onchidin (10) | 12 / C2 | 2,080 | 2,080 | 2,080 | 2,064 | 269 | 1,765 | 810 |

Stated in the paper: "none of the chiral fingerprints tested was able to cope with the C3 symmetrical
dodecadepsipeptide antibiotic valinomycin (12, 1,376 stereoisomers), the C4 symmetrical macrolide
ionophore antibiotic nonactin (13, 16,456 stereoisomers), or the C7 symmetrical hepta-arginine cyclic
peptide NP213 (14, 20 stereoisomers)" and "performance did not increase significantly when using much
larger bit sizes or without MinHashing or folding". MAP4C gets **714 / 1,376 = 51.9%** of valinomycin
stereoisomers and **13 / 20 = 65%** of NP213's — worse than MAP2C (17/20) on the latter. For NP213,
*every* MAP variant loses to the simpler MAP2C, and MAP6C is worst of all (7/20).

Note the paper also records a plain hash collision: "with 4,096 bits, only 135 different FP values are
obtained with 2,048 bits due to a bit collision".

Third-party summary of these limitations: LamaLab tool notes, "MAP4C: One chiral fingerprint to find
them all", 16 Oct 2024 — "MAPC performed well on sugars and peptides … performance decreased for
macrocyclic NPs". https://lamalab-org.github.io/toolminutes/papers/universal_chiral_fp_2024.html

### 4.2 Blog-level criticism

- **Pat Walters, Practical Cheminformatics (2020)** — the regression result in Section 1.3 above
  (Cohen's d < −0.8 in 20 of 24 ChEMBL datasets vs Morgan). Explicitly non-hostile in tone
  ("my intent here is not to be critical of the paper or the authors") but the numbers are negative.
  http://practicalcheminformatics.blogspot.com/2020/03/benchmarking-one-molecular-fingerprint.html
- **Iwatobipen (20 March 2020)** — found and fixed a bug preventing folded-fingerprint generation in
  the original repo (merged as GitHub PR; Iwatobipen is listed as a `reymond-group/map4` contributor).
  Also observed on his test set that "Morgan FP seems more sensitive to difference of compound
  structure. Because similarity is lower to MAP4FP."
  https://iwatobipen.wordpress.com/2020/03/20/new-molecular-fingerprint-for-chemoinformatics-map4-rdkit-memo-chemoinformatics/
- **Venkatraman et al. (2024)** — threshold criticism quoted in Sections 1.1 and 3.5.

### 4.3 Open GitHub issues (`reymond-group/map4`, 12 open of 27 total)

https://github.com/reymond-group/map4/issues

- **#21 "Totally different got similarity of 1"** (open since 2023-01-18). A user reports two
  structurally unrelated molecules receiving similarity 1.0. The evidence is a screenshot only, the
  maintainers have not responded, and the only comments are a third party (UnixJunkie) asking for
  more pairs. **Unverified**; recorded here as an open, unresolved bug report, not as a confirmed
  failure.
- **#12 "If a molecule has electric charge, some shingles have empty smiles"** (closed, confirmed by
  the author). Disconnected components (salts, counter-ions) yield shingles with an empty SMILES side
  and a sentinel topological distance of `100000000`, e.g. `b'|100000000|c(c)(c)[N+]'`. Author's
  explanation: single-atom fragments form no atom pairs. **Practical consequence: MAP4 must be fed
  desalted, single-component structures.**
- **#11 "List of shingles changes when I restart PC"** (closed). The shingle set was returned as a
  Python `set`, so ordering varied between runs; the reporter also saw SHA-1 values differ across
  restarts in the string-return path. Resolved by the reporter switching to a list; the author
  acknowledged the set was used for de-duplication. Affects the explainability/shingle-mapping path,
  not the folded bit vector.
- **#24 "AttributeError: module 'tmap' has no attribute 'Minhash'"** (open since 2024-05-02, 5
  comments). The reference implementation does not run against current `tmap`.
- **#20 "Looks like the search is linear and not using a 'forest'?"** (open since 2023-01-18, no
  reply) — about the hosted MAP4 similarity-search app.
- **#22 "facing issue while running MAP4"**, **#18 "Problem with installation"**,
  **#15 "Dependency on Python 3.6"** — all open installation/compatibility reports.
- **#26 "Publishing the package on Pypi"** (open since 2024-10-02), alongside **#27** (a closed
  external PR that refactored map4 into a PyPI-installable package). The maintainers never merged a
  PyPI release themselves.

### 4.4 Defect found and verified during this review (not previously published)

In `mapchiral` 0.0.7 (the MAP4C reference implementation, current on PyPI), the documented `seed`
parameter of `encode()` is **silently ignored on the default code path**:

```python
def encode(mol, max_radius=2, n_permutations=2048, mapping=False, seed=42):
    if mapping:
        return get_fingerprint_with_mapping(mol, max_radius, n_permutations, seed)
    else:
        return get_fingerprint(mol, max_radius, n_permutations)   # seed not passed
```

`get_fingerprint(mol, max_radius, n_permutations)` takes no `seed` argument. Verified by
`inspect.getsource` against the installed package. Consequence: `encode(m, seed=1)` and
`encode(m, seed=99)` return identical vectors, so a user cannot vary the permutation family to
estimate the MinHash error of their own results without bypassing the public API (Section 3.3 had to
seed `rdkit.Chem.rdMHFPFingerprint.MHFPEncoder` directly). Benign for reproducibility; misleading as
an API, and it blocks exactly the uncertainty check Section 3 calls for.

### 4.5 What was *not* found

No retraction, no published rebuttal, no review article arguing against MAP4, and no independent
paper reporting a reproduction failure of the 2020 benchmark.

---

## SECTION 5 — THE ALTERNATIVES, INDEPENDENTLY EVALUATED ON PEPTIDES / LARGE MOLECULES

### 5.1 Count-based ECFP, Topological Torsion and RDKit fingerprints — the strongest third-party evidence

Adamczyk, Ludynia, Czech, *"Molecular Fingerprints Are Strong Models for Peptide Function
Prediction"*, arXiv:2501.17901 (2025). AGH University of Krakow. 6 benchmarks, **126 datasets**, no
hyperparameter tuning, LightGBM head (500 trees), fingerprints from scikit-fingerprints at default
settings (ECFP radius 2, TT path length 4, RDKit path length 7). Code:
https://github.com/AGH-ML-and-Chemoinformatics-Group/peptides_molecular_fingerprints_classification

LRGB (Long Range Graph Benchmark) peptide tasks, their Table 1:

| Model | Peptides-func AUPRC ↑ | Peptides-struct MAE ↓ |
|---|---|---|
| GraphGPS | 65.35 ± 0.41 | 0.2500 ± 0.0005 |
| GCN | 68.60 ± 0.50 | 0.2460 ± 0.0007 |
| GRIT | 69.88 ± 0.82 | 0.2460 ± 0.0012 |
| DRew | 71.50 ± 0.44 | 0.2536 ± 0.0015 |
| HDSE | 71.56 ± 0.58 | 0.2457 ± 0.0013 |
| S²GCN (prior SOTA) | 73.11 ± 0.66 | 0.2447 ± 0.0032 |
| **RDKit (count)** | **73.11** | **0.2459** |
| **TT (count)** | **73.18** | **0.2438** |
| **ECFP (count)** | **74.60** | **0.2432** |

Binary → count makes a large difference on peptides (their Table 2): Peptides-func AUPRC ECFP
70.57 → 74.60 (+4.03), TT 66.18 → 73.18 (+7.00), RDKit 63.88 → 73.11 (+9.23); Peptides-struct MAE
ECFP 0.3049 → 0.2432, TT 0.3298 → 0.2438, RDKit 0.3331 → 0.2459. Cost: ECFP + LightGBM takes **19
seconds** on a 12-core CPU for Peptides-func, versus up to 60 GPU-hours for the graph transformers.
They also reach SOTA on three antimicrobial-peptide benchmarks (TT best average F1, up to +12.1 F1
over the best BERT-style model on DRAMP).

This is the best third-party evidence that exists for any representation on large flexible molecules,
and it is for **count-based ECFP / TT / RDKit**, which are exact, C++-implemented, and in RDKit.
MAP4 was not in the comparison.

### 5.2 MHFP6

Independently evaluated by Boldini et al. (2024) on natural products: MHFP gives the best **PR-AUC
(0.669)** of 20 fingerprints and excels on the antitumour dataset (0.89 ROC-AUC, 0.82 PR-AUC);
it is one of three "most promising" encodings named alongside ASP and LSTAR
(https://pmc.ncbi.nlm.nih.gov/articles/PMC10964529/). Independently evaluated by Venkatraman et al.
(2024) in VS: AUC 0.81 / 0.81 / 0.54 / 0.53 and DRF 0.17 / 0.10 / 0.97 / 0.94 — second only to MAP4
on DEKOIS/DUD-E (https://europepmc.org/articles/PMC11356940). MHFP6 carries the **same MinHash
estimation error** as MAP4 (Section 3) and the same sketch semantics; it is not an escape from that
issue, only from the C3-symmetry and speed issues.

### 5.3 RDKit atom pairs (with chirality) and Topological Torsion

- **TT** is independently strong: best MUV AUC (0.61) and best MUV DRF (0.72) of 32 fingerprints in
  Venkatraman et al. (2024), and near-best on DEKOIS/DUD-E (0.80/0.80, DRF 0.15/0.10). Best MAE on
  LRGB Peptides-struct among fingerprints in Adamczyk et al. (2025).
- **Atom pairs**: the picture is implementation-dependent and should not be over-read. Venkatraman's
  `AP2D` (CDK, 4096 bits) is the *worst* circular/path entry in their table (AUC 0.64/0.66/0.49/0.51,
  DRF 0.55/0.39/1.24/1.00), while their `AT2D` atom-triplet variant is strong (0.78/0.79/0.58/0.55).
  Boldini et al. (2024) find count-based AP captures natural-product repeat motifs better than binary
  AP. Landrum's collision study flags `hashap` as the fingerprint most damaged by folding
  (Spearman R 0.874 at 1,024 bits versus 0.990 for mfp2), i.e. atom pairs genuinely need more bits.
  Patterson et al. (1996) rate atom pairs mid-table for neighbourhood behaviour (mean enhancement
  1.15, significant in 33% of dataset/descriptor cases).
  The MAP4 paper's own Table 4 shows unhashed RDKit AP at 1.7% zero-distance NNs on HMDB.
- **Morgan with chirality**: ECFP is the best descriptor in MoleculeACE's 24-method activity-cliff
  benchmark ("ECFPs yielding the lowest average prediction error on average"; SVM + ECFP best overall;
  https://doi.org/10.1021/acs.jcim.2c01073), best in Walters' 24-dataset regression comparison against
  MAP4, and best on LRGB Peptides-func in Adamczyk et al. (2025). Its documented weakness is exactly
  what MAP4 was built to fix: 72.9% zero-distance nearest neighbours on HMDB.

### 5.4 MXFP

**No independent evaluation found.** Every evaluation located is from the Reymond group: Awale &
Reymond, *Mol. Inf.* (2019) https://doi.org/10.1002/minf.201900016; Capecchi, Zhang & Reymond,
*JCIM* (2019) PDGA; and the MAP4 2020 benchmark, in which MXFP "perform[s] significantly worse" than
MAP4 on small molecules. Code: https://github.com/reymond-group/mxfp_python (MIT) and a mirror at
https://github.com/markusorsi/mxfp-python. MXFP is a 217-D real-valued fuzzy pharmacophore atom-pair
vector compared by Manhattan distance, so it has **no MinHash sampling error** — but also no
third-party validation at all.

### 5.5 Monomer-level fingerprints

No third-party benchmark was found that pits monomer/residue-level (e.g. HELM-based) peptide
fingerprints against atom-level fingerprints with reported numbers. The relevant recent surveys
(e.g. *Drug Discovery Today* review of computational peptide tools,
https://www.sciencedirect.com/science/article/pii/S1359644626000176, 2026) catalogue tools rather
than benchmark them. Treat this as an evidence gap, not as a negative result.

### 5.6 Which alternative has the best third-party evidence for large flexible molecules?

**Count-based ECFP (Morgan) and count-based Topological Torsion**, on the strength of Adamczyk et al.
(2025): 126 datasets, 6 benchmarks, state of the art against tuned long-range GNNs and protein
language models, zero tuning, seconds of CPU time, and the representations are exact, deterministic
and in RDKit. MHFP6 is the best-evidenced *MinHashed* alternative but inherits MAP4's sampling error.
MXFP has no independent evidence whatever.

---

## SECTION 6 — PRACTICAL PROPERTIES

All repository facts verified against the GitHub REST API and raw file contents on 2 October 2026.
All speeds measured for this review (`bench.py`, same scratchpad): Python 3.12.3, RDKit 2026.03.5,
Intel Core i7-14700HX, **single-threaded**, best of 3 runs (2 for peptides); 20 drug-like molecules
(mean 15 heavy atoms) and 10 peptides built with `Chem.MolFromSequence` (mean 124 heavy atoms).

### 6.1 Measured speed

| Implementation | ms / drug-like molecule | ms / peptide | × slower than RDKit Morgan (peptides) |
|---|---|---|---|
| MAP4C — `mapchiral` 0.0.7, 2048 perm | 2.404 | 51.7 | ~490× |
| MAP4 — `skfp.fingerprints.MAPFingerprint`, 2048 | 0.618 | 23.1 | ~218× |
| MHFP6 — `skfp.fingerprints.MHFPFingerprint`, 2048 | 1.197 | 20.1 | ~190× |
| RDKit AtomPair, 2048, `includeChirality=True`, counts | 0.016 | 0.721 | ~6.8× |
| RDKit Morgan r=2, 2048, `includeChirality=True`, counts | 0.011 | 0.106 | 1× |

The authors' own assessment agrees in direction: "the current version of the MAP fingerprint is
implemented in Python and therefore it is relatively slow" (MAP4 paper, Conclusions). `mapchiral`
adds `multiprocessing` to amortise this across cores; the per-molecule cost is unchanged.

### 6.2 Implementation, licence, maintenance, packaging

| | MAP4 (reference) | MAP4C | MHFP6 | MXFP | RDKit AP / Morgan | scikit-fingerprints |
|---|---|---|---|---|---|---|
| Repo | reymond-group/map4 | reymond-group/mapchiral | reymond-group/mhfp | reymond-group/mxfp_python | rdkit/rdkit | scikit-fingerprints/scikit-fingerprints |
| Language | Python (repo is 44% CSS — it bundles a Flask search app) | Python | Python | Python | C++ with Python bindings | Python over RDKit C++ |
| **LICENSE file** | `LICENSE.md`, full **MIT** text, "Copyright (c) 2017 GDB / Reymond Research Group" | **`LICENSE` exists but is 0 bytes — empty.** GitHub reports `NOASSERTION`. README badge and a "## License" section claim MIT; PyPI metadata has **no** licence field. **No licence text is distributed anywhere.** | **MIT** | **MIT** | **BSD-3-Clause** | **MIT** (GitHub) |
| Last commit | **2021-05-10** (repo last push 2023-05-01) | **2024-12-20** (repo last push 2025-01-16) | 2023-02-16 | 2024-09-27 | active | active (PyPI 2026-08-27) |
| Stars / open issues | 125 / 12 open of 27 | 2 / 0 | 98 | 2 | — | — |
| PyPI | **`pip install map4` fails today.** `pypi.org/simple/map4/`, the JSON API, and the Tsinghua and Aliyun mirrors all return **HTTP 404** (checked 2 Oct 2026), although libraries.io records `map4` 1.1.3 published 2024-10-13 under MIT — the package appears to have been **removed**. Those 1.1.x releases were third-party repackagings (repo `setup.py` still says `version='1.0'`), opened as issue #26 and PR #27 in Oct 2024. Official install route remains `pip install git+https://github.com/reymond-group/map4@v1.0` | **yes** — `mapchiral` 0.0.7, 2024-11-25 | **yes** — `mhfp` 1.9.6, 2023-02-16 | **yes** — `mxfp` 1.1.3, 2024-09-27 | yes — `rdkit` 2026.3.6 | yes — 2.1.0 |
| conda-forge | **no** (404) | **no** (404) | **yes** | **no** (404) | **yes** | **no** (404) |
| Hard dependencies | `mhfp` **and `tmap`** (reymond-group/tmap, C++, MIT; the PyPI package literally named `tmap` is an unrelated 2019 project) | RDKit, NumPy only | RDKit | RDKit, NumPy | — | RDKit, scikit-learn |
| Default size | `dimensions=1024` in `MAP4Calculator.__init__` | `n_permutations=2048` | 2048 | 217 (fixed) | user-set | user-set |
| Similarity semantics | MinHash estimate of Jaccard; **standard Tanimoto/cosine on the raw vector is invalid** (README: "the similarity/dissimilarity between two MinHashed fingerprints cannot be assessed with 'standard' Jaccard, Manhattan, or Cosine functions … a custom kernel/loss function needs to be implemented for machine learning applications of MAP4") | same | same | Manhattan on real-valued vector | exact Tanimoto / Dice on the folded feature set | per-fingerprint |

The licence point is worth stating plainly because the task asked for it: **`reymond-group/mapchiral`,
the only reference implementation of MAP4C, ships an empty LICENSE file.** The README asserts MIT and
the PyPI package carries no licence metadata at all, so a downstream user has a claim of MIT in prose
and no licence grant in the distribution. For `reymond-group/map4` the MIT text is genuinely present
and complete.

A secondary friction: MAP4's reference implementation requires `tmap`, which is not the `tmap`
package on PyPI, and open issue #24 reports that current `tmap` no longer exposes `tmap.Minhash`.
The practical consequences are (a) MAP4 has no working one-command install, and (b) most people who
report using "MAP4" today are in fact using scikit-fingerprints' independent reimplementation, whose
numerical agreement with the reference implementation has not been published.

---

## VERDICT

**Independent evidence on MAP4 is thin, and on MAP4C it is absent.** Three third-party evaluations
exist — Venkatraman et al. 2024 (*Pharmaceuticals* 17:992), Boldini et al. 2024 (*J. Cheminform.*
16:35) and Pat Walters' 2020 blog post — and they do not agree. No independent group has evaluated
MAP4C at all; its only published numbers are the authors'.

What the independent evidence **does** support: MAP4 is a competitive, sometimes best-in-class
fingerprint for *retrieval of known actives on analogue-rich benchmarks*. On DUD-E and DEKOIS it
topped a 32-fingerprint field run by a group with no stake in it (AUC 0.83 / 0.81; DRF₀.₁ 0.07 /
0.14, both the best values in the table). Its zero-collision claim on HMDB is real and reproducible.
It produces fewer false positives than most encodings in natural-product classification.

What the independent evidence **does not** support:
1. **That high resolution implies meaningful similarity.** No source treats collision rate as evidence
   of a good similarity measure; the standing criterion since Patterson et al. (1996) is the
   descriptor-difference / activity-difference relationship, under which a random number scores zero
   despite perfect resolution. Direct measurements (Landrum's collision studies, 2022–2023;
   arXiv:2511.17078, 2025) show that eliminating collisions changes similarity values by ~0.01–0.04
   and model AUC by at most a few percent, in both directions. And MAP4's own Table 4 shows the plain
   RDKit atom-pair fingerprint already at 1.7% zero-distance neighbours — the headline 70-point gap
   is only against substructure fingerprints.
2. **That MAP4 is a good ML featurisation.** It lost to Morgan with a large effect size in 20 of 24
   ChEMBL regression datasets, and it won none of the three headline metrics in the 20-fingerprint
   natural-product QSAR benchmark.
3. **That MAP4 generalises past analogue-rich benchmarks.** On MUV and LIT-PCBA it is mid-pack, and
   the paper that ranked it first concluded that no fingerprint, MAP4 included, is a reliable proxy
   for binding activity, and that MAP4 similarity gives no basis for a meaningful threshold.

**On hard thresholds, the sceptic's instinct is correct and the quantity is now measured.** At
k = 2048 the MinHash Jaccard estimate has SE = √(J(1−J)/2048) — 0.011 at J = 0.5, confirmed here
empirically (observed SD 0.0076–0.0088 at J ≈ 0.13–0.15, matching prediction). Neither MAP4 paper
mentions standard error, variance, or thresholds anywhere in its text. In the regime MAP4 actually
occupies on large molecules — median Jaccard 0.011 across 50 M natural-product pairs — the 95%
interval on a single value is about ±41% of the value. A fixed permutation seed makes that error
reproducible, not small.

**Known failure modes are concrete and partly admitted by the authors:** C₃ and higher rotational
symmetry in macrocycles (valinomycin, 714 of 1,376 stereoisomers; NP213, 13 of 20 — worse than the
simpler MAP2C), with the authors noting that larger bit sizes do not fix it; disconnected/salt forms
producing empty shingles with a sentinel distance of 100,000,000; an unresolved open bug report of
unrelated molecules scoring similarity 1.0; and a verified defect in `mapchiral` 0.0.7 where the
documented `seed` argument is silently ignored on the default path.

**Practically, the costs are real.** MAP4C measured at 51.7 ms per peptide and 2.4 ms per drug-like
molecule single-threaded — roughly 490× and 150× slower than RDKit Morgan/atom-pair respectively.
`pip install map4` does not work today (PyPI 404 on three endpoints). The MAP4C repository's LICENSE
file is empty despite an MIT claim in the README. Neither is on conda-forge.

**If the decision is which representation to trust for large flexible molecules on third-party
evidence alone**, the answer is not MAP4: it is count-based ECFP and count-based Topological Torsion,
which reached state of the art across 126 peptide datasets and six benchmarks against tuned
long-range GNNs and protein language models (Adamczyk et al., arXiv:2501.17901, 2025) — exact,
deterministic, C++-fast, and with no sampling error to threshold around. MAP4 remains defensible for
the one job the independent evidence actually backs: nearest-neighbour *retrieval and visualisation*
across molecules of very different sizes, where ranking matters and the absolute similarity value
does not.
