# AMP oracle dependencies — availability check

Checked 2026-10-02. Every count below was obtained by downloading the file and counting
FASTA headers (`grep -c '^>'`) or by reading the licence file directly, except where marked
"not verified".

---

## PART 1 — AMP training data suitable for retraining

### Summary table

| Dataset | Positives | Negatives | How negatives built | Licence | Route | Redistributable in a benchmark package? |
|---|---|---|---|---|---|---|
| AMP Challenge 2027 `data/antibacterial.fasta` | 39,448 | **none** | n/a | BSD-3-Clause (repo); upstream MarLys CC0 | git clone, direct file | Yes — but positives only |
| BATTLE-AMP `data/amp_positive` / `amp_negative` | 20,313 | 20,181 | UniProt-sampled non-AMPs (per upstream models) | MIT (repo) | git clone, direct file | Yes |
| BATTLE-AMP `data/activity/broad_*` | 3,415 | 940 | **experimental**: DBAASP MIC above inactive threshold | MIT (repo), DBAASP attribution | git clone | Yes, with DBAASP acknowledgement |
| BATTLE-AMP `data/syntax/broad_positive_shuffled.fasta` | — | 37,565 | **composition-preserving shuffles of the positives** | MIT | git clone | Yes |
| DBAASP | ~17.6k in MarLys snapshot | activity-resolved | MIC measurements, not a label set | "access, download, copy, use, adapt, redistribute without restriction" + attribution | REST API / web | Yes, with acknowledgement |
| DRAMP v5.0 general | 12,784 | none | n/a | CC BY 4.0 | direct HTTP file | Yes |
| DRAMP v5.0 clinical | 97 | none | n/a | CC BY 4.0 | direct HTTP file | Yes |
| APD6 | 6,309 | none | n/a | "All Rights Reserved" | web browse; no bulk file found | **No** |
| AMPlify | 3,338 tr / 835 te | 3,338 tr bal / 835 te bal / 25,689 te imbal | UniProt non-AMP, length-matched | GPL-3.0 | git clone, direct file | Yes, but GPL-3 |
| AmPEP / amPEPpy | 3,268 | 3,268 (two samplings) + 166,791 full pool | UniProt non-AMP, length-proportional or random subsample | GPL-3.0 | git clone, direct file | Yes, but GPL-3 |
| Macrel | reuses AmPEP + iAMP-2L | same | same | MIT (code); data not shipped | script fetches, **plus manual PDF extraction** | No — must fetch, and partly manual |
| ampir | not shipped | not shipped | UniProt/SwissProt keyword split | GPL-2 (package) | reconstruct via `AMP_pub` Rmd pipeline | No |
| AMP Scanner v2 | 712 tr / 354 eval / 712 te | 712 / 354 / 712 | random UniProt draw ("DECOY") | GPL-3.0 | git clone, direct file | Yes, but GPL-3 |

### Detail

#### AMP Challenge 2027 — `szczurek-lab/amp-challenge-2027`
- https://github.com/szczurek-lab/amp-challenge-2027 (2026).
- `data/antibacterial.fasta`: **39,448 sequences** confirmed (4,318,298 bytes). Headers are
  `MLAMP#######` ids with `len=`, `charge=`, `disulfide=`, `dbs=` (AMPDB|APD|BaAMPs|CAMP|
  CancerPPD|DBAASP|DRAMP|InverPep|SATPdb|dbAMP|…) and `activity=` fields.
- **Positives only.** The repo README uses it purely as a novelty/exclusion reference:
  the 50,000-sequence submission library must contain no exact match to it, and the top-100
  list must stay at or below 80% Levenshtein ratio to every entry. It is explicitly *not*
  a labelled training set, and several participant repos state the organisers told them not
  to treat it as supervised training data
  (e.g. https://github.com/HarshMulodhia/amp-ctmc-2027/blob/main/DATA.md).
- Licence: repo `LICENSE` is **BSD 3-Clause, Copyright (c) 2026, Ewa Szczurek lab** (verified
  by reading the file). Upstream provenance: the `MLAMP` ids come from the **MarLys AMP
  database (MLAMP_db)**, Marczak, Bocian & Łyskowski, Mendeley Data v3, 2026-03-09,
  doi:10.17632/w4hb5grjwb.3, **CC0 1.0**, ~163,892 sequences aggregated from 13 primary AMP
  databases. So the 39,448 file is redistributable twice over.
- Format: plain FASTA. Not split. Route: direct file in a git clone (no registration).
- **Verdict for retraining: insufficient on its own — it has no negatives.**

#### BATTLE-AMP itself — `szczurek-lab/battleamp-snakemake`
- https://github.com/szczurek-lab/battleamp-snakemake (2026). Licence file verified:
  **MIT**, copyright Szymczak, Bukała, Zarzecki et al. 2026. Model submodules under `models/`
  carry their own upstream licences; the data and pipeline do not.
- This is the single most useful find in Part 1. Verified contents and counts:

| Path | Headers | What it is |
|---|---|---|
| `data/amp_positive.fasta` | 20,313 | binary-classification positives |
| `data/amp_negative.fasta` | 20,181 | binary-classification negatives |
| `data/activity/broad_positive.fasta` | 3,415 | MIC-derived actives, broad spectrum |
| `data/activity/broad_negative.fasta` | 940 | MIC-derived **experimentally inactive** peptides |
| `data/activity/gramplus_positive.fasta` | 2,064 | Gram-positive actives |
| `data/activity/gramplus_negative.fasta` | 1,136 | Gram-positive inactives |
| `data/activity/gramminus_positive.fasta` | 2,937 | Gram-negative actives |
| `data/activity/gramminus_negative.fasta` | 1,075 | Gram-negative inactives |
| `data/syntax/broad_positive_shuffled.fasta` | **37,565** | **composition-preserving shuffles of the broad positives** |
| `data/syntax/synthetic_random.fasta` | 10,000 | uniform-random decoys |
| `data/syntax/synthetic_realistic.fasta` | 10,000 | realistic non-AMP decoys |
| `data/activity/{species,strain}/`, `data/dbaasp/`, `data/paired/`, `data/slay/` | (not counted) | species- and strain-resolved MIC tasks, DBAASP export, activity-cliff pairs, SLAY |
| `datasets/battleamp-all/sequences.fasta` | (18.4 MB) | the union scoring set |
| `datasets/holdout-2026`, `paired-safe`, `pyelonephritis-borgards`, `slay`, `deep-amp` | — | evaluation splits incl. a post-training-cutoff holdout |

- Negatives in `data/activity/*` are **experimental inactives** (DBAASP MIC at or above the
  `inactive` threshold in `config/config.yaml`), not UniProt draws — this is the defining
  difference from every classic AMP dataset below and is exactly what the paper argues for.
- `data/activity/*.csv` carry the per-measurement MIC records, so a regression target is
  available as well as a binary one.
- Raw per-model predictions (309 MB) are on Zenodo, doi:10.5281/zenodo.22661234. Not needed
  for retraining.
- Route: `git clone` (no registration, no API key). Already split by task via
  `config/config.yaml`; `datasets/holdout-2026` is a time-based holdout.
- Caveat: DBAASP-derived rows carry DBAASP's attribution requirement (below), and at least
  one AMP-Challenge participant repo states DBAASP's terms stopped them redistributing their
  raw DBAASP export (https://github.com/Ragnk95/amp-challenge-2027_Francos_Group ABSTRACT.md).
  DBAASP's own terms page contradicts that — see next item.

#### DBAASP
- https://dbaasp.org. Terms (https://dbaasp.org/terms-and-conditions): data "may be accessed,
  downloaded, copied, used, adapted, and redistributed without restriction", conditional on
  acknowledging DBAASP as the source and citing the 2021 *Nucleic Acids Research* DBAASP v3
  paper. So **redistribution is permitted with attribution**; it is not an OSI licence but it
  is permissive in substance.
- Access: REST API advertised at https://dbaasp.org/api?page=rest. The site is JS-rendered;
  `GET /api/v1/peptides?limit=2` returned HTTP 200 with an empty body in this check, so the
  exact endpoint shape was **not verified** here. Practical route used by BATTLE-AMP: a
  pre-made export committed under `data/dbaasp/` in the benchmark repo.
- Content: MIC measurements per peptide/strain, not a pre-labelled positive/negative split.
  Negatives have to be defined by an activity threshold, as BATTLE-AMP does.
- Current live entry count **not verified** (statistics page returned HTTP 500). The MarLys
  aggregation records 17,628 DBAASP sequences as of early 2026.

#### DRAMP
- http://dramp.cpu-bioinfor.org. **v5.0, June 2026**, 32,296 entries site-wide. Licence stated
  on the site: **CC BY 4.0**. Patent data needs separate authorisation; general and clinical
  data require citing the original author articles.
- Download route: direct HTTP, no registration, e.g.
  `http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta`
  (the path still says `DRAMP3.0_new` even on the v5.0 site).
- Verified counts by download: `general_amps.fasta` **12,784**; `Antibacterial_amps.fasta`
  **28,702**; `patent_amps.fasta` 18,715; `natural_amps.fasta` 13,537;
  `synthetic_amps.fasta` 17,698. `clinical_amps.txt` is a 98-line TSV = **97 clinical
  entries**, with columns `DRAMP_ID, Sequence, Name, …, Activity, Stage_of_development,
  Company, Target_Organism, Pubmed_ID, clinicaltrial`.
- Formats: .fasta, .txt (TSV), .xlsx. v5.0 adds SMILES exports (`general_smiles.*`,
  `patent_smiles.*`) and a "FASTS" format covering non-natural residues.
- **Positives only, no negatives, no split.** Redistributable under CC BY 4.0 with attribution.

#### APD3 / APD6
- https://aps.unmc.edu. Current version is **APD6** (dated 2026-01-01): 6,309 peptides —
  3,379 natural AMPs, 2,290 synthetic, 373 predicted.
- Footer reads "Copyright 2003-present Dept of Pathology, Microbiology & Immunology, UNMC,
  All Rights Reserved". **No open licence.** No bulk download file was found from the landing
  page; the site is a browse/search interface.
- **Not redistributable inside a benchmark package.** Positives only in any case. Note that
  AMP Scanner v2's positives derive from APD3, and those 1,778 sequences *are* redistributed
  under GPL-3 in the AMP Scanner repo — a licensing inconsistency worth not relying on.

#### AMPlify — `bcgsc/AMPlify`
- https://github.com/bcgsc/AMPlify. Licence: **GPL-3.0** ("Copyright 2020 Chenkai Li … under
  the terms of the GNU General Public License … version 3", with a commercial-licensing
  contact). Data ships in-repo under `data/`.
- Verified counts: `AMPlify_AMP_train_common.fa` **3,338**; `AMPlify_AMP_test_common.fa` **835**;
  `AMPlify_non_AMP_train_balanced.fa` **3,338**; `AMPlify_non_AMP_test_balanced.fa` **835**;
  `AMPlify_non_AMP_test_imbalanced.fa` **25,689**; `AMPlify_non_AMP_train_imbalanced.fa`
  (14.9 MB, not counted — ~100k).
- Negatives: UniProt/SwissProt sequences without antimicrobial annotation, length-matched to
  the positives; the imbalanced variants preserve the real class ratio.
- **Already split** into train/test, balanced and imbalanced. FASTA. Direct file, no registration.
- Redistributable only under GPL-3, which is copyleft — fine for a GPL benchmark, a problem
  for a permissively licensed one.

#### AmPEP / amPEPpy — `tlawrence3/amPEPpy`
- https://github.com/tlawrence3/amPEPpy, licence **GPL-3.0**. Original AmPEP paper:
  Bhadra et al., *Sci Rep* 8:1697, 2018, https://www.nature.com/articles/s41598-018-19752-w —
  3,268 AMPs and 166,791 non-AMPs, 19 random-forest classifiers at different pos:neg ratios,
  10-fold CV.
- Verified in-repo counts: `M_model_train_AMP_sequence.numbered.fasta` **3,268**;
  `M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta` **3,268**;
  `M_model_train_nonAMP_sequence.numbered.randomsubsample.fasta` **3,268**. The full
  `M_model_train_nonAMP_sequence.numbered.fasta` is 28.8 MB (the 166,791 pool).
- Negatives: UniProt sequences lacking AMP annotation; two published subsampling schemes
  (length-proportional and random) at 1:1.
- Not pre-split into train/test (the paper used CV). FASTA, direct file.

#### Macrel — `BigDataBiology/macrel`
- https://github.com/BigDataBiology/macrel, code licence **MIT**. Training data is **not
  shipped**: `train/README.md` instructs you to fetch AmPEP's two zips from
  `ShirleyWISiu/AmPEP`, the HemoPI-1 FASTAs from
  `https://webs.iiitd.edu.in/raghava/hemopi/data/HemoPI_1_dataset/{main,validation}/{pos,neg}.fa`,
  and the iAMP-2L dataset **by extracting it from two supplementary PDFs**
  (`http://www.jci-bioinfo.cn/iAMP/Supp-S1.pdf`, `Supp-S2.pdf`), which the README itself says
  "still requires a bit of manual labour".
- So Macrel is a fetch-by-script dependency at best, and the iAMP-2L part is not automatable.
- Relevant to Part 2: Macrel trains **two** models, an AMP classifier and a **hemolysis
  classifier**, from `build-AMP-table.py` / `build-Hemo-table.py` / `train-models.py`.

#### ampir — `Legana/ampir`
- https://github.com/Legana/ampir, CRAN package, licence **GPL-2**. No training data in the
  package (`data-raw/` holds only `AAidx.rds` and model-build scripts). The dataset
  construction lives in https://github.com/Legana/AMP_pub as R Markdown
  (`01_collate_databases.Rmd`, `02_build_training_data.Rmd`) plus `build_data_package.sh`;
  the data itself is not committed. Positives from SwissProt antimicrobial-keyword entries,
  negatives from SwissProt non-AMP proteins, at the precursor-protein level rather than the
  mature-peptide level.
- **Must be fetched/rebuilt by script.** Two models, "precursor" and "mature".

#### AMP Scanner v2 — `dan-veltri/amp-scanner-v2`
- https://github.com/dan-veltri/amp-scanner-v2, licence **GPL-3.0**. Dataset also offered as a
  91.9 kB zip from https://www.dveltri.com/ascan/v2/about.html.
- Verified in-repo counts under `original-dataset/`: `AMP.tr.fa` **712**, `AMP.eval.fa` **354**,
  `AMP.te.fa` **712**; `DECOY.tr.fa` **712**, `DECOY.eval.fa` **354**, `DECOY.te.fa` **712**.
  Total 1,778 per class, matching the paper's "1,778 positive / 1,778 negative".
- Positives: APD3 entries with Gram-positive and/or Gram-negative activity. Negatives:
  1,778 sequences drawn at random from UniProt ("DECOY"), not composition-matched.
- **Already split** three ways. Smallest of the candidates, and the negative construction is
  precisely the one BATTLE-AMP shows is defeated by shuffles.

### Also found, not on the required list

- **AMPBenchmark** (Sidorczuk et al., *Brief Bioinform* 2022, bbac343,
  https://doi.org/10.1093/bib/bbac343; data at https://biogenies.info/NegativeDatasets/).
  Supplies **training sets generated by 11 different negative-sampling methods in five
  replications**, plus matched benchmark sets, expressly so that a new model can be trained
  once per sampling scheme. This is the closest existing thing to a controlled negative-set
  ablation and is directly relevant to a shuffle-robustness claim. Licence not verified.
- **ESCAPE** (arXiv 2511.04814, 2025): ~80,000 peptides, multilabel AMP benchmark with a
  transformer baseline.
- **AMPBench-MT** (arXiv 2607.25518, 2026): homology-controlled (MMseqs2 30% clusters)
  benchmark spanning binary AMP, MIC regression, spectrum, **low-toxicity, HC50 and
  selectivity** endpoints; balanced binary task of 60,946 peptides split 42,664 / 9,141 / 9,141.
  Relevant to both dependencies. Licence and download route not verified.
- **QMAP** (Lavertu, Corbeil & Germain, 2026, cited in AMPBench-MT): "a benchmark for
  standardized evaluation of antimicrobial peptide MIC and hemolytic activity regression".
  Not independently verified.
- **dbAMP 3.0** (2025), **MarLys MLAMP_db** (CC0, above) — positives-only aggregations.

### Part 1 answer

**Yes.** `szczurek-lab/battleamp-snakemake` is MIT-licensed, clones without registration, and
ships 20,313 / 20,181 binary positives and negatives, 3,415 / 940 experimentally-grounded
broad-spectrum actives and inactives (plus Gram-split variants), and — already built —
**37,565 composition-preserving shuffled negatives** in
`data/syntax/broad_positive_shuffled.fasta`, with per-measurement MIC CSVs and a 2026 holdout.
It can be redistributed inside a benchmark package provided DBAASP is acknowledged.

Redistributable as-is: BATTLE-AMP (MIT), AMP Challenge `antibacterial.fasta` (BSD-3/CC0,
positives only), DRAMP (CC BY 4.0, positives only), DBAASP (attribution).
Redistributable only under copyleft: AMPlify, AmPEP/amPEPpy, AMP Scanner v2 (all GPL).
Must be fetched by script, not shipped: Macrel (and its iAMP-2L part needs manual PDF
extraction), ampir.
Cannot be redistributed: APD6.

---

## PART 2 — Hemolysis and peptide-cytotoxicity predictors

Measured runtimes below are from this machine (CPU only, single core unless noted), in a clean
venv, 2026-10-02. They are evidence that the tool runs offline, not a benchmark of the hardware.

### Deciding question answered first

| Tool | Offline? | Evidence |
|---|---|---|
| **Macrel** | **Yes** | `pip install macrel` (1.6.1), ran 1,000 peptides in **0.68 s wall, 269 MB RSS**, output table has `Hemolytic` and `Hemolytic_probability` columns beside `AMP_probability` |
| **HemoPI2** | **Yes** | `pip install hemopi2` (1.3); the wheel ships **all weights, 224 MB** (`hemopi2_ml_clf.sav` 5.9 MB, `HemoPI2_reg.sav` 98 MB, ESM2-t6 `pytorch_model.bin` 31 MB, MERCI Perl). RF classifier: **14.5 s / 100 peptides**, **106 s / 1,000**. RF regressor (HC50 in µM): **11.2 s / 100** |
| ToxinPred2 / ToxinPred3 | Yes | `pip install toxinpred2` / `toxinpred3`, GPL-3, weights in repo (10–13 MB repos) — not run here |
| AMPDeep, Plisson et al., HEPAD, ConsAMPHemo, HemoNet | Code published | repos exist with weights or model dirs; HemoNet reported non-functional by third party |
| HemoPI (2016) | Yes (Java) | `HemoPI.zip`, JDK 7/8; superseded by HemoPI2 |
| HemoPImod | Web form | input is a **PDB tertiary structure**, not a sequence |
| **HemoPred** | **Dead** | `http://codes.bio/hemopred/` returned **HTTP 404** in this check |
| **HAPPENN** | **Web form; unreachable** | `https://research.timmons.eu/happenn` **timed out** in this check; no code or weights were ever released, only the dataset |
| HLPpred-Fuse, HemPepPred, BERT-HemoPep60 | Web form | see below |

### Two gotchas found by running them

1. **HemoPI2's pip package shells out to bare `python3`**, not `sys.executable`
   (`os.system(f'python3 {nf_path}/../Model/composition_calculate.py …')` in
   `hemopi2/python_scripts/hemopi2_classification.py` line 527). Inside a venv whose interpreter
   is not first on `PATH` this silently fails to write `{wd}/out2` and the run dies with
   `FileNotFoundError: 'wd/out2'`. Putting the venv's `bin` on `PATH` fixes it. It also declares
   neither `torch` nor `transformers` as dependencies — both must be installed by hand.
2. `-wd/--working` is **required** in HemoPI2, though the README's usage block shows it as optional.

### Per-tool detail

| Tool | Predicts | Scale | Training data | Test performance (split) | Input | Non-canonical | Code/weights | Licence | Last update |
|---|---|---|---|---|---|---|---|---|---|
| **HemoPI** (Chaudhary 2016, *Sci Rep* 6:22843) | hemolytic / not | binary + SVM score | HemoPI-1 552+552; HemoPI-2 462+462; HemoPI-3 885+738 | HemoPI-1 95.3% acc, MCC 0.91 (5-fold) / 96.4%, 0.93 (20% held out); HemoPI-2 78.0%, MCC 0.56 (CV) / 75.7%, 0.51 | single-letter or FASTA | no | **Yes** — Java jar `HemoPI.zip` + Android apk | article CC BY 4.0 | 2016 |
| **HemoPI2** (Rathore 2025, *Commun Biol* 8:176) | hemolytic / not **and HC50** | binary + score; HC50 in µM | 1,926 peptides (3,147 DBAASP v3 + 560 Hemolytik, dedup): 891 hemolytic / 1,035 non. Cut-off HC50 ≤ 100 µM. ≥6 residues, natural AA only | 80/20 split, 5-fold CV inside train. Independent test: RF+MERCI **AUROC 0.921**, acc 83.5%, MCC 0.670; ESM2-t6+MERCI AUROC 0.919; RF alone 0.888. Regression **R 0.739, R² 0.543** | FASTA or one-per-line | **no** | **Yes** — GitHub + `pip install hemopi2`, weights bundled; Zenodo 10.5281/zenodo.14676712 | **GPL-3.0** (repo `LICENSE.txt`) | repo pushed 2026-07-13 |
| **HemoPred** (Win 2017, *Future Med Chem* 9:275) | hemolytic / not | binary | HemoPI datasets | RF on amino-acid/dipeptide/physchem; AUC ≈0.73 as re-measured by HemoNet authors | web form | no | **No** | n/a | web server **404 in this check** |
| **HAPPENN** (Timmons & Hewage 2020, *Sci Rep* 10:10869) | hemolytic / not | probability 0–1 | **HAPPENN dataset: 3,738 peptides — 1,543 hemolytic, 2,195 non**; 3,408 from DBAASP + 1,174 from Hemolytik, 844 overlap. Also HAPPENN-RR90 and **HAPPENN-hard** (each positive paired with the most compositionally similar negative) | 10-fold CV **85.66% acc, MCC 0.71**; RR90 82.73% / 0.65; **HAPPENN-hard 77.54% / 0.55** | FASTA + flags for N-acetylation / C-amidation | only N-term acetylation and C-term amidation; "non-natural amino acids are not supported" | **No** — web server only; dataset released as supplementary | article CC BY | 2020; server did not respond here |
| **"HemolytikPred"** | — | — | — | — | — | — | — | **No tool of this name was found.** The nearest things are **Hemolytik** (Gautam 2014, *NAR* 42:D444), a database, now **Hemolytik 2.0** — `raghavagps/Hemolytik2`, **13,215 sequences** verified by download, GPL-3.0, includes chemically modified peptides and a REST API — and **HLPpred-Fuse** (Hasan 2020, *Bioinformatics*), a web-server classifier | | | | |
| **HemoPImod** (Kumar 2020, *Front Pharmacol* 11:54) | hemolytic potency of **chemically modified** peptides | binary | 583 modified hemolytic + 582 non-hemolytic | RF on fingerprints: 78.33% acc / AUC 0.86 (main), 78.29% / 0.85 (validation) | **PDB tertiary structure file** — not FASTA, not SMILES | yes, that is the point | repo `raghavagps/HemoPI-MOD` holds `hemopimod.tar.gz` (9.6 MB) + train/test zips, MIT | MIT (repo) | repo pushed 2026-05-13 |
| **HemoPI-MOD2** (`anandr88/HemoPI-MOD2`, 2026) | HC50 of chemically modified peptides | quantitative | — | — | — | yes | **repo is a stub: LICENSE + README only, 1 kB** | MIT | 2026-07-27 |
| **ToxinPred** (2013) | toxic / non-toxic peptide | binary | 1,805 toxic peptides ≤35 residues | — | FASTA | no | datasets in `raghavagps/ToxinPred` (MIT); original SVM via web | MIT (data repo) | 2026-05 re-upload |
| **ToxinPred2** (Sharma 2022, *Brief Bioinform* 23:bbac174) | toxic / non-toxic **protein** | binary | main 8,233 + 8,233; non-redundant 1,924 + 1,924 (CD-HIT 40%); realistic 1,924 + 19,240. Negatives: Swiss-Prot `NOT toxin NOT allergen AND reviewed:yes`, ≥35 aa, `BJOUXZ` discarded | hybrid BLAST+MERCI+RF: **AUC 0.99**, Sn 93.69%, Sp 97.39%, MCC 0.91 (validation) | FASTA or one-per-line | **no** (non-standard characters discarded) | **Yes** — GitHub + `pip install toxinpred2` | **GPL-3.0** | repo pushed 2026-08-05 |
| **ToxinPred3** (Rathore 2024, *Comput Biol Med* 179:108926) | toxic / non-toxic peptide | binary | ensemble ML + MERCI motifs | not extracted here | FASTA or one-per-line | no | **Yes** — GitHub (10 MB incl. `model/`) + `pip install toxinpred3` | **GPL-3.0** | repo pushed 2026-05-08 |
| **Macrel** (Santos-Júnior 2020, *PeerJ* 8:e10555) | AMP **and** hemolytic, in one pass | binary + probability | hemolysis RF trained on **HemoPI-1, 442 + 442**; 22 engineered descriptors | "performs similarly to the state of the art … with enhanced precision"; per-metric numbers in the paper | FASTA (peptides mode), also contigs/reads | no | **Yes** — `pip install macrel`, bioconda, GitHub | **MIT** | PyPI 1.6.1 |
| **HemoNet** (Yaseen 2021, *J Bioinform Comput Biol* 19:2150021) | hemolytic / not | binary | DBAASP + Hemolytik; SeqVec embedding + **circular fingerprint of the N/C-terminal modification's SMILES** | **AUC-ROC 0.88** vs 0.73 for HemoPI/HemoPred (5-fold and non-redundant CV) | sequence + N/C-terminal modification (name or SMILES) | **yes, partially** — N/C-terminal modifications, D- and L-amino acids | repo has `weights.hdf` + `options.json`, but **no licence file** and last pushed **2021-03-11**; the HemoPI2 authors report the GitHub code was **non-functional** when they tried to benchmark it | none declared | 2021 |
| **AMPDeep** (Salem 2022, *BMC Bioinformatics* 23:389) | hemolytic / not | binary | three hemolysis datasets; Prot-BERT-BFD fine-tuned via secretion-task transfer | "state of the art on three hemolysis datasets"; HemoPI2's independent benchmark put it at **AUC 0.602** | FASTA | no | **Yes** — `milad73s/AMPDeep`, training + preprocessing scripts | **MIT** | 2022-08-15 |
| **Plisson et al. 2020** (*Sci Rep* 10:16581) | non-hemolytic nature + activity | binary + activity | HemoPI datasets, XGBoost on physicochemical descriptors with an applicability domain | HemoPI2 benchmark: **AUC 0.740** (independent) | CSV/FASTA via notebooks | no | **Yes** — repo with `Models/`, `Data/`, `Scripts/` | **MIT** | 2026-05-16 |
| **HLPpred-Fuse** (Hasan 2020, *Bioinformatics* 36:3350) | hemolytic + activity | two-layer binary | HemoPI + PEPred-SUITE | — | web form | no | **No** | n/a | 2020 |
| **HemoDL** (2024, *Anal Biochem*) | hemolytic / not | binary | double-ensemble, sequence + transformer features | — | — | no | not verified | — | 2024 |
| **ConsAMPHemo** (Xie 2025, *Protein Sci* 34:e70087) | high/low hemolysis **and HC50** | binary + HC50 | ProtBert-BFD + Siamese contrastive learning, AAindex/BLOSUM62; 5–50 aa | — (two-stage: classifier then regressor) | FASTA | no | **Yes** — `Cpillar/ConsAMPHemo`, 120 MB incl. `model/` and `Dataset/`, **no licence file** | none declared | 2024-12-26 |
| **HEPAD** (Chen 2025, *BMC Bioinformatics* 26:234) | hemolytic / not | binary | 57 descriptors → 3,992 features, adaptive selection; HMP1 negatives are **random SwissProt peptides length-matched to the positives**, HMP2/HMP3 are harder | MCC **0.973 / 0.643 / 0.609** on three independent sets; +1.9–13.3% MCC over prior tools | FASTA/list via Python API | no | **Yes** — `csh07/HEPAD`, 11 MB with `MLProcess`, `userPackage`, models; **no licence file** | none declared | 2024-04-17 |
| **AmpLyze** (Qiu, Feng & Poczos, arXiv:2507.08162, 2025) | **HC50 regression** | pHC50 = −ln(HC50) | 1,926 peptides, Hemolytik + DBAASP v3 (same set as HemoPI2) | **PCC 0.756, MSE 0.987**, stratified 5-fold CV; mutation-direction PCC 0.844 | FASTA (ProtT5 + ESM2 embeddings, BiLSTM + cross-attention) | no | **No code link on the arXiv abstract page**; none found | paper CC BY 4.0 | v2 2025-08-13 |
| **HemPepPred** (Li 2025, *Foods* 14:4143) | **HC50 regression** | HC50 µM | 951 peptides with experimental-condition descriptors; PLM embeddings + handcrafted descriptors, SHAP + Calibrated_Explanation | not extracted | web form | no | **No** — web only, `http://hem.cqudfbp.net` (HTTP 200 in this check) | MDPI CC BY | 2025 |
| **BERT-HemoPep60** (Cai … Siu, *IEEE JBHI*, 2026, doi:10.1109/JBHI.2026.3717818) | **quantitative hemolysis**, HC5 / HC10 / HC50, six mammalian species via prefix prompts | regression | domain-adaptive pretraining; peptides **up to 60 aa**, human RBC focus | 5-fold CV **PCC 0.7431 (HC5), 0.8088 (HC10), 0.7606 (HC50)** | FASTA | no | **web app at `https://app.cbbio.online/hemopep60/home` (HTTP 200); downloadable code not confirmed** (page is JS-rendered) | not verified | 2026-07-28 |

### Benchmarks that supply a hemolysis test set rather than a predictor

- **QMAP** — Lavertu, Corbeil & Germain, *Sci Rep* 2026 (PMID 42236847; bioRxiv
  doi:10.64898/2026.02.03.703041; code https://github.com/anthol42/QMAP, pushed 2026-08-20,
  **no licence file**; dataset on Hugging Face `anthol42/qmap_benchmark_2025`;
  `pip install qmap-benchmark`). Two tasks: **MIC regression (~4,000 peptides, E. coli)** and
  **HC50 regression (~800 peptides)**, both from DBAASP v3 (18,033 sequences processed).
  Homology-aware splits at **60% identity** using a GraphPart-style similarity graph with Leiden
  community detection, five seeds. Its headline finding: **"low predictability for hemolytic
  activity"** — HC50 regression is substantially worse than MIC regression for every baseline,
  including an ESM2-650M linear probe. This is the single most important caveat for Part 2.
- **AMPBench-MT** (arXiv:2607.25518, 2026) — MMseqs2 30% cluster splits; includes
  19 low-toxicity, **17 HC50** and **15 selectivity** evaluations.

### Part 2 answer

A hemolysis predictor that **runs offline inside an optimisation loop exists and was verified
running here**. Two, in fact, at opposite ends of a speed/richness trade-off:

- **Macrel** (MIT, `pip install macrel`): binary hemolytic call plus probability, 1,000 peptides
  in **0.68 s**, and it emits the AMP probability in the same table — one call gives both terms of
  a selectivity objective. Trained on HemoPI-1's 442+442, so the label is coarse.
- **HemoPI2** (GPL-3, `pip install hemopi2`): binary **plus a quantitative HC50 in µM**, weights
  bundled in the wheel, AUROC 0.921 / R 0.739 on its own 20% independent test. Costs about
  **0.1 s per peptide**, so ~10 peptides/s, and needs `torch`/`transformers` installed by hand and
  the venv's `python3` first on `PATH`.

HemoPred and HAPPENN are web-only and neither responded in this check. HemoNet is the only tool
that models non-canonical residues from sequence, and its code is unlicensed, four years stale and
reported broken. Nothing in this inventory accepts SMILES or HELM for the whole peptide; HemoNet
takes SMILES only for terminal modifications, and HemoPImod takes a PDB structure.

---

## PART 3 — Prior art on the mitigation

### 3a. Has anyone trained an AMP classifier with composition-preserving shuffles as hard negatives?

**Yes — OmegAMP, from the same lab that wrote BATTLE-AMP, and it is the only clear instance found.**

**OmegAMP** — Soares, Hetzel, Szymczak, Der Torossian Torres, Sommer, de la Fuente-Nunez, Theis,
Günnemann & Szczurek, arXiv:2504.17247 (v1 2025-04, v2); code
https://github.com/szczurek-lab/OmegAMP, **MIT**, pushed 2026-10-01.

Its classifier training augmentation generates exactly three synthetic negative families, quoted
from the paper:

- **Random (R)**: "Generate sequences s=(a₁,…,aₗ), where each amino acid aᵢ∼U(𝒜)".
- **Shuffled (S)**: "Given a known AMP sequence s∈𝒮ₗ¹, generate a shuffled sequence π(s), where
  π∈𝒫ₗ is a random permutation." — composition, net charge and hydrophobicity are preserved by
  construction; only order is destroyed.
- **Mutated (M)**: "Starting from a known AMP sequence, randomly select 5 distinct positions and
  replace each amino acid with a new one sampled uniformly."

Reported misclassification rates for the OmegAMP classifier: **0.4% on random, 0.5% on shuffled,
0.7% on mutated**, against **20.4%–99.7%** for the baseline classifiers on the same decoys.
The wet-lab follow-up tested 25 designs and 24 (96%) showed antimicrobial activity.

The repository ships the trained classifiers directly: `models/broad-classifier.json` plus five
species classifiers (A. baumannii, E. coli, K. pneumoniae, P. aeruginosa, S. aureus) and seven
strain-level classifiers (ATCC-numbered). MIT licence, no registration.

Two caveats worth recording:

1. **The evaluation is partly circular.** OmegAMP reports its shuffled-decoy FPR on the same
   decoy family it trained against. BATTLE-AMP's result is the harder, external one, and
   **OmegAMP's classifier is not among BATTLE-AMP's 21 benchmarked variants** — the two papers do
   not cross-check each other.
2. OmegAMP's shuffle is a **full random permutation** (unigram-preserving). It does not preserve
   dipeptide/k-mer composition, which is the stricter control used in the genomics literature.

**Related but not the same thing:**

- **Porto et al., *J Theor Biol* 2018** (ScienceDirect S0022519317302205) —
  "Antimicrobial activity predictors benchmarking analysis using shuffled and designed synthetic
  peptides". Shuffles used **for evaluation only**: 78 sequences (40 designed, 38 shuffled) from
  Loose et al.'s linguistic model, against four web servers and one standalone. Finding: "the
  systems failed on predicting shuffled versions of designed peptides, as they are identical in
  AMPs composition, which implies in accuracies below 30%". So the shuffle attack was documented
  **eight years before BATTLE-AMP**, and nobody fixed it in the interim.
- **AMPBenchmark** / Sidorczuk et al., *Brief Bioinform* 2022, bbac343 (data at
  https://biogenies.info/NegativeDatasets/) — 11 negative-sampling schemes in five replications,
  from 5,190 DBAASP v3.0 AMPs. All 11 are **UniProt keyword/length/localisation filters**
  (Wang, CS-AMPpred, iAMP-2L, Gabere&Noble, AMAP, AMPScanner V2, dbAMP, Witten&Witten, AmpGram,
  ampir-mature, AMPlify). **None of them shuffles the positives.** This is the most systematic
  negative-sampling study in the field and the shuffle control is simply absent from it.
- **HAPPENN-hard** (Timmons & Hewage 2020, hemolysis not antimicrobial) — pairs each positive with
  "the most compositionally similar non-hemolytic peptide, as measured by the Euclidean distance
  between their amino acid composition vectors". Accuracy drops from 85.66% to **77.54%**, MCC
  0.71 → 0.55. This is composition-matched hard-negative *training and evaluation*, but by
  retrieval of real peptides rather than by shuffling.
- **"Learning peptide properties with positive examples only"** (*Digital Discovery* 2024,
  S2635098X24000767) — a PU-learning study whose negative pool of 13,585 sequences comes from
  "insoluble and hemolytic peptides, as well as, **the scrambled positives**". Scrambled negatives,
  but for solubility/hemolysis/non-fouling endpoints, not antibacterial activity.
- The genomics literature treats k-mer-preserving shuffled negatives as a known-double-edged
  control: a 2020 *PLOS One* study of regulatory-sequence prediction and a TF-binding study both
  report that models trained against dinucleotide-shuffled negatives learn to separate shuffled
  from real sequences rather than the biology (genomic-background negatives ~86% AUROC vs
  dinucleotide-shuffled ~76% and plain-shuffled ~63% on high-quality test sets). Expect the same
  failure mode here: shuffles must be *added* to realistic negatives, not substituted for them.

### 3b. What BATTLE-AMP itself recommends

BATTLE-AMP (Szymczak, Bukała, Zarzecki, Sala, Borišek, Fadavi, Olayo-Alarcon, Sroka,
Colomé-Tatché, Gambin, Müller, Setny, Szczurek; bioRxiv 2026, doi:10.64898/2026.06.19.733349;
pipeline https://github.com/szczurek-lab/battleamp-snakemake, MIT):

- Scope of the failure: on SyntheticShuffled decoys, **nine of 21 model variants exceed FPR 0.80**
  (MBC-Attention 0.87, AMPredictor 0.89; the APEX regressor variants are the exception at
  0.00–0.01), and **no model achieves MCC above 0.17** on that task.
- Its recommendation, verbatim: *"training on curated, species-resolved MIC data and using
  composition-preserving shuffled sequences as hard negatives during training, which would
  directly penalize the compositional shortcuts identified here."* It additionally suggests
  incorporating *"structural information through predicted conformations, dynamics-derived
  features, or structure-aware embeddings."*
- It also reports that **models trained on MIC data outperform binary classifiers regardless of
  architecture** — i.e. the recommended fix is regression-on-MIC *plus* shuffled hard negatives,
  not binary classification plus shuffled hard negatives.
- Survey context: of 48 published methods surveyed, **fewer than 25% were reproducible**.
- Explicit scope limit: **"the benchmark does not yet address toxicity prediction."** BATTLE-AMP
  says nothing about hemolysis or selectivity.

So the planned mitigation (a) is exactly what BATTLE-AMP asks for, and OmegAMP is an existing
implementation of it with MIT-licensed weights. Mitigation (b) is outside BATTLE-AMP's scope and
has no endorsement from it.

### 3c. Is selectivity an established objective, and in what formulation?

**Yes — it is the standard medicinal-chemistry endpoint for AMPs, under the names
"selectivity index" (SI) and "therapeutic index" (TI), and it is now a formal benchmark endpoint.**

Formulations in use, most to least common:

| Formulation | Where |
|---|---|
| **SI = HC50 / MIC** (ratio of 50%-hemolysis concentration to MIC; "sometimes also indicated as a 'therapeutic index'") | Juretić/Zoranić et al., *BBA Biomembranes* 2012 (adepantins) states the definition explicitly and builds a predictive model for SI itself from the AMPad database, using non-homologous sequences at <70% pairwise identity |
| **SI = HC50 / MIC or IC50 / MIC** (erythrocytes or a normal cell line) | *Sci Rep* 2023, s41598-023-43274-9, cell-selective AMP design from crocodile haemoglobin |
| **SI = HC10 / MIC** (the more conservative 10%-hemolysis threshold) | species-aware language-model AMP discovery, PMC12271573 |
| **minimise MIC / HC50** (the reciprocal, as a loss) | **QMAP** (Lavertu, Corbeil & Germain, *Sci Rep* 2026): "AMP development aims to maximize potency while minimizing toxicity, which can be expressed as minimizing the ratio MIC/HC50" |
| **SI_log = pMIC − pHC50** (log-scale difference, the numerically well-behaved version) | **AMPBench-MT** (arXiv:2607.25518, 2026) — a first-class benchmark endpoint with **9,774 paired rows**, alongside 1,595 HC50 rows and 7,143 low-toxicity rows, split by MMseqs2 at 30% identity / 0.8 coverage. Data at https://huggingface.co/datasets/ZihengZhou06/AMPBench-MT under a **"Research and Review Use License"** — *not* an open licence |
| **Pareto / multi-objective**: maximise antimicrobial score, minimise hemolysis, no scalarisation | MoFormer (arXiv:2406.02610), HMAMP (arXiv:2405.00753, hypervolume-driven, Kneedle knee-point selection), PepZOO (PMC11725395, multi-objective zeroth-order optimisation), *Int J Mol Sci* multi-objective de novo design for *S. aureus* (PMC11728188) |
| **constraint/filter**: generate, then drop anything a toxicity or hemolysis classifier flags | the common industrial pattern — Macrel's own pipeline (AMP classifier then hemolysis classifier), OmegAMP's discriminator-guided filtering, AIGCRS-AMP30 (classifier-guided diffusion + MIC regression + HC50 regression + multi-criterion screening) |

Three things to note about the log form. `SI_log = pMIC − pHC50` with `pX = −log₁₀(X [mol/L])` is
the same quantity as `log(HC50/MIC)`, it is finite and signed where the raw ratio blows up when
MIC → 0, and it is what the only benchmark that scores selectivity directly uses. If a single
scalar objective is wanted, that is the one with a published, homology-controlled evaluation
behind it.

The counterweight, from **QMAP** on homology-aware 60%-identity splits: **"low predictability for
hemolytic activity"**. HC50 regression is markedly worse than MIC regression for every baseline
they tried, including an ESM2-650M linear probe. A selectivity objective built on a predicted
HC50 inherits that error, and the error is in the denominator.

---

## VERDICT

### Dependency 1 — AMP training data suitable for retraining with shuffled hard negatives

**AVAILABLE.**

`git clone https://github.com/szczurek-lab/battleamp-snakemake` (MIT) gives, verified by download:
20,313 / 20,181 binary positives and negatives; 3,415 / 940 broad-spectrum MIC-derived actives and
**experimental** inactives, with Gram-positive (2,064 / 1,136) and Gram-negative (2,937 / 1,075)
variants; per-measurement MIC CSVs for the regression target BATTLE-AMP actually recommends; and
**37,565 composition-preserving shuffled negatives already generated** in
`data/syntax/broad_positive_shuffled.fasta`, with 10,000 random and 10,000 realistic decoys beside
them and a `holdout-2026` time split. The only condition on redistribution is acknowledging DBAASP.

Supplement with the AMP Challenge `antibacterial.fasta` (39,448, BSD-3 over CC0 MarLys) or DRAMP
v5.0 (12,784 general, CC BY 4.0) if more positives are wanted; both are positives-only.
AMPlify / AmPEP / AMP Scanner v2 are all GPL-3, which would make a permissively licensed benchmark
package copyleft. APD6 cannot be redistributed at all.

### Dependency 2 — Hemolysis predictor that runs offline in an optimisation loop

**AVAILABLE**, and verified running on this machine.

- **Macrel 1.6.1**, MIT, `pip install macrel`: **1,000 peptides in 0.68 s**, 269 MB RSS, CPU only,
  and the same output table carries `AMP_probability` and `Hemolytic_probability` — both terms of
  a selectivity objective from one call. Trained on HemoPI-1 (442 + 442), so a coarse binary label.
- **HemoPI2 1.3**, GPL-3, `pip install hemopi2`: binary **plus HC50 in µM**, all weights bundled in
  the 224 MB wheel, independent-test AUROC 0.921 and regression R 0.739. Costs **~0.1 s/peptide**
  (14.5 s per 100, 106 s per 1,000), needs `torch` and `transformers` installed by hand, requires
  `-wd`, and only works if the venv's `python3` is first on `PATH` (it shells out to bare `python3`).

The sensible pairing is Macrel inside the loop and HemoPI2 as the HC50 scorer on shortlists, or
HemoPI2 throughout if ~10 peptides/s is affordable.

**The caveat that matters more than availability:** QMAP (*Sci Rep* 2026) finds **low predictability
for hemolytic activity** under homology-aware 60%-identity splits — HC50 regression is far weaker
than MIC regression for every baseline tested. The predictor exists and runs; its accuracy on
peptides unlike the training set is the real risk, and a selectivity objective divides by it.

**Not usable:** HemoPred (web server returned HTTP 404), HAPPENN (web form only, server timed out,
no code or weights ever released), HLPpred-Fuse / HemPepPred / BERT-HemoPep60 (web forms),
HemoNet (unlicensed, last touched 2021, reported non-functional by the HemoPI2 authors).
**"HemolytikPred" does not exist** — the name most likely refers to the **Hemolytik** database
(now Hemolytik 2.0: 13,215 sequences, GPL-3, `raghavagps/Hemolytik2`, includes chemically modified
peptides and a REST API) or to **HLPpred-Fuse**, which is a web server.

### Prior art

Mitigation (a) is **already implemented and released**: OmegAMP (MIT, `szczurek-lab/OmegAMP`)
trains its AMP classifiers against random, **shuffled** and mutated synthetic negatives and reports
0.5% FPR on shuffled decoys versus 20.4–99.7% for baselines. Its evaluation is on the decoy family
it trained on, and BATTLE-AMP did not benchmark it, so the external check is still missing.
Mitigation (b) is standard practice under the name selectivity index; the formulation with a
published homology-controlled benchmark behind it is **SI_log = pMIC − pHC50** (AMPBench-MT), and
the equivalent ratio form **HC50/MIC** goes back at least to 2012.
