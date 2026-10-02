# Peptide / peptidomimetic similarity: what the field actually does

Scope: documented practice only (released code read where public, papers quoted where not). No recommendations.
Compiled 2026-10-02.

**Reference fact used throughout.** RDKit defaults chirality OFF. From the RDKit Python wrappers:
`GetMorganFingerprint`: `radius`, `useChirality = false`, `useBondTypes = true`, `useFeatures = false`, `useCounts = true`;
`GetHashedMorganFingerprint` / `GetMorganFingerprintAsBitVect`: `nBits = 2048`, `useChirality = false`, `useFeatures = false`
(https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp, lines ~1074–1103).
New-style generator: `rdFingerprintGenerator.GetMorganGenerator(radius = 3, countSimulation = false, includeChirality = false, …, fpSize = 2048)`
(https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Fingerprints/Wrap/MorganWrapper.cpp, lines 93–99).
So any call that does not pass `useChirality`/`includeChirality` silently makes L-Ala and D-Ala identical.

---

## SECTION 1 — GENERATIVE PEPTIDE DESIGN PAPERS: WHAT THEY ACTUALLY USED

### 1.0 Every fingerprint call in the seven codebases, verbatim

| Repo | file:line | exact call | chirality |
|---|---|---|---|
| MolecularAI/PepINVENT | `pepinvent/supervised_learning/utils/chem.py:28` | `AllChem.GetMorganFingerprint(mol, radius=4, useChirality=True, useCounts=True)` | **ON (explicit)** |
| MolecularAI/PepINVENT | `pepinvent/scoring_function/scoring_components/predictive_model.py:51` | `AllChem.GetMorganFingerprint(Chem.MolFromSmiles(query_smiles), radius=4, useChirality=True, useCounts=True)` | **ON (explicit)** |
| MolecularAI/PepINVENT, MSDLLCpapers/PepEvolve (via `reinvent_chemistry`) | `reinvent_chemistry/conversions.py:22-26` | `AllChem.GetMorganFingerprint(mol, radius, useCounts=use_counts, useFeatures=use_features)` with defaults `radius=3, use_counts=True, use_features=True` | **OFF (default)** |
| MSDLLCpapers/PepEvolve | `pepinvent/supervised_learning/utils/chem.py:28`, `predictive_model.py:52` | identical to PepINVENT's radius-4 chiral count FP | **ON (explicit)** |
| programmablebio/peptune, sophtang/PepTune | `src/scoring/functions/scoring_utils.py:19`, `permeability.py:40` | `AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=size)`, defaults `radius=3, size=2048` | **OFF (default)** — the repos contain no other fingerprint call and no metric code |
| charlesxu90/helm-gpt | `utils/metrics_utils.py:26` (and `smi_utils.py:21`, `agent/scoring/kras.py:37`, `kras_ic50.py:37`, `permeability.py:26`) | `AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=size)` / `AllChem.GetHashedMorganFingerprint(molecule, radius, nBits=size)`, defaults `radius=3, size=2048` | **OFF (default)** — this is the fingerprint behind the published diversity and SNN |
| MSDLLCpapers/PepThink-R1 | — | no fingerprint call in the public repo; the paper says only "Tanimoto similarity … between the Morgan fingerprints" | **unstated** |
| szczurek-lab/hydramp | — | none (sequence model) | n/a |
| BirolLab/AMPd-Up | — | none (sequence model) | n/a |
| novonordisk-research/pepfunn (comparator) | `src/pepfunn/similarity.py:172,195` | `AllChem.GetMorganFingerprintAsBitVect(mol, 4, 2048)` | **OFF (default)** |
| novonordisk-research/pepfunn | `src/pepfunn/clustering.py:115,196,232,289` | `AllChem.GetMorganFingerprintAsBitVect(x, 4, 1024)` | **OFF (default)** |
| novonordisk-research/pepfunn | `src/pepfunn/similarity.py:241`, `sequence.py:574` | `rdMolDescriptors.GetMorganFingerprint(mol, 2, useChirality=True)` (monomer matrix, Dice) | **ON (explicit)** |
| molecularsets/moses (template) | `moses/metrics/utils.py` | `Morgan(molecule, morgan__r, nBits=morgan__n)` = `AllChem.GetMorganFingerprintAsBitVect`, defaults `r=2, n=1024` | **OFF (default)** |
| BenevolentAI/guacamol (template) | `guacamol/utils/chemistry.py:289` | `AllChem.GetMorganFingerprintAsBitVect(m, radius, length)`, defaults `radius=2, length=4096` | **OFF (default)** |

**Summary of the chirality question.** The only place `useChirality=True` appears in any of the seven generative codebases is PepINVENT's (and its fork PepEVOLVE's) training-monitoring and QSAR featurisation — neither of which is a reported metric. Every fingerprint that produces a *published number* — HELM-GPT's diversity and SNN, PepTune's stated r3/2048 Morgan, and the MOSES/GuacaMol templates they cite — leaves RDKit's default in place and is therefore blind to D/L inversion, to epimerisation, and to any stereochemical edit the models are explicitly designed to make.

### 1.1 PepINVENT (AstraZeneca)
Paper: arXiv:2409.14040 (2024); Chem. Sci. 2025, DOI 10.1039/D4SC07642G.
Code: https://github.com/MolecularAI/PepINVENT (read at commit e976dd1, 2026-10).

- **Diversity:** not a scalar. The paper reports a t-SNE chemical-space plot of generated amino acids. Fingerprint stated verbatim: *"The 1024 bit Morgan fingerprints with radius=3, useChirality=True and useCounts=True computed with RDKit v.2024.03.5 was projected to 2-dimensional space with Scikit-learn v.0.24.2"* (arXiv:2409.14040, 2024). **Chirality explicitly ON.**
- **Novelty:** set membership at the amino-acid building-block level — generated AAs classified as natural / non-natural-in-training-set / novel (= canonical SMILES not present in training set). No fingerprint, no threshold (2024/2025).
- **Uniqueness:** reported at three levels deliberately separating stereochemistry — raw string, canonical **isomeric** SMILES, canonical SMILES with stereochemistry stripped. The paper notes uniqueness drops when chirality is removed. This is the only paper in this set that treats chirality as a first-class variable.
- **SNN:** not computed.
- **Similarity-based reward:** **none.** The released scoring components are exactly `maximum_ring_size`, `molecular_weight`, `substructure_match`, `predictive_model`, `lipophilicity`, `custom_alerts`
  (`pepinvent/scoring_function/component_names_enum.py`).
- **Similarity that IS in the code:**
  - Training-time monitoring, `pepinvent/supervised_learning/utils/chem.py` line 28:
    `AllChem.GetMorganFingerprint(mol, radius=4, useChirality=True, useCounts=True)` → `DataStructs.TanimotoSimilarity`. Unhashed count FP, radius 4, **chirality ON**. Used in `transformer_trainer.py:211` to score sampled-vs-source similarity each epoch.
  - QSAR featurisation, `scoring_components/predictive_model.py:51`: same radius-4, `useChirality=True`, `useCounts=True` fingerprint.
  - Inherited Mol2Mol paired-dataset Tanimoto, `reinvent_models/mol2mol/dataset/paired_dataset.py:143` → `reinvent_chemistry.Similarity.calculate_tanimoto` over `Conversions.smiles_to_fingerprints`, whose defaults are
    `AllChem.GetMorganFingerprint(mol, radius=3, useCounts=True, useFeatures=True)` — i.e. **FCFP6 counts with chirality OFF**
    (https://github.com/MolecularAI/reinvent-chemistry/blob/main/reinvent_chemistry/conversions.py lines 22–40; similarity.py uses `DataStructs.BulkTanimotoSimilarity`).
  → Three different Morgan parameterisations coexist in one repository (r3/1024/chiral for t-SNE, r4/counts/chiral for monitoring and QSAR, r3/counts/features/achiral inherited from REINVENT).
- **Diversity filter (RL):** REINVENT's, carried over unchanged. Shipped configs use
  `"name": "NoFilterWithPenalty"` or `"IdenticalMurckoScaffold"`, `"score_threshold": 0.4`, `"bucket_size": 25`, `"similarity_threshold": 0.4`, `"penalty": 0.5`
  (`data/experiment_configurations/config_scenario1_cyclic.json`, `config_crbp_peptide.json`).
  `IdenticalMurckoScaffold` computes `MurckoScaffold.GetScaffoldForMol` then `Chem.MolToSmiles(scaffold, isomericSmiles=False)` — **scaffold comparison with stereochemistry discarded**, and for a linear peptide with no ring the Murcko scaffold is the empty string. The `similarity_threshold` field is not read by that filter class.

### 1.2 PepEVOLVE (Merck & Co.)
Paper: arXiv:2511.16912 (2025). Code: https://github.com/MSDLLCpapers/PepEvolve.

- The paper contains **zero** occurrences of "Tanimoto", "similarity", "fingerprint" or "Morgan" (full-text scan of https://arxiv.org/html/2511.16912v1). Reported metrics are score trajectories, best/mean reward, convergence steps. No diversity metric, no novelty metric, no SNN.
- Code inherits PepINVENT/REINVENT wholesale: same `reinvent_models/mol2mol/dataset/paired_dataset.py` FCFP6-counts-achiral Tanimoto; same scoring components plus `synthetic_accessibility`, `hydrogen_bond_donor`, `logp` — **still no similarity component**.
- Every manuscript config (`data/manuscript/**/*.json`) carries `"name": "IdenticalMurckoScaffold"`, `"score_threshold": 0.4`, `"bucket_size": 25`, `"similarity_threshold": 0.4`, `"penalty": 0.5`. The README documents `similarity_threshold` default `0.4` as *"Tanimoto similarity threshold for grouping molecules into the same scaffold bucket"* — but no fingerprint is named and the filter class does not use the value.

### 1.3 PepTune / PepMDLM (Duke)
Paper: arXiv:2412.17780v4 (2025); ICML 2025 PMLR v267:59017–59065. Code: https://github.com/programmablebio/peptune.

- **Metric suite:** explicitly "the Moses metrics [Polykovskiy et al. 2020]" — validity, uniqueness, diversity, SNN, plus Shannon-entropy "randomness" and token-level KL.
- **Diversity** (Appendix F.2, eq. 38): `1 − mean_{i<j} Tanimoto(f(x_i), f(x_j))`, where *"f(x_i) and f(x_j) are the 2048-dimensional Morgan fingerprint with radius 3"*. No chirality statement.
- **SNN** (eq. 39): `max_j Tanimoto(f(x_i), f(x̃_j))` over the dataset — same fingerprint.
- **SNN/novelty reference set:** *"Due to the limit of memory and CPU time required to load all the training dataset of 11 million peptide SMILES, we chose to sample a subset of 1000 batches randomly (∼100k sequences) for novelty and SNN calculation."* SNN is therefore against ~0.9 % of the training set.
- **Validity** is their own `SMILES2PEPTIDE` filter, not RDKit-only.
- **Train/test split:** k-means into 1000 clusters on Morgan fingerprint, 0.8/0.2 by cluster. The training data is *"90 % canonical amino acids, 10 % unnatural amino acids from SwissSidechain, 10 % dextro-chiral alpha carbons, 20 % N-methylated amine backbone atoms, 10 % PEGylated peptides"* — i.e. D-centres are 10 % of residues while the similarity measure, as specified, has no chirality term.
- **Released code:** the repository contains **no** diversity/SNN/novelty implementation. The only Morgan calls are property-predictor featurisation:
  `AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=size)` with `radius=3, size=2048` in `src/scoring/functions/scoring_utils.py:19` and `src/scoring/functions/permeability.py:40`. `useChirality` not passed → **False**. Repo-wide grep for `tanimoto|Chirality|MACCS|diversity|novelty|FCD|MAP4` returns nothing else.
- **Similarity-based reward:** none. Objectives are binding affinity, permeability, solubility, hemolysis, non-fouling. One sibling-node similarity check exists as a flag only (`pareto_mcts.py:348 def expand(self, parentNode, eps=1e-5, checkSimilarity=True)`); the embedding-cosine helper in `src/utils/generate_utils.py` compares PeptideCLM encoder outputs, not fingerprints.

### 1.4 HELM-GPT (KAUST)
Paper: Bioinformatics 40(6) btae364 (2024). Code: https://github.com/charlesxu90/helm-gpt — `utils/metrics_utils.py`, read in full.

This is the most completely specified case, because the code is the metric.

```python
def fingerprints_from_mol(molecule, radius=3, size=2048, hashed=False):
    if hashed:  fp_bits = AllChem.GetHashedMorganFingerprint(molecule, radius, nBits=size)
    else:       fp_bits = AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=size)
```
- **Fingerprint: ECFP6 bit vector, radius 3, 2048 bits, `useChirality` not passed → False.** HELM peptides are converted to SMILES first (`get_cycpep_smi_from_helm`) via a 3104-monomer library built from ChEMBL + CycPeptMPDB + KRAS HELMs.
- **Diversity:** `diversity = 1 - (average_agg_tanimoto(fps, fps, agg='mean', p=1)).mean()` (line 156) — mean pairwise Tanimoto over *unique canonical SMILES only*.
- **SNN:** `snn = average_agg_tanimoto(self.ref_fps, fps, agg='max', p=1)` (line 158), reference = every SMILES in `data/prior/prior.csv` (the full ChEMBL prior). Labelled in the returned dict as `"snn": snn,  # "structural novelty"`.
- `average_agg_tanimoto` implements Jaccard on float tensors: `jac = tp / (x_stock.sum(1) + y_gen.sum(0) - tp)`, with `jac[np.isnan(jac)] = 1`.
- **Novelty:** pure set difference on RDKit canonical SMILES — `novelty = len(gen_smiles_set - train_set) / len(gen_smiles_set)` (line 163). `canonic_smiles` uses `Chem.MolToSmiles(mol)` with default `isomericSmiles=True`, so novelty *is* stereochemistry-aware while diversity and SNN are not.
- **Similarity-based reward:** none; the RL tasks are permeability and KRAS affinity (`agent/scoring/`), both ECFP6/2048 descriptor models with the same achiral fingerprint.
- Reported values: prior model SNN 0.750; PepTune's reproduction reports HELM-GPT diversity 0.595, SNN 0.975 (arXiv:2412.17780v4, 2025).
- The HELM-GPT authors themselves flag transfer failure of a small-molecule score: *"We suspect that the SAscore model, originally built for small molecules, may not be well suited to evaluate larger molecules such as cyclic peptides"* (Bioinformatics 2024). They do not extend that suspicion to the fingerprint.

### 1.5 PepThink-R1 (Merck & Co.)
Paper: arXiv:2508.14765 (2025), NeurIPS 2025. Code: https://github.com/MSDLLCpapers/PepThink-R1.

- **Similarity-based reward: yes, and it is the only one in this set.** Section 3.3:
  `R = dup_fac · (0.8 · prop_smooth + 0.2 · sim_fac)`, with
  *"we computed the Tanimoto similarity s between the Morgan fingerprints of the original and generated peptides"*, `sim_fac = σ(α·(s − s₀))`, **`s₀` = 0.6** ("target similarity center"). Radius, bit count and chirality are not stated anywhere in the paper, and the RL reward code is not in the public repo.
- **Duplication penalty:** `dup_fac = (1/max(1, n+1))^γ`, `n` = times the exact molecule has appeared, tracked in an LRU cache. Exact-identity, not similarity.
- **Novelty / Uniqueness (released evaluation notebook `cycpepMPDB.ipynb`, cell 32):** raw Python set arithmetic on the `SMILES` column with **no canonicalisation step**:
  ```python
  unique_smiles = set(smiles_series)
  uniqueness    = len(unique_smiles) / len(smiles_series)
  novel_smiles  = unique_smiles - orig_smiles_set
  novelty       = len(novel_smiles) / len(unique_smiles)
  ```
  No fingerprint, no Tanimoto, no SNN, no FCD anywhere in the repo (grep returns only the reward-free notebook strings).
- Reported novelty > 0.98 for nearly every model including GPT-4o — i.e. the metric is string-identity against a 1880-row CycPeptMPDB-derived training set.
- Stated limitation: *"reinforcement learning, while improving property control, reduces structural diversity"*, measured as uniqueness falling to 0.200–0.300.

### 1.6 HydrAMP (Univ. of Warsaw)
Paper: Nat. Commun. 14, 1453 (2023). Code: https://github.com/szczurek-lab/hydramp.

Sequence-space model; no molecular fingerprints anywhere in the repository.
- **Redundancy removal, positives:** exact-duplicate removal by MD5 hash of the sequence string (`scripts/dataset_preparation.ipynb`, cell 9).
- **Redundancy removal, negatives:** CD-HIT. Verbatim cell: *"## 3.2. CD-HIT - remove sequences sharing ≥40% sequence identity"*, command
  `cdhit -i ../data/raw/Uniprot_0_200_negative.fasta -o ../data/interim/uniprot_negatives_cdhit_40.fasta -c 0.4 -T 4 -n 2 -M 2000`.
  Applied **only to the UniProt negative set**, not to the AMP positives. No justification given for 0.4.
- **Similarity to prototype:** Levenshtein distance, implemented by hand (`scripts/z_sigma_improvement.ipynb`, `def levenshteinDistanceDP(token1, token2)`), reported as box plots of distance between generated analogues and Pexiganan / CAMEL (Fig. 2c,d). Unnormalised, in residues.
- **Novelty / diversity:** no scalar novelty metric; "creativity" is the count of unique generated analogues out of 10 000 attempts at temperature τ ∈ {1,2,5}. Acceptance criteria are classifier probabilities (P(AMP) > 0.8, P(low MIC) > 0.5), not similarity.
- Non-canonical residues: excluded by construction (20-letter alphabet, length ≤ 25).

### 1.7 AMPd-Up (BC Cancer / BCGSC)
Paper: Li et al., Protein Science 33:e5088 (2024). Code: https://github.com/BirolLab/AMPd-Up.

- **Similarity measure, verbatim:** *"The similarity between two sequences was calculated as 1 − d_{i,j} / max(l_i, l_j) × 100 %, where d_{i,j} is the edit distance and l_i, l_j are lengths of the sequences regarding the numbers of amino acid residues."*
- **Nearest-neighbour definition, verbatim:** *"The similarity of a sequence to a set of sequences was defined as the maximum of all similarity values calculated between that sequence and the sequences in the target set"* — i.e. a sequence-space SNN.
- **Reference sets:** the training set (APD3 antibacterial, 2019-03-20) and a "known AMP" set of 4538 distinct sequences from APD3 (2022-07-11) + DADP (2018-12-06). Both FASTA files ship in `data/`.
- **Reported result:** the generated-vs-training similarity distribution peaks at **50–55 %**; between-model-instance generated similarity 33.56 %, within-instance 39.14 %. No cutoff is defined for "too similar"; novelty is argued from the distribution.
- **Redundancy removal:** *"After removing duplicates and sequences with non-standard amino acids, we ended up with a non-redundant set of 2253 antibacterial sequences ≤50 aa"* — exact dedup plus hard deletion of non-canonical entries. No CD-HIT, no BLAST.
- No fingerprint, no Tanimoto, no FCD.

### 1.8 Section 1 cross-table

| Model (year) | diversity | novelty | SNN | similarity reward | chirality in the similarity measure |
|---|---|---|---|---|---|
| PepINVENT 2024/25 | t-SNE on Morgan r3/1024/**chiral**/counts | AA-level set membership | — | none | **yes** (t-SNE, monitoring, QSAR); **no** in inherited Mol2Mol FCFP6 and in the Murcko filter |
| PepEVOLVE 2025 | none reported | none reported | — | none | n/a (no similarity reported) |
| PepTune 2024/25 | 1 − mean pairwise Tanimoto, Morgan r3/2048 | token-KL + set diff | max Tanimoto vs 100k-sample of 11M train | none | not stated; code default = **no** |
| HELM-GPT 2024 | 1 − mean agg Tanimoto, ECFP6/2048 | canonical-SMILES set diff (isomeric) | max agg Tanimoto vs full ChEMBL prior | none | **no** (bit-vector FP), yes only in novelty's canonical SMILES |
| PepThink-R1 2025 | uniqueness only (raw SMILES strings) | raw-SMILES set diff | — | **Tanimoto/Morgan, centre s₀ = 0.6** | not stated |
| HydrAMP 2023 | count of unique analogues | — | — | none | n/a (canonical alphabet only) |
| AMPd-Up 2024 | pairwise sequence similarity distributions | max sequence similarity vs APD3+DADP | yes, sequence-space | none | n/a (canonical alphabet only) |

---

## SECTION 2 — DATABASES AND RESOURCES

### CycPeptMPDB (Li et al., J. Chem. Inf. Model. 63:2240–2250, 2023; http://cycpeptmpdb.com)
- **No similarity measure of any kind.** A full-text scan of the usage documentation returns zero hits for "similar", "Tanimoto", "fingerprint", "cluster", "t-SNE", "PCA", "UMAP" (fetched 2026-10-02).
- Search is metadata-only: publication year, permeability, assay type, original compound name, molecular weight (RDKit `MolWt`), monomer length, molecule shape (Circle / Lariat), and substring filtering of the result table.
- **Redundancy: deliberately not removed.** *"Some peptides overlapped in structure between different publications and had different membrane permeability measurements, they were recorded as separate data in CycPeptMPDB (there were 8,466 peptides including duplicated structures)"* for 7,991 unique structures (v1.2 browse page, 2026).
- Structure normalisation is at monomer level: peptides are cut at peptide/ester bonds into 385 monomer types, each assigned a "Natural Analog" class (20 AAs + unknown) by reference to PubChem and the ChEMBL monomer library. This monomer dictionary, not a fingerprint, is the database's notion of chemical relatedness.
- Downstream users supply their own split: the CycPeptMP implementation deduplicates to 6,889 PAMPA peptides and ships fixed index files (https://github.com/akiyamalab/cycpeptmp); the 13-method benchmark (J. Cheminform. 17, 2025, DOI 10.1186/s13321-025-01083-4) uses *"Murcko scaffolds … generated using the RDKit library, ignoring chirality differences"* for its scaffold split.

### DBAASP v3 (Pirtskhalava et al., Nucleic Acids Res. 49:D288–D297, 2021; https://dbaasp.org)
- **No fingerprint similarity search and no sequence-similarity search.** Retrieval is keyword/indexed-field search (ID, name, complexity, synthesis type, sequence substring, length, N/C-terminal modification, unusual amino acids, intrachain bonds, source, target, UniProt ID, activity ranges). Results are rows of *"chemically unique peptide"*.
- **Clustering** exists but is descriptor-based, not similarity-based in the fingerprint sense: a DBSCAN density clustering over physicochemical descriptors (Moon–Fleming hydrophobicity, normalised hydrophobic moment, charge, pI, penetration depth, tilt angle, propensity to disordering, linear moment, in-vitro aggregation), per Vishnepolsky & Pirtskhalava, J. Chem. Inf. Model. 58:1141–1151 (2018). The Statistics page additionally offers amino-acid composition and pairwise residue-distribution profiles.
- Its predictors require *"20 canonical amino acids … in FASTA format"* and *"length … should not exceed 30 amino acids"*, so DBAASP's own analytics exclude the non-canonical content the database stores.
- No published redundancy threshold; entries are per structurally unique peptide card including modifications.

### DRAMP (Shi et al., NAR 50:D488, 2022; Ma et al., NAR 53:D403–D410, 2025; v5.0 released 2026-09-11)
- **Similarity search = BLAST, FASTA and SSEARCH; alignment = Stretcher, Matcher, Clustal Omega.** All sequence-based, FASTA input only. No chemical fingerprint search.
- **Redundancy policy, verbatim (DRAMP 4.0):** *"whenever new data are collected, we will first perform sequence comparison, if there are duplicate sequences, we will first determine whether there are different sequence modifications, if there are, they will be included in the database as a new entry … If the sequence is different or the modification is different, the peptide will be included in the database as a new entry."* Exact-sequence matching plus modification matching; no identity threshold.
- Cross-database overlap counts (Fig. 2) explicitly *"excluded sequences that were redundant, predicted or longer than 100 amino acids"* and treat *"all peptides in the database with the same sequence but with different structural modifications … as one sequence"*.
- v5.0 (2026) added SMILES for a subset and a "FASTS" format (SMILES embedded in the FASTA annotation line) stated to *"enable the accurate representation of peptides containing non-natural residues and modifications"* — but the analysis tools remain BLAST/SSEARCH/FASTA/CD-Search. No structure-similarity search has been added.

### APD3 (Wang et al., NAR 44:D1087–D1093, 2016; APD6: NAR 54:D363–D374, 2026)
- **Redundancy rule, verbatim:** *"To reduce sequence redundancy, AMPs from different species that share the same sequence occupy a single entry in the database"* (55 such entries flagged "found in multiple species"). Exact identity only; synthetic fragments of a natural AMP are folded into the parent entry as derivatives. Inclusion criteria: natural, known sequence, demonstrated activity, ≤100 aa (later relaxed to 200).
- **Similarity:** the prediction interface returns *"five most similar sequences based on sequence alignment with all the peptides in the database"* (APD chapter, https://aps.unmc.edu/assets/pdf/apd2013w.pdf). The alignment scoring is described only as word-weighted; no algorithm name, no identity threshold, no E-value is published.
- APD3 documents a case where its own similarity search found neighbours that BLAST did not (bactofencin A, O'Shea et al.) — i.e. the database explicitly positions itself as an alternative to BLAST for short peptides.

### Pistoia Alliance HELM
- **The reference HELM toolkit has no similarity function and the Alliance says so.** HELM Challenges page: *"HELM adopters have a need to search through a repository of HELM strings: either to check for uniqueness, to find similarities… One of the limitations of the current toolkit is the lack of a search tool. Given a large set of HELM strings, it is not possible to find matches of a particular term using HELM notation or other descriptors."* The only prior art cited is a 2014 student proof-of-concept supporting exact match and substructure search.
  (https://pistoiaalliance.atlassian.net/wiki/spaces/PUB/pages/20512783/HELM+Challenges)
- What HELM does standardise is *identity*, not similarity: the canonicalisation challenge covers canonical HELM string, canonical monomer dictionary, and canonical monomer structures, aimed at *"registration of biomolecules … to create a database without redundancies"*. HELM2NotationToolkit carries some canonicalisation; chemistry is delegated to a plugin (CDK by default).
- `HELMCoreLibrary` (https://github.com/PistoiaHELM/HELMMonomerSets) is a monomer dictionary built by frequency analysis of public datasets with EBI and PubChem — again an identity/vocabulary resource.
- **Third-party fill-in:** `QuattroResearch/HELMSimilarityLibrary` (Java) — *"Based on the enumeration of monomer paths, a fingerprint for each notation can be generated. The similarity is calculated via Tanimoto coefficient of two fingerprints… there is a possible extension of taking the natural analogs of modified monomers into account."* Path length, bit size and default thresholds are not documented in the README (they are in the attached MSc thesis only). This is a monomer-level, not atom-level, fingerprint.

---

## SECTION 3 — THE MOSES / GUACAMOL METRICS AS APPLIED TO PEPTIDES

### 3.1 MOSES (Polykovskiy et al., Front. Pharmacol. 11:565644, 2020; https://github.com/molecularsets/moses)
Code read at master, 2026-10-02.

- **Fingerprint** (`moses/metrics/utils.py`): `def fingerprint(smiles_or_mol, fp_type='maccs', dtype=None, morgan__r=2, morgan__n=1024, …)`, and for Morgan: `np.asarray(Morgan(molecule, morgan__r, nBits=morgan__n), dtype='uint8')` where `Morgan = AllChem.GetMorganFingerprintAsBitVect`. **Only `mol`, `radius`, `nBits` are passed — `useChirality` is never set, so it is `False`.** Module default is MACCS (166 bits, key 0 dropped), but every metric overrides to `'morgan'`. The operative setting is therefore **ECFP4 bit vector, radius 2, 1024 bits, achiral**.
- **SNN** (`metrics.py`, `class SNNMetric`): `average_agg_tanimoto(pref['fps'], pgen['fps'], device=self.device)` with defaults `agg='max'`, `p=1`. Jaccard on float vectors: `jac = tp / (x_stock.sum(1,keepdim=True) + y_gen.sum(0,keepdim=True) - tp)`. Max over the **reference** set per generated molecule, then mean. Reference sets are `SNN/Test` and `SNN/TestSF` (held-out scaffold set), not the training set.
- **Internal diversity**: `return 1 - (average_agg_tanimoto(gen_fps, gen_fps, agg='mean', device=device, p=p)).mean()`; reported as `IntDiv` (p = 1) and `IntDiv2` (p = 2). Self-pairs are included, so the diagonal contributes T = 1.
- **Scaffold similarity (`Scaf`)**: `compute_scaffold(mol, min_rings=2)` → `MurckoScaffold.GetScaffoldForMol(mol)`, `Chem.MolToSmiles(scaffold)`, then
  `if scaffold_smiles == '' or n_rings < min_rings: return None`.
  Scores are a cosine similarity over `Counter`s of scaffold SMILES; `cos_similarity` returns `np.nan` if either Counter is empty.
  **Consequence for peptides:** a head-to-tail cyclic peptide (cyclo-(Ala)₆, cyclosporin A) has exactly one ring, so `compute_scaffold` returns `None`. For a set of such peptides the Counter is empty and **Scaf is NaN**. Linear peptides have no ring at all → empty scaffold string → also `None`.
- **Frag**: `AllChem.FragmentOnBRICSBonds`, cosine over BRICS-fragment counts. BRICS does not cut a macrolactam, so cyclo-(Ala)₆ yields 1 fragment (the whole molecule) and cyclosporin A yields 10, of which 9 are side chains and one is the entire 11-residue macrolactam — the metric degenerates to side-chain composition plus per-molecule singletons.
- **Novelty**: `len(set(canonic_smiles(gen)) - set(train)) / len(gen_set)`; `canonic_smiles` uses RDKit default `isomericSmiles=True` on the generated side only.
- **Dataset the metrics were tuned on (ZINC Clean Leads):** MW **250–350 Da**, ≤ 7 rotatable bonds, XlogP ≤ 3.5; then filters enforced in `mol_passes_filters`: allowed elements `{'C','N','S','O','F','Cl','Br','H'}`, no formal charges, and
  `if ring_info.NumRings() != 0 and any(len(x) >= 8 for x in ring_info.AtomRings()): return False`.
  4,591,276 → 1,936,962 molecules. MW 250–350 Da ≈ 17–25 heavy atoms. The same function is reported as the `Filters` metric, so **every macrocyclic peptide scores Filters = 0 by construction**, and the ≥8-atom-ring rule excludes macrocycles from the reference distribution entirely.

### 3.2 GuacaMol (Brown et al., JCIM 59:1096–1108, 2019; https://github.com/BenevolentAI/guacamol)
- **Fingerprint** (`guacamol/utils/chemistry.py`): `def get_fingerprints(mols, radius=2, length=4096): return [AllChem.GetMorganFingerprintAsBitVect(m, radius, length) for m in mols]` — **ECFP4, 4096 bits, chirality not passed → False**; similarity by `DataStructs.BulkTanimotoSimilarity`.
- **Distribution-learning benchmarks**: validity; uniqueness and novelty over `canonicalize_list(..., include_stereocenters=False)` — **stereochemistry stripped before the set comparison**; KL divergence over 9 RDKit descriptors (`BertzCT, MolLogP, MolWt, TPSA, NumHAcceptors, NumHDonors, NumRotatableBonds, NumAliphaticRings, NumAromaticRings`) plus max internal ECFP4 similarity, scored `sum(exp(-kl))/k`; FCD scored `exp(-0.2·FCD)` at `sample_size = 10000`.
- **Peptides are excluded by name.** `guacamol/data/get_data.py` docstring: *"Pre-filter molecules of 5 <= length <= 200, because processing larger molecules (e.g. peptides) takes very long."* Then in `filter_and_canonicalize`: `if len(smiles) > 200: return []`, element SMARTS `[!#1!#5!#6!#7!#8!#9!#14!#15!#16!#17!#34!#35!#53]`, charge neutralisation, ECFP4 Tanimoto cutoff to the holdout set, and
  `# Drop out if too long canonicalized:  if len(canon_smi) > 100: return []`.
  **There is no MW ceiling; the ceiling is a 100-character canonical SMILES**, roughly ≤550–600 Da, i.e. about cyclo-(Ala)₇. Cyclosporin A (canonical isomeric SMILES ≈ 227 characters) is excluded twice over.

### 3.3 FCD / ChemNet (Preuer et al., JCIM 58:1736–1741, 2018; https://github.com/bioinf-jku/FCD, insilicomedicine/fcd_torch)
- **ChemNet training distribution:** trained to predict bioactivities for *"about 6 000 assays"* drawn from **ChEMBL, ZINC and PubChem** (arXiv:1803.09518, 2018). Architecture: two 1-D convolutional layers (SELU), max-pool, two stacked LSTM layers, fully connected output. The training set itself has never been released (bioinf-jku/FCD issue #10, 2021).
- **Activation used:** the hidden state of the **second LSTM layer after the full input sequence**, dimension **512** (`fcd/fcd.py`, `fcd_torch/fcd.py`: `return np.zeros((0, 512))`).
- **Input:** character one-hot over a 35-symbol vocabulary (`fcd/utils.py`), normalised by sequence length.
- **Length limit:** `DEFAULT_PAD_LEN = 350`.
  - `fcd_torch` (the package MOSES imports): `__PAD_LEN = 350`; `get_one_hot` breaks at `dst == one_hot.shape[0] - 1`, i.e. **silently truncates at 349 tokens**, no warning.
  - `fcd` ≤ 1.2.0 raised `IndexError` on longer strings (issue #14, 2023).
  - `fcd` ≥ 1.2.1: `pad_len = max(DEFAULT_PAD_LEN, max_len)` plus `warnings.warn("Padding lengths differing from the default of 350 may affect FCD scores. See https://github.com/hogru/GuacaMolEval.")` — long SMILES no longer truncate, but the whole batch's padding changes and the score stops being comparable to published values.
- **Is 350 actually exceeded by peptides?** Measured canonical isomeric SMILES lengths: cyclosporin A **227**, octreotide **167**, cyclo-(Ala)ₙ ≈ 15 chars/residue (n = 24 → 357), cyclo-(Phe)ₙ ≈ 22 chars/residue (n = 15 → 329). **A typical 6–12-residue cyclic peptide is ~90–270 characters and does not exceed the pad length.** The limit bites around ≥16 bulky or ≥23 small residues, and for stapled / conjugated / glycosylated constructs. So the usual "FCD can't take peptide SMILES" claim is, as stated, wrong.
- **The real validity problem is the vocabulary and the training distribution.** The 35-symbol set omits `/`, `\`, `9`, `%`, `b`, `p`, `Se`; `@` is present, so tetrahedral stereocentres are resolved but E/Z double-bond stereoisomers produce bit-identical embeddings and FCD exactly 0 (bioinf-jku/FCD issue #24, 2026). No statement in the paper or README restricts FCD to drug-like space; the restriction is implicit in ChemNet's ChEMBL/ZINC/PubChem bioactivity training data and in the fact that the Frechet distance is a Gaussian fit to those activations.

### 3.4 Who applies these to peptides, and who objects
- **Applied, with substitutions.** HELM-GPT (2024) and PepTune (2024/25) both say "Moses metrics" and both change the fingerprint: HELM-GPT uses ECFP6/2048, PepTune uses Morgan r3/2048, against MOSES's r2/1024. **Neither reports FCD, Scaf or Frag.** PepTune additionally replaces RDKit validity with a bespoke peptide filter.
- **FCD explicitly refused.** ClaMP (Chemical Language Model Linker, arXiv:2410.20182, 2024/25), p. S18: *"We do not use this metric because our preliminary results and an independent evaluation found it is highly sensitive to the sample size and molecule padding length."*
- **The padding-length critique.** Holzgruber, *GuacaMolEval* (https://github.com/hogru/GuacaMolEval) — cited inside the `fcd` package's own warning string: FCD depends on both reference sample size and padding length; *"since the value of 350 is hard-coded in the `fcd` package, we could consider changing this value as 'cheating'"*.
- **Sample-size dependence.** Özçelik & Grisoni, *"How evaluation choices distort the outcome of generative drug discovery"* (2025, arXiv:2501.05457; PMC12613558): FCD and Fréchet Descriptor Distance both fall monotonically with library size and plateau only above ~10 000 designs; they recommend ≥10⁵, and report FCD failing to rank held-out actives closer to the fine-tuning set than inactives.
- **Novelty-as-set-difference is gameable.** Renz, Van Rompaey, Wegner, Hochreiter & Klambauer, *Drug Discov. Today Technol.* 32–33:55–63 (2019/2020), doi:10.1016/j.ddtec.2020.09.003: the "AddCarbon" model — add one carbon atom to a random training molecule — near-saturates GuacaMol's distribution-learning suite and beats every baseline but the LSTM on FCD. Directly relevant to a leakage audit: every peptide paper in Section 1 defines novelty as exact set difference, which an AddCarbon-style model passes at ~100 %.
- **Direct statement that small-molecule similarity fails on peptides.** Keller-Findeisen et al., *ACS Chem. Biol.* (2023), doi:10.1021/acschembio.3c00159, verbatim: *"previously developed similarity metrics for small molecules were not useful for comparing peptides, as the common backbone atoms and large size of peptides combine to reduce the dynamic range of these similarity metrics."*
- **The companion score transfers no better.** SAscore (Ertl & Schuffenhauer, *J. Cheminform.* 1:8, 2009) was validated on 40 drug-like PubChem molecules and carries an explicit penalty for rings > 8 atoms; RDKit's `sascorer.py` applies `macrocyclePenalty = math.log10(2)` whenever any macrocycle is present. HELM-GPT (2024) reports this as a visible failure in its own results.
- **No re-derivation exists.** No publication was found that recomputes MOSES/GuacaMol-style metrics with peptide-appropriate defaults and publishes them as a peptide benchmark. The field's practice is: borrow the metric names, drop FCD/Scaf/Frag, swap the fingerprint radius, and leave chirality off.

---

## SECTION 4 — SEQUENCE-BASED MEASURES IN THE AMP FIELD

### 4.1 What is actually used, per pipeline

| Tool / dataset (year) | positives | negatives | tool |
|---|---|---|---|
| AMPlify (2022), PMC8788131 | **no clustering at all**; duplicates only | keyword-filtered Swiss-Prot ≤200 aa, length-matched; no clustering | — |
| AmPEP (2018), PMC5785966 | duplicates + non-standard-AA removal only | same | — |
| amPEPpy 1.0 (2021), OSTI 1766413 | reuses AmPEP's files verbatim | same | — |
| AMPd-Up (2024), PMC11237553 | duplicates + non-standard removal only | n/a | — |
| AMP Scanner v2 (2018), PMC6084614 | **CD-HIT ≥0.9 removed** | **CD-HIT ≥0.4 removed**, plus BLAT v35 purge of anything matching a known AMP | CD-HIT + BLAT |
| Macrel (2020), PMC7751412 | **CD-HIT v4.8.1 -c 0.80 with 90 % coverage of the shorter sequence**, applied to positives *and* negatives jointly before splitting | same | CD-HIT |
| CAMPR4 (2023), PMC9825550 | CD-HIT web server, **90 %** | CD-HIT, **>90 %** | CD-HIT |
| ampir-mature | **0.9** | — | CD-HIT |
| iAMP-2L (2013), dbAMP, CS-AMPpred, AMAP, Witten & Witten | — | **0.4** | CD-HIT |
| Wang et al. 2011 | — | **0.7** | CD-HIT |
| HydrAMP (2023), PMC10017685 | **not clustered** | **CD-HIT ≥0.4 removed** (`-c 0.4 -n 2`) | CD-HIT |
| Sidorczuk et al. benchmark sets from DBAASP v3 (2022), PMC9487607 | **CD-HIT 4.8.1 >90 %** | varies | CD-HIT |
| AMPSphere / Santos-Júnior (Cell 2024), PMC10491242 | CD-HIT 4.8.1 hierarchically at **100 % → 85 % → 75 %**, on an **8-letter reduced alphabet** `[LVIMC][AG][ST][FYW][EDNQ][KR]` chosen to raise sensitivity for short peptides | same | CD-HIT + MMseqs2 |
| PepBenchmark (2026), arXiv:2604.10531 | **MMseqs2 at 90 %** for de-redundancy; **MMseqs2 at 30 %** for the split; ECFP-similarity connected components at **0.95** for non-canonical sets | same | MMseqs2 |
| iAMPpred (2017), PMC5304217; dbAMP 2.0 (2022), PMC8690246 | "removal of redundant sequences" — **no tool, no threshold stated** | — | — |

Sidorczuk et al. (*Brief. Bioinform.* 2022, PMC9487607) is the field's own audit; its Table 1 is the inventory the table above draws on. Their own choice of 90 % is justified purely by convention: *"This threshold was most frequently used for the reduction of positive data in the algorithms selected."*

**Modal pattern: 0.9 on positives, 0.4 on negatives, cited rather than derived.** No AMP paper was found that empirically justifies either number. A substantial minority (AMPlify, AmPEP, AMPd-Up, HydrAMP-positives) does no identity clustering at all, and AMPlify argues against it: *"AMPs that are highly similar to each other at the sequence level were kept as separate entries, since small changes in amino acid compositions may lead to large changes in AMP activity."*

### 4.2 Where 0.9 and 0.4 come from
- **0.9 is CD-HIT's compiled-in default.** User guide: *"`-c` sequence identity threshold, **default 0.9**"* (https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki).
- **0.4 is the algorithm's documented floor, not a biological cutoff.** Same guide, *Algorithm limitations*: *"A limitation of short word filter is that it can not be used below certain clustering thresholds… word size 5 is for thresholds 0.7 ~ 1.0 / word size 4 is for thresholds 0.6 ~ 0.7 / word size 3 is for thresholds 0.5 ~ 0.6 / word size 2 is for thresholds 0.4 ~ 0.5"*, and *"**Because of the algorithm, cd-hit may not be used for clustering proteins at <40 % identity.**"* This is enforced in the source: `cdhit-common.c++`, `Options::Validate()` line 360 —
  `if ((cluster_thd > 1.0) || (cluster_thd < 0.4)) bomb_error("invalid clstr");` (verified at master, 2026-10-02). HydrAMP's `-n 2` is exactly the word size the manual prescribes for 0.4–0.5.
- **Identity definition (verbatim from the manual):** the default is *"global sequence identity calculated as: number of identical amino acids in alignment divided by the full length of the shorter sequence"*. `-G 0` switches to *"local sequence identity… divided by the length of the alignment"* with the manual's own warning *"don't use -G 0 unless you use alignment coverage controls, see options -aL, -AL, -aS, -AS"*. `-aS`/`-aL` default `0.0`, `-s` length-difference cutoff default `0.0` — i.e. by default a short peptide fully contained in a longer one is 100 % identical to it.
- **Silent length filter:** `-l`, *"length of throw_away_sequences, default 10"* — CD-HIT discards sequences of ≤10 residues by default, which removes the shortest AMPs without comment.
- Neither Li & Godzik (*Bioinformatics* 22:1658, 2006) nor Fu et al. (*Bioinformatics* 28:3150, 2012) proposes or benchmarks 0.9 or 0.4 as meaningful similarity levels; both papers are about speed and scale. **The field's two canonical thresholds are a software default and a software limit.**

### 4.3 BLAST and alignment alternatives
- **Macrel (2020)** gives explicit parameters: blastp as a baseline classifier at *"maximum e-value of 1e−5, minimum identity of 50 %, word size of 5, 90 % query coverage, window size 10, subject besthit"*; for annotation against DRAMP/nr, *"maximum e-value of 1 × 10⁻⁵ and a word size of 3 … hits with a minimum of 70 % identity and 95 % query coverage"*. Result: *"Using blastp as a classification method was no better than random."*
- **AMP Scanner v2 (2018):** BLAT v35, default settings.
- **AMPSphere (2024):** MMseqs2 `easy-search` throughout (DRAMP matches at identity ≥75 %, E ≤ 1e−5; gene families at ≥30 % identity, ≥50 % coverage of the shorter sequence, E ≤ 1e−3; DIAMOND blastp at >50 % identity, >90 % query+target coverage). Notably they **did not trust default BLAST statistics on peptides** — Smith–Waterman with BLOSUM62 (gap −10/−0.5), E-values from Karlin–Altschul with *"the values of κ (0.132539) and λ (0.313667) constants adjusted to search for a short input sequence"*.
- **CAMPR3/R4, DRAMP, APD:** expose BLAST / PHI-BLAST / jackhmmer / SSEARCH / FASTA and global EMBOSS Stretcher–Matcher + Clustal Omega as user tools with no fixed cutoff.
- **"Analogue" convention:** ApexGO (*Nat. Mach. Intell.* 2026, PMC13201158) states the folk threshold explicitly — *"These similarity levels (typically 40 %–60 %) fall below the ≥70 % identity threshold commonly used to define analogues."*

### 4.4 Levenshtein / edit distance in AMP generative work
- HydrAMP (2023): raw, unnormalised Levenshtein distance to the prototype.
- AMPd-Up (2024): `1 − d/max(l_i, l_j) × 100 %`, max over the reference set.
- ApexGO (2026, PMC13201158): `1 − d / l_template`, candidates retained at **similarity ≥ 75 %**; the same constraint applied to HydrAMP and diffusion baselines for comparability.
- MOFormer (*Brief. Bioinform.* 2025, PMC12596111): novelty = fraction of generated sequences whose Levenshtein distance to **every** training sequence exceeds T; at **T = 3**, novelty 0.289; mean minimum Levenshtein distance 3.672.
- AMPCliff (*J. Adv. Res.* 2025, PMC12869253): compares Levenshtein, Smith–Waterman-aligned Levenshtein, sequence identity, mean Tanimoto and mean BLOSUM62, with activity-cliff thresholds *"1 for Levenshtein and Levenshtein aligned, and 0.9 for the other 3"*.
- Anti-diabetic peptide pipeline (*Sci. Rep.* 2026, DOI 10.1038/s41598-026-39985-4): generated candidates admitted only at *"global identity < 70 % to any Train positive and minimum edit distance ≥ 3–5 to the nearest Train neighbor"*, followed by a post-hoc nearest-neighbour identity audit across splits.

### 4.5 Non-canonical residues: stated plainly
- **CD-HIT does not reject them; it silently recodes them.** Source `cdhit-common.c++` lines 52–57:
  `int aa2idx[] = {0, 2, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 6};` over A…Z. Decoding: **B → 2 = N**, **Z → 6 = E**, and **J, O, U and X all → index 20**, one shared slot. Selenocysteine (U), pyrrolysine (O), and the ambiguity codes J and X are therefore the *same residue* internally and **count as identities with one another** in the identity numerator. No warning, no rejection. The embedded BLOSUM62 scores the X row −1 against everything except A/S/T (0).
- **BLASTP/BLOSUM62:** X scores −1 against all residues including itself, so X-rich peptides are penalised rather than matched. NCBI provides `-task blastp-short`, *"optimized for query sequences shorter than 30 residues"*, with different defaults (word_size 2, PAM30, gapopen 9, gapextend 1, threshold 16, vs blastp's 3 / BLOSUM62 / 11 / 1 / 11). Essentially no AMP paper reports using it.
- **D-amino acids, N-methylation, cyclisation, staples, terminal amidation and PEGylation have no representation in one-letter FASTA at all.** They are not approximated; they are absent. DRAMP 4.0 says so directly: *"all peptides in the database with the same sequence but with different structural modifications are considered as one sequence."* DBAASP v3 stores them as separate structured fields (terminal modifications, unusual amino acids, intrachain bonds), i.e. outside the sequence string, so any FASTA export silently discards them — and DBAASP's own predictors require *"20 canonical amino acids … in FASTA format"*.
- **Pipelines that simply delete non-canonical entries before training:** AmPEP 2018 (*"sequences with unnatural amino acids (B, J, O, U, X, and Z) were also eliminated"*), AMPd-Up 2024, AMPlify 2022 (*"non-standard amino acids are not taken into consideration in this study"*; such sequences are marked "Invalid"), Macrel 2020 (*"No peptides containing non-canonical amino acids were kept"*), iAMPpred 2017, CAMPR4 2023 (non-standard AAs **plus stapled and circular peptides** excluded), Sidorczuk 2022, HydrAMP 2023.
  **Net effect: the non-canonical peptides that generative models are built to produce are precisely the ones the field's redundancy and homology machinery removes from the reference sets.**

### 4.6 Published criticism of identity thresholds for peptides
- **AMPCliff (*J. Adv. Res.* 2025, PMC12869253)** is the most direct: *"our experimental results confirm that sequence identity is not a suitable metric for measuring similarity between short AMPs with canonical amino acids"*; identity *"averages similarity over the entire sequence length… This averaging operation biases sequence identity toward identifying longer sequences"* (one substitution in a 10-mer already reads as 0.9); and *"this also suggests that MMseqs2 is not suitable for clustering short peptides."*
- **NCBI BLAST FAQ:** *"virtually identical short alignments have relatively high E values. This is because the calculation of the E value takes into account the length of the query sequence… shorter sequences have a higher probability of occurring in the database purely by chance."*
- **amPEPpy (2021):** *"a major challenge in AMP discovery is the lack of sequence conservation of these peptides, limiting the effectiveness of traditional sequence homology-based search tools such as BLAST."*
- **Macrel (2020)**, empirically, under homology-controlled splits: blastp best-hit classification *"was no better than random, confirming that homology-based methods are not appropriate for this problem beyond very close homologs."*
- **Fernández-Díaz et al. (*J. Cheminform.* 2025, PMC12751563):** *"chemical fingerprint-based similarity measures outperform traditional sequence alignment-based metrics for partitioning standard peptide datasets, challenging conventional practice."*
- **PepBenchmark (2026):** sequence identity alone is insufficient even when applied — *k*-mer leakage survives an MMseqs2 split, because short peptides are dense in functional motifs.

---

## SECTION 5 — MACROCYCLES AND NON-CANONICAL PEPTIDES: WHAT WAS PICKED, AND WAS IT JUSTIFIED

### Measures designed for the problem (justification present)
- **MXFP** — macromolecule eXtended atom-pair FingerPrint, 217 dimensions: 7 pharmacophore categories × 31 exponentially spaced topological-distance bins (d₀ = 0 to d₃₀ = 317.8 bonds), each atom pair smeared as an 18 %-width Gaussian, each category normalised by N_c^1.5 *"to reduce the sensitivity of the fingerprint to molecule size"*. Compared with **Manhattan distance (City-Block Distance, CBD)**, not Tanimoto (https://github.com/reymond-group/mxfp_python, 2024). Explicit rationale: substructure fingerprints saturate on large molecules.
- **PDGA** (Capecchi, Zhang, Reymond, *J. Chem. Inf. Model.* 2019, DOI 10.1021/acs.jcim.9b01014; https://github.com/reymond-group/PeptideDesignGA): peptide genetic algorithm whose **fitness function is MXFP CBD to a query**, with a user-supplied `similarity-threshold` (example in the README: *"compounds with CBD smaller than 300 from Ala-Leu-Cys1-His-Gaba-Cys1-Ile will be annotated"*) and topology linear/cyclic. Handles cyclisation and non-canonical building blocks natively. Later version switched to MAP4 + RDKit fingerprints (`PDGA-MAP4_AP`).
- **MAP4** (Capecchi, Probst, Reymond, *J. Cheminform.* 12:43, 2020): MinHashed atom pairs of circular substructures up to diameter 4; **Jaccard** similarity (the repo warns standard Jaccard/Manhattan/cosine cannot be applied feature-wise to MinHashed vectors). Justified against an explicitly constructed **peptide benchmark** recovering BLAST analogues among scrambled and point-mutated peptides; the paper states substructure FPs win on small molecules, atom-pair FPs win on peptides, and MAP4 is the only one good at both.
- **MAP4C / MAPc** (Orsi & Reymond, *J. Cheminform.* 16:53, 2024): adds Cahn–Ingold–Prelog R/S/r/s annotation to shingles whose centre is a stereocentre, plus cis/trans. Directly motivated by peptides: it distinguishes all **2048 stereoisomers of antimicrobial undecapeptide ln65**, all 512 of nona-arginine, all 4096 of polymyxin B2, where *"ECFP6C only saw about half of them and ECFP4C and APC distinguished less than 10 %"*. It still **fails** on C3/C4/C7-symmetric macrocycles (valinomycin, nonactin, cyclic hepta-arginine NP213). Recommended by its authors as the default chiral fingerprint for peptides.
- **Hestia-GOOD / peptide partitioning study** (Fernández-Díaz, Ochoa, Hoang, López, Shields, *J. Cheminform.* 2025, PMC12751563). The only systematic comparison of similarity functions for peptide dataset partitioning; the paper states *"There is no prior work systematically evaluating the optimal similarity functions to use for similarity-based dataset partitioning for peptide datasets."* Compared MMseqs2 (± k-mer prefilter), EMBOSS `needleall` (Needleman-Wunsch), Jaccard on ECFP and on **MAPc** at diameters 4–20, plus Molformer-XL and ESM2-8M embedding similarities. Conclusion verbatim: *"The most versatile similarity function is Jaccard similarity with MAPc fingerprints, as it is the top choice in six of the eight scenarios… Attention has to be dedicated to the optimal fingerprint diameter, particularly with modified peptides, with the recommended options being diameters of at least 8."* And: *"chemical fingerprint-based similarity measures outperform traditional sequence alignment-based metrics for partitioning standard peptide datasets, challenging conventional practice."*

**Provenance caveat on the MXFP / MAP4 / MAP4C line.** Every performance claim quoted above for MXFP, MAP4 and MAP4C comes from the Reymond group's own papers and their own benchmarks (the "extended Riniker benchmark" they constructed, and the stereoisomer-enumeration test they designed). The one third-party peptide evaluation found here is Fernández-Díaz et al. 2025, which uses MAPc for dataset partitioning and reports it as the best of the functions *they* tested. Independent adoption in the generative peptide-design papers of Section 1 is **zero**.

### Monomer-level fingerprints
- **PepFuNN** (Ochoa & Deibler, Novo Nordisk, *J. Pept. Sci.* 2025; https://github.com/novonordisk-research/pepfunn). Two parallel notions:
  - monomer fingerprint `monomerFP(peptide, radius=2, nBits=1024, add_freq=False, prop_list=['heavy','nrot','hacc','hdon','nhet','tpsa','mw'])` — BILN peptide → igraph monomer graph → fragments of radius ≤ 2 → physchem-property token → hash → bit. Morgan-like but with monomers as atoms; handles cyclics, branches and non-canonicals.
  - atom-level fallback, `src/pepfunn/similarity.py:172` and `:195`: `AllChem.GetMorganFingerprintAsBitVect(mol, 4, 2048)` + `DataStructs.TanimotoSimilarity`. **Chirality not passed → False.** The paper justifies the radius: *"A radius of four in the fingerprints is used by default to capture larger fragments found in repetitive structures like peptides."* That is the clearest explicit radius justification in the peptide literature found here.
  - clustering, `src/pepfunn/clustering.py:115`: `GetMorganFingerprintAsBitVect(x, 4, 1024)` → `BulkTanimotoSimilarity` → `Butina.ClusterData(dists, nfps, cutoff, isDistData=True)` with **`run_clustering(self, cutoff=0.3)`** — i.e. default cluster at Tanimoto ≥ 0.7.
  - monomer–monomer substitution matrix, `similarity.py:241`: `rdMolDescriptors.GetMorganFingerprint(mol, 2, useChirality=True)` with **Dice** similarity and a default `threshold=60`. So PepFuNN turns chirality ON for monomer comparison and leaves it OFF for whole-peptide comparison.
- **HELMSimilarityLibrary** (above): monomer-path fingerprint + Tanimoto, with optional natural-analog back-off.
- **CycPeptMPDB**'s natural-analog monomer classes (above) are a coarse monomer-level similarity used for browsing.

### Graph edit distance
- **GLAMOUR** (Mohapatra, An, Gómez-Bombarelli, *Mach. Learn.: Sci. Technol.* 3:015028, 2022): macromolecule graph with monomers as nodes; **exact GED** where *"For insertion and deletion of node/edge, we add a fixed cost to the distance, while for substitution, we multiply a constant cost with the Tanimoto dissimilarity of the molecules being substituted"*, using **stereochemical fingerprints** for the monomer Tanimoto, node/edge insertion/deletion and substitution cost all set to 3 by grid search. Stated justification: BLOSUM62-type substitution matrices *"are based on evolutionary statistics and cannot be used for non-natural building blocks"*. Propagation-attribute graph kernels used as the scalable approximation (GED is NP-hard).
- **MacroSimGNN** (Shi et al., *Macromolecules* 2025; ChemRxiv 2024): a GNN trained to predict GLAMOUR-style pairwise GED, because *"Graph edit distance is accurate but computationally expensive, and graph kernel methods are computationally efficient but inaccurate."*
- **cyclicpeptide** Python package (Liu et al., *Brief. Bioinform.* 2024/25): `GraphAlignment` module built on NetworkX, offering **graph edit distance** and **maximum common subgraph** as the two similarity metrics, explicitly because *"traditional linear sequence alignment algorithms"* fail on cyclic structures (their benchmark: Biopython scores the cyclic permutation "PVFFAAGF" of "AAGFPVFF" at 0.54). A GCN surrogate (`GA GCN`) approximates the GED score for speed.

### Sequence/graph hybrids that avoid an explicit similarity function
PepLand (*Brief. Bioinform.* 26(4) bbaf367, 2025), PepFoundry (2026), SinCAA (2026) and PepMNet learn embeddings over atom/fragment/monomer heterogeneous graphs and report Euclidean distances in embedding space rather than a fingerprint similarity. PepFoundry reports L/D-analog embedding distances of 3.623–17.585 (mean 9.556 ± 3.097) as evidence that stereochemistry is retained — a representation-specific scale with no cross-paper comparability.

---

## SECTION 6 — CONVENTIONAL NUMERIC THRESHOLDS IN PEPTIDE WORK

| Threshold | Meaning as used | Measured under | Source and year | Re-derived for peptides? |
|---|---|---|---|---|
| **0.4** `similarity_threshold` | "same scaffold bucket" in the RL diversity filter | unspecified; the class that reads the config (`IdenticalMurckoScaffold`) never uses it, and compares exact non-isomeric Murcko scaffold SMILES instead | REINVENT default, copied verbatim into PepINVENT (2024/25) and PepEVOLVE (2025) config files | **No.** Small-molecule default, never re-derived; on linear peptides the Murcko scaffold is often empty |
| **0.4** `score_threshold`, **25** `bucket_size`, **0.5** `penalty` | memory admission / saturation / penalty in the same filter | n/a | same | **No** |
| **0.6** `s₀` | centre of the similarity reward — designs are pushed *towards* this similarity to the input peptide | "Tanimoto similarity between the Morgan fingerprints"; radius, nBits, chirality all unstated | PepThink-R1, arXiv:2508.14765 (2025) | **No justification given** |
| **0.3** Butina cutoff (= Tanimoto 0.7) | default peptide clustering | `GetMorganFingerprintAsBitVect(mol, 4, 1024)`, achiral | PepFuNN `clustering.py:107` (2025) | Radius 4 justified for peptides; the 0.3 cutoff is not |
| **60** (percent) | monomer "similar enough" in the substitution matrix | `GetMorganFingerprint(mol, 2, useChirality=True)` + **Dice** | PepFuNN `similarity.py:generate_matrix(threshold=60)` (2025) | Not justified |
| **> 0.9** | "same structure" when defining a permeability activity cliff | three alternatives, any one qualifying: Tanimoto of whole-molecule fingerprints; Tanimoto of Murcko scaffolds; **Levenshtein distance between SMILES strings** | MultiCycPermea, *BMC Biol.* 23, 2025 | **No.** Levenshtein on SMILES text as a peptide "sequence similarity" is a direct small-molecule-literature borrowing |
| **ECFP-12 at 50 %, MAPc-12 at 60–70 %, MAPc-8 at 60–70 %, MAPc-20 at 80 %, MMseqs2 at 80 %** | maximum train↔test similarity at which a partition still behaves monotonically | named per row; diameters, so ECFP-12 = radius 6 | Fernández-Díaz et al., *J. Cheminform.* 2025 | **Yes** — the only thresholds in this report derived on peptide data, per task |
| **0.4** CD-HIT `-c` (with `-n 2`) | redundancy removal, negatives only | sequence identity | HydrAMP (2022/2023), `dataset_preparation.ipynb` | Threshold stated, not justified |
| **50–55 %** | observed modal generated-vs-training sequence similarity, presented as evidence of novelty | `1 − editdistance/max(len)` | AMPd-Up, *Protein Sci.* 2024 | Descriptive, not a cutoff |
| **0.95** | "same cluster" for **non-canonical** peptide splitting; clusters = connected components of the similarity graph | 1024-bit **ECFP6** ("ECFP-split") | PepBenchmark, arXiv:2604.10531 (2026) | Chosen for peptides but not calibrated; sits far above the drug-like random baseline |
| **90 %** MMseqs2, **30 %** MMseqs2 | near-duplicate removal; train/test split for canonical peptides | sequence identity | PepBenchmark (2026) | 30 % is the protein-homology convention imported wholesale |
| **70 %** CD-HIT, plus "edit distance ≥ 3–5" | admission gate for generated candidates into the training pool; and the audit criterion | global sequence identity; Levenshtein in residues | anti-diabetic peptide generative pipeline, *Sci. Rep.* 2026, DOI 10.1038/s41598-026-39985-4 | Threshold stated as "operational"; a **post-hoc nearest-neighbour identity audit** is reported (max NN identity < 50 %, median ≈ 35 %, 0 % of pairs ≥ 70 %) — rare in this literature |
| **0.5** ECFP4 Tanimoto | train/holdout contamination cutoff for the whole benchmark corpus | ECFP4, 4096 bits, achiral | GuacaMol `filter_and_canonicalize(..., tanimoto_cutoff=0.5)` (2019) | Small-molecule; peptides are already excluded by the 100-char rule |
| **0.85** (MACCS/UNITY heritage), **0.3–0.42** (ECFP4/ECFP6) | "similar enough to share activity" | named | Patterson/Martin neighbourhood work; Jasial et al., *J. Cheminform.* 8:27 (2016) | Never re-derived on peptides; the 0.85 figure still circulates as a generic "same molecule" threshold |
| **≤ 100 aa** (DRAMP), **≤ 200 aa** (APD3), **≤ 50 aa** (AMPd-Up), **≤ 30 aa** (DBAASP predictors), **100-char canonical SMILES** (GuacaMol), **MW 250–350 Da** (MOSES) | inclusion/length limits that silently define the comparison population | n/a | respective papers | n/a |

**Calibration context that peptide papers generally omit.** Tanimoto is not comparable across fingerprints. Landrum's ChEMBL random-pair baselines (50 000 random pairs, molecules < 50 heavy atoms): the 90th-percentile "random" similarity is 0.167 for Morgan2 bits, 0.127 for Morgan3 bits, 0.226/0.174 for the count versions, 0.549 for MACCS, 0.344 for Atom-Pair bits (https://greglandrum.github.io/rdkit-blog/, notebook *Fingerprint Thresholds*, 2025 revision). The same author's 2025 post states the problem directly: *"the Tanimoto similarities calculated using different fingerprints can be very, very different… if you report a similarity value you should always mention the fingerprint used."* Separately, Jasial, Hu, Vogt & Bajorath (*J. Cheminform.* 8:27, 2016) show activity-relevant Tc is ~0.8 for MACCS but ~0.3 for ECFP4, and that *"generally applicable Tc threshold values as an indicator of activity similarity do not exist"*. None of these baselines were computed on peptide-sized molecules, and no peptide paper in this survey recomputes them.

---

## WHAT THE FIELD ACTUALLY DOES

**The modal choice, chemistry side:** Tanimoto on a folded Morgan/ECFP bit vector, radius 2 or 3, 1024 or 2048 bits, **with `useChirality` left at RDKit's `False`**, computed on a SMILES string obtained by expanding a HELM/CHUCKLES/BILN peptide. Diversity = 1 − mean pairwise Tanimoto; SNN = max Tanimoto to a reference set; novelty = exact set difference on canonical SMILES strings. HELM-GPT, PepTune and PepFuNN all land here; MOSES and GuacaMol supply the template. Radius 4 appears where someone thought about peptides specifically (PepINVENT's monitoring code, PepFuNN, which gives the only explicit reason — repetitive backbone fragments).

**The modal choice, sequence side:** exact-duplicate removal plus CD-HIT, 0.9 on positives and 0.4 on negatives; novelty against a database reported as max normalised edit-distance similarity or as a raw Levenshtein distribution. MMseqs2 is displacing CD-HIT in 2025–2026 work.

**Novelty is, almost everywhere, string identity.** Six of the seven generative papers in Section 1 define novelty as a set difference over SMILES or sequence strings. Only AMPd-Up and HELM-GPT compute anything graded against the reference corpus, and HELM-GPT's graded quantity (SNN) is reported as a model-quality statistic, not as a leakage check.

**How often is the choice justified?** Of the ~20 primary sources read here, five give a reason for the measure or the number: PepFuNN (radius 4, for repetitive peptide fragments), MAP4/MAP4C/MXFP (built and benchmarked for peptides), Fernández-Díaz et al. (thresholds derived per task on peptide data), AMPSphere (reduced alphabet and recalibrated Karlin–Altschul constants for short sequences), GLAMOUR (BLOSUM-free monomer substitution cost so non-natural monomers are representable). Everyone else states a parameter without a reason, or inherits one silently from REINVENT, MOSES, GuacaMol or CD-HIT's defaults. **No generative peptide-design paper in Section 1 justifies its fingerprint parameters, and none states its chirality setting except PepINVENT.**

**The two measures that cover non-canonical peptides and were validated on them — MAP4C and monomer-graph GED — are used by nobody in Section 1.**

## GAPS

- **Chirality is off by default and nobody says so.** Six of seven generative papers use or ship an achiral Morgan fingerprint while generating D-amino acids, N-methylations and epimers; PepTune's own training set is 10 % dextro-chiral α-carbons. Under the stated metric, an all-D peptide and its all-L enantiomer are the same molecule. MAP4C (2024) shows ECFP4C/APC resolve <10 % of a 2048-stereoisomer peptide set, so turning the flag on is not a fix either.
- **No peptide-calibrated noise floor exists.** Every reported Tanimoto is uninterpretable without a random-pair baseline on peptides; the only published baselines (Landrum, <50 heavy atoms, ChEMBL) are for drug-like molecules, and the ACS Chem. Biol. 2023 result says the shared backbone compresses the dynamic range in exactly the regime peptides occupy.
- **The MOSES suite is partly undefined on peptides, not merely untuned.** `Scaf` returns NaN (`min_rings=2`), `Filters` returns 0 (≥8-atom ring rule), `Frag` degenerates under BRICS. Papers report "Moses metrics" having quietly dropped the three that break, and nobody has published replacements.
- **Exact-match novelty is the standard and it is known to be gameable.** AddCarbon (Renz et al. 2019/2020) saturates it; PepThink-R1's released code does set arithmetic on *uncanonicalised* SMILES strings. A near-duplicate at Tanimoto 0.99 counts as fully novel in every Section 1 paper.
- **The redundancy machinery deletes the modality.** CD-HIT recodes J/O/U/X to one shared index and ignores stereochemistry, modification and cyclisation entirely; at least eight AMP pipelines explicitly delete non-canonical sequences, and CAMPR4 deletes stapled and circular peptides outright. Databases differ on whether a modification variant is a separate entry (DRAMP: merged; DBAASP: stored out-of-sequence; CycPeptMPDB: duplicates deliberately kept).
- **No published post-hoc leakage audit of a generative peptide model.** The only identity-controlled audits found are of *predictors* (BBB-peptide benchmarks, APD-derived AMP scoring, the anti-diabetic *Sci. Rep.* 2026 pipeline), not of generators, and none uses a chirality-aware or monomer-level measure.
