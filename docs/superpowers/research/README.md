# Research corpus — peptide representation benchmark

Background research gathered while designing a representation-controlled peptide
optimization benchmark. Twelve files, ~7,100 lines, every claim carrying a URL
and year.

## Provenance, and what that means for using these

**These were produced by AI research agents, not by a human reading the
literature.** They are a starting point and an audit trail, not a source of
record. Before any claim here goes into a paper:

- Re-read the cited source yourself. Several claims in this corpus were wrong on
  the first pass and corrected only when a second agent checked the primary
  source (see *Known disagreements* below).
- Treat a quoted number as a pointer to where to look, not as a verified value.
- Check that a cited paper says what the summary says it says. Venue
  attributions in particular were wrong more than once.

**Some files contain original measurements** rather than literature synthesis.
Those are marked below. They were produced by installing the software and
running it, which makes them more reliable than the citation work — but they are
single unreviewed runs, not validated experiments.

## Files

| File | What it settles | Original measurements |
|---|---|---|
| `peptide-similarity-and-holdout.md` | why ECFP fails on peptides; threshold conventions; how benchmarks define holdouts | no |
| `peptide-similarity-in-practice.md` | what the field actually uses: modal choice is achiral Morgan, novelty is string identity, MOSES metrics partly undefined on peptides | no |
| `map4-independent-evaluation.md` | independent evidence for MAP4 is thin and absent for MAP4C; count-based ECFP and topological torsion are SOTA on 126 peptide datasets | yes — MinHash error across 25 permutation seeds; `mapchiral` 0.0.7 ignores its `seed` argument; timings |
| `hard-realistic-peptide-tasks.md` | cyclic-peptide permeability is *not* size/lipophilicity-hackable; no usable solubility oracle exists; AMP classifiers are the hackable ones | yes — CycPeptMPDB v1.2 downloaded, adversarial probes run against PepINVENT's shipped XGBoost oracle |
| `peptide-docking-as-oracle.md` | docking fails validity, GPU determinism and objective validity for peptide ligands | no |
| `catalysis-oracles-and-cost.md` | xTB timings; ORD holds 1,430 organocatalysis reactions of 2.4M; no public superbase dataset | yes — `xtb` 6.7.1 installed and timed; ORD API queried directly |
| `conformer-free-descriptors.md` | conformer sensitivity ranks gap > LUMO > HOMO > IP/EA > ω,N > dipole; almost no NNP outputs frontier orbitals; phosphorus rules most out | no |
| `benchmark-statistics.md` | IQM is wrong for a proportion-valued metric; Wilson/Jeffreys for success rate; Multiple Comparison Matrix for ordering | yes — IQM-on-proportions behaviour derived and numerically verified |
| `benchmark-figure-conventions.md` | what benchmark and critique papers actually plot; PMO's parity-scatter form; variance reporting norms | no |
| `scale-sweep-protocols.md` | rankings do flip with scale, including for tokenizers; three sizes over ~10× is the credible minimum | yes — computed embedding-fraction table |
| `architecture-granularity-interaction.md` | tokenizer span ≈28 points vs architecture span ≈1.7 in chemistry; MambaByte's claim is memory, not inductive bias | no |
| `amp-oracle-dependencies.md` | BATTLE-AMP's repo ships 37,565 pre-generated shuffled negatives (MIT); Macrel returns AMP and hemolytic probability in one call; three licence traps | yes — Macrel 1.6.1 and HemoPI2 installed and timed; FASTA counts verified by download |

## Known disagreements and corrections

Places where this corpus contradicts itself or corrected an earlier claim.
Resolve these against the primary source before citing:

- **Ciepliński et al. docking vs molecular weight.** One agent reported
  `r = −0.79`; a later agent that read the paper reports **no numeric
  coefficient was published**, only the phrase "moderately strong correlation".
  Unresolved. Do not cite the number.
- **kraken conformer drift.** The "1% to >75%" range was initially applied to
  electronic descriptors. It is not: the >75% figure is for *steric* descriptors
  (octant volumes). The kraken paper states electronic properties are "generally
  less sensitive". Corrected in `conformer-free-descriptors.md`.
- **Tripp & Hernández-Lobato follow-up venue.** Cited early as ICML 2024 main
  track; it is the **AI-for-Science workshop**.
- **Renz et al. co-authors** are Van Rompaey, Wegner, Hochreiter and Klambauer.
- **The Langevin rebuttal** is J Cheminform 2022.
- **PepFoundry** is BilodeauGroup, AGPL-3.0 — not MolecularAI.
- **BILN** is Boehringer Ingelheim, not Pistoia Alliance.
- **SATURN's architecture comparison** on OB100 is RNN 10/10, transformer 5/10,
  Mamba 4/10. Mamba came third.
- **CHUCKLES is two different things.** The 1994 original (Siani, Weininger &
  Blaney) covers peptides and peptoids and is not α-restricted. What PepINVENT
  and PepEVOLVE call CHUCKLES is a different object — residue SMILES templates
  joined by `|` — and that one is α-restricted in deployment.
