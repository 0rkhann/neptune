# Benchmark statistics: what the literature actually says

Scope: facts and stated recommendations from the cited sources. Quoted fragments are short and marked. Where a paper states a precondition or scope limit for its own method, that limit is reported as close to verbatim as the source allows.

---

## SECTION 1 — Agarwal et al. 2021, "Deep RL at the Edge of the Statistical Precipice"

Sources: arXiv abs/2108.13264 (v4 HTML, 2022) <https://arxiv.org/html/2108.13264v4>; NeurIPS 2021 proceedings PDF <https://proceedings.neurips.cc/paper/2021/file/f514cec81cb148559cf475e7426eed5e-Paper.pdf>; Google Research blog, 2021 <https://research.google/blog/rliable-towards-reliable-evaluation-reporting-in-reinforcement-learning/>; rliable README <https://github.com/google-research/rliable/blob/master/README.md>.

### 1.1 The four recommendations and the exact reasoning for each

Their Table 1 ("Our recommendations for reliable evaluation, easily applicable with a handful of runs") pairs each current practice with a replacement:

| Desideratum | Their recommendation | Their stated reason |
|---|---|---|
| Uncertainty in aggregate performance | **Interval estimates via stratified bootstrap CIs** | Point estimates "ignore statistical uncertainty" and "hinder results reproducibility". Any estimate from a finite number of runs *is* a random variable, so "it should be treated as such" (2021/2022, §1). |
| Variability across tasks and runs | **Performance profiles (score distributions)** | Per-task tables are "overwhelming beyond a few tasks", standard deviations are "often omitted", and tables give an "incomplete picture for multimodal and heavy-tailed distributions". Profiles "show tail distribution of scores on combined runs across tasks" and let you "easily read any score percentile". |
| Aggregate metric | **Interquartile mean (IQM)** across all runs | Mean is "often dominated by performance on outlier tasks". Median "requires large number of runs to claim improvements" and is a "poor indicator of overall performance: zero scores on nearly half the tasks do not affect it". IQM is "robust to outlier scores but more statistically efficient than median". |
| Other aspects of the gain | **Probability of improvement** and **optimality gap** | Reported alongside IQM "to show other aspects of performance gains". |

Definitions as stated:

- **IQM** = "Also called 25% trimmed mean, IQM discards the bottom and top 25% of the runs and calculates the mean score of the remaining 50% runs (= ⌊NM/2⌋ for N runs each on M tasks)." It "interpolates between mean and median across runs, which are 0% and almost 50% trimmed means respectively."
- **Stratified bootstrap**: "Re-sample runs with replacement independently for each task to construct an empirical bootstrap sample with N runs each for M tasks from which we calculate a statistic and repeat this process many times." CIs are **percentile** bootstrap CIs.
- **Performance profile (run-score distribution)**: F̂_X(τ) = (1/M) Σ_m F̂_m(τ), the fraction of *runs* above threshold τ. They state this is an **unbiased** estimator of the underlying distribution, and that "an outlier run with extremely high score can change the output ... by at most 1/(MN)" — versus 1/M for the average-score (per-task-mean) distribution, which they call a **biased** estimate.
- **Optimality gap** = "the amount by which the algorithm fails to meet a minimum score of γ = 1.0". Explicit precondition: "This assumes that a score of 1.0 is a desirable target beyond which improvements are not very important" (e.g. human-level). γ may be chosen differently (their Appendix A.7).
- **Probability of improvement** = P(X>Y) = (1/M) Σ_m P(X_m > Y_m), "how likely it is for X to outperform Y on a randomly selected task", computed from the **Mann–Whitney U statistic** averaged across tasks (blog, 2021).

### 1.2 Stated preconditions and scope limits (the important part)

1. **The metric is a trimmed mean, so it inherits trimmed-mean semantics.** The paper is explicit that IQM *should* degrade under zero scores: "zero scores on nearly half of the tasks does not affect the median while IQM exhibits a severe degradation." This is sold as a feature on continuous normalized scores.
2. **Probability of improvement ignores effect size.** Verbatim: "unlike IQM and optimality gap, this metric does not account for the size of improvement." They note it "does not distinguish between two algorithms which uniformly improve on all tasks by 1% and 100%".
3. **Best aggregate metric is left open.** Verbatim: "finding the best aggregate metric is still an open question and is often dependent on underlying" (context: the score distribution / use case).
4. **Number of runs.** The headline claim is "easily applicable with 3–10 runs per task". But the coverage study (their Appendix A.5) says: percentile CIs "provide good interval estimates for as few as N=10 runs for both median and IQM scores"; **"With 3 runs, bootstrap CIs underestimate the true 95% CIs and might require a larger nominal coverage rate."** Per-task (non-aggregated) bootstrap CIs need roughly 20–30 runs per task for true 95% coverage. So: aggregate IQM CIs at 5–10 runs, not 3.
5. **Number of tasks.** No minimum M is stated anywhere. Their variance expression for the profile, σ²_X = (1/M²N) Σ_m F_m(τ)(1−F_m(τ)), decreases in M, and all their case studies use M = 6 to 55 tasks. **With few tasks and few runs, NM is small and the trim removes most of the data** (see §2).
6. **Commensurability: assumed, via per-task normalization, and not otherwise justified.** The formalism (their §2) *starts* from "a scalar, normalized score x_{m,n}", "obtained by linearly rescaling per-task scores based on two reference points" (Atari: random agent = 0, human = 1). Pooling runs across tasks into one IQM is only meaningful if that rescaling has made tasks comparable. They never test or defend the commensurability assumption; they treat the choice of reference points as a design decision and show it matters — for Procgen they explicitly **change the normalization** ("we recommend using normalization based on the estimated minimum and maximum scores on ProcGen") precisely because PPO-normalized scores are heavy-tailed and "mean scores highly dependent on performance on a small fraction of tasks". So: **their method does assume commensurable normalized scores, and they demonstrate that the conclusion changes with the normalization choice.**
7. **Does IQM require continuous scores?** Not stated anywhere in the paper — no assumption about continuity, ties, or boundedness is made explicit. The method is developed and validated entirely on continuous, unbounded-above normalized returns (Atari human-normalized, DM Control, Procgen). This silence is the single largest gap for a bounded/near-binary metric (see §2 and GAPS).
8. **Does not fix reproducibility by other means.** "the problem is not solved by fixing random seeds ... since it does not really address the question of whether an algorithm would perform well under similar conditions but with different seeds. Furthermore, fixed seeds might benefit certain algorithms more than others."
9. **They reject dichotomous significance testing**: "Nor can the problem be solved by the use of dichotomous statistical significance tests."
10. **Overlapping CIs are not a test.** Figure 2 caption: "when CIs overlap, properly accounting for uncertainty entails computing CIs for score differences."
11. **A later RL methodology paper adds a scope limit the original does not state.** Patterson, Neumann, White & White, "Empirical Design in Reinforcement Learning", arXiv 2304.01315 (2023), JMLR: "Generally, IQM is most usefully applied across collections of environments" — and they warn that on a *single* environment, trimming can delete exactly the low-probability catastrophic failures you need to see: removing "the agents whose performance do not conform to our pre-existing notions ... is not helping to create a clearer picture of our algorithm, as we are simply ignoring its shortcomings." They also caution against choosing IQM "because it provides a convenient way to reduce the number of runs" rather than because IQM is the estimand you want.

### 1.3 rliable: the actual API

Repo: <https://github.com/google-research/rliable> (README, accessed 2026).

```python
from rliable import library as rly, metrics, plot_utils
```

- `rly.get_interval_estimates(scores_dict, aggregate_func, reps=50000)` → `(aggregate_scores, aggregate_score_cis)`. `scores_dict` maps algorithm name → array of shape **(num_runs × num_tasks)** (or `(num_runs × num_tasks × num_frames)` for sample-efficiency curves). Stratified bootstrap, percentile CIs.
- `rly.create_performance_profile(scores_dict, thresholds)` → `(score_distributions, score_distributions_cis)`.
- Metrics: `metrics.aggregate_median`, `aggregate_iqm`, `aggregate_mean`, `aggregate_optimality_gap`, `metrics.probability_of_improvement` (the last takes a dict of pairs `{name: (scores_x, scores_y)}`).
- Plots: `plot_utils.plot_interval_estimates`, `plot_performance_profiles`, `plot_probability_of_improvement`, `plot_sample_efficiency_curve`.

**What it does not support:** no binomial/proportion intervals, no per-task random effects, no mixed models, no multiplicity control, no ordering/rank-stability statistic. `aggregate_optimality_gap` hardcodes the γ=1 convention. There is no API for "is this ordering stable"; the closest is probability-of-improvement with its own CI.

---

## SECTION 2 — When is IQM inappropriate?

### 2.1 Trimmed means: breakdown point and efficiency (classical results)

- The α-trimmed mean has **asymptotic breakdown point α** (CMU SDS 383C lecture notes, "Lecture 14 — Robust Statistics", <https://www.cs.cmu.edu/~psarkar/sds383c_16/lecture_scribe15.pdf>; also Wikipedia "Robust statistics", "The X% trimmed mean has a breakdown point of X%"). IQM is the 25% trimmed mean → **breakdown point 0.25**, versus 0.5 for the median and 0 for the mean.
- Breakdown point and efficiency trade off directly: "The breakdown point of each of these estimators increases as the trimming proportion increases, while the efficiency decreases" (Hubert & Vandervieren / comparison of trimming-based estimators, AStA 2008, <https://link.springer.com/article/10.1007/s10182-008-0099-5>).
- Trimmed means have "higher efficiency for mixed distributions and heavy-tailed distribution ... at the cost of lower efficiency for some other less heavily tailed distributions (such as the normal distribution)" (Wikipedia, "Truncated mean", <https://en.wikipedia.org/wiki/Trimmed_mean>).
- Standard practice is α ∈ [0.1, 0.2] (Zuev, Caltech Math 408 lecture notes, <https://www.its.caltech.edu/~zuev/teaching/2013Spring/Math408-Lecture-36.pdf>); IQM at α=0.25 sits at the robust end of the usual range.

All of this is derived for **location estimation on an unbounded, continuous-ish distribution**. None of it is a licence for bounded or discrete data.

### 2.2 IQM on a 0/1 variable: it is a saturating, biased transform of the proportion

I could find **no published source that applies or warns against IQM for proportions**. The behaviour is however elementary and worth stating because it is decisive here.

Let q = fraction of runs that succeeded, pooled over the NM runs. Sorted 0/1 data puts the zeros in [0, 1−q] and the ones in (1−q, 1]. The IQM averages the window [0.25, 0.75] of the order statistics, so

```
IQM(binary) = 0                 if q ≤ 0.25
            = (q − 0.25) / 0.5  if 0.25 < q < 0.75
            = 1                 if q ≥ 0.75
```

Verified numerically (N=100 pooled runs, exact match at every q tested):

| q (true success rate) | IQM |
|---|---|
| 0.00 / 0.10 / 0.25 | 0.000 |
| 0.26 | 0.020 |
| 0.40 | 0.300 |
| 0.50 | 0.500 |
| 0.60 | 0.700 |
| 0.74 | 0.980 |
| 0.75 / 0.90 / 1.00 | 1.000 |

Consequences for a success-rate benchmark:

1. **IQM is not an estimator of the success rate.** It is a piecewise-linear, 2×-amplified, clipped function of it. An arm at 70% and an arm at 80% success read as 0.90 and 1.00.
2. **It destroys ordering information in exactly the regime described in the task.** Any two arms both below 25% success read as 0.000; any two arms both above 75% read as 1.000. "Frequently exactly 0 or exactly 1" per cell is the worst case for this.
3. **It exaggerates differences in the middle**, so bootstrap CIs on IQM-of-binary will look tighter relative to the stretched scale and overstate separation.
4. **It is degenerate at small NM.** With 5 arms × 3 seeds × 3 sizes, a per-cell IQM over 3 runs is undefined/trim-free; rliable's own convention and common reimplementations fall back to the plain mean below 4 samples (e.g. the comment "with fewer than 4 samples the trim would empty the set, so it falls back to the plain arithmetic mean"), i.e. IQM silently becomes the mean.

**Verdict for metric (a):** IQM is the wrong aggregator for a proportion. Agarwal et al. never claimed otherwise; they simply never considered the case.

### 2.3 The standard treatment of a binary outcome aggregated over tasks and seeds

The established answer, from the meta-analysis-of-proportions literature, is a **one-step binomial model with random effects (GLMM / random-intercept logistic regression)**, not a transform-then-average two-step:

- Lin & Chu, "Meta-analysis of proportions using generalized linear mixed models", *Epidemiology* 31(5), 2020, <https://pmc.ncbi.nlm.nih.gov/articles/PMC7398826/>: two-step methods "impractically treat within-study variances as fixed, known values and require ad hoc corrections for zero counts"; "In general, GLMMs led to smaller biases and mean squared errors, and higher coverage probabilities than two-step methods." Model: event counts with a binomial likelihood plus θ_i ~ N(0, τ²) random effects on the link scale.
- Lin & Xu, "Arcsine-based transformations for meta-analysis of proportions: Pros, cons, and alternatives", *Health Science Reports* 3(3), 2020, <https://pmc.ncbi.nlm.nih.gov/articles/PMC7384291/>: arcsine transforms have "a bounded domain, implying truncations for the assumed normal distribution"; "we highly recommend the use of GLMMs or Bayesian models for synthesizing proportions".
- Schwarzer, Chemaitelly, Abu-Raddad & Rücker, *Research Synthesis Methods* 10(3), 2019, <https://pmc.ncbi.nlm.nih.gov/articles/PMC6767151/>: the Freeman–Tukey double-arcsine back-transformation can be "seriously misleading"; "Generalized linear mixed models seem to be a promising alternative."
- Trikalinos, Trow & Schmid (AHRQ, NCBI Bookshelf NBK179162), <https://www.ncbi.nlm.nih.gov/books/NBK179162/>: "Discrete likelihood methods are preferable for the meta-analyses of proportions and rates. We discourage the use of approximate methods that require continuity corrections."
- Lin, *J Gen Intern Med* 2021, <https://link.springer.com/content/pdf/10.1007/s11606-021-07098-5.pdf>: across 43,644 real Cochrane datasets methods mostly agree, **but** two-step methods diverge badly from the logit GLMM "for small total sample sizes (< 50) and crude event rates within 10–20% and 90–95%", and in fold-change for "small total event counts (< 10) and low crude event rates (< 20%)". That is precisely the 3–5-seeds, p-near-0-or-1 regime.
- **Dissent exists.** Doi et al., *BMC Med Res Methodol* 25, 2025, <https://link.springer.com/article/10.1186/s12874-025-02527-z>, argue the Freeman–Tukey transform beats logit for *extreme* proportions (coverage 94.5% vs 86% for very small p) and "seems to be the transformation of choice". This is a live disagreement; it is about *which two-step transform*, not about GLMM vs two-step.

Also relevant: with a cloglog or logit link and random intercepts for task, the GLMM estimate is interpretable as a median-across-tasks prevalence (Jackson et al., *Res Synth Methods*, 2024, <https://pubmed.ncbi.nlm.nih.gov/39668973/>), and the link choice can be selected by AIC; they warn that "misspecification of the link function can introduce bias" and that a sandwich SE "may not sufficiently avoid undercoverage due to link function misspecification".

---

## SECTION 3 — Proportion / binomial intervals at small n

Primary reference: Brown, Cai & DasGupta, "Interval Estimation for a Binomial Proportion", *Statistical Science* 16(2):101–117, 2001. PDF: <http://stat.wharton.upenn.edu/~tcai/paper/Binomial-StatSci.pdf>; JSTOR <https://www.jstor.org/stable/2676784>.

### 3.1 The intervals and their citations

| Interval | Original citation | Form |
|---|---|---|
| **Wald** (standard) | textbook / Laplace | p̂ ± z·√(p̂(1−p̂)/n) |
| **Wilson score** | Wilson (1927), *JASA* 22:209–212 | score interval; BCD call it "the Wilson interval" |
| **Clopper–Pearson "exact"** | Clopper & Pearson (1934), *Biometrika* 26:404–413 | inversion of the exact binomial test |
| **Agresti–Coull** | Agresti & Coull (1998), *The American Statistician* 52:119–126 | p̃ ± z_{0.025}·√(p̃(1−p̃)/ñ) with ñ = n+4, p̃ = (X+2)/(n+4) ("add two successes and two failures") |
| **Jeffreys (equal-tailed)** | Jeffreys prior Beta(1/2,1/2); equal-tailed posterior interval | Beta(X+½, n−X+½) quantiles |

### 3.2 What Brown, Cai & DasGupta conclude (near-verbatim)

- On Wald: "the chaotic coverage properties of the Wald interval are far more persistent than is appreciated"; "common textbook prescriptions regarding its safety are misleading and defective in several respects and cannot be trusted"; "The performance is so erratic and the qualifications given in the influential texts are so defective that **the standard interval should not be used**"; it "deserves not to be used at all."
- Concrete coverage pathologies they give: at n=100, nominal 95% coverage is 0.952 at p=0.106 but 0.911 at p=0.107; at p=0.5 coverage is 0.953 at n=17 but 0.919 at n=40; at p=0.005 coverage climbs to 0.945 at n=591 then drops to **0.792 at n=592**. "the coverage of the standard interval can be significantly lower at quite large sample sizes, and this happens in an unpredictable and rather random way."
- **The rule.** "we recommend the Wilson or the equal-tailed Jeffreys prior interval for small n (n ≤ 40). These two intervals are comparable in both absolute error and length for n ≤ 40, and we believe that either could be used, depending on taste." "For larger n (n > 40), the Wilson, the Jeffreys and the Agresti–Coull interval are all very similar, and so for such n, due to its simplest form ... the Agresti–Coull interval should be recommended." They summarise: "we recommend the Agresti–Coull interval for practical use when n ≥ 40."
- Practical note: "The Wilson interval has a closed-form formula. The Jeffreys interval does not" (they tabulate limits for n ≤ 30 and give accurate closed-form approximations).
- Even at small n, "the Agresti–Coull interval is strongly preferable to the standard one and so might be the choice where simplicity is a paramount objective."
- **Behaviour at p near 0 or 1.** They note both Wilson and Jeffreys have "disturbing downward spikes in the coverages ... very close to the two boundaries", and that **modified** versions of Wilson and Jeffreys correct these spikes. If your p̂ is routinely at the boundary, use the boundary-modified Wilson/Jeffreys, or accept Clopper–Pearson's conservatism.
- **Clopper–Pearson** was examined and deliberately **not** recommended — it is in the "additional intervals" group, rejected for the usual reason (it is guaranteed ≥ nominal coverage, hence conservative and unnecessarily wide). BCD do not recommend it for general use.

### 3.3 Applying this to 3–5 seeds

With n = 3–5 runs per cell, every interval is wide and no interval has good coverage uniformly in p. Practical implications:

- Wald is unusable: at X=0 or X=n it returns a zero-width interval [p̂, p̂]. This is the failure mode exactly matching "frequently exactly 0 or exactly 1".
- Wilson and Jeffreys never collapse at the boundary and are BCD's recommendation at n ≤ 40.
- Clopper–Pearson is the only one that guarantees ≥95% coverage; at X=0, n=5 it gives [0, 0.522], at X=5, n=5 it gives [0.478, 1]. This is honest and is the number you should quote if you want to be unattackable; BCD would call it conservative.
- **Do not report a per-cell CI as the headline.** With n=3–5 the per-cell interval spans most of [0,1] regardless of method. The informative quantity is the arm-level rate pooled over tasks and seeds with a task random effect (§2.3), or the paired per-task difference between arms (§7).
- Anthropic/Miller, "Adding Error Bars to Evals", arXiv 2411.00640 (2024), <https://arxiv.org/abs/2411.00640>, is the closest ML-side statement of the same discipline: report SEM; **cluster standard errors on the unit of randomisation** when items arrive in related groups (they find clustered SEs "up to 3X" larger than naive ones, and warn naive CIs are "likely anti-conservative (too narrow)"); and **analyse paired per-item differences** between two systems, since "paired differences represent a 'free' reduction in estimator variance".

---

## SECTION 4 — Comparing many methods over many tasks

### 4.1 Demšar, JMLR 7:1–30, 2006 — the recommended procedure and its assumptions

PDF: <https://www.jmlr.org/papers/volume7/demsar06a/demsar06a.pdf>

**Procedure.**
1. Two methods, many data sets: **Wilcoxon signed-ranks test**. Not the paired t-test.
2. k > 2 methods, N data sets: **Friedman test** on per-data-set ranks (average ranks for ties), with the **Iman–Davenport (1980)** correction, because "Friedman's χ²_F is undesirably conservative". F_F = (N−1)χ²_F / (N(k−1) − χ²_F), F-distributed with k−1 and (k−1)(N−1) df.
3. If Friedman rejects: **Nemenyi** post-hoc for all-pairs, CD = q_α·√(k(k+1)/(6N)); or **Bonferroni–Dunn** (and Holm/Hochberg/Hommel) when comparing all methods to one control, which is more powerful because it makes k−1 rather than k(k−1)/2 comparisons.
4. **CD diagram**: an axis of average ranks (best to the right), with bars joining groups whose average ranks differ by less than CD. For the control case, mark ±CD around the control's rank.

**Stated assumptions and limits.**
- *Why not averaging or the t-test:* "If the results on different data sets are not comparable, their averages are meaningless." The t-test "suffers from three weaknesses. The first is commensurability". Averages are "susceptible to outliers".
- *Wilcoxon:* "assumes commensurability of differences, but only qualitatively: greater differences still count more ... but the absolute magnitudes are ignored." And: "The Wilcoxon test assumes continuous differences d_i, therefore they should not be rounded to, say, one or two decimals since this would decrease the power of the test due to a high number of ties." **← this bites directly on a 0/1 success rate, which is maximally tied.**
- *Sign test:* "does not assume any commensurability of scores or differences" but "is much weaker than the Wilcoxon signed-ranks test".
- *Friedman:* the χ² approximation holds "when N and k are big enough (as a rule of a thumb, **N > 10 and k > 5**)". For smaller N and k, use exact critical values (Zar 1998; Sheskin 2000). **A 5-arm study with a handful of tasks is below this rule of thumb.**
- *Repeated-measures ANOVA* (the parametric alternative) assumes normality and **sphericity**, which Demšar argues is implausible here.
- *Independence:* "running the algorithms on multiple data sets naturally gives a sample of independent measurements". The whole framework rests on **independent data sets** as the sampling unit.

**How seeds / multiple runs are meant to be handled — this is explicit and restrictive.**
- "We do not record the variance of these results over multiple samples, and therefore assume nothing about the sampling scheme. The only requirement is that the compiled results provide reliable estimates of the algorithms' performance on each data set."
- "In our task, multiple resampling from each data set is used **only to assess the performance score and not its variance**. The sources of the variance are the differences in performance over (independent) data sets and not on (usually dependent) samples."
- §3.2.3, verbatim: "Could we also consider the variance, or even the results of individual folds? There are variations of the ANOVA and the Friedman test which can consider multiple observations per cell provided that the observations are independent (Zar, 1998). This is not the case here, since training data in multiple random samples overlaps. **We are not aware of any statistical test that could take this into account.**"

So Demšar's framework says: collapse your seeds to one number per (arm, task) cell, and use that cell value as the datum. It gives you **no way to report seed variance**, and explicitly disclaims one. That is the binding limitation for a 3–5-seed design where seed variance is part of the claim.

### 4.2 García & Herrera, JMLR 9:2677–2694, 2008 — post-hoc procedures

<https://www.jmlr.org/papers/v9/garcia08a.html>

- Diagnosis: "the Nemenyi test is **very conservative** and it may not find any difference in most of the experimentations."
- Recommendation: for all-pairs (n×n) comparisons, exploit the **logical relations among hypotheses**. Use **Shaffer's static** procedure (Shaffer 1986) or **Bergmann–Hommel** (1988). "As a general rule, the Bergmann–Hommel procedure is the most powerful one but it requires intensive computation in comparisons involving numerous classifiers. The second one, Shaffer's procedure, can be used instead ... in these cases."
- Also: **report adjusted p-values**, not just reject/accept at a fixed α — "valid p-values associated to each comparison useful to be compared with any level of significance without restrictions".
- Note Holm vs Hochberg/Hommel/Rom: "they have more power than Holm's procedure, but the difference between them is not very notable".

### 4.3 Benavoli, Corani & Mangili, JMLR 17, 2016 — "Should we really use post-hoc tests based on mean-ranks?"

<https://jmlr.org/papers/v17/benavoli16a.html>

- The defect: the mean-ranks (Nemenyi) comparison of A and B depends on the *other* algorithms in the pool. "the difference between A and B could be declared significant if the pool comprises algorithms C, D, E and not significant if the pool comprises algorithms F, G, H." This "can increase the type I error when comparing two equivalent algorithms and conversely decrease the power when comparing algorithms whose performance is truly different."
- Their recommendation: "we recommend to avoid the mean-ranks test for the post-hoc analysis. One should instead perform the multiple comparison using tests whose decision depend only on the two algorithms being compared, such as the **sign test or the Wilcoxon signed rank test**." Sign test is "more robust, as it only assumes the observations to be identically distributed", but low power; Wilcoxon is more powerful but "makes the additional assumption of a symmetric distribution of the differences".
- "Regardless the adopted test, the multiple comparisons should be performed adjusting the significance level to control the family-wise Type-I error" via Demšar (2006) / García & Herrera (2008) procedures.
- They point the reader to the Bayesian counterparts.

**This invalidates the CD diagram as Demšar drew it**, because the CD diagram's pairwise bars are mean-rank comparisons.

### 4.4 Benavoli, Corani, Demšar & Zaffalon, JMLR 18, 2017 — "Time for a change" (the Bayesian alternative)

PDF: <https://www.jmlr.org/papers/volume18/16-305/16-305.pdf>

**Criticisms of NHST, as they state them:**
- "it does not answer the question we ask." A 90% confidence level "is not the probability of one classifier outperforming another."
- "NHST does not estimate probabilities of hypotheses."
- "Point-wise null hypotheses are practically always false" — "by adding enough data points it is possible to claim significance".
- "The p-value does not separate between the effect size and the sample size."
- "NHST ignores magnitude and uncertainty."
- p ≤ 0.05 "leads to 'black and white thinking'".
- On multiplicity specifically: "the correction factor depends on the way the analyst intends to conduct the comparison ... two analysts can draw different conclusions from the same data because of the variety of comparisons that they made."

**What they propose instead:**
1. **rope** (region of practical equivalence, after Kruschke & Liddell): declare an interval of differences that are practically equivalent. "In classification 1% seems to be a reasonable choice. However, in other domains a different value could" be appropriate — i.e. **the rope is domain-specific and the analyst must justify it**. Report three posterior probabilities: P(A ≻ B), P(A ≈ B) (mass inside the rope), P(B ≻ A). Decision rule used in the paper: declare A ≻ B if P(A ≻ B) > 0.95.
2. **Bayesian correlated t-test** — one data set, cross-validation; corrects for the fold correlation ρ = n_test/n_train.
3. **Bayesian signed-rank test** — multiple data sets, uses only the per-data-set mean differences x̄_i.
4. **Bayesian hierarchical correlated t-test** — multiple data sets *and* multiple runs:
   - x_i ~ MVN(1μ_i, Σ_i) with diag σ_i² and off-diagonal ρσ_i² (models the correlation among runs/folds within a data set);
   - μ_1…μ_q ~ t(μ_0, σ_0, ν) — a **Student** high-level distribution, chosen so the model can "robustly deal with data sets whose μ_i's are far away from the others", and because "the heavy tails of the Student make more cautious the conclusions drawn by the model";
   - σ_1…σ_q ~ unif(0, σ̄), with σ̄ = 1000·s̄ following Gelman (2006).
   - Prior on μ_0: uniform on [−1, 1] — **"This choice works for all the measures bounded within ±1, such as accuracy, AUC, precision and recall."** (A precondition: the metric must be bounded in ±1.)
   - Benefit: joint estimation shrinks the x̄_i, so "the hierarchical model thus estimates the μ_i's more accurately than the x̄_i's adopted by the other tests".
   - Implemented in **Stan**; "Computations in hierarchical models are obtained by Markov-Chain Monte Carlo sampling"; it is "slower than the Bayesian signed rank".
   - Hyper-prior sensitivity: they ran a sensitivity analysis and report "inferences and decisions are stable w.r.t. the choice of the hyper-parameters".
5. **Which of the two for multiple data sets?** "In our opinion, the hierarchical model is preferable because it takes as inputs the m runs of the k-fold cross-validation results for each dataset and so it makes inference ... by exploiting all available information: the sample mean (x̄_i), the variability of the data (sample standard deviation σ̂_i) and the correlation due to the overlapping training set (ρ). Conversely, the Bayesian signed rank only considers x̄_i."
6. **On multiplicity in the Bayesian frame**, quoting Gelman et al. (2012): "in Bayesian analysis we usually do not have to worry about multiple comparisons. The reason is that we do not worry about Type I error, because the null hypothesis is hardly believable to be true." Gelman's own remedy is multilevel partial pooling; Benavoli et al. say that for a hierarchical model over *multiple classifiers* "This may be a direction to pursue in future research" — i.e. **they did not implement the multi-classifier hierarchical model**; in the paper they "instead mitigate false alarms through the rope."
7. Conclusion: "We discourage the use of frequentist null hypothesis significance tests (NHST) in machine learning".
   Code: <https://github.com/BayesianTestsML/tutorial/> (R and Python).

### 4.5 CD diagrams are misleading — Ismail-Fawaz et al. 2023

"An Approach to Multiple Comparison Benchmark Evaluations that is Stable Under Manipulation of the Comparate Set", arXiv 2305.11921 (2023), Ismail-Fawaz, Dempster, Tan, Herrmann, Miller, Schmidt, Berretti, Weber, Devanne, Forestier, Webb. <https://arxiv.org/abs/2305.11921>; software <https://github.com/MSD-IRIMAS/Multi_Comparison_Matrix>.

Three stated problems with the CD diagram's mean rank:
1. "rank ignores the magnitude of differences between comparates, and so can both hide large differences and exaggerate insignificant differences";
2. "rank is sensitive to which comparates are included in the comparison" (citing Benavoli 2016, Berrar 2022);
3. by a no-free-lunch argument, "the rank of all algorithms will converge as the number of datasets increases" (Berrar 2022).

They show CD diagrams are "open to both inadvertent and intentional manipulation" — e.g. adding a weaker variant of one comparate changes which pairs are declared different.

**What they recommend instead: the Multiple Comparison Matrix (MCM).** Design requirements they state: "give primacy to pairwise comparisons"; "emphasize descriptive statistics over statistical hypothesis testing"; each pairwise statistic δ(c_i,c_j) must be **invariant to the rest of the comparate set**; comparates are ordered by the **average of the performance measure itself**, not by mean rank, because that order is invariant too. Each cell reports the mean of γ(c_i,t) − γ(c_j,t) over tasks, the win/tie/loss count, and a Wilcoxon signed-rank p-value (bolded below a threshold). They note the Bayesian signed-rank probabilities could be substituted for the p-value. They are explicit that this is "not ... any new statistical method ... but rather a change in mindset".

---

## SECTION 5 — Mixed-effects / hierarchical models in benchmarking

### 5.1 Is there precedent? Yes, going back twenty years, but it is a minority practice.

**Statistics / classic ML:**
- Hothorn, Leisch, Zeileis & Hornik, "The Design and Analysis of Benchmark Experiments", *JCGS* 14(3):675–699, 2005. <https://www.zeileis.org/papers/Hothorn+Leisch+Zeileis-2005.pdf>. They give a framework in which performance is a random variable induced by the data-generating process, so "results of the experiment do not require specialized methods for the analysis: the full standard statistical tool box can be applied directly". They explicitly note prior art: "**Mixed models are applied for the comparison of algorithms across benchmark problems** (for example Lim, Loh, and Shih 2000; Kim and Loh 2003)." Their own stated caveat for the pre-existing approaches: "A basic problem common to these approaches is that the correlation between internal performance estimates, such as those calculated for each fold in k-fold cross-validation, violates the assumption of independence."
- Corani, Benavoli, Demšar, Mangili & Zaffalon, "Statistical comparison of classifiers through Bayesian hierarchical modelling", *Machine Learning*, 2017 — the hierarchical model behind §4.4.

**Modern ML / NLP — concrete published examples:**
- Hagmann, Meier & Riezler, "Towards Inferential Reproducibility of Machine Learning Research", arXiv 2302.04054 (2023, ICLR), <https://arxiv.org/abs/2302.04054>. Verbatim: "We show how to use **linear mixed effects models (LMEMs)** to analyze performance evaluation scores, and to conduct statistical inference with a **generalized likelihood ratio test (GLRT)**. This allows us to incorporate arbitrary sources of noise like meta-parameter variations into statistical significance testing, and to assess performance differences conditional on data properties. Furthermore, a **variance component analysis (VCA)** enables the analysis of the contribution of noise sources to overall variance and the computation of a reliability coefficient by the ratio of substantial to total variance." Their framing: "Instead of removing noise, we propose to incorporate several sources of variance ... into an analysis of significance and reliability".
- Benito-Santos, Ghajari & Fresno, "Robust Estimation of Population-Level Effects in Repeated-Measures NLP Experimental Designs", ACL 2025 Long Papers, pp. 33076–33089, <https://aclanthology.org/2025.acl-long.1586.pdf>. They "demonstrate that LMMs can uncover significant population-level effects — **even under low-resource (small-N) experimental designs** — while mitigating confounds and random noise", and position this as "a transparent blueprint for repeated-measures experimentation". Their justification is Clark's (1973) "language-as-fixed-effect fallacy": "ignoring nested structure — such as items within a dataset, annotator-specific biases, or model parameterizations — can inflate Type I error rates."
- Sanchez Carmona, Jiang & Dong, "Towards Robust Comparisons of NLP Models: A Case Study", COLING 2025, <https://aclanthology.org/2025.coling-main.332.pdf>. They fit a **cross-classified mixed-effects model (CCMEM)** — a design that is the exact analogue of arms × tasks × seeds — "to isolate the effect of both nuisance factors (such as random seeds) and datasets from the effects of the models' capabilities". Result: "after isolating nuisance factors and datasets, our results show that the difference between BioLinkBERT and MSR BiomedBERT is, actually, **7 times smaller than previously reported**." They cite the CCMEM as "widely used in the Social Sciences ..., Psychology ..., and Health".
- Related: "A Multilevel Analysis of PubMed-only BERT-based Biomedical Models", ClinicalNLP 2024, <https://aclanthology.org/2024.clinicalnlp-1.10.pdf>.

### 5.2 When is it preferred over rank tests or bootstrap CIs?

The literature's stated grounds:
- **When seed variance is part of the claim.** Demšar (2006) explicitly cannot handle multiple correlated observations per cell (§4.1). A mixed model can, by giving seed (or run) its own variance component. Benavoli et al. (2017) make the same argument for their hierarchical model over the signed-rank test: the signed-rank test "only considers x̄_i".
- **When you want to decompose where the variance lives.** VCA (Hagmann et al. 2023); Bouthillier et al. (2021) show empirically that in deep learning, weight initialisation "contribute[s] a small part of the variance ... much smaller than the variance due to perturbing the split of the data" — so seed-only error bars systematically understate uncertainty.
- **When n per cell is small.** Partial pooling shrinks noisy cell estimates toward the grand mean, which Benavoli et al. (2017) state "estimates the μ_i's more accurately than the x̄_i's", and Benito-Santos et al. (2025) claim works "even under low-resource (small-N) experimental designs".
- **When the outcome is binary** — then it is a GLMM, and the meta-analysis literature (§2.3) prefers it precisely at small n and extreme p.
- Against: a mixed model imposes distributional assumptions (normal random effects, often normal residuals) that rank tests and the bootstrap avoid; Demšar's whole argument for non-parametrics is that these assumptions are implausible across heterogeneous tasks, and he notes we "cannot" verify them with the number of data sets typically available ("the number of data sets is usually much less than 30").
- Against (specific): with **3 tasks, 3 sizes, 3–5 seeds**, a random-effects variance over tasks is estimated from very few levels and will frequently be estimated at or near zero (a known small-#groups pathology; not specific to ML). Nobody in the ML sources above reports a minimum number of groups.

### 5.3 Do ML reviewers accept it?

Partial evidence, not a clean answer:
- Acceptance exists: ACL 2025 long paper, COLING 2025 main, ClinicalNLP 2024, ICLR-track arXiv 2302.04054 all got through peer review with LMM/GLMM/CCMEM as the central method. So it is reviewable.
- It is clearly a minority. Dror, Baumer, Shlomov & Reichart, "The Hitchhiker's Guide to Testing Statistical Significance in NLP", ACL 2018, <https://aclanthology.org/P18-1128.pdf>: of 180 experimental long papers at ACL 2017, **only 63 ran any significance test**, 21 of those didn't name the test, and 6 of the named ones used the wrong test; the most common test was the plain t-test. Sadeqi Azer, Khashabi, Sabharwal & Roth (2020) annotated 439 ACL-2018 papers and found **6 papers using confidence intervals and 0 using Bayesian tools**.
- Practical reading: a mixed model will not be *rejected* for being a mixed model at an NLP/ML venue, but it is far enough off the modal reporting style that you should also show the conventional view (per-arm estimates with CIs, and a rank/ordering summary) alongside it.

---

## SECTION 6 — Multiple comparisons in a factorial ML benchmark

### 6.1 What is recommended

- **Demšar 2006:** control family-wise error. Omnibus Friedman first, then Bonferroni–Dunn/Holm/Hochberg/Hommel against a control; Nemenyi for all-pairs. Note his own framing: the usual research goal is "to test whether a newly proposed method is better than the existing ones", which is a **comparison with a control**, not all-pairs — and the control case is more powerful.
- **García & Herrera 2008:** for all-pairs, Shaffer static or Bergmann–Hommel (exploiting logical relations between hypotheses), and report adjusted p-values.
- **Benavoli et al. 2016:** whatever pairwise test you use, "it is necessary to control the family-wise type I error".
- **Dror, Baumer, Shlomov & Reichart (TACL 2017), "Replicability Analysis for NLP"**, <https://aclanthology.org/anthology-files/pdf/Q/Q17/Q17-1033.pdf>: when the question is "on how many datasets does A beat B, and which ones", per-dataset p-values thresholded at α "does not guarantee to bound the probability to make at least one erroneous claim" and "is error-prone when the number of participating datasets is large". They import Benjamini–Heller partial-conjunction replicability analysis, which bounds the probability of overestimating the number of datasets with a true effect, plus a multiple-testing procedure bounding the probability of any false superiority claim.
- **Benavoli et al. 2017:** reject the whole frame; use rope + posterior probabilities, and note that Bayesian multilevel partial pooling is Gelman's answer to multiplicity, though they did not implement it across classifiers.
- **Ismail-Fawaz et al. 2023:** de-emphasise inference entirely; report a stable descriptive pairwise matrix.
- **Miller 2024 (arXiv 2411.00640):** plan the experiment — they give a sample-size formula "so that model evaluators can determine in advance the size of difference that may be reliably detected", which is the pre-registration-adjacent recommendation in the eval literature.

### 6.2 What is actually done

- ACL 2017: "out of 110 papers that used multiple datasets **only 3 corrected for multiplicity** (all using the Bonferroni correction)." TACL 2017: 4 of 19. (Dror et al., ACL 2018.)
- ACL 2018: 73 of 439 papers used p-values, **6 used CIs, 0 used Bayes factors or HDIs**. (Sadeqi Azer et al. 2020, <https://exa.ai/library/publication/hcc7ndhc99d>.)
- A 2024 survey of learning-to-rank-in-NLP papers found "69 of 108 papers (64%) that do not report statistical significance", with the paired t-test the most common when anything is reported.
- In deep RL, Agarwal et al. 2021 describe the *status quo* as point estimates with 3–5 runs and no uncertainty at all; their proposal does not include any multiplicity correction either — rliable has no multiplicity API.

**So the honest summary: the recommended practice is FWER control (Holm/Shaffer/Bergmann–Hommel) or a hierarchical model; the actual practice in ML benchmark papers is overwhelmingly *nothing*, and the current reform wave (rliable, MCM, Bayesian tests) deliberately replaces the correction question with interval estimates / descriptive pairwise statistics rather than answering it.** Pre-registration of a primary comparison is not an established norm in ML benchmarking; the nearest thing in the literature is Miller's (2024) advance sample-size planning.

---

## SECTION 7 — Practical verdict

Design under discussion: 5 arms × several tasks × 3 model sizes × 3–5 seeds; non-commensurable raw task scores; claim is about **ordering of arms and its stability**.

### (a) Success rate — a proportion, often exactly 0 or 1

| Option | Supported by | Caveat / who disagrees |
|---|---|---|
| **Per-(arm,task,size) Wilson or equal-tailed Jeffreys interval** on k/n | Brown, Cai & DasGupta 2001 — "Wilson or the equal-tailed Jeffreys prior interval for small n (n ≤ 40)" | At n=3–5 the interval spans most of [0,1]; both have coverage spikes very near p=0 and p=1, for which BCD's *modified* Wilson/Jeffreys exist |
| **Clopper–Pearson** when you need a guaranteed-conservative number | Clopper & Pearson 1934 | BCD examined it and declined to recommend it (conservative, over-wide) |
| **Agresti–Coull** once pooled n ≥ 40 (e.g. an arm pooled over tasks × sizes × seeds) | BCD 2001 — "the Agresti–Coull interval for practical use when n ≥ 40" | Only after pooling; pooling across heterogeneous tasks needs the random effect below |
| **Random-intercept logistic GLMM** (success ~ arm + size + (1\|task) + (1\|task:seed)), report marginal success probabilities per arm with CIs and arm-vs-arm log-odds contrasts | Lin & Chu 2020; Lin & Xu 2020; Schwarzer et al. 2019; Trikalinos et al.; and, in ML, Hagmann et al. 2023, Benito-Santos et al. ACL 2025, Sanchez Carmona et al. COLING 2025 | Link misspecification biases estimates and sandwich SEs may not fix undercoverage (Jackson et al. 2024); few task levels → variance components poorly identified; complete separation when a cell is all-0 or all-1 (use a penalised/Bayesian fit) |
| **Per-task paired win/tie/loss counts between arms + sign test** | Demšar 2006 ("does not assume any commensurability of scores or differences"); Benavoli et al. 2016 (prefer a test depending only on the two arms compared) | Low power; Demšar: the sign test "is much weaker than the Wilcoxon" |
| **Descriptive pairwise matrix (MCM): mean per-task difference in success rate, win/tie/loss, Wilcoxon p** | Ismail-Fawaz et al. 2023 | They themselves keep the Wilcoxon p despite reservations; Wilcoxon is weak on heavily tied data |
| **Do NOT use IQM** | Derived in §2.2 (saturates to 0 below 25% and 1 above 75%, doubles the scale in between); no source endorses IQM for proportions | Agarwal et al. 2021 neither recommend nor forbid it — they are silent on bounded/discrete scores. Disagreement here is between "rliable is the RL norm, just use it" and the statistics of binary data; the literature does not resolve it, the arithmetic does. |

Transform-based pooling of proportions (logit / arcsine / Freeman–Tukey) is a genuinely contested middle ground: Lin & Xu 2020 and Schwarzer et al. 2019 advise against Freeman–Tukey and for GLMM; Doi et al. 2025 argue the Freeman–Tukey transform is better than logit at extreme p. **Both camps agree the naive Wald/arithmetic average of per-cell proportions with normal SEs is the wrong thing.**

### (b) Regret — continuous, distance to known optimum

| Option | Supported by | Caveat |
|---|---|---|
| **Per-task normalisation to [0,1], then IQM with stratified-bootstrap percentile CIs** | Agarwal et al. 2021; rliable | Requires the normalisation to make tasks commensurable — their own Procgen example shows conclusions move with the normalisation choice. 3 runs/cell is below their own coverage finding ("With 3 runs, bootstrap CIs underestimate the true 95% CIs"); aim for ≥5, ideally 10 |
| **Performance profiles (run-score distributions) over all arms** | Agarwal et al. 2021 — "statistically unbiased, more robust to outliers, and require fewer runs" | Profiles usually intersect; they themselves say finer comparisons then need an aggregate metric |
| **Optimality gap** | Agarwal et al. 2021 | Assumes a meaningful target γ; default γ=1 "assumes that a score of 1.0 is a desirable target beyond which improvements are not very important" |
| **Probability of improvement (Mann–Whitney U, averaged over tasks), with CI** | Agarwal et al. 2021; Bouthillier et al. 2021 propose P(A>B) with threshold 0.75 as a decision criterion | "does not account for the size of improvement" (Agarwal et al.); it is a pairwise statistic so it does not by itself give a total order |
| **Linear mixed model: regret ~ arm × size + (1\|task) + (1\|task:arm) + seed nesting; GLRT + variance components** | Hagmann et al. 2023; Benito-Santos et al. 2025; Sanchez Carmona et al. 2025; precedent back to Hothorn et al. 2005 | Normality/variance assumptions; few task levels; a minority reporting style at ML venues |
| **Bayesian hierarchical correlated t-test with a rope, per arm pair** | Benavoli et al. 2017 | Pairwise only (no multi-arm version implemented); the rope must be justified ("in classification 1% seems to be a reasonable choice ... in other domains a different value could" apply); prior on μ₀ assumes the measure is bounded within ±1 |
| **Wilcoxon signed-rank on per-task differences, per arm pair, with FWER control** | Demšar 2006; Benavoli et al. 2016 | Demšar: "assumes commensurability of differences, but only qualitatively"; assumes continuous differences and symmetric differences |

### (c) The ordering claim and whether the ordering is stable

| Option | Supported by | Caveat |
|---|---|---|
| **Friedman + Iman–Davenport, then CD diagram** | Demšar 2006 | Rule of thumb **N > 10 data sets and k > 5 algorithms** is violated by this design; and three separate papers say the CD diagram's pairwise claims are not trustworthy |
| **Replace Nemenyi with pairwise Wilcoxon/sign tests + FWER control (Holm, Shaffer static, Bergmann–Hommel); report adjusted p-values** | Benavoli et al. 2016; García & Herrera 2008 | Nemenyi is "very conservative"; Bergmann–Hommel is most powerful but computationally heavy |
| **Multiple Comparison Matrix: order arms by the mean metric (not mean rank); all-pairs cells with mean difference, win/tie/loss, p** | Ismail-Fawaz et al. 2023 | Explicitly descriptive, not inferential. This is the only cited proposal that directly targets *stability of the ordering*: pairwise outcomes are "invariant to C\{c_i,c_j}" and "cannot be gamed by manipulating the set of other comparates" |
| **Report the ordering separately per model size and per task; show where it flips** | Implied by Agarwal et al.'s performance profiles and by Hagmann et al.'s "performance differences conditional on data properties" | No single source prescribes a stability statistic |
| **Bootstrap the ordering: resample seeds within (arm, task) and report P(arm j ranks 1st), P(arm j > arm k)** | Stratified bootstrap machinery from Agarwal et al. 2021; probability-of-improvement is its pairwise special case; Bouthillier et al. 2021 use P(A>B) > 0.75 as the decision rule | No cited source reports a full rank-probability table for ML benchmarks; this is an extension, not an established convention |
| **Mixed model with arm × size interaction; the "ordering is stable" claim is the interaction being small** | Hagmann et al. 2023 ("assess performance differences conditional on data properties"); Sanchez Carmona et al. 2025 | Testing an absence of interaction is not the same as establishing equivalence; needs a rope/equivalence margin |

**Where the literature disagrees, explicitly:**
1. *Is NHST appropriate at all?* Demšar 2006 and García & Herrera 2008 say yes, done correctly. Benavoli et al. 2017 say no ("We discourage the use of frequentist null hypothesis significance tests"). Agarwal et al. 2021 say no for RL ("Nor can the problem be solved by the use of dichotomous statistical significance tests"). Ismail-Fawaz et al. 2023 say de-emphasise it and report descriptives.
2. *Mean rank vs mean metric as the ordering statistic.* Demšar: mean rank (no commensurability needed). Ismail-Fawaz et al.: mean metric, because mean rank is unstable under changes to the comparate set. Both cannot be satisfied: mean rank is commensurability-free but unstable; mean metric is stable but requires commensurable normalised scores.
3. *Robust vs efficient aggregator.* Agarwal et al. prefer IQM over median for efficiency. Patterson et al. 2023 warn that the trimming can hide real catastrophic failures and that IQM is for "collections of environments", not single ones.
4. *Proportion pooling.* GLMM (Lin & Chu 2020, Schwarzer et al. 2019, Trikalinos et al.) vs Freeman–Tukey (Doi et al. 2025) at extreme p.

---

## GAPS

- **No published source treats IQM applied to a proportion.** Agarwal et al. 2021 state no continuity, boundedness or tie assumption; the saturation result in §2.2 is derived here, not cited. If this matters to a reviewer, it needs a one-paragraph derivation in the paper, not a citation.
- **No source gives a minimum number of tasks** for IQM, stratified-bootstrap CIs, or performance profiles. Agarwal et al. quantify runs (3 too few, 5–10 adequate, 10 good) but are silent on M.
- **No established statistic for "is the ordering stable".** The MCM paper (2023) is the only one that targets stability, and it does so by construction (invariance to the comparate set), not by reporting a stability measure. Rank-probability tables from a bootstrap are an obvious construction with no cited precedent in ML benchmarking.
- **No source covers a factorial arms × sizes × tasks design.** Every cited framework is arms × tasks. How to treat model size (a third fixed factor, or a stratification variable with the whole analysis repeated) is unaddressed.
- **Benavoli et al.'s hierarchical model is pairwise only**; they state the multi-classifier hierarchical model "may be a direction to pursue in future research". So there is no off-the-shelf Bayesian hierarchical multi-arm ordering analysis in this literature.
- **No guidance on minimum number of random-effect levels** for a mixed model in an ML benchmark; with 3 model sizes and a handful of tasks, variance components will often be estimated near zero, and none of the ML mixed-model papers cited report a lower bound.
