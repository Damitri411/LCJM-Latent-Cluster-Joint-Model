# LCJM: Latent Cluster Joint Model

Bayesian latent cluster joint models for structurally complex multivariate longitudinal and time-to-event data, implemented in R and JAGS.

The model identifies latent subgroups of subjects. Within each subgroup, the joint behaviour of several longitudinal biomarkers and a time-to-event outcome is allowed to differ. Cluster membership is estimated jointly with every other model parameter, so subjects are grouped by their full longitudinal–survival profile and not by a single outcome.

This repository contains the code accompanying:

> Kundu, D., Basu, S., Gogoi, M. P., & Das, K. (2025). A joint latent-class Bayesian model with application to ALL maintenance studies. *Journal of Applied Statistics*. [https://doi.org/10.1080/02664763.2025.2511935](https://doi.org/10.1080/02664763.2025.2511935)

---

## Model overview

For subject *i* belonging to latent cluster *G<sub>i</sub> = c*:

**Longitudinal sub-model (mixed effects)**

- A **cluster-specific** marker *y₁* has fixed effects β₁<sup>(c)</sup>, random intercept and slope **u**<sub>1i</sub><sup>(c)</sup> ~ N(0, D₁<sup>(c)</sup>), and error variance σ²₁<sup>(c)</sup>, all of which vary by cluster.
- **Global** markers *y₂, y₃* have parameters shared across clusters, with a joint 4-dimensional random effect **u**<sub>23,i</sub> ~ N(0, D₂₃).
- Mean trajectories are quadratic in time (polynomial degree is selectable). They include time-varying covariates (e.g. drug doses) and baseline covariates.

**Survival sub-model (proportional hazards)**

```
h_i(t | G_i = c) = h0_c(t) · exp( ψ1_c·m1_i(t) + ψ2_c·m2_i(t) + ψ3_c·m3_i(t) + θ_c' z_i )
```

- *m<sub>k,i</sub>(t)* is the true (error-free) longitudinal trajectory, so the association is current-value.
- The association parameters ψ<sub>c</sub> and baseline covariate effects θ<sub>c</sub> are cluster-specific.
- The baseline hazard *h0<sub>c</sub>(t)* is either a constant (simulations) or a **penalized cubic B-spline (P-spline)** with 17 basis coefficients and a second-order difference penalty (real-data analysis).
- The cumulative hazard is approximated with **15-point Gauss–Kronrod quadrature**.

**Mixture / cluster layer**

- Mixing weights have either a Dirichlet prior, π ~ Dir(α), or a **covariate-dependent multinomial logit**, π<sub>ic</sub> ∝ exp(φ<sub>c</sub>' x<sub>i</sub>).
- The marginal likelihood Σ<sub>c</sub> π<sub>c</sub> L<sub>i</sub>(c) is entered with the JAGS **ones trick**: `ones[i] ~ dbern(L_i / Const)`.
- Posterior membership probabilities π\*<sub>ic</sub> ∝ π<sub>c</sub> L<sub>i</sub>(c) are monitored, and subjects are assigned by MAP or by the posterior mode of `clus[i]`.

**Priors**

- Vague normal priors on regression and association coefficients.
- Gamma(0.001, 0.001) priors on precisions.
- Wishart or factor-type (ττ' + εI) priors on random-effect precision matrices.
- Gamma prior on the P-spline smoothing parameter.

---

## Model selection

`List_of_criterion()` computes the following from the posterior draws, for choosing the number of clusters (K = 2, 3, 4) and the polynomial degree:

- log-likelihood at the posterior mean
- DIC
- AIC
- BIC
- LPML (via CPO)
- WAIC₁ and WAIC₂
- minimum cluster size (%)
- misclassification rate (simulations only)

---

## Repository contents

| File | Purpose |
|---|---|
| `Simulation_JM.R` | Main simulation study: 200 replicates, n = 250, 10 visits each, 2 true clusters (π = 0.6/0.4). Fits K = 2, 3, 4 and reports information criteria, bias, 95% CrI length and coverage, with label-switching alignment. |
| `Simulation.R` | Single-replicate version of the simulation with extra null covariates. |
| `JM_comp_Lymph_Cluster.R` | **Main real-data model.** Lymphocyte is cluster-specific; ANC and platelets are global. Uses a P-spline baseline hazard and a staged initialisation (pilot run → warm-started 2-chain runs). Includes diagnostics. |
| `JM_comp_Lymph_Dependence_Cluster.R` | Variant where the cluster-specific lymphocyte random effects are centred on subject-level effects linked to ANC/platelets. |
| `JM_comp_Lymph_Cluster_Polynomial_Selection.R` | Fits and compares polynomial degrees 1–3 for the longitudinal trajectories. |
| `Cluster_JM_Model_Selection.R` | Computes DIC/AIC/BIC/LPML/WAIC for real-data fits with different K. |
| `Cluster_JM_viz.R` | Post-processing and figures: cluster-wise mean trajectories, model-based and Kaplan–Meier survival curves, π\* membership plot, forest plots of β, θ and ψ, and MCMC trace, density and cumulative-mean diagnostics. |
| [`Plots/`](Plots) | EPS figures from the two-cluster real-data fit: cluster-wise trajectories, model-based and Kaplan–Meier survival curves, posterior membership probabilities, and forest plots of ψ, θ and β. See [`Plots/Readme.md`](Plots/Readme.md) for a description of each figure. |

---

## Real-data application

The methods were developed for a **paediatric acute lymphoblastic leukaemia (ALL) maintenance-therapy** cohort.

**Longitudinal outcomes** (log scale, median-centred):

- lymphocyte count (WBC − ANC)
- absolute neutrophil count (ANC)
- platelet count

**Time-varying covariates:** 6-mercaptopurine (6-MP) and methotrexate (MTX) doses, median-normalised.

**Baseline covariates:** age, gender, lineage, presenting WBC, NCI risk, bulky disease, CNS status, cytogenetic risk, risk group at presentation, morphological CR1.

**Event:** relapse.

**Inclusion:** patients with at least 5 valid visits. Time is aligned to each patient's first visit and scaled to [0, 1].

> The clinical data are **not** included in this repository.

---

## Requirements

- R (≥ 4.0)
- [JAGS](https://mcmc-jags.sourceforge.io/) (≥ 4.3)
- R packages:

```r
install.packages(c("rjags", "R2jags", "jagsUI", "MASS", "dplyr", "data.table",
                   "survival", "splines", "lme4", "monomvn", "emdbook",
                   "JMbayes", "ggplot2", "gridExtra", "ggpubr", "survminer",
                   "latex2exp"))
```

---

## Getting started

1. Install JAGS and the R packages listed above.
2. Edit the file-path variable `fp` at the top of each script. Paths are currently hard-coded to local Windows directories.
3. **Simulation:** run `Simulation_JM.R`. Each replicate saves its data list (`dat_list_<i>.rds`) and accumulates the `diagnostics` (model selection) and `Basic_diagnostic` (bias / coverage) data frames.
4. **Real data:** place the prepared dataset (`data_0.rds`) in the `DATA` folder. Run `JM_comp_Lymph_Cluster.R` with `nCluster` set as required, then `Cluster_JM_Model_Selection.R` and `Cluster_JM_viz.R`.

**Computation notes**

- Each fit runs one or two chains with 20k–25k iterations.
- The cluster-wise likelihood is evaluated for every subject, cluster and quadrature node, so runs on real data can take several hours.
- Use `parallel = TRUE` in `jagsUI::jags()` to run chains on separate cores.
- Label switching is handled after fitting by aligning estimated clusters to the true labels (simulation) or by warm-start initial values that fix the label ordering (real data).

---

## Authors

Damitri Kundu, Sevantee Basu, Manash Pratim Gogoi, and Kiranmoy Das.

Code maintained by Damitri Kundu.

## Citation

If you use this code, please cite:

> Kundu, D., Basu, S., Gogoi, M. P., & Das, K. (2025). A joint latent-class Bayesian model with application to ALL maintenance studies. *Journal of Applied Statistics*. https://doi.org/10.1080/02664763.2025.2511935

```bibtex
@article{kundu2025jointlatentclass,
  title   = {A joint latent-class Bayesian model with application to {ALL} maintenance studies},
  author  = {Kundu, Damitri and Basu, Sevantee and Gogoi, Manash Pratim and Das, Kiranmoy},
  journal = {Journal of Applied Statistics},
  year    = {2025},
  doi     = {10.1080/02664763.2025.2511935}
}
```

