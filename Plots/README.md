# Plots

Figures from the real-data analysis in:

> Kundu, D., Basu, S., Gogoi, M. P., & Das, K. (2025). A joint latent-class Bayesian model with application to ALL maintenance studies. *Journal of Applied Statistics*. https://doi.org/10.1080/02664763.2025.2511935

All figures come from the selected **two-cluster** joint model, fitted with `JM_comp_Lymph_Cluster.R` and plotted with `Cluster_JM_viz.R`.

- Patients are assigned to clusters by the maximum posterior membership probability π\*.
- Cluster 1 is shown in **red** and Cluster 2 in **blue**.
- Time is in weeks from each patient's first maintenance visit.

The files are EPS (vector) and can be included directly in LaTeX. GitHub does not preview EPS, so download a file to view it.

---

## Cluster profiles

| File | What it shows |
|---|---|
| `Lymph_Long.eps` | Mean fitted log(lymphocyte) trajectory in each cluster over weeks 1–110. The dashed horizontal line marks a reference level. |
| `ANC_PLT_post_analysis.eps` | Mean fitted log(ANC) (left) and log(platelet) (right) trajectories, grouped by the estimated clusters. ANC and platelets share parameters across clusters, so any difference between the curves comes from the patients in each group. Dashed lines mark reference levels. |
| `Survival.eps` | Model-based relapse-free probability, 1 − P(relapse), averaged within each cluster over weeks 1–210 and smoothed with loess. |
| `KaplanMeierPlot.eps` | Kaplan–Meier relapse-free curves for the two estimated clusters, as an empirical check on the model-based curves. |
| `Posterior_prob_clus1.eps` | Each patient's posterior probability of belonging to Cluster 1, P(G<sub>i</sub> = 1 \| rest), coloured by assigned cluster. Dashed and dotted lines are at 0.50 and 0.45/0.55. The legend gives the percentage of patients in each cluster. |

## Parameter estimates

Each forest plot shows the posterior mean (point) and 95% credible interval (bar). A vertical dashed line marks zero.

| File | What it shows |
|---|---|
| `Association.eps` | Cluster-specific association parameters ψ linking the current values of lymphocyte, ANC and platelet count to the relapse hazard. |
| `Baseline_cov.eps` | Cluster-specific effects θ of baseline covariates on the relapse hazard. |
| `Beta_LYM.eps` | Cluster-specific effects of baseline covariates on the lymphocyte trajectory. |
| `Medicine_LYM.eps` | Cluster-specific effects of 6-MP and MTX doses on the lymphocyte trajectory. |
| `Beta_ANC_PLT.eps` | Effects of baseline covariates on the ANC and platelet trajectories (shared across clusters). |
| `Medicine_ANC_PLT.eps` | Effects of 6-MP and MTX doses on the ANC and platelet trajectories (shared across clusters). |

---

## Reproducing the figures

1. Fit the two-cluster model with `JM_comp_Lymph_Cluster.R` (`nCluster = 2`). This saves `JM_comp_Lymph_Cl_2.rds`.
2. Run `Cluster_JM_viz.R` with `nCluster = 2`. It rebuilds each plot from that fit.
3. Save each plot as EPS, e.g. `ggsave("Association.eps", p_asso, device = cairo_ps)`.

The clinical data are not included in this repository, so the figures cannot be regenerated without access to the original dataset.
