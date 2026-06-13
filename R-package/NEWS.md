# spatiocompo 0.1.0

- Initial public release.
- Bayesian spatial compositional data modeling with a Gaussian Markov Random
  Field (GMRF) prior and Dirichlet likelihood.
- MALA sampler with Fisher information preconditioning (`mala_step()`,
  `mcmc_sampling()`).
- Gibbs sampling for spatial dependence parameters
  (`gibbs_kappa_rho()`).
- Covariate support with horseshoe prior for regression coefficients.
- Kronecker-structured sparse precision matrix operations
  (`create_Q()`, `Q_rhoxQ()`, `log_det_Q()`) for scalability on large grids.
- Simulation tools: `simulate_compo_data()`, `simulate_landcover()`.
- European land cover dataset (`europe_landcover`) and convenience wrapper
  `run_europe_model()`.
- Diagnostics and visualization: `plot_trace()`, `plot_field()`,
  `plot_europe_map()`, `plot_results()`, `summary_mcmc()`.
- Convergence and model checking: `check_convergence()`, `ess()`,
  `posterior_predictive()`.
- Vignettes: introduction, simulation study, land cover example, model
  specification.
