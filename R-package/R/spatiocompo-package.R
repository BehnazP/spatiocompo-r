#' @title spatiocompo: Bayesian Spatial Compositional Data Modeling
#'
#' @description
#' Bayesian spatial compositional data modeling using Gaussian Markov
#' Random Fields (GMRF) and Dirichlet likelihood. The package implements
#' Metropolis-Adjusted Langevin Algorithm (MALA) for the latent field and
#' precision parameter, combined with Gibbs sampling for spatial range and
#' covariance parameters.
#'
#' @details
#' The main function is [mcmc_sampling()] which performs posterior inference
#' for spatial compositional data. The model assumes:
#'
#' \itemize{
#'   \item Observations \eqn{y_i} follow a Dirichlet distribution with
#'         concentration parameter \eqn{\alpha \cdot z_i}
#'   \item The compositional mean \eqn{z_i} is the inverse ALR transform
#'         of a latent Gaussian field \eqn{x_i}

#'   \item The latent field follows a GMRF prior with Matern-like covariance
#'   \item Covariates can be included through a linear predictor
#' }
#'
#' For simulated data examples, see [simulate_compo_data()] and
#' [simulate_landcover()].
#'
#' @section Key Functions:
#' \describe{
#'   \item{Compositional transforms}{[alr()], [invalr()], [dalr()]}
#'   \item{Precision matrices}{[create_Q()], [Q_rhoxQ()], [log_det_Q()]}
#'   \item{MCMC sampling}{[mcmc_sampling()], [mala_step()], [gibbs_kappa_rho()]}
#'   \item{Simulation}{[simulate_compo_data()], [simulate_landcover()]}
#'   \item{Visualization}{[plot_results()], [plot_trace()]}
#' }
#'
#' @references
#' Pirzamanbein, B., Lindstrom, J., Poska, A., Gaillard, M.-J. (2018).
#' Modelling spatial compositional data: Reconstructions of past land
#' cover and uncertainties. \emph{Spatial Statistics} 24: 14--31.
#' \doi{10.1016/j.spasta.2018.03.005}
#'
#' @author Behnaz Pirzamanbin, Johan Lindstrom
#'
#' @docType package
#' @name spatiocompo-package
#' @aliases spatiocompo
#'
#' @importFrom Matrix Diagonal bandSparse bdiag kronecker chol solve
#' @importFrom Matrix sparseMatrix t colSums rowSums symmpart
#' @importFrom stats rgamma rnorm runif dgamma optim rchisq acf density
#'   optimize quantile sd var
#' @importFrom grDevices hcl.colors
#' @importFrom graphics abline image par rect segments
#' @importFrom utils data flush.console
"_PACKAGE"
