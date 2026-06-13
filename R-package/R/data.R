#' European Land Cover Dataset
#'
#' Pollen-based land cover reconstructions for Europe at 5 time periods,
#' along with LPJ-GUESS vegetation model predictions used as covariates.
#'
#' @format A list with 7 components:
#' \describe{
#'   \item{observations}{List of 5 data frames (one per time period), each with:
#'     \describe{
#'       \item{lon, lat}{Observation coordinates}
#'       \item{coniferous, broadleaved, unforested}{Pollen-based compositions}
#'       \item{elevation}{Elevation at observation location}
#'       \item{adj_coniferous, adj_broadleaved, adj_unforested}{LPJ-GUESS
#'         predictions adjusted with KK10 land use}
#'       \item{lpj_coniferous, lpj_broadleaved, lpj_unforested}{Raw LPJ-GUESS
#'         predictions}
#'     }
#'   }
#'   \item{grid}{List of 5 data frames with predictions on the full grid:
#'     \describe{
#'       \item{lon, lat}{Grid coordinates}
#'       \item{adj_coniferous, adj_broadleaved, adj_unforested}{LPJ-GUESS
#'         adjusted predictions}
#'       \item{lpj_coniferous, lpj_broadleaved, lpj_unforested}{Raw LPJ-GUESS
#'         predictions}
#'       \item{elevation}{Grid elevation}
#'     }
#'   }
#'   \item{time_periods}{Character vector: "AD1900", "AD1725", "AD1425",
#'     "BC1000", "BC4000"}
#'   \item{components}{Character vector: "coniferous", "broadleaved",
#'     "unforested"}
#'   \item{description}{Dataset description}
#'   \item{source}{Data source reference}
#'   \item{citation}{Full citation for the paper}
#' }
#'
#' @details
#' This dataset contains:
#' \itemize{
#'   \item Pollen-based land cover observations at 175-204 locations per
#'     time period
#'   \item LPJ-GUESS vegetation model predictions on a 40x27 grid (679 cells)
#'   \item Elevation data for both observations and grid
#' }
#'
#' The three land cover components (coniferous forest, broadleaved forest,
#' unforested) form a compositional dataset that sums to 1 at each location.
#'
#' @source
#' Pirzamanbein, B., Lindstrom, J., Poska, A., Gaillard, M.-J. (2018).
#' Modelling spatial compositional data: Reconstructions of past land
#' cover and uncertainties. \emph{Spatial Statistics} 24: 14--31.
#' \doi{10.1016/j.spasta.2018.03.005}
#'
#' @examples
#' # Load the dataset
#' data(europe_landcover)
#'
#' # See available time periods
#' europe_landcover$time_periods
#'
#' # Get AD1900 observations
#' obs <- europe_landcover$observations$AD1900
#' head(obs)
#'
#' @seealso [load_europe_data()], [run_europe_model()]
#'
"europe_landcover"


#' Simulated Land Cover Data
#'
#' A simulated dataset of spatial compositional land cover data on a 20x20 grid.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{y}{Matrix (400 x 3) of observed compositional data. Columns are
#'     Forest, Open, and Other land cover types.}
#'   \item{z_true}{Matrix (400 x 3) of true compositional means.}
#'   \item{x_true}{Matrix (400 x 2) of true latent field in ALR space.}
#'   \item{G}{Sparse precision matrix (400 x 400) for the GMRF.}
#'   \item{coords}{Matrix (400 x 2) of grid coordinates.}
#'   \item{params}{List of true parameter values:
#'     \describe{
#'       \item{kappa}{Spatial range parameter (0.4)}
#'       \item{alpha}{Dirichlet precision (80)}
#'       \item{rho}{2x2 component covariance matrix}
#'     }
#'   }
#' }
#'
#' @details
#' The data were simulated using `simulate_landcover()` with:
#' - Grid size: 20 x 20 = 400 cells
#' - Pattern: North-south gradient (more forest in north)
#' - Dirichlet precision: 80 (moderate observation noise)
#' - Spatial range: kappa = 0.4
#' - Cross-component correlation: -0.3 between forest and open land
#'
#' @examples
#' data(simulated_landcover)
#'
#' # View dimensions
#' dim(simulated_landcover$y)
#'
#' # Plot forest cover
#' plot_field(simulated_landcover$y[, 1], n_side = 20, main = "Forest")
#'
#' # True parameters
#' simulated_landcover$params
#'
#' @source Generated with `simulate_landcover(n_side = 20, pattern = "gradient",
#'   kappa_true = 0.4, alpha_true = 80, seed = 2024)`
"simulated_landcover"
