#' Plot MCMC Results
#'
#' Creates diagnostic plots for MCMC output from [mcmc_sampling()].
#'
#' @param mcmc Output from [mcmc_sampling()].
#' @param true_params Optional list of true parameter values for comparison
#'   (e.g., from simulation). Should have components alpha, kappa, rho.
#' @param burnin Number of initial samples to discard. Default is iter/4.
#' @param which Which plots to create:
#'   \itemize{
#'     \item "trace": Trace plots for main parameters
#'     \item "density": Posterior density plots
#'     \item "acf": Autocorrelation plots
#'     \item "all": All of the above
#'   }
#'
#' @return Invisibly returns NULL. Called for side effect of plotting.
#'
#' @details
#' The function creates a multi-panel plot showing:
#' \itemize{
#'   \item Trace plots showing the sampled values over iterations
#'   \item Posterior density estimates with optional true values marked
#'   \item Autocorrelation to assess mixing
#' }
#'
#' If true_params is provided, vertical lines are added to density plots
#' and horizontal lines to trace plots showing the true values.
#'
#' @seealso [plot_trace()] for trace plots only, [mcmc_sampling()]
#'
#' @examples
#' \dontrun{
#' # After running MCMC
#' sim <- simulate_compo_data(n_side = 10, seed = 123)
#' priors <- create_priors(sim$G, d = 2)
#' theta0 <- find_initial_values(sim$y, sim$G)
#' fit <- mcmc_sampling(theta0, priors, sim$y, iter = 2000)
#'
#' # Plot diagnostics
#' plot_results(fit, true_params = sim$params)
#' }
#'
#' @export
plot_results <- function(mcmc, true_params = NULL, burnin = NULL,
                         which = c("all", "trace", "density", "acf")) {

  which <- match.arg(which)

  iter <- length(mcmc$alpha)

  if (is.null(burnin)) {
    burnin <- floor(iter / 4)
  }

  post_idx <- (burnin + 1):iter

  # Extract samples
  alpha <- mcmc$alpha[post_idx]
  kappa <- mcmc$kappa[post_idx]

  # Get rho diagonal elements
  d <- sqrt(nrow(mcmc$rho))
  rho_diag <- mcmc$rho[seq(1, d^2, by = d + 1), post_idx]

  if (which %in% c("all", "trace")) {
    plot_trace(mcmc, true_params, burnin)
  }

  if (which %in% c("all", "density")) {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    on.exit(par(old_par), add = TRUE)

    # Alpha density
    plot(density(alpha), main = expression(alpha), xlab = "")
    if (!is.null(true_params)) {
      abline(v = true_params$alpha, col = "red", lwd = 2)
    }

    # Kappa density
    plot(density(kappa), main = expression(kappa), xlab = "")
    if (!is.null(true_params)) {
      abline(v = true_params$kappa, col = "red", lwd = 2)
    }

    # Rho diagonals
    for (i in seq_len(min(2, d))) {
      plot(density(rho_diag[i, ]),
           main = bquote(rho[.(i) * "," * .(i)]), xlab = "")
      if (!is.null(true_params)) {
        abline(v = true_params$rho[i, i], col = "red", lwd = 2)
      }
    }
  }

  if (which %in% c("all", "acf")) {
    old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
    on.exit(par(old_par), add = TRUE)

    acf(alpha, main = expression(alpha), lag.max = 50)
    acf(kappa, main = expression(kappa), lag.max = 50)

    for (i in seq_len(min(2, d))) {
      acf(rho_diag[i, ], main = bquote(rho[.(i) * "," * .(i)]), lag.max = 50)
    }
  }

  invisible(NULL)
}


#' Plot MCMC Trace
#'
#' Creates trace plots for monitoring MCMC convergence.
#'
#' @param mcmc Output from [mcmc_sampling()].
#' @param true_params Optional list of true parameter values.
#' @param burnin Number of burnin samples (shown in different color).
#'
#' @return Invisibly returns NULL.
#'
#' @details
#' Trace plots show the sampled parameter values over iterations. Good mixing
#' appears as a "fuzzy caterpillar" pattern with no trends or long excursions.
#'
#' The burnin period is shown in gray, post-burnin in black. True values
#' (if provided) are shown as horizontal red lines.
#'
#' @seealso [plot_results()] for comprehensive diagnostics.
#'
#' @examples
#' \dontrun{
#' fit <- mcmc_sampling(theta0, priors, y, iter = 2000)
#' plot_trace(fit)
#' }
#'
#' @export
plot_trace <- function(mcmc, true_params = NULL, burnin = NULL) {

  iter <- length(mcmc$alpha)

  if (is.null(burnin)) {
    burnin <- floor(iter / 4)
  }

  d <- sqrt(nrow(mcmc$rho))

  # Set up plotting
  n_plots <- 2 + d  # alpha, kappa, rho diagonals
  old_par <- par(mfrow = c(min(4, n_plots), 1), mar = c(2, 4, 2, 1))
  on.exit(par(old_par), add = TRUE)

  # Colors: gray for burnin, black for post-burnin
  cols <- ifelse(seq_len(iter) <= burnin, "gray70", "black")

  # Alpha trace
  plot(mcmc$alpha, type = "l", col = cols[1],
       ylab = expression(alpha), xlab = "", main = "")
  for (i in 2:iter) {
    segments(i - 1, mcmc$alpha[i - 1], i, mcmc$alpha[i], col = cols[i])
  }
  if (!is.null(true_params)) {
    abline(h = true_params$alpha, col = "red", lwd = 2)
  }
  abline(v = burnin, lty = 2, col = "blue")

  # Kappa trace
  plot(mcmc$kappa, type = "n",
       ylab = expression(kappa), xlab = "", main = "")
  for (i in 2:iter) {
    segments(i - 1, mcmc$kappa[i - 1], i, mcmc$kappa[i], col = cols[i])
  }
  if (!is.null(true_params)) {
    abline(h = true_params$kappa, col = "red", lwd = 2)
  }
  abline(v = burnin, lty = 2, col = "blue")

  # Rho diagonals
  rho_idx <- seq(1, d^2, by = d + 1)
  for (j in seq_len(min(2, d))) {
    rho_j <- mcmc$rho[rho_idx[j], ]
    plot(rho_j, type = "n",
         ylab = bquote(rho[.(j) * "," * .(j)]), xlab = "", main = "")
    for (i in 2:iter) {
      segments(i - 1, rho_j[i - 1], i, rho_j[i], col = cols[i])
    }
    if (!is.null(true_params)) {
      abline(h = true_params$rho[j, j], col = "red", lwd = 2)
    }
    abline(v = burnin, lty = 2, col = "blue")
  }

  invisible(NULL)
}


#' Plot Spatial Field
#'
#' Visualizes a spatial field on a regular grid using base R graphics.
#'
#' @param z Vector or matrix of values to plot. If matrix, columns are
#'   different components.
#' @param n_side Grid side length (total grid is n_side x n_side). Used for
#'   square grids. Ignored if \code{grid_size} is provided.
#' @param grid_size Integer vector c(n_lon, n_lat) for rectangular grids.
#' @param lon,lat Vectors of unique longitude and latitude values for axis
#'   labels. If NULL, integer indices are used.
#' @param component Which component to plot if z is a matrix. Default is 1.
#' @param main Plot title.
#' @param ... Additional arguments passed to image().
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' sim <- simulate_compo_data(n_side = 15, seed = 42)
#' plot_field(sim$z_true[, 1], n_side = 15, main = "Component 1")
#'
#' @export
plot_field <- function(z, n_side = NULL, grid_size = NULL, lon = NULL,
                       lat = NULL, component = 1, main = "", ...) {

  if (is.matrix(z)) {
    z <- z[, component]
  }

  if (!is.null(grid_size)) {
    n_lon <- grid_size[1]
    n_lat <- grid_size[2]
  } else if (!is.null(n_side)) {
    n_lon <- n_side
    n_lat <- n_side
  } else {
    stop("Provide either n_side (square grid) or grid_size (rectangular grid)")
  }

  z_mat <- matrix(z, nrow = n_lon, ncol = n_lat)

  x_ax <- if (!is.null(lon)) sort(unique(lon)) else seq_len(n_lon)
  y_ax <- if (!is.null(lat)) sort(unique(lat)) else seq_len(n_lat)

  xlab <- if (!is.null(lon)) "Longitude" else "x"
  ylab <- if (!is.null(lat)) "Latitude" else "y"

  image(x_ax, y_ax, z_mat,
        col = hcl.colors(100, "viridis"),
        xlab = xlab, ylab = ylab, main = main,
        asp = 1, ...)

  invisible(NULL)
}


#' Plot Europe Model Results
#'
#' Maps the posterior mean compositions from [run_europe_model()] on the
#' European grid with country borders and ocean masked.
#'
#' @param result Output from [run_europe_model()].
#' @param burnin Number of initial samples to discard. Default is iter/4.
#' @param components Which components to plot: "all" (default) or an integer
#'   vector (e.g., 1:2).
#'
#' @return Invisibly returns the posterior mean compositions (N x D matrix).
#'
#' @examples
#' \dontrun{
#' result <- run_europe_model(time = "AD1900", model = "intercept", iter = 1000)
#' plot_europe_map(result)
#' }
#'
#' @seealso [run_europe_model()], [plot_field()]
#'
#' @export
plot_europe_map <- function(result, burnin = NULL, components = "all") {

  fit  <- result$fit
  data <- result$data
  grid <- data$grid

  iter <- length(fit$alpha)
  if (is.null(burnin)) burnin <- floor(iter / 4)
  post_idx <- (burnin + 1):iter

  N <- grid$N
  d <- data$d
  D <- data$D

  # Posterior mean of x (latent field in ALR space)
  # fit$x is (N*d) x (n_saved); thin index maps to post_idx
  n_x <- ncol(fit$x)
  x_burnin <- max(1, floor(n_x * burnin / iter))
  x_post <- fit$x[, (x_burnin + 1):n_x, drop = FALSE]
  x_mean <- rowMeans(x_post)

  # Reshape to N x d and add covariate contribution if present
  x_mat <- matrix(x_mean, nrow = N, ncol = d)

  if (!is.null(fit$beta)) {
    beta_post <- fit$beta[, post_idx, drop = FALSE]
    beta_mean <- rowMeans(beta_post)
    p <- data$p
    beta_mat <- matrix(beta_mean, nrow = p, ncol = d)

    # Grid covariates (same as used in model)
    # For intercept model: B_grid is just a column of 1s
    B_grid_single <- matrix(1, nrow = N, ncol = 1)

    if (data$model == "elevation") {
      # The bundled grid is a land-masked subset; align its elevation to the
      # model's full rectangular grid (length N) by matching (lon, lat) so it
      # conforms with the latent field x_mat.
      gd <- europe_landcover$grid[[data$time]]
      ge <- gd$elevation
      ge[ge < 0] <- 0
      grid_elev <- rep(0, N)
      idx <- match(paste(gd$lon, gd$lat), paste(grid$lon, grid$lat))
      grid_elev[idx[!is.na(idx)]] <- ge[!is.na(idx)]
      # Same standardization as the observation covariates
      obs_elev <- europe_landcover$observations[[data$time]]$elevation
      obs_elev[obs_elev < 0] <- 0
      grid_elev_std <- (grid_elev - mean(obs_elev)) / sd(obs_elev)
      B_grid_single <- cbind(1, grid_elev_std)

    } else if (data$model == "LPJ_KK10_ESM") {
      # LPJ grid covariates: intercept + alr(LPJ-adjusted composition) +
      # standardized elevation, aligned to the model grid by (lon, lat).
      gd  <- europe_landcover$grid[[data$time]]
      adj <- cbind(gd$adj_coniferous, gd$adj_broadleaved, gd$adj_unforested)
      adj <- adj / rowSums(adj)
      adj <- no_zero_one(adj, eps = 1e-5)
      alr_adj <- alr(adj)
      ge <- gd$elevation; ge[ge < 0] <- 0
      obs_elev <- europe_landcover$observations[[data$time]]$elevation
      obs_elev[obs_elev < 0] <- 0
      elev_std <- (ge - mean(obs_elev)) / sd(obs_elev)
      idx <- match(paste(gd$lon, gd$lat), paste(grid$lon, grid$lat))
      ok  <- !is.na(idx)
      B_grid_single <- matrix(0, nrow = N, ncol = p)
      B_grid_single[idx[ok], ] <- cbind(1, alr_adj[ok, , drop = FALSE], elev_std[ok])
    }

    eta <- x_mat + B_grid_single %*% beta_mat
  } else {
    eta <- x_mat
  }

  # Convert from ALR to compositions
  z_mean <- invalr(eta)

  comp_names <- data$data$components
  if (is.null(comp_names)) {
    comp_names <- c("Coniferous", "Broadleaved", "Unforested")
  }

  if (identical(components, "all")) {
    components <- seq_len(D)
  }

  lon_unique <- sort(unique(grid$lon))
  lat_unique <- sort(unique(grid$lat))

  # Mask ocean cells using maps package
  has_maps <- requireNamespace("maps", quietly = TRUE)
  ocean <- rep(FALSE, N)
  if (has_maps) {
    land <- maps::map.where("world", grid$lon, grid$lat)
    ocean <- is.na(land)
  }

  old_par <- par(mfrow = c(1, length(components)), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  # Tight axis limits = data extent (half a cell beyond the outer cell
  # centres). Using xaxs/yaxs = "i" and no fixed aspect keeps the latitude
  # axis pinned to the data instead of being expanded to fill the panel.
  lon_range <- range(lon_unique) + c(-0.5, 0.5)
  lat_range <- range(lat_unique) + c(-0.5, 0.5)

  for (k in components) {
    z_k <- z_mean[, k]
    z_k[ocean] <- NA

    z_mat <- matrix(z_k, nrow = grid$size[1], ncol = grid$size[2])

    # Light blue background (ocean)
    plot(NA, xlim = lon_range, ylim = lat_range, xaxs = "i", yaxs = "i",
         xlab = "Longitude", ylab = "Latitude", main = comp_names[k])
    rect(lon_range[1], lat_range[1], lon_range[2], lat_range[2],
         col = "lightblue", border = NA)

    # Data on land (ocean cells are NA = transparent)
    image(lon_unique, lat_unique, z_mat,
          col = hcl.colors(100, "viridis"),
          add = TRUE)
  }

  invisible(z_mean)
}


#' Plot Observed Compositions at the Sites
#'
#' Maps the observed (pollen-based) compositions from the bundled
#' \code{europe_landcover} data at the observation locations, one panel per
#' component. A site-level counterpart to [plot_europe_map()].
#'
#' @param time Time period: "AD1900", "AD1725", "AD1425", "BC1000", or "BC4000".
#' @param components Which components to plot: "all" (default) or an integer
#'   vector (e.g., 1:2).
#'
#' @return Invisibly returns the observed compositions (n x D matrix).
#'
#' @examples
#' \dontrun{
#' plot_europe_obs("AD1900")
#' }
#'
#' @seealso [plot_europe_map()]
#'
#' @export
plot_europe_obs <- function(time = "AD1900", components = "all") {

  utils::data("europe_landcover", package = "spatiocompo", envir = environment())
  if (!time %in% europe_landcover$time_periods) {
    stop("Time '", time, "' not available. Choose from: ",
         paste(europe_landcover$time_periods, collapse = ", "))
  }

  obs   <- europe_landcover$observations[[time]]
  comps <- europe_landcover$components
  y     <- as.matrix(obs[, comps])
  y     <- y / rowSums(y)
  D     <- ncol(y)

  comp_names <- c("Coniferous", "Broadleaved", "Unforested")
  if (length(comp_names) != D) comp_names <- comps

  if (identical(components, "all")) components <- seq_len(D)

  pal <- hcl.colors(100, "viridis")
  old_par <- par(mfrow = c(1, length(components)), mar = c(4, 4, 2, 1))
  on.exit(par(old_par))

  for (k in components) {
    cols <- pal[pmax(1L, ceiling(y[, k] * 100))]
    plot(obs$lon, obs$lat, col = cols, pch = 19, cex = 1.2,
         xlab = "Longitude", ylab = "Latitude",
         main = paste("Observed:", comp_names[k]))
  }

  invisible(y)
}
