# HELPERS ======================================================================
#' Get Sample Code from File Names
#'
#' @param file A [`character`] vector of path names.
#' @return A [`character`] vector.
#' @noRd
get_bdx <- function(file) {
  regmatches(x = file, m = regexpr(pattern = "BDX[0-9]{5}", text = file))
}

#' Geometric Mean
#'
#' @param x A [`numeric`] vector.
#' @param trim A length-one [`numeric`] vector specifying the fraction (0 to 0.5)
#'  of observations to be trimmed from each end of `x` before the mean is
#'  computed.
#' @param na.rm A [`logical`] scalar: should `NA` values be stripped before the
#'  computation proceeds?
#' @return A [`numeric`] vector.
#' @noRd
gmean <- function(x, trim = 0, na.rm = FALSE) {
  index <- is.finite(x) & x > 0
  exp(mean(log(unclass(x)[index]), trim = trim, na.rm = na.rm))
}

# XRD ==========================================================================
## READ ------------------------------------------------------------------------
## List all diffractograms (Brucker RAW)
xrd_files <- here::here("analysis/data/raw_data") |>
  list.files(pattern = ".raw", full.names = TRUE)

## List all samples
xrd_samples <- get_bdx(xrd_files)

## Read files
xrd_data <- lapply(
  X = xrd_files,
  FUN = function(x) {
    diffractogram <- rxylib::read_xyData(x, verbose = FALSE, metaData = FALSE)
    xy <- diffractogram$dataset[[1]]$data_block
    colnames(xy) <- c("theta", get_bdx(x))
    as.data.frame(xy)
  }
)

## Check that all diffractograms share the same scale (2 theta)
theta_scale <- xrd_data[[1]]$theta
for (i in 2:length(xrd_data)) {
  if (isFALSE(all.equal(theta_scale, xrd_data[[i]])))
     warning(sprintf("The theta positions of sample %d is different.", i))
}

## Visual inspection
plot(
  x = NA, y = NA,
  xlab = "", ylab = "",
  xlim = c(3, 70),
  ylim = c(0, 120000)
)
for (i in xrd_data) {
  lines(x = i[[1]], y = i[[2]])
}

## PRE-PROCESS -----------------------------------------------------------------
## Penalized likelihood smoothing
lambda <- seq(from = 1, to = 8, length.out = 40)
lambda <- 10^lambda

xrd_clean <- vector(mode = "list", length = length(xrd_data))
for (i in seq_along(xrd_data)) {
  ## Save plot for visual inspection
  sprintf("analysis/figures/process_%s.png", xrd_samples[i]) |>
    here::here() |>
    png(width = 1024, height = 768, res = 100)

  selected <- xrd_data[[i]]
  names(selected) <- c("x", "y")

  ## Subset from 5° to 55° (2theta)
  selected <- alkahest::signal_select(selected, from = 3.5, to = 70)
  plot(selected, type = "l", ylim = c(0, max(selected$y)),
       xlab = expression(2*theta), ylab = "Count", main = xrd_samples[i])

  ## Strip ka2 (J. J. de Rooi et al. 2014)
  selected <- alkahest::ka2_strip_penalized(selected, lambda = lambda,
                                            wave = c(1.5406, 1.54443),
                                            tau = 0.5, nseg = 1)
  lines(selected, col = "red")

  ## Correct sample displacement
  ### Find closest peak to 26.64 (2 theta)
  pks <- alkahest::peaks_find(selected, SNR = 1, m = 3)
  qtz <- which.min(abs(pks$x - 26.64))
  delta <- 26.64 - pks$x[qtz]
  ### Shift
  selected <- alkahest::signal_shift(selected, lag = delta)
  cat(sprintf("Sample displacement: %g", delta), sep = "\n")

  # lines(selected, col = "green")

  ## Linearly Interpolate
  selected <- alkahest::resample_interpolate(selected, from = 3.6, to = 69.6, by = 0.02)

  ## 4S Peak Filling baseline estimation (Liland, K. H. 2015)
  baseline <- alkahest::baseline_peakfilling(selected, n = 10, m = 5, by = 10,
                                             sparse = TRUE)
  lines(baseline, col = "grey")

  ## Remove baseline
  selected$y <- selected$y - baseline$y
  lines(selected, col = "blue")

  ## Subset from 4° to 54° (2theta)
  selected <- alkahest::signal_select(selected, from = 4, to = 54)

  ## Replace negative values by 0
  selected <- alkahest::replace_negative(selected, value = 0)

  ## Rescale intensities to 0-1
  # selected <- alkahest::rescale_area(selected)
  # selected <- alkahest::rescale_range(selected, min = 0, max = 1)

  names(selected) <- names(xrd_data[[i]])
  xrd_clean[[i]] <- selected

  legend("topright", legend = c("raw", "ka2 stripped", "baseline", "corrected"),
         col = c("black", "red", "grey", "blue"), lty = 1)
  dev.off()
}

## WRITE -----------------------------------------------------------------------
## Build a data.frame
xrd_data_table <- Reduce(
  f = function(x, y) merge(x, y, by = "theta", all = TRUE, sort = FALSE),
  x = xrd_data
)
xrd_clean_table <- Reduce(
  f = function(x, y) merge(x, y, by = "theta", all = TRUE, sort = FALSE),
  x = xrd_clean
)

## Export table
utils::write.table(
  x = xrd_data_table,
  file = here::here("analysis/data/derived_data/xrd_raw.csv"),
  sep = ",", dec = ".", row.names = FALSE,
)
utils::write.table(
  x = xrd_clean_table,
  file = here::here("analysis/data/derived_data/xrd_clean.csv"),
  sep = ",", dec = ".", row.names = FALSE,
)

# WD-XRF =======================================================================
## List all files
here::here("analysis/data/raw_data") |>
  list.files(pattern = "WDXRF", full.names = TRUE) |>
  lapply(
    FUN = function(x) {
      xrf <- read.table(x, header = TRUE, sep = ",", dec = ".")
      xrf$site <- gsub("\\s*\\([^\\)]+\\)", "", xrf$site)
      xrf[, c("sample", "site", "CaO", "SiO2", "Al2O3")]
    }
  ) |>
  do.call(rbind, args = _) |>
  subset(sample %in% xrd_samples) |>
  (function(x) split(x, f = x$sample))() |>
  lapply(
    FUN = function(x) {
      if (nrow(x) == 1) return(x)
      x$sample = unique(x$sample)
      x$site = unique(x$site)
      x$CaO = gmean(x$CaO)
      x$SiO2 = gmean(x$SiO2)
      x$Al2O3 = gmean(x$Al2O3)
    }
  ) |>
  do.call(rbind, args = _) |>
  utils::write.table(
    file = here::here("analysis/data/derived_data/xrf_cas.csv"),
    sep = ",", dec = ".", row.names = FALSE,
  )
