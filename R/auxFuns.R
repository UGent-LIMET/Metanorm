#' Calculate observation weights
#'
#' @param residuals (numeric) Residuals from a (generalized additive) model fit
#'
#' @details
#' This function calculates weights based on the residuals of a model fit, here
#' used for generalized additive model fits. Residuals that are large are given
#' lower or zero weight, suggesting that these are outliers and should
#' contribute less to the model fit. The intention is to obtain fits robust to
#' outliers and overfitting.
#'
#' The weighting approach is Tukey's biweight approach.
#'
#' Residuals may include NA values if observations were not used for fitting
#' the initial model. A zero weight will be returned for these cases.
#'
#' @returns Estimated observation weights, giving no or reduced weight to
#' outlying observations.
#'
#' @import stats
#' @export
calcWeights <- function(residuals){
  est.weights <- numeric(length(residuals))
  res.mad <- mad(residuals, na.rm = TRUE)
  c.est <- 4.6851*res.mad*1.4826
  for(weighti in 1:length(est.weights)){
    if(is.na(abs(residuals[weighti]))){
      est.weights[weighti] <- 0
    } else if (abs(residuals[weighti])>c.est){
      est.weights[weighti] <- 0
    } else {
      est.weights[weighti] <- (1-(residuals[weighti]/c.est)^2)^2
    }
  }
  # rarely, all (but one) weights are zero (two distinct groups), in that case,
  #   set all to one as otherwise no observations remain for model fitting
  if(sum(est.weights == 0) >= (length(est.weights) - 1)){
    est.weights <- rep(1, length(est.weights))
  }
  return(est.weights)
}

#' Generate a PCA plot of a data matrix, coloring observations by a factor
#' variable 'type'
#'
#' @param x (matrix) Numerical matrix, columns containing samples, and rows
#' containing compounds
#' @param type (character) Type of sample, e.g., QCs vs. non QCs, or different
#' batches.
#'
#' @returns A ggplot
#' @import ggplot2
#' @importFrom rlang .data
#' @export
plotPCA <- function(x, type) {
  pca.result <- prcomp(t(x), scale = TRUE)
  pca.data <- data.frame("PC1" = pca.result$x[, 1],
                         "PC2" = pca.result$x[, 2],
                         "Type" = type)
  pca.data$Type <- as.factor(pca.data$Type)
  plot <- ggplot(pca.data, aes(x = .data$PC1, y = .data$PC2,
                               color = .data$Type,
                               shape = .data$Type)) +
                 geom_point(size = 3) +
                 labs(title = element_blank(), x = "PC1", y = "PC2") +
                 theme_minimal() +
                 scale_color_grey() +
                 scale_shape_manual(values = seq(1:length(levels(pca.data$Type))))
  return(plot)
}

#' LOESS span selection with LOOCV / MSE as criterion
#'
#' @param x (numeric) Injection order
#' @param y (numeric) Observed intensities
#' @param method (character) Either "QC-RLSC" or "rLOESS"
#'
#' @returns The optimal span for fitting a (robust) LOESS model.
#' @import stats
selectSpan <- function(x, y, method = c("QC-RLSC", "rLOESS")) {
  # Initialize a vector of spans to evaluate
  # Dunn et al. suggest lambda + 1 as start of this sequence,
  #  i.e., 3 for a second degree polynomial,
  #  however, this often gives warnings (span too small)
  #  so set this higher, in practice anything above +- 0.075
  #  appears to work fine
  stepSize <- 1 / length(x)

  # if many QCs, set the step size a bit higher, but choose a multiple
  #   of 1/length(x)
  if(stepSize < 0.025){
    stepSize <- c(1:length(x))[which((1:length(x))/length(x) > 0.025)[1]]/
      length(x)
  }
  spanVals <- seq(3 / length(x), 1, stepSize)
  spanVals <- spanVals[spanVals > 0.075]

  # initialize a vector to store the MSE for each span
  mse <- rep(NA_real_, length(spanVals))

  # loop over the spans
  for (s in 1:length(spanVals)) {
    span <- spanVals[s]
    msei <- vector("numeric", (length(x)))
    try({
      for (i in 1:length(x)) {
        # leave one out, then fit and predict the left out
        xtrain <- x[-i]
        ytrain <- y[-i]
        if(method == "QC-RLSC"){
          fit <- loess(ytrain ~ xtrain, span = span)
          ypred <- predict(fit, newdata = data.frame(xtrain = x[i]))
  
        } else if(method == "rLOESS"){
          fit <- loess(ytrain ~ xtrain, span = span, family = "symmetric")
          ypred <- predict(fit, newdata = data.frame(xtrain = x[i]))
  
        }
  
        # store the mse_i
        msei[i] <- (y[i] - ypred)^2
      }
  
      # calculate the mse
      #  remove NAs: cannot predict out of range values
      #  (i.e., lowest and highest)
      mse[s] <- mean(msei, na.rm = TRUE)
    }, silent = TRUE)
  }
  # get optimal span, and return
  finalSpan <- spanVals[which.min(mse)]
  return(finalSpan)
}
