#' Perform normalization using a given compound's intensities.
#'
#' @param raw (numeric) Vector of length the number of samples, containing the
#'   compounds' intensities.
#' @param order (numeric) Injection order of the samples. May be NULL if the
#'   column order is chronological and no samples are omitted. If not
#'   chronological, the argument's length should equal the number of samples.
#'   Default NULL.
#' @param keepScale (logical) Whether the original intensity scale should be
#'   retained, or whether rescaling can be performed. In the former case, fitted
#'   values are first subtracted from observed values, and the median observed
#'   value is re-added. In the latter case, normalization is achieved by dividing
#'   by the fitted value, and values are centered around one. Default TRUE.
#' @param QConly (logical) Whether only QC samples should be used for fitting
#'   the normalization model (FALSE), or whether also other samples contribute to
#'   the model fit (TRUE). Sample contributions can be modified by the parameter
#'   weights. Default FALSE.
#' @param QCcheck (logical) Whether to check for significantly different GAM
#'   fits between QCs and samples (TRUE) or not (FALSE). In the case of
#'   significantly different fits, QC samples are omitted and data normalized
#'   using the sample data only. Default FALSE.
#' @param QCcheckp (numeric) Cutoff p-value for determining significantly 
#'   different fits between QCs and samples separately, versus QCs and samples
#'   combined. Default 0.1.
#' @param changepoints (logical) Whether to check for significant changepoints
#'   within batches (TRUE) or not (FALSE). Only for tGAM/rGAM. Default FALSE.
#' @param type (character) Whether the sample is a QC (value 'QC') or not (any
#'  other value). May be NULL, if QCOnly is set to FALSE. Default NULL.
#' @param batch (character) Vector of length the number of samples, indicating
#'   for each sample to which batch it belongs. May be NULL, in which case samples
#'   are presumed to belong to a single batch. Default NULL.
#' @param batchwise (logical) Whether to fit batchwise GAMs. Batchwise fitting
#'   may result in significant speed gains for complex experiments with many
#'   batches, especially when QCcheck is set to TRUE. Note that this changes the
#'   QCcheck hypothesis test from a per-batch test (TRUE) to cross-batch tests
#'   (FALSE). Default TRUE.
#' @param weights (numeric) Vector of length the number of samples, indicating
#'   what weight should be attributed to each sample when fitting the
#'   normalization model. Can be NULL, in which case all samples will be assigned
#'   equal weight. Note: when fitting only to QC samples, the QCOnly approach
#'   will be faster. Weights are only used for the rGAM and tGAM models, and
#'   ignored otherwise. Default NULL.
#' @param model (character) One of 'rGAM', a robustified Gaussian GAM; 'tGAM',
#'   a scaled-t GAM; 'QC-RLSC', as in Dunn et al. (2011); 'QC-RSC', as in Kirwan
#'   et al. (2013); or 'rLOESS', a robustified version of QC-RLSC. Default: tGAM.
#' @param k (numeric) The maximum basis dimension to use when fitting GAMs.
#'   Lower values result in faster analysis time. Calculated in the wrapper
#'   function, as a function of two parameters.
#' @param cv (character) Type of cross-validation to use for QC-RLSC, QC-RSC and
#'   rLOESS approaches. Either 'GCV' (faster) or 'LOOCV' (as suggested in the
#'   original publications, but slower). Default: GCV.
#' @param plotdir (character) File location of plots (or NULL for no plots).
#' @param plottype (character) File type for writing plots, either "pdf" or 
#'   "png", Default: pdf.
#' @param i (numeric) Compound (row) index, used for plot file naming.
#'
#' @importFrom fANCOVA loess.as
#' @importFrom mgcv gam
#' @importFrom mgcv scat
#' @importFrom changepoint cpt.mean
#' @importFrom changepoint cpts
#' @import stats
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
#' @importFrom grDevices pdf
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom methods is
#' @importFrom callr r
#' @export
#'
#' @returns Normalized intensities for the compound.
metanormWorker <- function(raw, order, keepScale, QConly, QCcheck, QCcheckp,
                           changepoints, type, batch, batchwise, weights, model,
                           k, cv, plotdir, plottype, i){

  ##########################
  # Some data preprocessing
  ##########################

  dat <- data.frame(x = order,
                    y = raw,
                    batch = batch,
                    weights = weights,
                    type = type,
                    isQC = (type == "QC"))
  dat$isQC <- as.factor(dat$isQC)

  # check for changepoints, if needed, create sub-batches
  if(model %in% c("tGAM", "rGAM") & changepoints){
    for(batchId in 1:length(unique(dat$batch))) {
      batchName <- unique(dat$batch)[batchId]
      batchLoc <- which(dat$batch == batchName)
      datSub <- dat[batchLoc, ]
      y <- datSub$y
      # cpt.mean cannot deal with NAs, need a workaround
      notNA <- which(!is.na(y))
      if(length(notNA) <= k + 1) next

      # try, because an error is thrown if minseglen cannot be met
      chPoints <- numeric(0)
      try(
        chPoints <- changepoint::cpts(changepoint::cpt.mean(y[notNA],
                                                            method = "BinSeg",
                                                            minseglen = k + 3)),
        silent = TRUE
      )

      if(length(chPoints) != 0) {
        batchSuffix <- letters[1:(length(chPoints) + 1)]
        mapChPoints <- notNA[chPoints]
        segStart <- c(1, mapChPoints + 1)
        segStop <- c(mapChPoints, nrow(datSub))
        for(segs in 1:length(segStart)) {
          datSub$batch[segStart[segs]:segStop[segs]] <- paste0(batchName,
                                                               batchSuffix[segs])
        }
        dat[batchLoc, ] <- datSub
      }
    }
  }
  dat$batch <- as.factor(dat$batch)

  # heuristic check for (too) many NAs, if needed, adjust k
  #   could be improved, because we'll be fitting for batches, potentially
  if(sum(is.na(dat$y)) > k*0.9 & model %in% c("tGAM", "rGAM")){
    warning(paste0("Too many NAs, reduced basis dimension to k = ", k*0.9,
                   " (for compound ", i, " only!)"))
    k <- k*0.9
  }

  # remove nonQC samples to speed up analysis in case of QConly
  if(QConly){
    datfit <- dat[type == "QC",]
  } else {
    datfit <- dat
  }

  # only fit using complete cases, because mgcv::gam removes NA values, leading
  #   to issues with reweighting
  # ! note: prediction is done on "dat", so that original NAs in dat remain NAs
  #   and no length mismatches occur
  datfit <- datfit[complete.cases(datfit),]
  
  ######################################################
  # Start of normalization section, different methods:
  # rGAM, tGAM, rLOESS, QC-RLSC, QC-RSC
  ######################################################
  predVals <- numeric(length(raw))

  if(model == "rGAM"){
    if(length(unique(batch)) > 1){
      if(batchwise){
        normVal <- numeric(length(dat$y))
        for(batchid in 1:length(levels(dat$batch))){
          # fit the model on fitting data only (depends on QCOnly, see above)
          batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
          batchDat <- datfit[batchIds, ]
          # but predict on all batch data
          allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
          allBatchDat <- dat[allBatchIds, ]

          # fit/refit
          rGAM <- mgcv::gam(y ~ s(x, k = k), data = batchDat,
                            weights = batchDat$weights, method = "REML")
          est.weights <- calcWeights(rGAM$residuals)
          
          # weights may be mostly 0 for very messy data, causing the fit to fail,
          #   use a try
          try({
            rGAM <- mgcv::gam(y ~ s(x, k = k), data = batchDat,
                            method = "REML", weights = est.weights*batchDat$weights)
          }, silent = TRUE)
          
          if(QCcheck & !QConly){
            rGAM2 <- mgcv::gam(y ~ isQC + s(x, k = k, by = isQC), data = batchDat,
                               weights = batchDat$weights, method = "REML")
            est.weights <- calcWeights(rGAM2$residuals)
            try({
              rGAM2 <- mgcv::gam(y ~ isQC + s(x, k = k, by = isQC), data = batchDat,
                               method = "REML", weights = est.weights*batchDat$weights)
            }, silent = TRUE)
            # if QCs differ significantly from samples, refit using samples only
            # if NA (fitting issue), also refit using samples only
            QCp <- anova(rGAM, rGAM2, test = "Chisq")$`Pr(>Chi)`[2]
            if(!is.na(QCp)){
              if(QCp < QCcheckp){
                batchDat <- batchDat[batchDat$type != "QC",]
                new.k <- ifelse(k > nrow(batchDat), nrow(batchDat), k)
                rGAM <- mgcv::gam(y ~ s(x, k = new.k), data = batchDat,
                                  weights = batchDat$weights, method = "REML")
                est.weights <- calcWeights(rGAM$residuals)
                
                try({
                  rGAM <- mgcv::gam(y ~ s(x, k = new.k), data = batchDat,
                                    method = "REML", weights = est.weights*batchDat$weights)
                }, silent = TRUE)  
              }
            } else {
              batchDat <- batchDat[batchDat$type != "QC",]
              new.k <- ifelse(k > nrow(batchDat), nrow(batchDat), k)
              rGAM <- mgcv::gam(y ~ s(x, k = new.k), data = batchDat,
                                weights = batchDat$weights, method = "REML")
              est.weights <- calcWeights(rGAM$residuals)
              try({
                rGAM <- mgcv::gam(y ~ s(x, k = new.k), data = batchDat,
                                  method = "REML", weights = est.weights*batchDat$weights)
              }, silent = TRUE)
            }
          }

          preds <- predict(rGAM, newdata = allBatchDat)
          predVals[allBatchIds] <- preds

          if(keepScale){
            normVal[allBatchIds] <- allBatchDat$y  - preds
          } else {
            normVal[allBatchIds] <- allBatchDat$y / preds
          }

        }
        if(keepScale){
          normVal <- normVal + mean(dat$y, na.rm = TRUE)
        }
      } else {
        rGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                          weights = datfit$weights, method = "REML")
        est.weights <- calcWeights(rGAM$residuals)
        try({
          rGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                            method = "REML", weights = est.weights*datfit$weights)
        }, silent = TRUE)
        if(QCcheck & !QConly){
          rGAM2 <- mgcv::gam(y ~ batch*isQC + s(x, k = k, by = interaction(batch, isQC)), data = datfit,
                             weights = datfit$weights, method = "REML")
          est.weights <- calcWeights(rGAM2$residuals)
          try({
            rGAM2 <- mgcv::gam(y ~ batch*isQC + s(x, k = k, by = interaction(batch, isQC)), data = datfit,
                               method = "REML", weights = est.weights*datfit$weights)
          }, silent = TRUE)
          # if QCs differ significantly from samples, refit using samples only
          # if NA (fitting issue), also refit using samples only
          QCp <- anova(rGAM, rGAM2, test = "Chisq")$`Pr(>Chi)`[2]
          if(!is.na(QCp)){
            if(QCp < QCcheckp){
              datfit <- dat[dat$type != "QC",]
              rGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                                weights = datfit$weights, method = "REML")
              est.weights <- calcWeights(rGAM$residuals)
              try({
                rGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                                  method = "REML", weights = est.weights*datfit$weights)
              }, silent = TRUE)
            }
          } else {
            datfit <- dat[dat$type != "QC",]
            rGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                              weights = datfit$weights, method = "REML")
            est.weights <- calcWeights(rGAM$residuals)
            try({
            rGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                              method = "REML", weights = est.weights*datfit$weights)
            }, silent = TRUE)
          }
        }
        predVals <- predict(rGAM, newdata = dat)
        if(keepScale){
          normVal <- dat$y  - predVals + mean(dat$y, na.rm = TRUE)
        } else {
          normVal <- dat$y/predVals
        }
      }
    } else {
      rGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                        weights = datfit$weights, method = "REML")
      est.weights <- calcWeights(rGAM$residuals)
      try({
        rGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                          method = "REML", weights = est.weights*datfit$weights)
      }, silent = TRUE)
      if(QCcheck & !QConly){
        rGAM2 <- mgcv::gam(y ~ isQC + s(x, k = k, by = isQC), data = datfit,
                           weights = datfit$weights, method = "REML")
        est.weights <- calcWeights(rGAM2$residuals)
        try({
          rGAM2 <- mgcv::gam(y ~ isQC + s(x, k = k, by = isQC), data = datfit,
                             method = "REML", weights = est.weights*datfit$weights)
        }, silent = TRUE)
        # if QCs differ significantly from samples, refit using samples only
        # if NA (fitting issue), also refit using samples only
        QCp <- anova(rGAM, rGAM2, test = "Chisq")$`Pr(>Chi)`[2]
        if(!is.na(QCp)){
          if(QCp < QCcheckp){
            datfit <- dat[dat$type != "QC",]
            rGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                              weights = datfit$weights, method = "REML")
            est.weights <- calcWeights(rGAM$residuals)
            try({
              rGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                                method = "REML", weights = est.weights*datfit$weights)
            }, silent = TRUE)
          }
        } else {
          datfit <- dat[dat$type != "QC",]
          rGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                            weights = datfit$weights, method = "REML")
          est.weights <- calcWeights(rGAM$residuals)
          try({
            rGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                              method = "REML", weights = est.weights*datfit$weights)
          }, silent = TRUE)
        }
      }
      predVals <- predict(rGAM, newdata = dat)
      if(keepScale){
        normVal <- dat$y  - predVals + mean(dat$y, na.rm = TRUE)
      } else {
        normVal <- dat$y/predVals
      }
    }

  } else if(model == "tGAM"){
    # multiple batches
    if(length(unique(batch)) > 1){
      if(batchwise){
        normVal <- numeric(length(dat$y))
        for(batchid in 1:length(levels(dat$batch))){
          # fit the model on fitting data only (depends on QCOnly, see above)
          batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
          batchDat <- datfit[batchIds, ]
          # but predict on all batch data
          allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
          allBatchDat <- dat[allBatchIds, ]

          # fit/refit
          tGAM <- mgcv::gam(y ~ s(x, k = k), data = batchDat,
                            method = "REML", family=mgcv::scat(),
                            weights = batchDat$weights)
          if(QCcheck & !QConly){
            # complex model rarely crashes R when few QCs/samples in a batch
            #  circumvent this by calling an external R process within a try
            #  if fails, fit to samples only to salvage
            tGAMw <- mgcv::gam(y ~ s(x, k = k), data = batchDat,
                              method = "REML", family=mgcv::scat())
            tGAM2 <- NULL
            tGAM2 <- tryCatch({
              callr::r(function(batchDat, k) {
                tryCatch({
                  mgcv::gam(y ~ isQC + s(x, k = k, by = isQC),
                            data = batchDat,
                            method = "REML",
                            family = mgcv::scat())
                }, error = function(e) NULL)
              }, args = list(batchDat = batchDat, k = k))
            }, error = function(e) NULL)
            if(is(tGAM2)[1] == "gam"){
              QCp <- anova(tGAMw, tGAM2, test = "Chisq")$`Pr(>Chi)`[2]
            }

            if(!is.na(QCp)){
              if(QCp < QCcheckp){
                batchDat <- batchDat[batchDat$type != "QC",]
                new.k <- ifelse(k > nrow(batchDat), nrow(batchDat), k)
                tGAM <- mgcv::gam(y ~ s(x, k = new.k), data = batchDat,
                                  method = "REML", family=mgcv::scat(),
                                  weights = batchDat$weights)
              }
            } else {
              batchDat <- batchDat[batchDat$type != "QC",]
              new.k <- ifelse(k > nrow(batchDat), nrow(batchDat), k)
              tGAM <- mgcv::gam(y ~ s(x, k = new.k), data = batchDat,
                                method = "REML", family=mgcv::scat(),
                                weights = batchDat$weights)
            }
          }

          preds <- predict(tGAM, newdata = allBatchDat)
          predVals[allBatchIds] <- preds

          if(keepScale){
            normVal[allBatchIds] <- allBatchDat$y  - preds
          } else {
            normVal[allBatchIds] <- allBatchDat$y / preds
          }

        }
        if(keepScale){
          normVal <- normVal + mean(dat$y, na.rm = TRUE)
        }
      } else {
        tGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                          method = "REML", family = mgcv::scat(),
                          weights = datfit$weights)
        if(QCcheck & !QConly){
          tGAMw <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                            method = "REML", family = mgcv::scat())
          tGAM2 <- mgcv::gam(y ~ batch*isQC +
                               s(x, k = k, by = interaction(batch, isQC)),
                             data = datfit,
                             method = "REML", family = mgcv::scat(),
                             weights = datfit$weights)
          # if QCs differ significantly from samples, refit using samples only
          # if NA (fitting issue), also refit using samples only
          QCp <- anova(tGAMw, tGAM2, test = "Chisq")$`Pr(>Chi)`[2]
          if(!is.na(QCp)){
            if(QCp < QCcheckp){
              datfit <- dat[dat$type != "QC",]
              tGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                                method = "REML", family = mgcv::scat(),
                                weights = datfit$weights)
            }
          } else {
            datfit <- dat[dat$type != "QC",]
            tGAM <- mgcv::gam(y ~ batch + s(x, k = k, by = batch), data = datfit,
                              method = "REML", family = mgcv::scat(),
                              weights = datfit$weights)
          }
        }
        predVals <- predict(tGAM, newdata = dat)
        if(keepScale){
          normVal <- dat$y  - predVals + mean(dat$y, na.rm = TRUE)
        } else {
          normVal <- dat$y / predVals
        }
      }

    # single batch
    } else {
      tGAM <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                        method = "REML", family=mgcv::scat(),
                        weights = datfit$weights)
      if(QCcheck & !QConly){
        tGAMw <- mgcv::gam(y ~ s(x, k = k), data = datfit,
                          method = "REML", family=mgcv::scat())
        tGAM2 <- mgcv::gam(y ~ isQC + s(x, k = k, by = isQC),
                           data = datfit,
                           method = "REML", family = mgcv::scat())
        # if QCs differ significantly from samples, refit using samples only
        # if NA (fitting issue), also refit using samples only
        QCp <- anova(tGAMw, tGAM2, test = "Chisq")$`Pr(>Chi)`[2]
        if(!is.na(QCp)){
          if(QCp < QCcheckp){
            datfit <- dat[dat$type != "QC",]
            new.k <- ifelse(k > nrow(datfit), nrow(datfit), k)
            tGAM <- mgcv::gam(y ~ s(x, k = new.k), data = datfit,
                              method = "REML", family=mgcv::scat(),
                              weights = datfit$weights)
          }
        } else {
          datfit <- dat[dat$type != "QC",]
          new.k <- ifelse(k > nrow(datfit), nrow(datfit), k)
          tGAM <- mgcv::gam(y ~ s(x, k = new.k), data = datfit,
                            method = "REML", family=mgcv::scat(),
                            weights = datfit$weights)
        }
      }
      predVals <- predict(tGAM, newdata = dat)
      if(keepScale){
        normVal <- dat$y  - predVals + mean(dat$y, na.rm = TRUE)
      } else {
        normVal <- dat$y / predVals
      }
    }

  } else if(model == "rLOESS"){
    normVal <- numeric(length(dat$y))
    if(cv == "LOOCV"){
      for(batchid in 1:length(levels(dat$batch))){
        # fit the model on fitting data only (depends on QCOnly, see above)
        batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
        batchDat <- datfit[batchIds, ]
        # but predict on all batch data
        allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
        allBatchDat <- dat[allBatchIds, ]

        spanval <- selectSpan(batchDat$x, batchDat$y, method = "rLOESS")
        fit <- loess(y ~ x, data = batchDat, span = spanval, family = "symmetric")
        preds <- predict(fit, newdata = allBatchDat)
        predVals[allBatchIds] <- preds

        if(keepScale){
          normVal[allBatchIds] <- allBatchDat$y  - preds
        } else {
          normVal[allBatchIds] <- allBatchDat$y / preds
        }

      }
      if(keepScale){
        normVal <- normVal + mean(dat$y, na.rm = TRUE)
      }
    }
    if(cv == "GCV"){
      for(batchid in 1:length(levels(dat$batch))){
        # fit the model on fitting data only (depends on QCOnly, see above)
        batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
        batchDat <- datfit[batchIds, ]
        # but predict on all batch data
        allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
        allBatchDat <- dat[allBatchIds, ]
        # loess.as does not deal with NAs, so remove incomplete cases
        batchDat <- batchDat[complete.cases(batchDat),]
        fit <- fANCOVA::loess.as(batchDat$x, batchDat$y, degree = 2,
                                 criterion = "gcv", family = "symmetric")
        preds <- predict(fit, newdata = allBatchDat)
        predVals[allBatchIds] <- preds

        if(keepScale){
          normVal[allBatchIds] <- allBatchDat$y  - preds
        } else {
          normVal[allBatchIds] <- allBatchDat$y / preds
        }
      }
      if(keepScale){
        normVal <- normVal + mean(dat$y, na.rm = TRUE)
      }
    }
  } else if(model == "QC-RLSC"){
    normVal <- numeric(length(dat$y))
    if(cv == "LOOCV"){
      for(batchid in 1:length(levels(dat$batch))){
        # fit the model on fitting data only (depends on QCOnly, see above)
        batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
        batchDat <- datfit[batchIds, ]
        # but predict on all batch data
        allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
        allBatchDat <- dat[allBatchIds, ]

        spanval <- selectSpan(batchDat$x, batchDat$y, method = "QC-RLSC")
        fit <- loess(y ~ x, data = batchDat, span = spanval)
        preds <- predict(fit, newdata = allBatchDat)
        predVals[allBatchIds] <- preds
        if(keepScale){
          normVal[allBatchIds] <- allBatchDat$y  - preds
        } else {
          normVal[allBatchIds] <- allBatchDat$y / preds
        }
      }
      if(keepScale){
        normVal <- normVal + mean(dat$y, na.rm = TRUE)
      }
    }
    if(cv == "GCV"){
      for(batchid in 1:length(levels(dat$batch))){
        # fit the model on fitting data only (depends on QCOnly, see above)
        batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
        batchDat <- datfit[batchIds, ]
        # but predict on all batch data
        allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
        allBatchDat <- dat[allBatchIds, ]
        # loess.as does not deal with NAs, so remove incomplete cases
        batchDat <- batchDat[complete.cases(batchDat),]
        fit <- fANCOVA::loess.as(batchDat$x, batchDat$y, degree = 2,
                                 criterion = "gcv")
        preds <- predict(fit, newdata = allBatchDat)
        predVals[allBatchIds] <- preds

        if(keepScale){
          normVal[allBatchIds] <- allBatchDat$y  - preds
        } else {
          normVal[allBatchIds] <- allBatchDat$y / preds
        }
      }
      if(keepScale){
        normVal <- normVal + mean(dat$y, na.rm = TRUE)
      }
    }
  } else if(model == "QC-RSC"){
    normVal <- numeric(length(dat$y))
    if(cv == "LOOCV"){
      for(batchid in 1:length(levels(dat$batch))){
        # fit the model on fitting data only (depends on QCOnly, see above)
        batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
        batchDat <- datfit[batchIds, ]
        # but predict on all batch data
        allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
        allBatchDat <- dat[allBatchIds, ]
        # smooth.spline does not deal with NAs, so remove incomplete cases
        batchDat <- batchDat[complete.cases(batchDat),]
        fit <- smooth.spline(batchDat$y ~ batchDat$x, cv = TRUE)
        preds <- predict(fit, x = allBatchDat$x)$y
        predVals[allBatchIds] <- preds

        if(keepScale){
          normVal[allBatchIds] <- allBatchDat$y  - preds
        } else {
          normVal[allBatchIds] <- allBatchDat$y / preds
        }
      }
      if(keepScale){
        normVal <- normVal + mean(dat$y, na.rm = TRUE)
      }
    }
    if(cv == "GCV"){
      for(batchid in 1:length(levels(dat$batch))){
        # fit the model on fitting data only (depends on QCOnly, see above)
        batchIds <- which(datfit$batch == levels(datfit$batch)[batchid])
        batchDat <- datfit[batchIds, ]
        # but predict on all batch data
        allBatchIds <- which(dat$batch == levels(dat$batch)[batchid])
        allBatchDat <- dat[allBatchIds, ]
        # smooth.spline does not deal with NAs, so remove incomplete cases
        batchDat <- batchDat[complete.cases(batchDat),]
        fit <- smooth.spline(batchDat$y ~ batchDat$x, cv = FALSE)
        preds <- predict(fit, x = allBatchDat$x)$y
        predVals[allBatchIds] <- preds

        if(keepScale){
          normVal[allBatchIds] <- allBatchDat$y  - preds
        } else {
          normVal[allBatchIds] <- allBatchDat$y / preds
        }
      }
      if(keepScale){
        normVal <- normVal + mean(dat$y, na.rm = TRUE)
      }
    }
  } else {
    stop(paste0("unsupported method (", model, ")"))
  }


  if(!is.null(plotdir)){

    batchlabels <- seq(length(levels(dat$batch)))
    if(length(batchlabels < 6)){
      batchlabels <- 15:(15+length(batchlabels)-1)
    }
    if(plottype == "pdf"){
      dat$pred <- predVals
      preNorm <- ggplot(dat, aes(x = .data$x,
                                 y = .data$y,
                                 shape = batch,
                                 col = type))+
        geom_point()+
        scale_shape_manual(values = batchlabels)+
        geom_line(data = dat, aes(x = .data$x, y = .data$pred), col = "darkgrey")+
        theme_minimal()+
        xlab("injection order")+
        ylab("intensity")+
        theme(legend.text = element_text(size = 5),
              legend.title = element_text(size = 5),
              legend.key.height = unit(5, "points"),
              legend.spacing.y = unit(0, "cm"))
      datn <- dat
      datn$y <- normVal
      dat
      postNorm <- ggplot(datn, aes(x = .data$x,
                                   y = .data$y,
                                   shape = batch,
                                   col = type))+
        geom_point()+
        scale_shape_manual(values = batchlabels)+
        theme_minimal()+
        xlab("injection order")+
        ylab("normalized ntensity")+
        theme(legend.text = element_text(size = 5),
              legend.title = element_text(size = 5),
              legend.key.height = unit(5, "points"),
              legend.spacing.y = unit(0, "cm"))
      p <- plot_grid(preNorm, postNorm, nrow = 2)
      pdf(paste0(plotdir, "/", i, ".pdf"), width = 8, height = 5)
      print(p)
      dev.off()  
    } else if(plottype == "png"){
      dat$pred <- predVals
      preNorm <- ggplot(dat, aes(x = .data$x,
                                 y = .data$y,
                                 shape = batch,
                                 col = type))+
        geom_point(size=2)+
        scale_shape_manual(values = batchlabels)+
        geom_line(data = dat, aes(x = .data$x, y = .data$pred), col = "darkgrey")+
        theme_minimal()+
        xlab("injection order")+
        ylab("intensity")+
        theme(legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.key.height = unit(12, "points"),
              legend.spacing.y = unit(0, "cm"))
      datn <- dat
      datn$y <- normVal
      dat
      postNorm <- ggplot(datn, aes(x = .data$x,
                                   y = .data$y,
                                   shape = batch,
                                   col = type))+
        geom_point(size=2)+
        scale_shape_manual(values = batchlabels)+
        theme_minimal()+
        xlab("injection order")+
        ylab("normalized ntensity")+
        theme(legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              legend.key.height = unit(12, "points"),
              legend.spacing.y = unit(0, "cm"))
      p <- plot_grid(preNorm, postNorm, nrow = 2)
      png(paste0(plotdir, "/", i, ".png"), width = 800, height = 500)
      print(p)
      dev.off()
    }
    
  }

  return(normVal)
}

#' Main metanorm function
#'
#' @param mat (numeric) Numerical matrix, columns containing samples, and rows
#'   containing compounds
#' @param order (numeric) Injection order of the samples. May be NULL if the
#'   column order is chronological and no samples are omitted. If not
#'   chronological, or if samples were omitted, the argument's length should
#'   equal the number of samples. Default NULL.
#' @param keepScale (logical) Whether the original intensity scale should be
#'   retained, or whether rescaling can be performed. In the former case, fitted
#'   values are first subtracted from observed values, and the median observed
#'   value is re-added. In the latter case, normalization is achieved by dividing
#'   by the fitted value, and values are centered around one. Default TRUE.
#' @param QConly (logical) Whether only QC samples should be used for fitting
#'   the normalization model (FALSE), or whether also other samples contribute to
#'   the model fit (TRUE). Sample contributions can be modified by the parameter
#'   'weights'. Default FALSE.
#' @param QCcheck (logical) Whether to check for significantly different GAM
#'   fits between QCs and samples (TRUE) or not (FALSE). In the case of
#'   significantly different fits, QC samples are omitted and data normalized
#'   using the sample data only. Note: tGAM models do not allow for hypothesis
#'   testing with weights, therefor, with model = "tGAM", weights are ignored
#'   when checking for differences between QCs and samples. Default FALSE.
#' @param QCcheckp (numeric) Cutoff p-value for determining significantly 
#'   different fits between QCs and samples separately, versus QCs and samples
#'   combined. Default 0.1.
#' @param changepoints (logical) Whether to check for significant changepoints
#'   within batches (TRUE) or not (FALSE). Only for tGAM/rGAM. Default FALSE.
#' @param type (character) Whether the sample is a QC (value 'QC') or not (any
#'   other value). May be NULL, if 'QCOnly' is set to FALSE. Default NULL.
#' @param batch (character) Vector of length the number of samples, indicating
#'   for each sample to which batch it belongs. May be NULL, in which case samples
#'   are presumed to belong to a single batch. Default NULL.
#' @param batchwise (logical) Whether to fit batchwise GAMs. Batchwise fitting
#'   may result in significant speed gains for complex experiments with many
#'   batches, especially when QCcheck is set to TRUE. Note that this changes the
#'   QCcheck hypothesis test from a per-batch test (TRUE) to cross-batch tests
#'   (FALSE). Default TRUE.
#' @param weights (numeric) Vector of length the number of samples, indicating
#'   what weight should be attributed to each sample when fitting the
#'   normalization model. Can be NULL, in which case all samples will be assigned
#'   equal weight. Note: when fitting only to QC samples, the 'QCOnly' approach
#'   will be faster. Weights are only used for the rGAM and tGAM models, and
#'   ignored otherwise. Default NULL.
#' @param model (character) One of 'rGAM', a robustified Gaussian GAM; 'tGAM',
#'   a scaled-t GAM; 'QC-RLSC', as in Dunn et al. (2011); 'QC-RSC', as in Kirwan
#'   et al. (2013); or 'rLOESS', a robustified version of QC-RLSC. Default: tGAM.
#' @param gam.k (numeric) The maximum basis dimension to use when fitting GAMs.
#'   Lower values result in faster analysis time. Default 10.
#' @param gam.frac (numeric) 0 < gam.frac <= 1. If not exceeding 'gam.k', the
#'   fraction of samples to use as basis dimension. Lower values result in faster
#'   analysis time. Note: in presence of batches, should be lower than one.
#'   Default: 0.9. Note: most often use of this feature will require setting a
#'   higher value for 'gam.k'.
#' @param pb (logical) Whether to display a progress bar. While allowing
#'   progress tracking, there is no load balancing when set to TRUE, this often
#'   results in significantly longer analysis times. Default: FALSE.
#' @param ncpu (numeric/character) Either 'auto', in which case all but two CPUs
#'   will be used, or a number specifying the number of CPUs to use. Default:
#'   'auto'.
#' @param cv (character) Type of cross-validation to use for QC-RLSC, QC-RSC and
#'   rLOESS approaches. Either 'GCV' (faster) or 'LOOCV' (as suggested in the
#'   original publications, but slower). Default: GCV.
#' @param plotdir (character) When not set to NULL, a plot of intensities before
#'   and after normalization will be generated for each compound. Recommended for
#'   detailed, compound-level evaluation of normalization performance. The plots
#'   are stored in the given directory. Default: NULL.
#' @param plottype (character) File type for writing plots, either "pdf" or 
#'   "png", Default: pdf.
#'
#' @returns A matrix of normalized intensities.
#' @export
#'
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' @import parallel
metanorm <- function(mat, order = NULL, keepScale = TRUE,
                     QConly = FALSE, QCcheck = FALSE, QCcheckp = 1e-1,
                     changepoints = FALSE, type = NULL, batch = NULL,
                     batchwise = TRUE, weights = NULL,
                     model = c("tGAM","rGAM", "QC-RLSC", "QC-RSC", "rLOESS"),
                     gam.k = 10, gam.frac = 0.9,
                     pb = FALSE, ncpu = "auto", cv = c("LOOCV", "GCV"),
                     plottype = "pdf", plotdir = NULL) {

  model <- match.arg(model)
  cv <- match.arg(cv)


  ### Input validation

  if(!is.numeric(gam.k)){stop("'gam.k' should be numeric")}
  if(!is.numeric(gam.frac)){stop("'gam.frac' should be numeric")}
  if(!is.numeric(QCcheckp)){stop("'QCcheckp' should be numeric")}
  if(!is.logical(keepScale)){stop("'keepScale' should be TRUE/FALSE")}
  if(!is.logical(QConly)){stop("'QConly' should be TRUE/FALSE")}
  if(!is.logical(QCcheck)){stop("'QCcheck' should be TRUE/FALSE")}
  if(!is.logical(changepoints)){stop("'changepoints' should be TRUE/FALSE")}
  if(!is.logical(batchwise)){stop("'batchwise' should be TRUE/FALSE")}
  if(!is.logical(pb)){stop("'pb' should be TRUE/FALSE")}
  if(is.numeric(gam.frac)){
    if(gam.frac<= 0 | gam.frac > 1){stop("'gam.frac' should be strictly positive and not larger than 1 ([0 < gam.frac <= 1])")}
  }
  if(is.null(order)){
    order <- 1:ncol(mat)
  } else if(is.numeric(order)){
    if (ncol(mat) != length(order)){
      stop("Length of intensities ('mat' number of columns) should match the injection order ('order') length")
    }
  } else if(!is.numeric(order)){
    stop("'order' should be numeric")
  }
  if(is.null(type)){
    if(QConly){
      stop("'type' cannot be NULL when QConly is TRUE; its length should equal the number of samples ('mat' number of columns)")
    } else {
      type <- rep("t1", ncol(mat))
    }
  } else if(length(type) != ncol(mat)){
    stop("'type' should be NULL or its length should equal the number of samples ('mat' number of columns)")
  }

  if(QConly & (sum(type == "QC") < 4)){
    stop("At least 4 QC samples are needed to perform normalization with QCs only")
  }

  if(is.null(batch)){
    batch <- rep("b1", ncol(mat))
  } else if(length(batch) != ncol(mat)){
    stop("'batch' should be NULL or its length should equal the number of samples ('mat' number of columns)")
  }

  if(is.null(weights)){
    weights <- rep(1, ncol(mat))
  } else if(length(weights) != ncol(mat)){
    stop("'weights' should be NULL or its length should equal the number of samples ('mat' number of columns)")
  }

  if(!(model %in% c("rGAM", "tGAM", "QC-RLSC", "QC-RSC", "rLOESS"))){
    stop("Model should be one of 'rGAM', 'tGAM', 'QC-RLSC', 'QC-RSC', 'rLOESS'")
  }

  if(model %in% c("QC-RLSC", "QC-RSC", "rLOESS")){
    if(!(cv %in% c("GCV", "LOOCV"))){
      stop("cv should be one of 'GCV', 'LOOCV'")
    }
  }

  if(!is.null(plotdir)){
    if (!dir.exists(plotdir)) {
      dir.create(plotdir, recursive = TRUE)
      warning("The specified 'plotdir' did not exist so was created")
    }
  }

  if(!(plottype %in% c("pdf", "png"))){
    stop("plotyype should be one of 'pdf', 'png'")
  }

  # Set number of cores to use
  if(ncpu == "auto"){
    ncpus <- detectCores() - 2
  } else if(is.numeric(ncpu)) {
    ncpus <- ncpu
  } else {
    stop("ncpu argument should be either 'auto' or of type numeric.")
  }

  # set the basis dimension for GAM, minimum of number of observations
  #   times gam.frac,
  #   or cap at gam.k for large sample series
  k <- min(ncol(mat)*gam.frac, gam.k)

  # check that k is not too large
  if(QConly){
    if(gam.k < ncol(mat)*gam.frac){
      if(sum(type == "QC") < (gam.k  + length(levels(batch)) - 1)){
        warning("gam.k and/or gam.frac may be too large, reduce in case of errors")
      }
    } else {
      if(sum(type == "QC") < (ncol(mat)*gam.frac + length(levels(batch)) - 1)){
        warning("gam.k and/or gam.frac may be too large, reduce in case of errors")
      }
    }
  } else {
    if(gam.k < ncol(mat)*gam.frac){
      if(ncol(mat) < (gam.k  + length(levels(batch)) - 1)){
        warning("gam.k and/or gam.frac may be too large, reduce in case of errors")
      }
    } else {
      if(ncol(mat) < (ncol(mat)*gam.frac + length(levels(batch)) - 1)){
        warning("gam.k and/or gam.frac may be too large, reduce in case of errors")
      }
    }
  }

  cl <- makeCluster(ncpus)


  if(pb){
    pboptions(use_lb = TRUE)
    result <- pblapply(seq_len(nrow(mat)), function(i, mat, order = order,
                                                    keepScale = keepScale,
                                                    QConly = QConly,
                                                    QCcheck = QCcheck,
                                                    QCcheckp = QCcheckp,
                                                    changepoints = changepoints,
                                                    type = type, batch = batch,
                                                    batchwise = batchwise,
                                                    weights = weights,
                                                    model = model, k = k, cv = cv,
                                                    plotyype = plottype, 
                                                    plotdir = plotdir) {
      row <- unname(unlist(mat[i, ]))
      metanormWorker(row, order = order, keepScale = keepScale, QConly = QConly,
                     QCcheck = QCcheck, QCcheckp = QCcheckp,
                     changepoints = changepoints, type = type,
                   batch = batch, batchwise = batchwise, weights = weights,
                   model = model, k = k, cv = cv, plottype = plottype,
                   plotdir = plotdir, i = i)
    }, cl = cl, mat,  order = order, keepScale = keepScale, QConly = QConly,
    QCcheck = QCcheck, QCcheckp = QCcheckp, changepoints = changepoints,
    type = type, batch = batch, batchwise = batchwise, weights = weights,
    model = model, k = k, cv = cv, plottype = plottype, plotdir = plotdir)

  } else {
    result <- clusterApplyLB(cl, seq_len(nrow(mat)),
                             function(i, mat, order = order,
                                        keepScale = keepScale, QConly = QConly,
                                        QCcheck = QCcheck, QCcheckp = QCcheckp,
                                        changepoints = changepoints, type = type,
                                        batch = batch, batchwise = batchwise,
                                        weights = weights,
                                        model = model, k = k, cv = cv,
                                        plottype = plottype, plotdir = plotdir) {
      row <- unname(unlist(mat[i, ]))
      metanormWorker(row,  order = order, keepScale = keepScale, QConly = QConly,
                     QCcheck = QCcheck, QCcheckp = QCcheckp,
                     changepoints = changepoints, type = type,
                   batch = batch, batchwise = batchwise, weights = weights,
                   model = model, k = k, cv = cv, plottype = plottype,
                   plotdir = plotdir, i = i)
    }, mat, order = order, keepScale = keepScale, QConly = QConly,
    QCcheck = QCcheck, QCcheckp = QCcheckp, changepoints = changepoints,
    type = type, batch = batch, batchwise = batchwise,
    weights = weights, model = model, k = k, cv = cv, plottype = plottype, 
    plotdir = plotdir)
  }
  stopCluster(cl)
  return(do.call(rbind, result))
}
