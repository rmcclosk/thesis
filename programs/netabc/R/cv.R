#' Generic cross-validation.
#'
#' Given data x and labels y, perform replicate 2-fold cross-validations to
#' determine the accuracy of a classifier. The parameter x can by anything; y
#' must be either a factor or a numeric vector.
#'
#' The other mandatory parameters are functions for doing CV-related tasks.

#' get.train will be passed two parameters, x and holdout, where x was input to
#' the function and holdout is an integer vector of indices of y which will be
#' held out and should be removed from x. 

#' get.test will be passed x, holdout, and m, where m is the fitted model. It
#' is necessary to pass the fitted model so that we can for example get the
#' indices of the support vectors.
#'
#' fit.model will be passed training data (as returned by get.train) and a
#' corresponding subset of y with the test labels held out. It should return a
#' model trained on the training data.
#'
#' predict will be passed the fitted model (as returned by fit.model) and the 
#' test data (as returned by get.test). It should return predicted labels for
#' the test set.
#'
#' Some summary statistics are calculated for each cross-validation, and these
#' are returned as a data.frame. Currently the only available statistics are
#' "rsquared" for numeric labels or "accuracy" for categorical labels.
#'
#' @param x data, can by anything
#' @param y labels, either a factor or a numeric
#' @param get.train function to create training subset of x
#' @param get.test function to create testing subset of x
#' @param n.cv number of cross-validations to perform
#' @param stats statistics to calculate
#' @param show.progress whether to display a progress bar
#' @param nthread number of threads to use
do.cv <- function (x, y, get.train, get.test, fit.model, predict, n.cv=1000,
                   stats=c("accuracy"), show.progress=TRUE, nthread=1)
{
    if (is.factor(y)) {
        allowed.stats <- c("accuracy")
    } else {
        allowed.stats <- c("rsquared")
    }

    if (!all(stats %in% allowed.stats)) {
        bad.stats <- paste0(stats[!stats %in% allowed.stats], collapse=", ")
        stop(sprintf("Unrecognized statistics: %s", bad.stats))
    }

    stat.functions <- list(
        accuracy = function (pred, y) sum(as.integer(pred) == as.integer(y)) / length(y),
        rsquared = function (pred, y) cor.test(pred, y)$estimate^2
    )

    if (show.progress) {
        pb <- txtProgressBar(max=n.cv-1, style=3, file=stderr())
    }

    result <- do.call(rbind, mclapply(1:n.cv, function (i) {
        holdout <- sample.int(nrow(x), nrow(x)/2)
        
        train <- get.train(x, holdout)
        m <- fit.model(train, y[-holdout])
        
        test <- get.test(x, holdout, m)
        pred <- predict(m, test)
        if (show.progress) {
            setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
        }
        sapply(stat.functions[stats], do.call, list(pred, y[holdout]))
    }, mc.cores=nthread))
    colnames(result) <- stats

    if (show.progress) {
        close(pb)
    }
    result
}

#' Perform replicate 2-fold cross-validations on a kSVM classifier.
#'
#' Returns a data.frame with the indicated statistics for each individual
#' cross-validation. Currently the only supported statistics are "accuracy" for
#' categorical labels or "rsquared" for numeric labels.
#'
#' @param kmat kernel matrix
#' @param y labels for data
#' @param stats statistics to calculate
#' @param show.progress whether to display a progress bar
#' @param nthread number of threads to use
#' @param ksvm.type type of kernel-SVM classifier
#' @return A data.frame with the indicated statistics for each CV
#' @seealso \url{http://stackoverflow.com/questions/1753299}
#' @seealso \code{\link{kernlab::ksvm}}
#' @export
ksvm.cv <- function (kmat, y, n.cv=1000, 
                     stats=ifelse(is.factor(y), "accuracy", "rsquared"),
                     show.progress=TRUE, nthread=1,
                     ksvm.type=ifelse(is.factor(y), "C-svc", "eps-svr"))
{
    get.train <- function(k, holdout) as.kernelMatrix(k[-holdout,-holdout])
    get.test <- function(k, holdout, m) as.kernelMatrix(k[holdout, -holdout][,SVindex(m), drop=F])
    fit.model <- function(k, train.y) ksvm(k, train.y, kernel="matrix", type=ksvm.type)
    predict <- kernlab::predict
    do.cv(kmat, y, get.train, get.test, fit.model, predict, n.cv=n.cv,
          stats=stats, show.progress=show.progress, nthread=nthread)
}

#' Perform replicate 2-fold cross-validations on an SVM classifier.
#'
#' Returns a data.frame with the indicated statistics for each individual
#' cross-validation. Currently the only supported statistics are "accuracy" for
#' categorical labels or "rsquared" for numeric labels.
#'
#' @param dmat distance matrix
#' @param y labels for data
#' @param stats statistics to calculate
#' @param show.progress whether to display a progress bar
#' @param nthread number of threads to use
#' @return A data.frame with the indicated statistics for each CV
#' @seealso \url{http://stackoverflow.com/questions/1753299}
#' @seealso \code{\link{kernlab::ksvm}}
#' @export
dsvm.cv <- function (dmat, y, n.cv=1000, 
                    stats=ifelse(is.factor(y), "accuracy", "rsquared"),
                    show.progress=TRUE, nthread=1)
{
    get.train <- function(d, holdout) d[-holdout,-holdout]
    get.test <- function(d, holdout, m) d[holdout, -holdout]
    fit.model <- function(d, train.y) svm(d, train.y)
    predict <- predict
    do.cv(dmat, y, get.train, get.test, fit.model, predict, n.cv=n.cv,
          stats=stats, show.progress=show.progress, nthread=nthread)
}

#' Perform replicate 2-fold cross-validations on a linear model classifier.
#'
#' Returns a data.frame with the indicated statistics for each individual
#' cross-validation. Currently the only supported statistic is "rsquared".
#'
#' @param x predictors
#' @param y labels for data
#' @param stats statistics to calculate
#' @param show.progress whether to display a progress bar
#' @param nthread number of threads to use
#' @return A data.frame with the indicated statistics for each CV
#' @seealso \url{http://stackoverflow.com/questions/1753299}
#' @seealso \code{\link{kernlab::ksvm}}
#' @export
lm.cv <- function (x, y, n.cv=1000, stats="rsquared", show.progress=TRUE,
                   nthread=1)
{
    get.train <- function(x, holdout) x[-holdout,,drop=FALSE]
    get.test <- function(x, holdout, m) x[holdout,,drop=FALSE]
    fit.model <- function(x, train.y) lm(train.y~pred, data=x)
    predict <- function(m, x.test) predict.lm(m, newdata=x.test)
    x.df <- data.frame(pred=x)
    do.cv(x.df, y, get.train, get.test, fit.model, predict, n.cv=n.cv,
          stats=stats, show.progress=show.progress, nthread=nthread)
}
