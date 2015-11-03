# http://stackoverflow.com/questions/1753299/help-using-predict-for-kernlabs-svm-in-r
ksvm.cv <- function (kmat, y, n.cv=1000, stats=c("accuracy"), show.progress=TRUE, nthread=1,
                     ksvm.type=ifelse(is.factor(y), "C-svc", "eps-svr"))
{
    get.train <- function(k, holdout) as.kernelMatrix(k[-holdout,-holdout])
    get.test <- function(k, holdout, m) as.kernelMatrix(k[holdout, -holdout][,SVindex(m), drop=F])
    fit.model <- function(k, train.y) ksvm(k, train.y, kernel="matrix", type=ksvm.type)
    predict <- kernlab::predict
    do.cv(kmat, y, get.train, get.test, fit.model, predict, n.cv=n.cv,
          stats=stats, show.progress=show.progress, nthread=nthread)
}

lm.cv <- function (data, y, n.cv=1000, stats=c("rsquared"), show.progress=TRUE, nthread=1)
{
    get.train <- function (d, holdout) d[-holdout,]
    get.test <- function (d, holdout, m) d[holdout,]
    fit.model <- function (d, train.y) {
        train.data <- cbind(d, y=train.y)
        frm <- as.formula(paste0("y~", paste(colnames(d), collapse="+")))
        lm(frm, data=train.data)
    }
    predict <- stats::predict
    do.cv(data, y, get.train, get.test, fit.model, predict, n.cv=n.cv,
          stats=stats, show.progress=show.progress, nthread=nthread)
}

rpart.cv <- function (data, y, n.cv=1000, stats=c("accuracy"), show.progress=TRUE, nthread=1)
{
    get.train <- function (d, holdout) d[-holdout,]
    get.test <- function (d, holdout, m) d[holdout,]
    fit.model <- function (d, train.y) {
        train.data <- cbind(d, y=train.y)
        frm <- as.formula(paste0("y~", paste(colnames(d), collapse="+")))
        model <- rpart(frm, data=train.data, method="class")
        # http://www.statmethods.net/advstats/cart.html
        prune(model, cp = model$cptable[which.min(model$cptable[,"xerror"]),"CP"])
    }
    predict <- function (m, d) stats::predict(m, data=d, type="vector")
    do.cv(data, y, get.train, get.test, fit.model, predict, n.cv=n.cv,
          stats=stats, show.progress=show.progress, nthread=nthread)
}

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
