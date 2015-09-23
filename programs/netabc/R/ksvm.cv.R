# http://stackoverflow.com/questions/1753299/help-using-predict-for-kernlabs-svm-in-r
ksvm.cv <- function (kmat, y, n.cv=1000, stats=c("accuracy"), show.progress=TRUE)
{
    if (is.factor(y))
    {
        allowed.stats <- c("accuracy")
        ksvm.type <- "C-svc"
    } 
    else 
    {
        allowed.stats <- c("rsquared")
        ksvm.type <- "eps-bsvr"
    }

    if (!all(stats %in% allowed.stats))
    {
        bad.stats <- paste0(stats[!stats %in% allowed.stats], collapse=", ")
        stop(sprintf("Unrecognized statistics: %s", bad.stats))
    }

    stat.functions <- list(
        accuracy = function (pred, y) sum(pred == y) / length(y),
        rsquared = function (pred, y) cor.test(pred, y)$estimate^2
    )

    if (show.progress)
        pb <- txtProgressBar(max=n.cv-1, style=3, file=stderr())

    result <- do.call(rbind, replicate(n.cv, {
        holdout <- sample.int(nrow(kmat), nrow(kmat)/2)
        
        train <- as.kernelMatrix(kmat[-holdout,-holdout])
        m <- ksvm(train, y[-holdout], kernel="matrix", type=ksvm.type)
        
        test <- as.kernelMatrix(kmat[holdout, -holdout][,SVindex(m), drop=F])
        pred <- kernlab::predict(m, test)
        if (show.progress)
            setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
        sapply(stat.functions[stats], do.call, list(pred, y[holdout]))
    }, simplify=FALSE))
    colnames(result) <- stats

    if (show.progress)
        close(pb)
    result
}
