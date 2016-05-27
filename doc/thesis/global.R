#!/usr/bin/env Rscript

library(netabc)
library(Hmisc)
library(xtable)

pp <- function (p, eq=FALSE) {
    if (p < 1e-5) {
        "{<}10^{-5}"
    } else {
        n <- latexSN(round(p, -floor(log10(p))))
        n <- sub("^1\\\\!\\\\times\\\\!", "", n)
        if (eq) {
            paste0("=", n)
        } else {
            n
        }
    }
}

f <- Sys.glob("../../simulations/abc-pa-free-m/point-estimate/*")
d <- fread(f)
d <- d[,m := round(m)]
d <- d[,alpha_error := abs(true_alpha - alpha)]
d <- d[,N_error := abs(true_N - N)]
d <- d[,I_error := abs(true_I - I)]
d <- d[,m_error := abs(true_m - m)]

make.glm <- function (param, family, link) {
    frm <- paste0(param, "_error~true_alpha+true_I+true_m")
    m <- glm(as.formula(frm), family(link=link), d)
    df <- as.data.frame(coef(summary(m)))
    df$Parameter <- sub(").*", "", sub("true_", "", rownames(df), fixed=TRUE))
    df$Parameter[1] <- "(Intercept)"
    setDT(df)
    setnames(df, "Std. Error", "Standard error")
    if ("Pr(>|t|)" %in% colnames(df)) {
        df$`t value` <- NULL
        setnames(df, "Pr(>|t|)", "p")
    } else {
        df$`z value` <- NULL
        setnames(df, "Pr(>|z|)", "p")
    }
    setcolorder(df, c("Parameter", "Estimate", "Standard error", "p"))
    df[,p := p.adjust(p, method="holm", n=12)]
    df
}
 
# make GLMs
alpha.glm <- make.glm("alpha", gaussian, "inverse")
I.glm <- make.glm("I", gaussian, "inverse")
N.glm <- make.glm("N", gaussian, "inverse")
m.glm <- make.glm("m", poisson, "log")
