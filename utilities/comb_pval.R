###
# Usage:
###
# This script is used to combine p-values using different methods.
# R codes for sumlog (Fisher's) and sumz (Liptak's) functions are modified from metap package (Michael Dewey, January 22, 2017), for liptak function is modified from Dr. Blanca Himes' codes.

# R functions
## weight and statistics
pcomp <- function(dt) {    
    ref_study <- dt[order(-dt$Total,dt$adj.P.Val),"Unique_ID"][1] # the study with largest sample size (if the studies have equal sample size, take the one with smallest p-values) is used for reference effect direction, and is used as a discovery study with two-sided p-values. Sort sample size by decreasing order and p-values by increasing order
    ref_sign <- sign(dt$logFC[which(dt$Unique_ID==ref_study)]) # use direction in reference study as reference
    dt$stats <- rep(NA,length(dt$Total))
    dt$stats[which(dt$Unique_ID==ref_study)] <- dt$adj.P.Val[which(dt$Unique_ID==ref_study)] # use two-sided p-value for reference study
    dt$stats[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)==ref_sign))] <- dt$adj.P.Val[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)==ref_sign))]/2 # if direction is same as the reference direction, convert two-sided p-value to one-sided p-value
    dt$stats[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)!=ref_sign))] <- 1-dt$adj.P.Val[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)!=ref_sign))]/2 # if direction is opposite to the reference direction, convert two-sided p-value to one-sided p-value and invert using 1-(p-value/2)
    return(dt)
}

# Liptak-combined p-values. Liptak: P=sum(qnorm(1 - x)) (formular: PMID: 21605215)
liptak <- function(t, w){
    d <- 1/sqrt(sum(w^2))
    pnorm( d*sum(w*qnorm(t)) )
}

# Combine p-values by the sum of logs method, also known as Fisher’s method, and sometimes as the chi-square (2) method (taken from metaP function sumlog)
sumlog <- function(p) {
    keep <- (p > 0) & (p <= 1)
    lnp <- log(p[keep])
    chisq <- (-2) * sum(lnp)
    df <- 2 * length(lnp)
    if(sum(1L * keep) < 2)
        stop("Must have at least two valid p values")
    if(length(lnp) != length(p)) {
        warning("Some studies omitted")
    }
    res <- list(chisq = chisq, df = df,
                p = pchisq(chisq, df, lower.tail = FALSE), validp = p[keep])
    return(res)
}

# Combine p-values using the sum of p method also known as Edgington’s method
sumz <- function(p, weights = NULL, data = NULL, subset = NULL, na.action = na.fail)  {
    if(is.null(data)) data <- sys.frame(sys.parent())
    mf <- match.call()
    mf$data <- NULL
    mf$subset <- NULL
    mf$na.action <- NULL
    mf[[1]] <- as.name("data.frame")
    mf <- eval(mf, data)
    if(!is.null(subset)) mf <- mf[subset,]
    mf <- na.action(mf)
    p <- as.numeric(mf$p)
    weights <- mf$weights
    if(is.null(weights)) weights <- rep(1, length(p))
    if(length(p) != length(weights)) warning("Length of p and weights differ")
    keep <- (p > 0) & (p < 1)
    if(sum(1L * keep) < 2)
        stop("Must have at least two valid p values")
    if(sum(1L * keep) != length(p)) {
        warning("Some studies omitted")
        omitw <- weights[!keep]
        if(sum(1L * omitw) > 0) warning("Weights omitted too")
    }
    zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep]) /
        sqrt(sum(weights[keep]^2))
    res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE),
                validp = p[keep], weights = weights)
    return(res)
}

#Liptak-combined p-values, weighted by sample size:
liptak_stat <- function(dt){
    dt <- pcomp(dt)
    weights <- dt$Total # sample size
    stats <- dt$stats # modified p-value
    pcomb <- liptak(stats,weights)
    pcomb <- ifelse(is.numeric(pcomb) & length(pcomb)>0, round(pcomb,3), pcomb)
    return(pcomb)
}

#Edgington-combined p-values, weighted by sample size:
sumz_stat <- function(dt){
    dt <- pcomp(dt)
    weights <- dt$Total # sample size
    stats <- dt$stats # modified p-value
    pcomb <- sumz(stats,weights)$p[1,1]
    pcomb <- ifelse(is.numeric(pcomb) & length(pcomb)>0, round(pcomb,3), pcomb)
    return(pcomb)
}

#Sum of logs method, known as Fisher/chisqr-combined p-values
sumlog_stat <- function(dt){
    dt <- pcomp(dt)
    stats <- dt$stats # modified p-value
    pcomb <- sumlog(stats)$p
    pcomb <- ifelse(is.numeric(pcomb) & length(pcomb)>0, round(pcomb,3), pcomb)
    return(pcomb)
}

