###
# Usage:
###
# This script is used for meta-analysis. The function meta.summaries is modified from rmeta package
# It estimates the variance of a weighted average using the squares of the weights (standard error). Random model is applied.

meta.summaries<-function(d,se,method=c("fixed","random"),
                         weights=NULL,logscale=FALSE,names=NULL,data=NULL,
                         conf.level=.95,subset=NULL,na.action=na.fail)
{
    if (conf.level>1 & conf.level<100)
        conf.level<-conf.level/100
    if ( is.null( data ) ) 
        data <- sys.frame( sys.parent() )
    mf <- match.call()
    mf$data <- NULL 
    mf$method<-NULL  
    mf$logscale<-NULL
    mf$subset <- NULL
    mf[[1]] <- as.name( "data.frame" )
    mf <- eval( mf,data )
    if ( !is.null( subset ) ) 
        mf <- mf[subset,]
    mf <- na.action( mf )
    
    if (is.null(mf$names)){
        if (is.null(mf$d) || is.null(names(mf$d)))
            mf$names<-seq(along=mf$d)
        else
            mf$names<-names(mf$d)
    }
    mf$names<-as.character(mf$names)
    method<-match.arg(method)
    vars<-mf$se^2
    
    vwts<-1/vars
    fixedsumm<-sum(vwts*mf$d)/sum(vwts)
    Q<-sum( ( ( mf$d - fixedsumm )^2 ) / vars ) 
    df<-NROW(mf)-1
    
    tau2<-max(0, (Q-df)/(sum(vwts)-sum(vwts^2)/sum(vwts)))
    
    if(is.null(mf$weights)){
        if (method=="fixed"){
            wt<-1/vars
        } else {
            wt<-1/(vars+tau2)
        }
    } else
        wt<-mf$weights
    
    summ<-sum(wt*mf$d)/sum(wt)
    if (method=="fixed")
        varsum<-sum(wt*wt*vars)/(sum(wt)^2)
    else
        varsum<-sum(wt*wt*(vars+tau2))/(sum(wt)^2)
    
    summtest<-summ/sqrt(varsum)
    
    df<-length(vars)-1
    rval<-list(effects=mf$d, stderrs=mf$se, summary=summ,se.summary=sqrt(varsum),
               test=c(summtest,1-pchisq(summtest^2,1)),
               het=c(Q,df,1-pchisq(Q,df)),
               call=match.call(), names=mf$names,tau2=tau2,
               variance.method=method, weights=wt, 
               weight.method=if(is.null(mf$weights)) method else "user",
               conf.level=conf.level,logscale=logscale)
    class(rval)<-"meta.summaries"
    rval
}

meta_stat <- function(dat) {
    d <- dat$logFC
    se <- dat$SD
    names=dat$Unique_ID
    method="random"
    conf.level=0.95
    x <- meta.summaries(d=d, se=se, names=names, method=method, conf.level=conf.level, logscale=FALSE)
    meta_fc=x$summary
    ci.value<- -qnorm((1-conf.level)/2)
    ci<-x$summary+c(-ci.value,0,ci.value)*x$se.summary
    meta_lower=ci[1]
    meta_upper=ci[3]
    meta_pval=x$test[2]
    res=list(meta_pval=meta_pval,meta_fc=2^meta_fc,meta_lower=2^meta_lower,meta_upper=2^meta_upper)
    return(res)
}