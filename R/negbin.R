#' Dynamic Negative Binomial model
#'
#' Construct NegBin model on time series data
#'
#' The model constructs Negative Binomial model parametrised over the mean and
#' dispersion parameter using Simple Exponential Smoothing.
#'
#' @param data the vector of values (time series vector).
#' @param h forecasting horizon.
#' @param intervals binary, defining whether to construct prediction intervals or not.
#' @param levels confidence levels for prediction intervals.
#' @param holdout whether to use holdout of h observations or not.
#' @param cumulative if \code{TRUE}, cumulative values are produced.
#' @param side defines, whether to provide \code{"both"} sides of prediction
#' interval or only \code{"upper"}, or \code{"lower"}.
#' @param nsim the amount of sample paths for the estimation of the prediction intervals.
#'
#' @return Function returns a model of a class "counter", which contains:
#' \itemize{
#' \item model - the name of the constructed model,
#' \item occurrence - ETS(A,N,N) model for demand intervals (from es() function),
#' \item sizes - ETS(A,N,N) model for demand sizes (from es() function),
#' \item fitted - fitted values,
#' \item forecast - forecasts of the function,
#' \item lower - lower bound of the prediction interval,
#' \item upper - upper bound of the prediction interval,
#' \item actuals - the provided actual values,
#' \item holdout - the actual values from the holdout (if holdout was set to TRUE),
#' \item levels - confidence levels used,
#' \item probability - the probability of success,
#' \item dispersion - the dispersion variable,
#' \item alpha - smoothing parameter of the model,
#' \item initial - the initial value for the time varying mean,
#' \item accuracy - the error measures for the data if the holdout was TRUE.
#' }
#'
#' @references 
#' \itemize{
#' \item Snyder, R. D., Ord, J. K., & Beaumont, A. (2012). Forecasting the
#' intermittent demand for slow-moving inventories: A modelling approach.
#' International Journal of Forecasting, 28(2), 485–496.
#' https://doi.org/10.1016/j.ijforecast.2011.03.009
#' }
#' 
#' @author Ivan Svetunkov, \email{ivan@svetunkov.ru}
#' 
#' @keywords ts models
#' 
#' @seealso \code{\link[counter]{hsp}, \link[smooth]{adam}}
#'
#' @examples
#' y <- c(rpois(50,0.3),rpois(50,0.8))
#' test <- negbin(y)
#'
#' @importFrom nloptr nloptr
#' @importFrom stats deltat dnbinom frequency quantile rnbinom rpois start time ts
#' @export negbin
negbin <- function(data, h=10, intervals=TRUE, levels=0.95, holdout=FALSE,
                   cumulative=FALSE, side=c("both","upper","lower"), nsim=100000){
    
    side <- match.arg(side);
    
    # This function constructs Negative Binomial model proposed in Snyder et al.(2012).
    #
    ### Values:
    # probability - probability of success,
    # size - target number of successful trials,
    # alpha - smoothing parameter used in the model,
    # initial - initial value of the model,
    #
    ### Author
    # Ivan Svetunkov
    #
    ### Examples:
    # y <- c(rpois(50,0.3),rpois(50,0.8))
    # test <- negbin(y)
    
    # Define obs, the number of observations of in-sample
    obsInsample <- length(data) - holdout*h;
    
    # Define obsAll, the overal number of observations (in-sample + holdout)
    obsAll <- length(data) + (1 - holdout)*h;
    
    # If obsInsample is negative, this means that we can't do anything...
    if(obsInsample<=0){
        stop("Not enough observations in sample.",call.=FALSE);
    }
    # Define the actual values
    y <- matrix(data[1:obsInsample],obsInsample,1);
    datafreq <- frequency(data);
    
    at <- yFitted <- rep(NA,obsInsample);
    
    if(cumulative){
        hFinal <- 1;
    }
    else{
        hFinal <- h;
    }
    yForecast <- rep(NA,hFinal);
    yUpper <- yLower <- array(dim=c(hFinal, length(levels)))
    atForecast <- rep(NA,h);
    
    muFitter <- function(mu0, alpha){
        yFitted[1] <- mu0;
        for(i in 2:obsInsample){
            yFitted[i] <- alpha * y[i-1] + (1-alpha)*yFitted[i-1];
        }
        return(yFitted);
    }
    
    CF <- function(A){
        yFitted <- muFitter(A[1], A[2]);
        p <- A[3];
        at <- yFitted;
        if(p!=1){
            at[] <- at * p / (1-p);
        }
        
        CFValue <- -sum(dnbinom(y, at, p, log=TRUE));
        return(CFValue);
    }
    
    # Initial value, smoothing parameter and probability
    A <- c(mean(y), 0.1, sum(y!=0)/obsInsample);
    ALower <- rep(0, 3);
    AUpper<- c(Inf, 1, 1);
    
    res <- nloptr::nloptr(A, CF, lb=ALower, ub=AUpper,
                          opts=list("algorithm"="NLOPT_LN_BOBYQA", "xtol_rel"=1e-8, "maxeval"=500));
    
    A <- res$solution;
    
    p <- A[3];
    yFitted <- muFitter(A[1], A[2]);
    at <- yFitted;
    if(p!=1){
        at[] <- at * p / (1-p);
    }
    atForecast[] <- at[length(at)];
    
    levelsLow <- levelsUp <- levelsNew <- levels;
    if(side=="both"){
        levelsLow[] <- (1-levelsNew)/2;
        levelsUp[] <- (1+levelsNew)/2;
    }
    else if(side=="upper"){
        levelsLow[] <- 0;
        levelsUp[] <- levelsNew;
    }
    else{
        levelsLow[] <- 1-levelsNew;
        levelsUp[] <- 1;
    }
    levelsLow[levelsLow<0] <- 0;
    levelsUp[levelsUp<0] <- 0;
    
    ySimulated <- matrix(NA, nsim, h);
    
    for(i in 1:h){
        ySimulated[,i] <- rnbinom(nsim,atForecast[i],p);
    }
    
    if(cumulative){
        if(intervals){
            yUpper[] <- quantile(rowSums(ySimulated),probs=levelsUp);
            yLower[] <- quantile(rowSums(ySimulated),probs=levelsLow);
        }
        yForecast[] <- mean(rowSums(ySimulated));
    }
    else{
        if(intervals){
            yUpper[] <- t(apply(ySimulated,2,quantile,probs=levelsUp));
            yLower[] <- t(apply(ySimulated,2,quantile,probs=levelsLow));
        }
        yForecast[] <- apply(ySimulated,2,mean);
    }
    
    if(side=="upper"){
        yLower[] <- 0;
    }
    
    yHoldoutStart <- time(data)[obsInsample]+deltat(data);
    
    if(holdout){
        yHoldout <- ts(data[(obsInsample+1):obsAll],start=yHoldoutStart,frequency=datafreq);
        if(cumulative){
            errormeasures <- measures(sum(yHoldout),yForecast,y,digits=5);
        }
        else{
            errormeasures <- measures(yHoldout,yForecast,y,digits=5);
        }
    }
    else{
        yHoldout <- NA;
        errormeasures <- NA;
    }
    
    yForecast <- ts(yForecast,start=yHoldoutStart,frequency=datafreq);
    yUpper <- ts(yUpper,start=yHoldoutStart,frequency=datafreq,names=as.character(levelsUp));
    yLower <- ts(yLower,start=yHoldoutStart,frequency=datafreq,names=as.character(levelsLow));
    yFitted <- ts(yFitted,start=start(y),frequency=datafreq);
    at <- ts(c(at,atForecast),start=start(y),frequency=datafreq);
    
    model <- list(model="NegBin", fitted=yFitted,forecast=yForecast,lower=yLower,upper=yUpper,
                  probabilty=p,dispersion=at,alpha=A[2],initial=A[1],
                  actuals=data,holdout=yHoldout,accuracy=errormeasures,levels=levels);

    return(structure(model, class="counter"));
    
}