#' Hurdle Shifted Poisson model
#'
#' Construct HSP model on time series data
#'
#' HSP model consists of two parts: Simple Exponential Smoothing for demand intervals
#' and the same method for the mean demand sizes (lambda parameter in Poisson distribution).
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
#' @seealso \code{\link[counter]{negbin}, \link[smooth]{adam}}
#'
#' @examples
#' y <- c(rpois(50,0.3),rpois(50,0.8))
#' test <- hsp(y)
#'
#' @importFrom greybox measures
#' @importFrom smooth es
#' @export hsp
hsp <- function(data, h=10, intervals=TRUE, levels=0.95, holdout=FALSE,
                cumulative=FALSE, side=c("both","upper","lower"), nsim=100000){
    
    side <- match.arg(side);
    
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
    
    # Define occurrence variable
    ot <- (y!=0)*1;
    obsNonzero <- sum(ot);
    
    # Define non-zero demand
    yot <- matrix(y[y!=0],obsNonzero,1);
    
    #### Fit the model to the demand intervals ####
    # Define the matrix of demand intervals
    zeroes <- c(0,which(y!=0));
    zeroes <- diff(zeroes);
    # Number of intervals in Croston
    iyt <- matrix(zeroes,length(zeroes),1);
    newh <- which(y!=0);
    newh <- newh[length(newh)];
    newh <- obsInsample - newh + h;
    
    modelOccurrence <- smooth::es(iyt,"ANN",h=newh,silent=TRUE);
    
    pFitted <- rep((modelOccurrence$fitted),zeroes);
    tailNumber <- obsInsample - length(pFitted);
    if(tailNumber>0){
        pFitted <- c(pFitted,modelOccurrence$forecast[1:tailNumber]);
    }
    pForecast <- modelOccurrence$forecast[(tailNumber+1):newh];
    
    pFitted <- ts(1/pFitted,start=start(y),frequency=frequency(y));
    pForecast <- ts(1/pForecast, start=time(y)[obsInsample]+deltat(y),frequency=frequency(y));
    states <- 1/modelOccurrence$states;
    
    #### Fit the model to the demand sizes ####
    modelSizes <- smooth::es(yot,"ANN",h=h,silent=TRUE);
    
    lambdaFitted <- rep(0,obsInsample);
    m <- 0;
    for(i in 1:obsInsample){
        if(y[i]!=0){
            m <- m + 1;
        }
        
        if(m==0){
            lambdaFitted[i] <- modelSizes$fitted[1];
        }
        else{
            lambdaFitted[i] <- modelSizes$fitted[m];
        }
    }
    
    lambdaForecast <- modelSizes$forecast;
    
    yFitted <- lambdaFitted*pFitted;
    
    if(cumulative){
        hFinal <- 1;
    }
    else{
        hFinal <- h;
    }
    
    yForecast <- rep(NA,hFinal);
    yUpper <- yLower <- array(dim=c(hFinal, length(levels)));
    
    levelsLow <- levelsUp <- levelsNew <- levels;
    levelsNew[levelsNew<0] <- 0;
    
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
        ySimulated[,i] <- (rpois(nsim,lambdaForecast[i]-1)+1)*(as.integer(runif(nsim)<pForecast[i]));
    }
    
    if(cumulative){
        pForecast <- pForecast[1];
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
        otHoldout <- yHoldout!=0;
        if(cumulative){
            errormeasures <- measures(sum(yHoldout),sum(pForecast*yForecast/h),y,digits=5);
        }
        else{
            errormeasures <- measures(yHoldout,yForecast,y,digits=5);
        }
    }
    else{
        yHoldout <- NA;
        errormeasures <- NA;
    }
    
    yForecast <- ts(pForecast*yForecast,start=yHoldoutStart,frequency=datafreq);
    yUpper <- ts(yUpper,start=yHoldoutStart,frequency=datafreq,names=as.character(levelsUp));
    yLower <- ts(yLower,start=yHoldoutStart,frequency=datafreq,names=as.character(levelsLow));
    yFitted <- ts(yFitted,start=start(y),frequency=datafreq)
    
    model <- list(model="HSP",occurrence=modelOccurrence,sizes=modelSizes,
                  fitted=yFitted,forecast=yForecast,lower=yLower,upper=yUpper,
                  actuals=data,holdout=yHoldout,accuracy=errormeasures,levels=levels);

    return(structure(model, class="counter"));
}