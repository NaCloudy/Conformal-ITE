fit_and_predict_band <- function(object, X_test, cfun, quantiles=c(0.4,0.6),outparams=list()){
  ########################################################################################
  ##  fit the regression for predicting interval, and predict the interval given data.  ##
  ########################################################################################

  chat <- object$chat
  upper <- chat$upper
  lower <- chat$lower



  Xval <- object$X_val
  Yval <- object$Y_val
  Tval <- object$T_val
  ite <- upper-lower

  print(paste0("t = 1 interval length ", mean(ite[which(Tval==1)])))
  print(paste0("t = 0 interval length ", mean(ite[which(Tval==0)])))



  ##################################################################
  ##           train the interval using validation data           ##
  ##################################################################


  low_outparams <- c(list(Y = lower, X = Xval, quantiles = quantiles[1]), outparams)
  up_outparams <- c(list(Y = upper, X = Xval, quantiles = quantiles[2]), outparams)


  low_Cmodel <-  function(X){
    do.call(cfun, c(low_outparams, list(Xtest=X)))
  }
  up_Cmodel <-function(X){
    do.call(cfun, c(up_outparams, list(Xtest=X)))
  }



  #################################################################
  ##                    predict the interval                     ##
  #################################################################


  interval <- data.frame(lower=low_Cmodel(X_test), upper= up_Cmodel(X_test), y1_mean= mean(ite[which(Tval==1)]), y0_mean = mean(ite[which(Tval==0)]))


  return(interval)
}
