

##' Stage (III) estimation
#'
#' Fit nonparametric regression based on the repreducing kernel approach.
#' @param  input.x (matrix): the latent factor loadings for the training reduced matrix based on the output of 2nd stage.
#' @param  output.y (vector):  the outcome Y
#' @param  keropt (scalar): means the reproducing kernel. Specifically, there are options as 1=rbf, 2=linear, and 3= polynomial (default: 1 )
#' @param  penalty (scalar): means regularized constant (default: 0:001 X n^(-2/3) )
#' @param  test.x (matrix): the latent factor loadings for the test reduced matrix based on the output of 2nd stage  (default: NULL ) .
#' @param  test.y (vector): the outcome Y test
#' @param  type (scalar): means the fitting model type. Specifically, there are options as 1=kernel ridge regression and  2=support vector regression (required field )
#' @return provides the fallowing outputs.
#'	fx (vector): represents the fitted value for the training.
#'	alpha (vector: represents the primal solution in the reproducing kernel space when the kernel ridge was used.
#'	ybar (scalar): represents the average of the output.
#'	labmda (scalar): number used as the regularized constant when the kernel ridge was used.
#'	residual (vector): produces the residual for the training.
#'	sigma2 (list[[2]): estimates the conditional variance of the outcome.
#'	T_fx (vector): predicted values for the test.
#'	T_residual (vector): produces the prediction error for the test dataset.
#' @export
#' @examples
#' ###  generate dataset
#' ####  stage 3 ##### run regression model
#' stage3_1 = fit_stage3(input.x=stage2$factors_train, output.y= Y.train, test.x=stage2$factors_test,test.y=Y.test, type= 1 )
#' ## test error from kernel ridge estimator
#' test.error_1 =    sum( stage3_1$T_residual^2)/n.test

#' stage3_2 = fit_stage3(input.x=stage2$factors_train, output.y= Y.train, test.x=stage2$factors_test,test.y=Y.test, type= 2 )
#' ## test error from support vector regression estimator
#' test.error_2 =    sum( stage3_2$T_residual^2)/n.test









fit_stage3 = function( input.x=stage2$factors_train, output.y= Y.train , keropt=1 , penalty=NULL,
                       test.x=stage2$factors_test, test.y=Y.test ,  type= 2   ){

  kerns  = c("radial", "linear", "polynomial")

  library(kernlab)
  library(e1071)
  n = dim(input.x )[1]
  p = dim(input.x )[2]

  if(type==2){
    TRAIN_1 <- data.frame( input.x  ,  Y =output.y,  row.names = NULL)
    FIT_svm_1= svm(  Y ~ .,  TRAIN_1, probability =TRUE, type="eps-regression", kernel=kerns[keropt] ,cost=1 )
    Hat_fx =as.matrix( predict(FIT_svm_1, TRAIN_1) )
    residual = output.y - Hat_fx
    sigma2 = sum( residual^2 )/length( n)

    if( is.null(test.x)==TRUE ){   TEST_1= NULL; PRE_SVR_1 = NULL;   }
    if( is.null(test.x)!=TRUE ){
      n.test = length(test.y)
      TEST_1  <- data.frame(  test.x ,  Y =rep("NA",length(n.test) ), row.names = NULL )
      T_fx  =as.matrix(predict(FIT_svm_1, TEST_1) )
      T_residual = test.y- T_fx
    }   ## if null close
    alpha=lambda = NULL
    output = list(   Hat_fx ,  alpha ,  mean(output.y), lambda ,residual, sigma2 ,   T_fx ,   T_residual      )
    names(output)=list("fx", "alpha", "ybar", "labmda","residual", "sigma2" ,  "T_fx" ,   "T_residual"    )
  }    ## if type  close

  if(type==1){

    rk =  keropt
    liner.ker <- polydot(degree = 1, scale = 1, offset = 0)
    quadr.ker <- polydot(degree = 2, scale = 1, offset = 0)
    rbf.ker <- rbfdot(sigma = 1/p )
    rk =  keropt
    RKkerS =c("rbf.ker", "liner.ker", "quadr.ker")
    RKker_text= RKkerS[rk]
    RKker = eval(parse( text = RKker_text) )


    if( is.null(penalty) ==TRUE){   lambda=0.001*(n)^(-0.1) }
    if( is.null(penalty) !=TRUE){   lambda = penalty    }

    KMP.train = kernelMatrix(RKker ,  input.x,  input.x)
    I=diag( rep(1,n) , nrow=n, ncol=n)
    int=I-(1/n)*rep(1,n)%*%t(rep(1,n))

    GKMP= int%*%KMP.train %*%int
    I_GKMP=solve( GKMP + lambda*I)

    Hat_fx =  GKMP%*%I_GKMP%*%output.y  + mean(output.y)*rep(1,n)
    alpha =   I_GKMP%*%output.y
    residual =  output.y -  Hat_fx
    sigma2 = sum( residual^2 )/length( n)

    if( (is.null(test.x)!=TRUE ) & (is.null(test.y) !=TRUE)  ) {
      n.test = length(test.y)
      if( n.test>1 ){         KMP.test  = kernelMatrix(RKker ,  test.x,   input.x)   }
      if( n.test==1 ){        KMP.test  = kernelMatrix(RKker ,  matrix( test.x, nrow=1 ),   input.x)   }
      cenRK =  (KMP.test%*%int - (1/n)*rep(1,n.test)%*%t(rep(1,n))%*%KMP.train%*%int)
      T_fx = cenRK%*%I_GKMP%*%output.y  + mean(output.y )*rep(1,n.test)
      T_residual =  test.y -T_fx
    }

    if( (is.null(test.x)==TRUE ) | (is.null(test.y) ==TRUE) ) {
      T_fx = NULL
      T_residual = NULL
    }
    output= list(   Hat_fx ,  alpha ,  mean(output.y), lambda ,residual, sigma2 , T_fx ,   T_residual      )
    names(output)=list("fx", "alpha", "ybar", "labmda","residual", "sigma2" ,  "T_fx" ,   "T_residual"  )
  }    ## if type  close

  return(output) }
