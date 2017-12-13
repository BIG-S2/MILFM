dataset = function(   v.nonsparse=c( 300, 300),  n.sparse=2500 ,  level= 0.02, N= 200   ){

  library(MASS)

  #v.nonsparse=c( 300, 300)
  n.nonsparse=sum( v.nonsparse )
  #n.sparse=2500

  CR1 =block_corr( v.nonsparse, c(0.9,0.8) , 0.5)
  X_1 = mvrnorm(n = N, rep(0, sum(v.nonsparse) ), CR1)


  mu_T = apply( X_1 , 1, mean)
  mu_G1= apply( X_1[, (1:v.nonsparse[1])] , 1, mean)
  mu_G2= apply( X_1[,(  (v.nonsparse[1]+1):sum(v.nonsparse)  )] , 1, mean)
  POP_mu =cbind( mu_T,  mu_G1,  mu_G2 )

  X_0 = matrix( rnorm(n.sparse*N,0,1),ncol=n.sparse, nrow=N)
  X=cbind(X_1,X_0)

  Y =   sqrt(  mu_T^2+mu_G1^2+mu_G2^2)   +  log ( sqrt( mu_T^2+mu_G1^2+mu_G2^2)  )   +level*rnorm(  N,0,1)
  outputs = list(X,Y)
  return(outputs )
}




block_corr=function( n.groups, corr.groups , corr.betweens ){
  CORR=matrix(corr.betweens, ncol= sum(n.groups) , nrow=sum(n.groups))

  for( i in 1:length(n.groups)){ ## i=1
    end=sum( n.groups[1:i] )
    if(i==1){start=1}
    if(i!=1){start= sum( n.groups[(1:(i-1))] ) +1 }

    CORR[(start:end),(start:end)]    =matrix(corr.groups[i], ncol= n.groups[i] , nrow=n.groups[i])
  }
  diag(CORR)=1
  return(CORR) }
