##' Stage (I) estimation
#'
#' Use a sure Independence  Screening (SIS) procedure based on a Hilbert-Schmidt Independence Criterion (HSIC) in the reproducing Hilbert kernel space to select
#'  a set of important features.
#' @param  x is matrix n times p containing  all candidate features (p) in columns and all observations (n) in rows.
#' @param  y is outcome variable.
#' @return provides the HSIC statistic as between marginal feature and the outcome as two-variable independence test
#'  in Reproducing Kernel Hilbert Spaces (RKHS).
#' @export
#' @examples
#' x=matrix( runif(100,0,1), ncol=10)
#' y=matrix( runif(10,0,1), ncol=1)
#' fit_stage1( x, y)

fit_stage1 = function( x, y) {


  sigma=1

  raw.X=standard_nonlinear( x )

  raw.Y=standard_nonlinear( y )

  n=dim(raw.X)[1]
  p=dim(raw.X)[2]

  eps=0.001*(n)^(-0.1)

  int=rep(1 ,n )
  H=(diag(int)-int%*%t(int)/n)
  I=diag(int)


  HSIC=rep(0,p)

  int=rep(1,n )
  H=(diag(int)-int%*%t(int)/n)

  status=unique(raw.Y)
  if(length(status)==2){
    Class=unique(raw.Y)
    Class1=(raw.Y==Class[1])+0
    Ly1=Class1%*%t(Class1)
    Class2=(raw.Y==Class[2])+0
    Ly2=Class2%*%t(Class2)
    Ly=Ly1+Ly2
  }
  if(length(status)==3){
    Class=unique(raw.Y)
    Class1=(raw.Y==Class[1])+0
    Ly1=Class1%*%t(Class1)
    Class2=(raw.Y==Class[2])+0
    Ly2=Class2%*%t(Class2)
    Class3=(raw.Y==Class[3])+0
    Ly3=Class3%*%t(Class3)
    Ly=Ly1+Ly2+Ly3
  }


  if(length(status)>3){
    Ly=Rbfkernel(raw.Y,raw.Y, sigma)
  }
  CLy= Ly%*%H
  Ry=  solve( CLy+ n*eps*I )


  for(i in 1:p ){
    Kx=Rbfkernel(raw.X[,i],raw.X[,i], sigma)

    CKx= Kx%*%H
    Rx=  solve( CKx+ n*eps*I )

    HSIC[i]= round( sum(  diag( Ly%*%H%*%Ry%*%Kx%*%H%*%Rx) ), 6 )
  }
  return(HSIC) }

###############################################################################

Rbfkernel =function( z, x, sigma){

  x=as.matrix(x)
  z=as.matrix(z)

  n.sb=dim(z)[1]
  n.ob=dim(x)[1]
  p.sb=dim(z)[2]
  p.ob=dim(x)[2]

  if(p.sb!=p.ob){ stop( "Check p")}
  p=p.sb

  Kernel_Matrix=matrix(0,n.sb,n.ob)
  distance=array(0,c(n.sb, p ,n.ob))
  for(i in 1:n.ob){
    distance[,,i]= ( z -  matrix(rep( x[i,],n.sb),ncol=p, byrow = TRUE) )^2
  }

  Kernel_Matrix=exp(-sigma*apply( distance,c(1,3),sum) )
  return(Kernel_Matrix)
}
###############################################################################

standard_nonlinear=function(x){
  x=as.matrix(x)
  n=dim(x)[1]
  p=dim(x)[2]
  output=matrix(0,ncol=p,nrow=n)
  for(i in 1:p){
    output[,i]=(x[,i]-min(x[,i]) )/(max(x[,i])-min(x[,i]) )
  }
  return(output)}



