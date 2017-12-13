
##' Stage (II) estimation
#'
#' Conduct aggregation procedure based on the output of 1st stage.
#' @param  train.x (matrix): the reduced feature training matrix (X) represents  n x q matrix
#' @param  test.x (matrix):  the reduced feature test matrix (X) represents  n x q matrix  (default: NULL )
#' @param  thetas (vector): represents thresholdings for the correlation matrix (default:  c(0, 0.2, 0.4, 0.6, 0.8  )  )
#' @param  groups (vector): represents the number of the group correponding to the thetas (default:  c( 1, 2, 4, 8, 10)  )
#' @return provides output, where
#'  factors_train (list[[1]): represents the latent factor loadings for the training reduced matrix, and
#'  factors_train (list[[2]): represents the latent factor loadings for the test reduced matrix.
#' @export
#' @examples
#' ###  generate dataset
#' data1 = dataset(v.nonsparse=c( 300, 300),  n.sparse=2500 ,  level= 0.02, N= 200   )
#' X = data1[[1]]
#' Y = data1[[2]]
#' ###  spit dataset into train and test
#' N= dim(X)[1]
#' T.row=seq( N )
#' Test=sort( sample(N)[1:(N*0.5)] )
#' Train=sort( T.row[-Test] )

#' X.train=X[Train, ]
#' X.test=X[Test,  ]
#' Y.train=Y[Train ]
#' Y.test=Y[Test ]

#' n.train=length(Y.train)
#' n.test=length(Y.test)
#' p=dim(X.train)[2]

#' SET_HS =  order(stage1,decreasing=TRUE)
#' plot(  stage1, xlab="X", ylab="signals")


### suppose cut-off to be the first 450 features! and focus on reduced matrix

#' G.NUM = 450
#' RX.train =X.train[, sort(SET_HS[1:G.NUM])  ]
#' RX.test = X.test [, sort(SET_HS[1:G.NUM])  ]



####  stage 2 ##### set thresholds and groups (length of these two should be the same )
#' stage2 =  fit_stage2( RX.train , RX.test,  thetas = c(0, 0.2, 0.4, 0.6, 0.8  ) , groups = c( 1, 2, 4, 8, 10) )




#########################################################################
fit_stage2 = function(  train.x=RX.train, test.x= NULL ,  thetas = c(0, 0.2, 0.4, 0.6, 0.8  ) ,
                        groups = c( 1, 2, 4, 8, 10)   ){

  outputs= vector("list", 2)
  RCov = cor(train.x )
  #thetas = c(0, 0.2, 0.4, 0.6, 0.8  )
  #groups = c( 1, 2, 4, 8, 10)

  aggrigate.PC =NULL
  aggrigate.EG =NULL


  for(i in 1:length(thetas) ){   ##  i = 1
    r0 = ComputeCluster(corrmat=RCov, cutoff= thetas[i] , groups= groups[i], "cen" )
    bpcs  = BLOCKPCAs(  RCov*(abs(RCov)> thetas[i] )  , r0[2,] , 1  )
    aggrigate.PC = cbind( aggrigate.PC , bpcs[[1]][,,1 ] )
    aggrigate.EG  = c( aggrigate.EG, bpcs[[2]][,1]      )
    aggrigate.EG  = ifelse( aggrigate.EG==0, 1, aggrigate.EG  )
  }

  n.agg = length(aggrigate.EG)
  outputs[[1]] =  train.x%*%aggrigate.PC%*%(solve( diag(aggrigate.EG^(0.5), ncol= n.agg,nrow= n.agg ) ) )
  if( is.null(test.x)!=TRUE ){
    outputs[[2]] =  test.x%*%aggrigate.PC%*%(solve( diag(aggrigate.EG^(0.5), ncol= n.agg,nrow= n.agg ) ) )   }
  if( is.null(test.x)==TRUE ){    outputs[[2]] = NULL  }

  names(outputs)=list("factors_train", "factors_test" )
  return(outputs) }
#########################################################################



##' estimates the pricipal components and corresponding eigen values for each cluster.
#'
#' @param  corrmat correlation matrix
#' @param  set the clustering indicator
#' @param  num The number of principal components
#' @return returns the block pcs and corresponding eigenvalues.
#' @export



BLOCKPCAs = function( corrmat  , set ,  num ){
  p =dim( corrmat)[2]
  n.b = length( unique(set) )
  varnames =seq(p)
  PCs =array(0, c(p, n.b, num) )
  Eigen = matrix(0, nrow=n.b, ncol=num)

  for( i in 1:n.b){   ## i=3
    if(num==1){
      idx = ( set == i )
      svd.obj = svd(  corrmat[ idx , idx ]  )
      local.pc=svd.obj$v[,num]
      if(  abs( sum(local.pc) )/sum(local.pc)  == -1 ){ local.pc=-1*local.pc }
      var.idx = varnames[ idx ]
      PCs[ var.idx  ,i, num]= round( local.pc, 7)
      Eigen[i,num]=round( svd.obj$d[num] ,7 )
    }

    if(num!=1){
      for( j in 1:num){   ## j=1  ##i=1
        idx = ( set == i )
        if(  sum(idx)>= num ){
          svd.obj = svd(  corrmat[ idx , idx ]  )
          local.pc =svd.obj$v[,j]
          if(  abs( sum(local.pc) )/sum(local.pc)  == -1 ){ local.pc=-1*local.pc }
          var.idx = varnames[ idx ]
          PCs[ var.idx  ,i, j]= round( local.pc, 7)
          Eigen[i,j]= round( svd.obj$d[j] ,7  )
        }  ## if close
        if(  sum(idx)< num ){
          svd.obj = svd(  corrmat[ idx , idx ]  )
          local.pc =svd.obj$v[,1]
          if(  abs( sum(local.pc) )/sum(local.pc)  == -1 ){ local.pc=-1*local.pc }
          var.idx = varnames[ idx ]
          if(j==1){   PCs[ var.idx  ,i, 1:j]= round( local.pc, 7); Eigen[i,j]= round( svd.obj$d[j] ,7 )  ;    }
          if(j!=1){   PCs[ var.idx  ,i, 1:j]=  0                 ; Eigen[i,j]= 0 ;    }
        }    ## if close
      }    ## for close
    }   ##if close
  } ## for n.b close
  return( list(PCs,Eigen) )       }

#########################################################################




##' Hierarchical Clustering based on Correlation distance
#'
#' Use the covariance thresholding method introduced by Bickel and Levina (2008) and the hierarchical clustering
#' method to partition the features into the multiple disjoint clusters.
#'
#' @param  corrmat correlation matrix
#' @param  cutoff thresholding value
#' @param  groups the number of the groups to be cluster
#' @param  hclustopt  the hierarchical cluster option (ex: "cen", "ward.D" )
#' @return matrix containing the features in the 1st row and the clustering indicator corresponding to the feature in the 2nd row
#' @export
#' @examples


ComputeCluster= function(  corrmat, cutoff, groups , hclustopt="ward.D"  ){
  distance0 = as.dist( round( matrix(1, ncol= dim(corrmat)[2] , nrow= dim(corrmat)[2])  -corrmat   ,5 ) )
  hc0 <- hclust( distance0 , hclustopt )
  r0  <- cutree(hc0, k = groups )
  features = seq(dim(corrmat)[2])
  clusters = r0
  outputs = rbind( features, clusters)
  return( outputs) }

#########################################################################

