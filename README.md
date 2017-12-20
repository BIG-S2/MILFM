# MILFM


 
## Overview  <a name="Overview"></a>
The aim of this paper is to develop a multiple-index latent factor modeling (MILFM) framework to  build an accurate prediction model   for clinical outcomes based on a massive number of features.   We develop a three-stage estimation procedure to build the prediction model. MILFM uses an independent screening method to select a set of  informative features, which may have a complex  nonlinear relationship with the outcome variables. Moreover, we develop   a  latent factor model  to project all  informative predictors onto  a small number of local subspaces, which lead to a few key features  that capture reliable and informative covariate information.  Finally, we fit the regularized empirical estimate to those key features in order to accurately  predict clinical outcomes. 


# Data analysis <a name="Dataanalysis"></a>
###  generate dataset and  split dataset into train and test
data1 = dataset(v.nonsparse=c( 300, 300),  n.sparse=2500 ,  level= 0.02, N= 200   )

X = data1[[1]];

Y = data1[[2]];

N= dim(X)[1];

T.row=seq( N );

Test=sort( sample(N)[1:(N*0.5)] );

Train=sort( T.row[-Test] );

X.train=X[Train, ]   ;

X.test=X[Test,  ] ;   

Y.train=Y[Train ];

Y.test=Y[Test ];

n.train=length(Y.train);

n.test=length(Y.test);

p=dim(X.train)[2];

###  stage 1 
stage1 =  fit_stage1( X.train , Y.train);

SET_HS =  order(stage1,decreasing=TRUE);

plot(  stage1, xlab="X", ylab="signals");

(suppose cut-off to be the first 450 features! and focus on reduced matrix)

G.NUM = 450

RX.train =X.train[, sort(SET_HS[1:G.NUM])  ]

RX.test = X.test [, sort(SET_HS[1:G.NUM])  ]


###  stage 2
stage2 =  fit_stage2( RX.train , RX.test,  thetas = c(0, 0.2, 0.4, 0.6, 0.8  ) , groups = c( 1, 2, 4, 8, 10) )

f2_out_train= stage2$factors_train

f2_out_test= stage2$factors_test

###  stage 3
stage3_1 = fit_stage3(input.x=f2_out_train, output.y= Y.train, test.x=f2_out_test,test.y=Y.test, type= 1 )

stage3_2 = fit_stage3(input.x=f2_out_train, output.y= Y.train, test.x=f2_out_test,test.y=Y.test, type= 2 )






##Installing <a name="Installing"></a>
```{r, eval=FALSE, warning=FALSE,message=FALSE}
install.packages("MILFM")
```

## Summary statistics <a name="Summarystatistics"></a>

###  stage 1 
  fit_stage1(x, y) --  how to compute non-parametric correlation measure based on the reproducing kernel. 
  Input x (matrix):  X represents n x p matrix. y (vector):  Y represents n x 1 vector.
  Output vector:  represents marginal HSIC norm based on the RBF kernel.



###  stage 2
fit_stage2( train.x, test.x,  thetas , groups)--  how to conduct aggregation procedure based on the output of 1st stage.
 Input
	 train.x (matrix):  reduced feature training matrix (X) represents  n x q matrix.
 	 test.x (matrix):  reduced feature test matrix (X) represents  n x q matrix  (default: NULL ). 
 	 thetas (vector): represents thresholdings for the correlation matrix (default:  c(0, 0.2, 0.4, 0.6, 0.8  )  ).
	 groups (vector): represents the number of the group correponding to the thetas (default:  c( 1, 2, 4, 8, 10)  ).
 Output
	factors_train (list[[1]): represents the latent factor loadings for the training reduced matrix. 
	factors_train (list[[2]): represents the latent factor loadings for the test reduced matrix.


###  stage 3
fit_stage3(input.x, output.y , keropt=1 , penalty,test.x, test.y,  type= 2)--  fit nonparametric regression based on the repreducing kernel approach.
Input
 input.x (matrix):  the latent factor loadings for the training reduced matrix based on the output of 2nd stage.
output.y (vector):  the outcome Y. 
keropt (scalar): means the reproducing kernel. Specifically, there are options as 1=rbf, 2=linear, and 3= polynomial (default: 1 ). 
enalty (scalar): means regularized constant (default: 0:001 X n^(-2/3) ). 
test.x (matrix): the latent factor loadings for the test reduced matrix based on the output of 2nd stage  (default: NULL ).
test.y (vector): the outcome Y test.
type (scalar): means the fitting model type. Specifically, there are options as 1=kernel ridge regression and  2=support vector regression (required field ). 
Output
	fx (vector): represents the fitted value for the training. 
	alpha (vector: represents the primal solution in the reproducing kernel space when the kernel ridge was used.
	ybar (scalar): represents the average of the output. 
	labmda (scalar): number used as the regularized constant when the kernel ridge was used.
	residual (vector): produces the residual for the training.
	sigma2 (list[[2]): estimates the conditional variance of the outcome.
	T_fx (vector): predicted values for the test. 
	T_residual (vector): produces the prediction error for the test dataset.












