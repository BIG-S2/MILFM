# MILFM
multiple-index latent factor modeling

The aim of this package is to develop a multiple-index latent factor modeling (MILFM) framework to build an accurate 
prediction model for clinical outcomes based on a massive number of features. 
We develop a three-stage estimation procedure to build the prediction model. 
MILFM uses an independent screening method to select a set of informative features, 
which may have a complex nonlinear relationship with the outcome variables. Moreover, 
we develop a latent factor model to project all informative predictors onto a small number of local subspaces, 
which lead to a few key features that capture reliable and informative covariate information. 
Finally, we fit the regularized empirical estimate to those key features in order to accurately predict clinical outcomes.

##  Please read the followings and run example.R. 


dataset.R  -- Generate simulation dataset (simulation scenario)

fit_stage1.R -- source code: how to compute non-parametric correlation measure based on the reproducing kernel 


fit_stage2.R -- source code: how to conduct aggregation procedure based on the output of 1st stage


fit_stage3.R -- fit nonparametric regression based on the repreducing kernel approach 



## Main- source code 


### dataset.R  -- produces dataset in the paper 


input --
 
	 v.nonsparse (vector):  represents the length of nonzero signals for each group to be generated. 
	 n.sparse (scalar):  indiates the lenght of zero signals  to be generated. 
	 level (scalar): noise level.   N (scalar): sample size 
 
output --
 
	 list[[1]]: covariate X.  list[[2]]: outcome  Y 



### fit_stage1.R -- source code: how to compute non-parametric correlation measure based on the reproducing kernel. 


input --

	 x (matrix):  X represents n x p matrix. y (vector):  Y represents n x 1 vector
   
output --

       vector:  represents marginal HSIC norm based on the RBF kernel.




### fit_stage2.R -- source code: how to conduct aggregation procedure based on the output of 1st stage
 

input --

	 train.x (matrix):  reduced feature training matrix (X) represents  n x q matrix
	 test.x (matrix):  reduced feature test matrix (X) represents  n x q matrix  (default: NULL ). thetas (vector): represents thresholdings for the correlation matrix (default:  c(0, 0.2, 0.4, 0.6, 0.8  )  ). groups (vector): represents the number of the group correponding to the thetas (default:  c( 1, 2, 4, 8, 10)  )

output --
 
	factors_train (list[[1]): represents the latent factor loadings for the training reduced matrix. 
	factors_train (list[[2]): represents the latent factor loadings for the test reduced matrix.




### fit_stage3.R -- fit nonparametric regression based on the repreducing kernel approach 
 

input --
  
	 input.x (matrix):  the latent factor loadings for the training reduced matrix based on the output of 2nd stage.
	 output.y (vector):  the outcome Y 
	 keropt (scalar): means the reproducing kernel. Specifically, there are options as 1=rbf, 2=linear, and 3= polynomial (default: 1 ). penalty (scalar): means regularized constant (default: 0:001 X n^(-2/3) ) 
	 test.x (matrix): the latent factor loadings for the test reduced matrix based on the output of 2nd stage  (default: NULL ) .
	 test.y (vector): the outcome Y test
	 type (scalar): means the fitting model type. Specifically, there are options as 1=kernel ridge regression and  2=support vector regression (required field ) 


output --
 
	fx (vector): represents the fitted value for the training. 
	alpha (vector: represents the primal solution in the reproducing kernel space when the kernel ridge was used.
	ybar (scalar): represents the average of the output. 
	labmda (scalar): number used as the regularized constant when the kernel ridge was used.
	residual (vector): produces the residual for the training. sigma2 (list[[2]): estimates the conditional variance of the outcome.
	T_fx (vector): predicted values for the test. 
	T_residual (vector): produces the prediction error for the test dataset.
