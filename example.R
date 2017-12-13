
rm(list=ls())


source("dataset.R")

source("fit_stage1.R")
source("fit_stage2.R")
source("fit_stage3.R")




###  generate dataset
data1 = dataset(v.nonsparse=c( 300, 300),  n.sparse=2500 ,  level= 0.02, N= 200   )
X = data1[[1]]
Y = data1[[2]]

###  spit dataset into train and test
N= dim(X)[1]

T.row=seq( N )
Test=sort( sample(N)[1:(N*0.5)] )
Train=sort( T.row[-Test] )

X.train=X[Train, ]   
X.test=X[Test,  ]    
Y.train=Y[Train ]
Y.test=Y[Test ]

n.train=length(Y.train)
n.test=length(Y.test)
p=dim(X.train)[2]




####  stage 1 #####
stage1 =  fit_stage1( X.train , Y.train)
SET_HS =  order(stage1,decreasing=TRUE)
plot(  stage1, xlab="X", ylab="signals")


### suppose cut-off to be the first 450 features! and focus on reduced matrix

G.NUM = 450
RX.train =X.train[, sort(SET_HS[1:G.NUM])  ]
RX.test = X.test [, sort(SET_HS[1:G.NUM])  ]



####  stage 2 ##### set thresholds and groups (length of these two should be the same )
stage2 =  fit_stage2( RX.train , RX.test,  thetas = c(0, 0.2, 0.4, 0.6, 0.8  ) , groups = c( 1, 2, 4, 8, 10) )

##stage2$factors_train
##stage2$factors_test



####  stage 3 ##### run regression model 
stage3_1 = fit_stage3(input.x=stage2$factors_train, output.y= Y.train, test.x=stage2$factors_test,test.y=Y.test, type= 1 )
## test error from kernel ridge estimator 
test.error_1 =    sum( stage3_1$T_residual^2)/n.test

stage3_2 = fit_stage3(input.x=stage2$factors_train, output.y= Y.train, test.x=stage2$factors_test,test.y=Y.test, type= 2 )
## test error from support vector regression estimator 
test.error_2 =    sum( stage3_2$T_residual^2)/n.test


## our performance

test.error_1

test.error_2

















































