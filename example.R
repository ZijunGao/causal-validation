# example of heterogeneous treatment effect (HTE) estimator validation based on proximity score and fact match. Cross-validation based.
# HTE is estimated via LASSO
# input: all parameters are set inside the function
# output:
# oracle: MSE of the estimators from LASSO corresponding to a sequence of tuning parameters
# MSE: the MSE of the estimator from LASSO chosen by the validation procedure

# library("randomForest")
# library("glmnet")
# library("gurobi")
example = function(){
  n = 200 # number of samples
  d = 10 # dimension of predictors
  sigma = 1 # noise magnitude
  alpha = beta = theta = rep(0, (d+1)); beta[(d+1)] = c(1); 
  alpha[1:5] = c(2,1.5,0,0,0); beta[1:5] = -c(2,1.5,-1,-0.8,-0.5); theta[1] = 1
  delta = 2 # magnitude of misspecification of the linear model for mu0
  lambdaSeq = (sqrt(2))^seq(-6, -16, -1); nfold = 10; nTree = 1000; nodeSize = 10
  max.controls = 5
  
  # data generation
  X = matrix(runif(d * n, -1, 1), ncol = d)
  prop =  exp(cbind(X,1) %*% theta)/(1 + exp(cbind(X,1) %*% theta))
  W = rbinom(n, 1, prop)
  tau = cbind(X, 1) %*% beta
  mu0 = cbind(X, 1) %*% alpha + delta * (X[,1])^2
  noise = rnorm(n, 0, sigma)
  Y = W * tau + mu0 + noise
  
  # preprocessing
  indexTreatment = which(W == 1); indexControl = which(W == 0); nTreatment = sum(W == 1); nControl = sum(W == 0)
  if(nTreatment <= nControl) {reverse = 1} else {reverse = -1; W = 1 - W; mu0 = mu0 + tau; tau = -tau; prop = 1 - prop; indexTreatment = which(W == 1); indexControl = which(W == 0); nTreatment = sum(W == 1); nControl = sum(W == 0)}
  XSML = cbind(X, diag(as.vector(W)) %*% cbind(X, 1)) 
  distMatrix = proximityScore(X = X, Y = Y, W = W, nTree = nTree, nodeSize = nodeSize)
  # fact match with pruning
  pairs = factMatch(distMatrix = distMatrix, maxDegreeTreatment = max.controls, maxDegreeControl = max.controls)$solution
  pairs = factMatchPrune(matchMatrix = pairs, distMatrix = distMatrix, indexTreatment = indexTreatment, indexControl = indexControl, Y = Y) 
  # split matched groups into folds
  folds = splitSubgroup(groupMatrix = pairs$groupMatrix, tauMatrix = pairs$tauMatrix, fold = nfold, indexControl = indexControl, indexTreatment = indexTreatment, treatmentMatrix = pairs$treatmentMatrix, controlMatrix = pairs$controlMatrix)
  foldIndicator = folds$foldIndicator; foldTau =  folds$foldTau; foldTreatment = folds$foldTreatment; foldControl = folds$foldControl
  
  # record
  record = rep(0, length(lambdaSeq))
  # cross-validation
  for(i in 1:nfold) {
    modelLASSO = glmnet(x = XSML[-foldIndicator[i,],], y = Y[-foldIndicator[i,]], family = c("gaussian"), lambda = lambdaSeq) 
    betaSML = modelLASSO$beta[(d+1) : (2*d+1),] 
    nonInfLen = sum(foldTau[i, ] != Inf)
    record = record + apply((foldTau[i,1:nonInfLen] -  cbind(X[foldTreatment[i, 1:nonInfLen],], 1) %*% betaSML)^2, 2, sum)
  }
  
  # train on the whole data
  modelLASSOFinal = glmnet(x = XSML, y = Y, family = "gaussian", lambda = lambdaSeq)
  betaSMLFinal = modelLASSOFinal$beta[(d+1) : (2*d+1),]
  # oracle MSE of beta
  oracle = apply((betaSMLFinal - beta * reverse)^2, 2, mean)
  indexMin = which.min(record); MSE = oracle[indexMin]
  
  # result
  result = list(); result$oracle = oracle; result$MSE = MSE
  return(result)
}
