# the function constructs proximity score distance 
# input:
# X: matrix of covariates
# Y: vector of responses
# W: vector of treatment assignment indicators
# nTree: number of trees in the random forest
# nodeSize: the minimal size of terminal nodes
# inputForest: whether to input a forest
# forest: the input forest if inputForest = TRUE
# output:
# distMatrix: matrix of distances. Each row represents a treated unit, and each column represents a control unit, the ij-th element represents the distance between the i-th treated unit and the j-th control unit.
# library("randomForest")
proximityScore = function(X, Y, W, nTree = 100, nodeSize = 10, inputForest = FALSE, forest = NULL) {
  # preprocessing
  indexControl = which(W == 0); indexTreatment = which(W == 1)
  # compute proximity score
  if(!inputForest) {forestControl = randomForest(x = X[indexControl,], y = Y[indexControl], ntree = nTree, nodesize = nodeSize, proximity = TRUE)} else {forestControl = forest}
  predictMu0 = predict(forestControl, newdata = X, predict.all = TRUE, proximity = TRUE)
  distMatrix = 1 - predictMu0$proximity[indexTreatment, indexControl]
  return(distMatrix)
}
