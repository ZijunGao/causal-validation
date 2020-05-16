# the function randomly equally splits subgroups into folds
# input:
# groupMatrix: matrix of units in each group. Each row represents a group, the indices therein denote the treated and control units in the group
# fold: number of folds
# indexTreatment: indices of treated units
# indexControl: indices of control units
# tauMatrix: matrix of treatment effects within a group. Each row represents a group, the values denote the differences between the treat unit(s) and control unit(s) within the group
# treatmentMatrix: matrix of treated units in each group. Each row represents a group, the indices therein denote the treated units in the group
# controlMatrix: matrix of control units in each group. Each row represents a group, the indices therein denote the control units in the group
# output:
# foldIndicator: matrix of units in each fold. Each row stands for a fold, the indices therein denote the units in the fold
# foldTau: matrix of group treatment effects in each fold. Each row stands for a fold, the indices therein denote the units in the fold
# foldTreatment: matrix of units in each fold. Each row stands for a fold, the indices therein denote the treated units in the fold
# foldControl: matrix of units in each fold. Each row stands for a fold, the indices therein denote the control units in the fold
splitSubgroup = function(groupMatrix, fold = 5, indexTreatment, indexControl, tauMatrix, treatmentMatrix, controlMatrix){
  groupSize = apply((groupMatrix > 0), 1, sum)
  edgeSize = apply((groupMatrix != Inf), 1, sum)
  foldSize = rep(0, fold)
  foldIndicatorTemp = matrix(0, nrow = 2 * fold, ncol = ceiling(length(groupSize)/2/fold))
  foldIndicatorTemp[1:length(groupSize)] = order(groupSize, decreasing = TRUE)
  foldIndicator = matrix(0, nrow = fold, ncol = 2 * ceiling(length(groupSize)/2/fold) * max(groupSize))
  foldTau = foldTreatment = foldControl = matrix(Inf, nrow = fold, ncol = 2 * ceiling(length(groupSize)/2/fold) * max(edgeSize))
  for (i in 1:fold){
    temp = groupMatrix[as.vector(foldIndicatorTemp[c(i, 2*fold + 1 - i),]),]
    foldSize[i] = sum(temp > 0)
    foldIndicator[i,1:foldSize[i]] = temp[which(temp>0)]
    edgeTemp = tauMatrix[as.vector(foldIndicatorTemp[c(i, 2*fold + 1 - i),]),]
    foldTau[i,1:sum(edgeTemp!=Inf)] = edgeTemp[(edgeTemp!=Inf)]
    foldTreatment[i,1:sum(edgeTemp!=Inf)] = treatmentMatrix[as.vector(foldIndicatorTemp[c(i, 2*fold + 1 - i),]),][(edgeTemp!=Inf)]
    foldControl[i,1:sum(edgeTemp!=Inf)] = controlMatrix[as.vector(foldIndicatorTemp[c(i, 2*fold + 1 - i),]),][(edgeTemp!=Inf)]
  }
  foldIndicator = foldIndicator[, 1:max(foldSize)]
  maxEdgeSize = max(apply(foldTau != Inf, 1, sum))
  foldTau = foldTau[, 1:maxEdgeSize]
  foldTreatment = foldTreatment[,1:maxEdgeSize]
  foldControl = foldControl[,1:maxEdgeSize]
  
  result = list(); result$foldIndicator = foldIndicator; result$foldTau = foldTau; result$foldTreatment = foldTreatment; result$foldControl = foldControl
  return(result)
} 