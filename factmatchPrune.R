# the function prunes a match output from the function "factmatch" to obey the multiple-to-one or one-to-multiple constraint
# input:
# matchMatrix: matrix of pairs: the ij-th (row i column j) element is 1 if treated unit i and control unit j are matched, otherwise is 0
# distMatrix: matrix of distances. Each row represents a treated unit, and each column represents a control unit, the ij-th element represents the distance between the i-th treated unit and the j-th control unit.
# indexTreatment: indices of treated units
# indexControl: indices of control units
# Y: responses (treated and control)
# output: 
# numberGroup: number of multiple-to-one or one-to-multiple groups
# groupSize: vector of numbers of samples in each group
# tauMatrix: matrix of treatment effects within a group. Each row represents a group, the values denote the differences between the treat unit(s) and control unit(s) within the group
# groupMatrix: matrix of units in each group. Each row represents a group, the indices therein denote the treated and control units in the group
# treatmentMatrix: matrix of treated units in each group. Each row represents a group, the indices therein denote the treated units in the group
# controlMatrix: matrix of control units in each group. Each row represents a group, the indices therein denote the control units in the group
# matchMatrixPrune: the matrix of pairs after pruning
factMatchPrune = function(matchMatrix, distMatrix, indexTreatment, indexControl, Y){
  # preprocessing
  nTreatment = dim(matchMatrix)[1]; nControl = dim(matchMatrix)[2]; 
  degreeTreatment = apply(matchMatrix, 1, sum); degreeControl = apply(matchMatrix, 2, sum)
  edges = matrix(0, nrow = sum(matchMatrix), ncol = 5)
  edgeIndex = as.vector(which(matchMatrix > 0))
  # edge value
  edges[, 1] = as.vector(distMatrix)[edgeIndex]
  # edge start vertex
  edges[, 2] = rep(seq(1,nTreatment), nControl)[edgeIndex]
  # edge end vertex
  edges[, 3] = rep(seq(1,nControl), rep(nTreatment, nControl))[edgeIndex]
  # degree of start vertex
  edges[, 4] = rep(degreeTreatment, nControl)[edgeIndex]
  # degree of end vertex
  edges[, 5] = rep(degreeControl, rep(nTreatment, nControl))[edgeIndex]
  
  # removable edges 
  removableIndex = which((edges[,4] > 1) & (edges[,5] > 1))
  edges = matrix(edges[removableIndex,], ncol = 5)
  
  if (length(removableIndex) != 0) {
    while(max(edges[,1]) >= 0){
      removeIndex = which.max(edges[,1])
      removeStartNode = edges[removeIndex, 2]; removeEndNode = edges[removeIndex, 3]
      degreeTreatment[removeStartNode] = degreeTreatment[removeStartNode] - 1; degreeControl[removeEndNode] = degreeControl[removeEndNode] - 1
      removeEdge = which((edges[, 2] == removeStartNode[(degreeTreatment[removeStartNode] == 1)]) & (edges[, 1] >= 0))
      removeEdge = c(removeEdge, which((edges[, 3] == removeEndNode[(degreeControl[removeEndNode] == 1)]) & (edges[, 1] >= 0)))
      edges[removeEdge, 1] = -1
      edges[removeIndex, 1] = -10
    } 
  }
  matchMatrixPrune = matchMatrix
  matchMatrixPrune[which(matchMatrix > 0)[removableIndex[(edges[,1] == -10)]]] = 0
  
  temp = matchToGroup(matchMatrix = matchMatrixPrune, indexControl = indexControl, indexTreatment = indexTreatment, Y = Y)
  result = list(); result$numberGroup = dim(temp$groupMatrix)[1]; result$groupSize = apply((temp$groupMatrix > 0),1,sum); result$groupMatrix = temp$groupMatrix; result$tauMatrix = temp$tauMatrix; result$treatmentMatrix= temp$treatmentMatrix; result$controlMatrix = temp$controlMatrix; result$matchMatrixPrune = matchMatrixPrune
  return (result)
}

# helper functions for function "factMatchPrune"

# convert a multiple-to-one/one-to-multiple match matrix to a group matrix
# input:
# matchMatrix: matrix of pairs: the ij-th (row i column j) element is 1 if treated unit i and control unit j are matched, otherwise is 0
# indexTreatment: indices of treated units
# indexControl: indices of control units
# Y: responses (treated and control)
# output: 
# tauMatrix: matrix of treatment effects within a group. Each row represents a group, the values denote the differences between the treat unit(s) and control unit(s) within the group
# groupMatrix: matrix of units in each group. Each row represents a group, the indices therein denote the treated and control units in the group
# treatmentMatrix: matrix of treated units in each group. Each row represents a group, the indices therein denote the treated units in the group
# controlMatrix: matrix of control units in each group. Each row represents a group, the indices therein denote the control units in the group
matchToGroup = function(matchMatrix, indexControl, indexTreatment, Y) {
  # dominant class in each group, i.e. multiple-to-one or one-to-multiple
  leaderControl = which(apply(matchMatrix, 2, sum) > 1)
  if (length(leaderControl) == 0){
    leaderTreatment = seq(1, dim(matchMatrix)[1])
    numberGroup = length(leaderTreatment); groupSize = c(apply(matchMatrix[leaderTreatment, ], 1, sum)) + 1
  } else {
    leaderTreatment = which(apply(matrix(matchMatrix[, leaderControl], ncol = length(leaderControl)), 1, sum) == 0) 
    numberGroup = length(leaderControl) + length(leaderTreatment); groupSize = c(apply(matrix(matchMatrix[, leaderControl], ncol = length(leaderControl)), 2, sum), apply(matrix(matchMatrix[leaderTreatment, ], nrow = length(leaderTreatment)), 1, sum)) + 1
  }
  groupMatrix = matrix(-1, nrow = numberGroup, ncol = max(groupSize))
  tauMatrix = treatmentMatrix = controlMatrix = matrix(Inf, nrow = numberGroup, ncol = max(groupSize)-1)
  if(length(leaderControl) != 0){
    for(i in 1:length(leaderControl)){
      groupMatrix[i, 1:groupSize[i]] = c(indexControl[leaderControl[i]], indexTreatment[which(matchMatrix[, leaderControl[i]] > 0)])
      tauMatrix[i, 1:(groupSize[i]-1)] = Y[groupMatrix[i,2:groupSize[i]]] - Y[groupMatrix[i,1]]
      treatmentMatrix[i, 1:(groupSize[i]-1)] = groupMatrix[i,2:groupSize[i]] 
      controlMatrix[i, 1:(groupSize[i]-1)] = groupMatrix[i,1]  
    }
  }
  lensLeaderControl = length(leaderControl)
  if(length(leaderTreatment) != 0){
    for(i in (1+lensLeaderControl):(length(leaderTreatment)+lensLeaderControl)){
      groupMatrix[i, 1:groupSize[i]] = c(indexTreatment[leaderTreatment[i-lensLeaderControl]], indexControl[which(matchMatrix[leaderTreatment[i-lensLeaderControl],] > 0)])
      tauMatrix[i, 1:(groupSize[i]-1)] = Y[groupMatrix[i,1]] - Y[groupMatrix[i,2:groupSize[i]]]
      treatmentMatrix[i, 1:(groupSize[i]-1)] = groupMatrix[i,1]  
      controlMatrix[i, 1:(groupSize[i]-1)] = groupMatrix[i, 2:groupSize[i]]
    }
  }
  groupMatrix[which(groupMatrix == 0)] = -1
  
  result = list(); result$groupMatrix = groupMatrix; result$tauMatrix = tauMatrix; result$treatmentMatrix = treatmentMatrix; result$controlMatrix = controlMatrix
  return(result)
}


