# the function computes the match that minimizes the average distance of the treatment-control pairs. (A treated/control unit may be used multiple times.)
# input:
# distMatrix: matrix of distances. Each row represents a treated unit, and each column represents a control unit, the ij-th element represents the distance between the i-th treated unit and the j-th control unit.
# maxDegreeTreatment: the maximal number of control units matched to a treated unit
# maxDegreeControl: the maximal number of treated units matched to a control unit
# output:
# valueF: the total number of treatment-control pairs constructed. A treated/control unit may be used multiple times
# totalCost: the total distance between all treatment-control pairs constructed
# solution: matrix of treatment-control pairs: the ij-th element is 1 if treated unit i and control unit j are matched, otherwise 0
# library("gurobi")
factMatch = function(distMatrix, maxDegreeTreatment = 5, maxDegreeControl = 5) {
  # pre-processing
  nTreatment = dim(distMatrix)[1]; nControl = dim(distMatrix)[2]; n = nTreatment + nControl
  
  # create causal network
  edgelist = data.frame(
    upperCapacity = c(rep(maxDegreeTreatment, nTreatment), rep(1, nTreatment * nControl), rep(maxDegreeControl, nControl)),
    lowerCapacity = c(rep(1, nTreatment), rep(0, nTreatment * nControl), rep(1, nControl)),
    cost = c(rep(0, nTreatment), as.vector(t(distMatrix)), rep(0, nControl))
  )
  
  # create capacity constraints for gurobi solver
  constraintsMatrix = createConstraintGurobi(n = n, nTreatment = nTreatment, nControl = nControl, upperCapacity = edgelist$upperCapacity, lowerCapacity = edgelist$lowerCapacity, valueF = 0)
  
  # solve the average minimal cost problem with binary search
  # initial values
  lowF = max(nTreatment, nControl); highF = min(nTreatment * maxDegreeTreatment, nControl * maxDegreeControl, nTreatment * nControl); epsilon = 1
  
  # lower boundary
  valueF = lowF
  tempFlow = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = valueF); totalCost = tempFlow$totalCost
  totalCostPlus = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = (valueF+epsilon))$totalCost
  deltaPlus = totalCostPlus - totalCost
  comparePlus = sign((deltaPlus/epsilon) - (totalCost/valueF))
  if(comparePlus >= 0){
    result = list(); result$valueF = valueF; result$totalCost = totalCost; result$solution = tempFlow$solution; print(c("lowF", lowF, highF, valueF, "comparePlus:", comparePlus)); return(result)
  } else {lowF = lowF + 1}
  
  # upper boundary
  valueF = highF
  tempFlow = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = valueF); totalCost = tempFlow$totalCost
  totalCostMinus = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = (valueF-epsilon))$totalCost
  deltaMinus = totalCost - totalCostMinus
  compareMinus = sign((deltaMinus/epsilon) - (totalCostMinus/(valueF - epsilon)))
  if(compareMinus <= 0){
    result = list(); result$valueF = valueF; result$totalCost = totalCost; result$solution = tempFlow$solution; print(c("highF", lowF, highF, valueF, "compareMinus:", compareMinus)); return(result)
  } else {highF = highF - 1}
  
  while(TRUE) {
    # solve min cost flow problem with |f| = midF
    midF = floor((lowF + highF)/2)
    valueF = midF
    tempFlow = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = valueF)
    totalCost = tempFlow$totalCost
    totalCostMinus = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = (valueF-epsilon))$totalCost
    totalCostPlus = minCostFlow(n = n, nTreatment = nTreatment, nControl = nControl, cost = edgelist$cost, constraintsMatrix = constraintsMatrix, valueF = (valueF+epsilon))$totalCost
    deltaPlus = totalCostPlus - totalCost; deltaMinus = totalCost - totalCostMinus
    comparePlus = sign((deltaPlus/epsilon) - (totalCost/valueF))
    compareMinus = sign((deltaMinus/epsilon) - (totalCostMinus/(valueF - epsilon)))
    
    # highF == lowF + 1: candidates are exhausted
    if(highF == (lowF + 1)){
      if(comparePlus == 0 | (comparePlus * compareMinus < 0)){print(c("highF = lowF + 1", lowF, highF, valueF, "comparePlus:", comparePlus, "compareMinus:", compareMinus)); break
      } else {valueF = highF; print(c("highF = lowF + 1", lowF, highF, valueF, "comparePlus:", comparePlus, "compareMinus:", compareMinus)); break}
    }
    
    # binary search
    if(comparePlus == 0 | (comparePlus *  compareMinus < 0)) {print(c(lowF, highF, valueF, "comparePlus:", comparePlus, "compareMinus:", compareMinus)); break
    } else if(comparePlus > 0) {highF = midF; midF = floor((lowF + highF)/2)
    } else {lowF = midF; midF = floor((lowF + highF)/2)}
  }
  result = list(); result$valueF = valueF; result$totalCost = totalCost; result$solution = tempFlow$solution
  return(result)
}

# helper functions for function "factMatch"

# create constraints for an average distance minimization matching problem (causal network)
# input:
# n: sample size
# nTreatment: number of treated units
# nControl: number of control units
# upperCapacity: vector of maximal edge capacities
# lowerCapacity: vector of minimal edge capacities
# valueF: the total number of treatment-control pairs constructed. A treated/control unit may be used multiple times
# output: a system of linear equations
# lhs: coefficient matrix 
# dir: symbols of equalities/inequalities
# rhs: constant terms
createConstraintGurobi = function(n, nTreatment, nControl, upperCapacity, lowerCapacity, valueF) {
  # preprocessing
  numberEdges = nTreatment * nControl + n
  constraints = list(lhs = NA, dir = NA, rhs = NA)
  
  # capacity constraints:
  # flow through each edge should not be larger than capacity upper limit
  iList = c(seq(1, numberEdges))
  jList = c(seq(1, numberEdges))
  xList = rep(1, numberEdges)
  constraints$dir = rep('<=', times = numberEdges)
  constraints$rhs = upperCapacity
  
  # flow through each edge should not be smaller than capacity lower limit
  iList = c(iList, seq((numberEdges + 1), 2 * numberEdges))
  jList = c(jList,  seq(1, numberEdges))
  xList = c(xList, rep(1, numberEdges))
  constraints$dir = c(constraints$dir, rep('>=', times = numberEdges))
  constraints$rhs = c(constraints$rhs, lowerCapacity)
  
  # node flow constraints:
  # for each node (source and sink edges excluded), the sum of all inputs and all outputs should be zero 
  iList = c(iList, seq((numberEdges * 2 + 1), (numberEdges * 2 + nTreatment)), rep(seq((numberEdges * 2 + 1), (numberEdges * 2 + nTreatment)), rep(nControl, nTreatment)))
  iList = c(iList, rep(seq((numberEdges * 2 + nTreatment + 1), (numberEdges * 2 + nTreatment + nControl)), (nTreatment + 1)))
  jList = c(jList,  seq(1, (numberEdges - nControl)))
  jList = c(jList, seq((nTreatment + 1), numberEdges))
  xList = c(xList, rep(1, nTreatment), rep(-1, nTreatment * nControl), rep(1, nTreatment * nControl), rep(-1, nControl))
  constraints$dir = c(constraints$dir, rep('=', times = n))
  constraints$rhs = c(constraints$rhs, rep(0, times = n))
  
  # source and the sink flow constraints:
  # for the source and the sink, all outbound nodes and all inbound nodes should be equal to the sum of flow through the network
  iList = c(iList, rep((numberEdges * 2 + n + 1), nTreatment), rep((numberEdges * 2 + n + 2), nControl))
  jList = c(jList,  seq(1, nTreatment), seq((numberEdges - nControl + 1), numberEdges))
  xList = c(xList, rep(-1, nTreatment), rep(1, nControl))
  mySparse = sparseMatrix(i = iList, j = jList, x = xList)
  constraints$lhs = mySparse
  constraints$dir = c(constraints$dir, rep('=', times = 2))
  constraints$rhs = c(constraints$rhs, -valueF, valueF)
  
  return(constraints)
}

# solve a minimal cost flow problem
# input:
# n: sample size
# nTreatment: number of treated units
# nControl: number of control units
# cost: vector of egde costs
# constraintsMatrix: capacity constraints for a causal network created by function "createConstraintGurobi"
# valueF: the total number of treatment-control pairs constructed. A treated/control unit may be used multiple times
# output:
# solution: matrix of treatment-control pairs: the ij-th element is 1 if treated unit i and control unit j are matched, otherwise 0
# totalCost: the total distance between all treatment-control pairs constructed
minCostFlow = function(n, nTreatment, nControl, cost, constraintsMatrix, valueF) {
  # update the constraints
  constraintsMatrix$rhs[c(length(constraintsMatrix$rhs) - 1, length(constraintsMatrix$rhs))] = c(-valueF, valueF)
  
  model <- list()
  model$A = constraintsMatrix$lhs
  model$obj = cost
  model$modelsense = 'min'
  model$rhs = constraintsMatrix$rhs
  model$sense = constraintsMatrix$dir
  params <- list(OutputFlag = 0)
  resultGurobi <- gurobi(model, params) # solver from gurobi
  
  pairMatrix = matrix(0, nrow = nControl, ncol = nTreatment)
  pairMatrix[which(resultGurobi$x[-c(seq(1, nTreatment),seq((nTreatment*nControl+nTreatment+1),(nTreatment*nControl+nTreatment+nControl)))] > 0)] = 1
  pairMatrix = t(pairMatrix)
  
  result = list(); result$solution = pairMatrix; result$totalCost = resultGurobi$objval
  return(result)
}




