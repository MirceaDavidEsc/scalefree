
#' scalefree: measure spatial correlations in collective systems
#'
#' The scalefree package provides functions that can take movement data
#' (velocities) from a collective system, decompose the full velocities into
#' collective and fluctuation components, measure the spatial correlations among
#' the fluctuations, and quantify statistical properties like correlation length
#' and susceptibility.
#'
#' @docType package
#' @name scalefree
NULL


#' Get correlation lengths from profiles
#'
#' This function takes a data frame of correlation profiles (grouped by a certain variable) and determines the zero-intercept (correlation length) for each ccorrelation profile.
#' It is dependent on nesting ahd purrr map functions to iterate over all the profiles. You should install those two packages before using this.
#'
#' @param correlationProfiles A data frame of correlation profiles, grouped by a certain column.
#'
#' @return
#' @export
#'
calcCorrelationLengths = function(groupedProfiles) {

  # Nest the profiles appropriately so that it can work with my functions.
  nestedDist = groupedProfiles %>% select(domain) %>% nest(.key = "domain")
  nestedV = groupedProfiles %>% select(vCorr) %>% nest(.key = "mVCorr")
  nestedD = groupedProfiles %>% select(dCorr) %>% nest(.key = "mDCorr")
  nestedS = groupedProfiles %>% select(sCorr) %>% nest(.key = "mSCorr")
  nestedProfiles = inner_join(nestedV, nestedD) %>% inner_join(nestedS) %>% inner_join(nestedDist)

  # Get the zero-crossing of all the correlation profiles
  nestedProfiles = nestedProfiles %>% mutate(vZero = map2(mVCorr, .y = domain, zerosFromDF),
                                             dZero = map2(mDCorr, .y = domain, zerosFromDF),
                                             sZero = map2(mSCorr, .y = domain, zerosFromDF))

  corrLengths = nestedProfiles %>% select(-domain, -mVCorr, -mDCorr, -mSCorr) %>% unnest()
  return(corrLengths)
}



#' Calculate fluctuation correlations in a full velocity field.
#'
#' @param velocityField A data frame containing all the vectors in a velocity field (full velocities)
#' @param sampleSize How many vectors to use from the original field for pair-wise correlation comparisons.
#'
#' @return The correlation functions/profiles for the vector field.
#' @export
calcCorrelationFunction = function(velocityField, fluctuationField = NULL, sampleSize = 3600) {
  if (dim(velocityField)[1] > sampleSize) {
    sampleRows = sample(1:dim(velocityField)[1], sampleSize)
  } else {
    sampleRows = 1:dim(velocityField)[1]
  }


  # Get pairwise indices for a manageable subset of vectors, and get their distances and speed correlations
  if (is.null(fluctuationField))
    fluctuationField = calculateFluctuationField(velocityField)

  forCalculation = velocityField %>% slice(sampleRows)
  vectorFluctuations = fluctuationField %>% slice(sampleRows)


  distances = getPairwiseDistances(forCalculation[, 1:2])
  speedCorrelations = getSpeedCorrelation(forCalculation[, 3:4])
  velocityCorrelations = getVelocityCorrelation(vectorFluctuations[, 3:4])
  directionCorrelations = getDirectionCorrelation(vectorFluctuations[, 3:4])

  pairwisecorrelations = data.frame(distances, velocityCorrelations, directionCorrelations, speedCorrelations)

  # From the pair-wise correlations, calculate the correlation profiles
  vCorrOut = normCorrFunction(pairwisecorrelations$distances, pairwisecorrelations$velocityCorrelations)
  dCorrOut = normCorrFunction(pairwisecorrelations$distances, pairwisecorrelations$directionCorrelations)
  sCorrOut = normCorrFunction(pairwisecorrelations$distances, pairwisecorrelations$speedCorrelations)

  corrOut = vCorrOut %>% mutate(directionCorr = dCorrOut[, -1], speedCorr = sCorrOut[, -1])
  colnames(corrOut) = c("domain", "vCorr", "dCorr", "sCorr")
  return(corrOut)
}



#' Calculate fluctuation field from a regular velocity field
#'
#' @param frameVectorField A data frame of 4 columns describing points in space with velocities (2D data)
#'
#' @return The fluctuation field once affine movements are removed.
#' @export
calculateFluctuationField = function(frameVectorField) {
  colnames(frameVectorField) = c("X", "Y", "vX", "vY")
  # From the initial positions and instantenous velocity, determine the future position.
  positionsInterp = frameVectorField %>% mutate(X2 = X + vX, Y2 = Y + vY) %>% select(X, Y, X2, Y2)

  # Identify the ideal affine transform that turns the old positions into the new positions
  transformParams = calculateOptimalAffine(positionsInterp)
  paramX = transformParams[[1]]
  paramY = transformParams[[2]]
  paramR = transformParams[[3]]
  paramD = transformParams[[4]]

  # Apply optimal affine transform on future points to get the closest superposition of points.
  realignedPoints = affineTransform(positionsInterp[,3:4], paramX, paramY, paramR, paramD)

  # Subtract past positions from future positions post-transformation to get fluctuation vectors.
  alignedPositions = cbind(positionsInterp, realignedPoints)
  colnames(alignedPositions) = c("X", "Y", "vX", "vY", "newx","newy")
  alignedPositions = alignedPositions %>% mutate(fluctX = newx - X, fluctY = newy - Y) %>%
    select(X, Y, fluctX, fluctY)
  return(alignedPositions)
}




#' Get the collective motion in a vector field.
#'
#' The vectors that define the movement of components in a collective system can be thought of as being composed of two parts: a fluctuation and a collective
#' component.
#'
#' @inheritParams calculateFluctuationField
#'
#' @return A data frame specifying the collective component of the collective movement; i.e. what can be explained by idealized affine transformations.
#' @export
#'
#' @examples
#'
#' @family vector decomposition
#' @seealso \code{\link{calculateFluctuationField} for fluctuation component}
calculateCollectiveField = function(frameVectorField) {
  colnames(frameVectorField) = c("X", "Y", "vX", "vY")
  # From the initial positions and instantenous velocity, determine the future position.
  positionsInterp = frameVectorField %>% mutate(X2 = X + vX, Y2 = Y + vY) %>% select(X, Y, X2, Y2)
  transformParams = calculateOptimalAffine(positionsInterp)
  paramX = -transformParams[[1]]
  paramY = -transformParams[[2]]
  paramR = -transformParams[[3]]
  paramD = 1/transformParams[[4]]
  print(paste(paramX, paramY, paramR, paramD))


  oldPoints = select(frameVectorField, X, Y)
  newPoints = affineTransform(oldPoints, paramX, paramY, paramR, paramD)
  collectiveVs = newPoints - oldPoints
  collectiveField = cbind(oldPoints, collectiveVs)
  colnames(collectiveField) = c("X", "Y", "collectiveX", "collectiveY")
  return(collectiveField)
}


#' Determine affine transformation parameters for best point-wise fit
#'
#' @param positions A N-by-4 data frame containing X1, Y1, X2, Y2, data for 2 sets of 2D Cartesian points.
#'
#' @return Affine x translation, y translation, rotation angle, and dilatation to match points.
#' @export
calculateOptimalAffine = function(positions) {
  # Calculate the optimal affine transformation parameters that minimize the difference between the initial and final positions.
  ## Do a coarse-grained hyperparameter optimization at first. Run 10000 iterations.
  rangeVX = range(positions[,3] - positions[,1])
  rangeVY = range(positions[,4] - positions[,2])
  randX = runif(1000, rangeVX[1], rangeVX[2])
  randY = runif(1000, rangeVY[1], rangeVY[2])
  randR = runif(1000, -pi/2, pi/2)
  randD = runif(1000, 0.5, 2) # A halving or doubling in size

  i = 1
  bestParams = c(0,0,0,1) # Start with no affine transform at all (no translation, rotation, dilatation of 1 [same size])
  bestDeviation = measureDeviation(bestParams, positions)
  while (i < 1000) {
    tempParams = c(randX[i], randY[i], randR[i], randD[i])
    tempDeviation = measureDeviation(tempParams, positions)
    if (tempDeviation < bestDeviation) {
      bestParams = tempParams
      bestDeviation = tempDeviation
    }
    i = i + 1
  }

  ## Use identified hyperparameters as initial condition for optimization, extract the actual optimal parameters.
  transformParams = optim(bestParams, measureDeviation, data = positions)
  return(transformParams$par)
}


#' Measure sum of rms distance between two sets of points after an affine transform.
#'
#' @param transformParams The parameters to use when applying an affine transformation
#' @param data A N-by-4 data frame with columns X and Y of first positions and X and Y of comparison positions.
#'
#' @return The sum of pair-wise distances.
#' @export
measureDeviation = function(transformParams, data) {
  transPoints = affineTransform(data[,3:4], transformParams[1], transformParams[2], transformParams[3], transformParams[4])
  sum(sqrt((data[,1] - transPoints[,1])^2 + (data[,2] - transPoints[,2])^2))
}



#' Perform affine transformation on two-dimensional points.
#'
#'
#'
#' @param points A N-by-2 data frame representing a set of two-dimensional Cartesian points.
#' @param rotate A scalar representing the angle (in radians) by which to rotate the points.
#' @param dilate A scalar representing the amount by which to magnify (or shrink) the point positions relative to the center of mass.
#' @param translateX A scalar representing the amount by which to move the points in the X direction.
#' @param tranlsateY A scalar representing the amount by which to move the points in the Y direction.
#'
#' @return
#' @export
affineTransform = function(points, translateX = 0, translateY = 0, rotate = 0, dilate = 1) {
  colnames(points) = c("X", "Y")
  center = colMeans(points)
  centeredPoints = points %>% mutate(X = X - center[[1]], Y = Y - center[[2]])

  centeredMat = as.matrix(centeredPoints, ncol = 2)
  rotationMatrix = t(matrix(c(cos(rotate), sin(rotate), -sin(rotate), cos(rotate)), byrow=T, nrow=2))
  rotatedMatrix = centeredMat %*% rotationMatrix
  rotatedPoints = as.data.frame(rotatedMatrix)
  dilatedPoints = rotatedPoints*dilate
  colnames(dilatedPoints) = c("newx", "newy")

  translatedPoints = dilatedPoints %>% mutate(newx = newx + center[[1]] + translateX, newy = newy + center[[2]] + translateY)
  colnames(translatedPoints) = c("X", "Y")
  return(translatedPoints)
}

#' Title
#'
#' @param thisFrame
#'
#' @return
#' @export
removeTranslationalComponent = function(thisFrame) {
  ###FUNCTION: Taking a vector field as input, this function removes any mean translational motion captured by the vector field (u,v) and produces (relu,relv) and places the coordinates of the vector field (R,C) in a center-of-mass reference frame (relR,relC).

  require(dplyr)
  returnFrame = mutate(thisFrame,relR=R-mean(R),relC=C-mean(C),relu=u-mean(u),relv=v-mean(v))
  return(returnFrame)
}

#' Title
#'
#' @param thisFrame
#'
#' @return
#' @export
removeRotationalComponent = function(thisFrame) {
  ###FUNCTION: This function removes the rotational movement (from the center-of-mass reference frame) of a vector field. This function assumes the object(s) represented by the vector field rotates as a solid object, such that the angular speed of rotation is linearly proportional to the radial distance from the center of mass.
  require(dplyr)
  #For each position, find the projected velocity that is perpendicular to radial vector
  #Calculate the speed/magnitude of this perpendicular vector.
  returnedFrame = mutate(thisFrame,mag = sqrt(relv^2+relu^2)) %>%
    mutate(magFromCenter = sqrt(relR^2+relC^2)) %>%
    mutate(unitrelR = relR/magFromCenter,unitrelC = relC/magFromCenter) %>%
    mutate(compS = (relR*relu+(-relC)*relv)/magFromCenter) %>%
    #Using the speed and radial distance, calculate the angular velocity for that point.
    mutate(angularS = compS/magFromCenter)

  #Determine the average angular velocity over all vectors.
  averageAngularS = summarise(returnedFrame,meanS=mean(angularS,na.rm = T))
  #Using the average angular velocity and radial distance, calculate the orthogonal vector to be subtracted from each vector.
  finalFrame = mutate(returnedFrame,tosubtu = averageAngularS$meanS*relR,tosubtv = -relC*averageAngularS$meanS) %>%

    mutate(fluctu = relu - tosubtu, fluctv = relv - tosubtv)
  return(finalFrame)
}



#' Measure correlation decay by slope at correlation length
#'
#' This function takes in two vectors that define the correlation profile of any
#' collective system and measures the slope at the zero-intercept (correlation
#' length) by using the 5-point stencil method on the five data points closest
#' to the x-intercept of the correlation profile. Based off of advice from
#' http://www.theanalysisfactor.com/r-tutorial-4/
#'
#' @param distances A vector representing distance measures
#' @param correlations A vector of the corresponding average pair-wise
#'   correlation at that distance.
#'
#' @return
#' @export
correlationSlope = function(distances, correlations) {

  quadraticFit = lm(correlations ~ distances + I(distances^2))
  definePoly = as.polynomial(quadraticFit$coefficients)
  derivativePoly = deriv(definePoly)
  slopeAtZero = predict(derivativePoly, 1)
  return(slopeAtZero)
}




#' Calculate Susceptibility (Maximum Cumulative Correlation)
#'
#' Uses the method
#' developed in Attanasi et al. 2014 to calculate the finite-size susceptibility
#' or maximal cumulative correlation of a collectively moving system. The
#' cumulative correlation is calculated using the trapezoidal (Simpson's) rule
#' for numerical integration. This function filters the correlation profile to
#' only consider values before the first zero-crossing if it detects any
#' negative correlation values.
#'
#' @param distance The domain over which to integrate. If not provided, daata
#'   points are assumed to be evenly spaced with distance 1.
#' @param correlation The correlation at the corresponding distance.
#'
#' @return
#' @export
calculateSusceptibility = function(distance = 1:length(correlation), correlation) {
  # Calculate susceptibility using trapezoidal rule of the curve.
  # Filter at the first zero crossing if not done already
  if (min(correlation) < 0) {
    zeroCrossing = getFirstZeroCrossing(correlation)
    distance = distance[1:zeroCrossing]
    correlation = correlation[1:zeroCrossing]
  }

  susceptibility = trapz(distance, correlation)
  return(susceptibility)
}



#' Title
#'
#' Calculates the pair-wise correlations between all vectors in a velocity field.
#' A vector can be provided to allow within- and between-group comparisons (not implemented).
#' Vectors can be subsampled to deal with quadratic memory complexity of calculating pair-wise correlations.
#'
#' @param frameData a N-by-5 data frame storing N velocity vectors with column names R, C, u, v, and Frame.
#' @param memoryMaxVectors the number of vectors to subset the data frame by in order to avoid memory issues.
#'
#' @return The pairwise velocity, direction, and speed correlation for all pairs of vectors found in frameData
#' @export
getPairwiseCorrelations = function(frameData, memoryMaxVectors = 3000) {
  # From positions and full velocities, get the positions of each vector relative to the center and
  # velocities relative to the mean.
  colnames(frameData) = c("X", "Y", "vX", "vY")
  velocityFrame = frameData %>% mutate(speed = sqrt(vX^2+vY^2))
  velocityFrame = cbind(velocityFrame, calculateFluctuationField(velocityFrame[, 1:4]))

  velocityFrame = velocityFrame %>% select(relR,relC,speed,fluctu,fluctv) %>%
    mutate(unitfluctu = fluctu/sqrt(fluctu^2+fluctv^2), unitfluctv = fluctv/sqrt(fluctu^2+fluctv^2)) %>%
    mutate(speedFluctuation = speed - mean(speed))

  # Produce all index pairs, excluding repeat entries by only getting upper-triangle.
  numVectors = dim(velocityFrame)[1]
  pairwiseInds = expand.grid(1:numVectors, 1:numVectors)
  inds1 = matrix(pairwiseInds[,1], nrow = numVectors)
  inds2 = matrix(pairwiseInds[,2], nrow = numVectors)
  indsi = inds1[upper.tri(inds1,T)]
  indsj = inds2[upper.tri(inds2,T)]
  pairwiseInds = cbind(indsi,indsj)

  # Calculate Euclidean distances
  coordinates = cbind(velocityFrame[pairwiseInds[,1],1:2], velocityFrame[pairwiseInds[,2],1:2])
  distances = sqrt((coordinates[,1] - coordinates[,3])^2 + (coordinates[,2] - coordinates[,4])^2)

  velocityCorrelation = velocityFrame[pairwiseInds[,1], 4]*velocityFrame[pairwiseInds[,2], 4] + velocityFrame[pairwiseInds[,1], 5]*velocityFrame[pairwiseInds[,2], 5]

  directionCorrelation = velocityFrame[pairwiseInds[,1], 6]*velocityFrame[pairwiseInds[,2], 6] + velocityFrame[pairwiseInds[,1], 7]*velocityFrame[pairwiseInds[,2], 7]

  speedCorrelation = velocityFrame[pairwiseInds[,1], 8]*velocityFrame[pairwiseInds[,2], 8]

  pairwiseCorrelations = data.frame(distance = distances, velocityCorrelation = velocityCorrelation, directionalCorrelation = directionCorrelation, speedCorrelation = speedCorrelation)
  return(pairwiseCorrelations)
}

getPairwiseDistances = function(coordinates) {
  coordinates = as.matrix(coordinates)
  distances = as.matrix(dist(coordinates, diag=T))
  distances = distances[upper.tri(distances, diag=T)]
  return(distances)
}


getVelocityCorrelation = function(fluctuationVectors) {
  fluctuationVectors = as.matrix(fluctuationVectors)
  tFluctuationVectors = t(fluctuationVectors)
  velocityCorrelations = fluctuationVectors %*% tFluctuationVectors
  velocityCorrelations = velocityCorrelations[upper.tri(velocityCorrelations, diag=T)]
  return(velocityCorrelations)
}

getDirectionCorrelation = function(fluctuationVectors) {
  colnames(fluctuationVectors) = c("vX", "vY")
  unitVectors = fluctuationVectors %>% mutate(modulus = sqrt(vX^2 + vY^2)) %>%
    mutate(uvX = vX/modulus, uvY = vY/modulus) %>%
    select(uvX, uvY) %>% as.matrix()
  tUnitVectors = t(unitVectors)
  directionCorrelations = unitVectors %*% tUnitVectors
  directionCorrelations = directionCorrelations[upper.tri(directionCorrelations, diag=T)]
  return(directionCorrelations)
}

getSpeedCorrelation = function(velocityVectors) {
  colnames(velocityVectors) = c("vX", "vY")
  speedFlucts = velocityVectors %>% mutate(speed = sqrt(vX^2 + vY^2)) %>%
    mutate(speedFluct = speed - mean(speed)) %>% select(speedFluct) %>% as.matrix
  tSpeedFlucts = t(speedFlucts)
  speedCorrelations = speedFlucts %*% tSpeedFlucts
  speedCorrelations = speedCorrelations[upper.tri(speedCorrelations, diag=T)]
  return(speedCorrelations)
}



#' Title
#'
#' @param numVectors
#'
#' @return
#' @export
#'
#' @keywords internal
getPairwiseIndices = function(numVectors) {
  # Produce all index pairs, excluding repeat entries by only getting upper-triangle.
  pairwiseInds = expand.grid(1:numVectors, 1:numVectors)
  inds1 = matrix(pairwiseInds[,1], nrow = numVectors)
  inds2 = matrix(pairwiseInds[,2], nrow = numVectors)
  indsi = inds1[upper.tri(inds1,T)]
  indsj = inds2[upper.tri(inds2,T)]
  pairwiseInds = cbind(indsi,indsj)
}


#' Title
#'
#' @param frameData
#'
#' @return
#' @export
normCorrFunction = function(domain,correlations, resolution=1, nknots = 50) {
  splineVal = bigspline(domain, correlations, type="cub", nknots)

  # Define prediction range and intervals
  maxDist = max(domain)
  predictDomain = seq(from = 0, to = maxDist, by = resolution)
  predictOut = predict(splineVal, predictDomain)
  normalizedOut = predictOut/predictOut[[1]]
  # plot(predictDomain, normalizedOut)
  # head(predictDomain[normalizedOut < 0], 1)
  # head(roughPredict$domain[roughPredict$meanCorr < 0], 1)
  corrDF = data.frame(domain = predictDomain, out = normalizedOut)
  return(corrDF)
}



#' Get first zero crossing of a function
#'
#' @param thisVector A numeric vector
#' @param distMapping A numeric vector
#'
#' @return
#' @export
getFirstZeroCrossing = function(thisVector,distMapping=1:length(thisVector)) {
  signDiff = diff(sign(thisVector))
  crossing = first(which(signDiff != 0))
  if (is.na(crossing))
    crossing=NA
  else {
    distances = distMapping[(crossing):(crossing+1)]
    values = thisVector[crossing:(crossing+1)]
    crossing = predict(lm(distances ~ values),data.frame(values=c(0)))
  }
  #Fit a linear model between the crossing point data line and interpolate the crossing point.
  return(crossing)
}


#' Identify the points of a collective that are part of the boundary
#'
#' This function takes a data frame of points and identifies the points that form the boundary through alpha-hull methods.
#'
#'
#' @param X X coordinate
#' @param Y Y coordinate
#' @param marginPercent The percent of the marginal points (from center radius) to be considered boundary
#' @param ballSize Parameter specifying ball size that represents alpha-hull.
#'
#' @return
#' @export
identifyBoundary = function(X, Y, marginPercent, ballSize) {
  centered = cbind(X - mean(X), Y - mean(Y))
  jitterPoints = as.data.frame(apply(centered,2, jitter))

  hull = ahull(jitterPoints$V1, jitterPoints$V2, alpha=ballSize)

  inahull(hull, cbind(jitterPoints$V1/(1-marginPercent), jitterPoints$V2/(1-marginPercent)))
}


#' Remove boundary points using alpha-hull
#'
#' This function takes a data frame of points and identifies the points that form the boundary through alpha-hull methods, removing the points.
#'
#'
#' @param pointsDF 4-column data frame that specifies points in space and their velocities.
#'
#' @return
#' @export
removeBoundary = function (pointsDF)
{
  pointsDF$bulk = scalefree::identifyBoundary(pointsDF$X, pointsDF$Y, 0.1, 5)
  pointsDF = pointsDF[pointsDF$bulk,]
  pointsDF = pointsDF[,1:4]
  return(pointsDF)
}


#' Standardize Placozoa data frames.
#'
#' @param dataFrame
#'
#' @return
#' @export
standardizePlacozoaVectors = function(dataFrame) {
  returnDF = select(dataFrame, X = V2, Y = V1, vX = V4, vY = V5, Frame = V3)
  return(returnDF)
}
