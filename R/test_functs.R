# Test the affine transformation capabilities.
test_affineNoChange = function(initialPoints) {
  testTransform = affineTransform(initialPoints, 0, 0, 0, 1)
  offset = testTransform - initialPoints
  finalPoints = cbind(initialPoints, offset)
  colnames(finalPoints) = c("X", "Y", "vX", "vY")
  return(max(abs(c(range(finalPoints$vX), range(finalPoints$vY)))) == 0)
}
test_affineNoChange(initialPoints)


# Try a translation in either direction
test_affineTranslate = function(initialPoints) {
  testTransform = affineTransform(initialPoints, 3, -2, 0, 1)
  offset = testTransform - initialPoints
  finalPoints = cbind(initialPoints, offset)
  colnames(finalPoints) = c("X", "Y", "vX", "vY")
  quiverPlot(finalPoints, 1)
  write_rds(finalPoints, "D:/Mircea/ScriptLibraries/R/Collective/translate.rds")
}

test_affineRotate = function(initialPoints) {
  # Try a rotation
  testTransform = affineTransform(initialPoints, 0, 0, pi/16, 1)
  offset = testTransform - initialPoints
  finalPoints = cbind(initialPoints, offset)
  colnames(finalPoints) = c("X", "Y", "vX", "vY")
  quiverPlot(finalPoints, 1)
  write_rds(finalPoints, "D:/Mircea/ScriptLibraries/R/Collective/rotate.rds")
}


test_affineDilate = function(initialPoints) {
  # Try a dilatation
  testTransform = affineTransform(initialPoints, 0, 0, 0, 1.3)
  offset = testTransform - initialPoints
  finalPoints = cbind(initialPoints, offset)
  colnames(finalPoints) = c("X", "Y", "vX", "vY")
  quiverPlot(finalPoints, 1)
  write_rds(finalPoints, "D:/Mircea/ScriptLibraries/R/Collective/dilate.rds")
}

# Can you reverse an affine transformation to reproduce the original point positions?
test_reverseAffine = function(initialPoints) {
  testTransform = affineTransform(initialPoints, 5, 5, 0, 1.05)
  offset = testTransform - initialPoints
  finalPoints = cbind(initialPoints, offset)
  colnames(finalPoints) = c("X", "Y", "vX", "vY")
  quiverPlot(finalPoints, 1)
  reverseTransform = affineTransform(testTransform, -5, -5, 0, 1/1.05)
  return(max(abs(initialPoints - reverseTransform)) < 10^-10)
}
test_reverseAffine(initialPoints)

# Test whether measurement of deviation is a reasonable criterion for matching points.
test_deviationReversed = function(initialPoints) {
  testTransform = affineTransform(initialPoints, 5, 5, pi/100, 1.005)
  transformedPoints = cbind(initialPoints, testTransform)
  colnames(transformedPoints) = c("X1", "Y1", "X2", "Y2")
  initialDeviation = measureDeviation(c(0,0,0,1), transformedPoints)
  reverseDeviation = measureDeviation(c(-5,-5,-pi/100, 1/1.005), transformedPoints)
  return(initialDeviation > reverseDeviation & reverseDeviation < 10^-4)
}
test_deviationReversed(initialPoints)

# Try to find the optimal reverse affine transform.
test_optimalReverseAffine = function(initialPoints) {
  affineParams = c(5, 5, pi/100, 1.005)
  testTransform = affineTransform(initialPoints, affineParams[1], affineParams[[2]], affineParams[[3]], affineParams[[4]])
  velocity = testTransform - initialPoints
  testVelocity = cbind(initialPoints, velocity)
  colnames(testVelocity) = c("X", "Y", "vX", "vY")
  quiverPlot(testVelocity, 1)
  positionsInterp = testVelocity %>% mutate(X2 = X + vX, Y2 = Y + vY) %>% select(X1 = X, Y1 = Y, X2, Y2)
  optimalAffine = calculateOptimalAffine(positionsInterp)

  reversal = c(affineParams[[1]] + optimalAffine[[1]], affineParams[[2]] + optimalAffine[[2]], affineParams[[3]] + optimalAffine[[3]], affineParams[[4]]*affineParams[[4]])
  return(reversal - c(0,0,0,1) < 0.05)
}
test_optimalReverseAffine(initialPoints)



test_retrieveFluctuations = function(testField) {
  testFlucts = calculateFluctuationField(testField)
  realCollective = testField %>% inner_join(testFlucts) %>% mutate(colX = vX - fluctX, colY = vY - fluctY) %>%
    select(X, Y, colX, colY)

  (fullVF = quiverPlot(testField, 1, colormapped = F))
  (fullFF = quiverPlot(testFlucts, 1, colormapped = F))
  (fullCF = quiverPlot(realCollective, 1, colormapped = F))
  max(sqrt(testFlucts$fluctX^2 + testFlucts$fluctY^2))
  max(sqrt(testField$vX^2 + testField$vY^2))
  testField %>% mutate(speed = sqrt(vX^2 + vY^2)) %>% ggplot(aes(speed)) + geom_histogram()
  testFlucts %>% mutate(speed = sqrt(fluctX^2 + fluctY^2)) %>% ggplot(aes(speed)) + geom_histogram()
  realCollective %>% mutate(speed = sqrt(colX^2 + colY^2)) %>% ggplot(aes(speed)) + geom_histogram()
}


#' Title
#'
#' @param testField
#'
#' @return
#' @export
#'
#' @examples
test_retrieveCollective = function(testField) {
  fluctField = calculateFluctuationField(testField)
  collectiveField = calculateCollectiveField(testField)
  remainderField = testField %>% mutate(remainderX = vX - fluctField$fluctX, remainderY = vY - fluctField$fluctY) %>%
    mutate(consistX = remainderX - collectiveField$collectiveX, consistY = remainderY - collectiveField$collectiveY)
}

#testField = read_rds("exampleField.rds")
#test_optimalReverseReal(testField)


#ahullVectors = read_rds("D:/Mircea/Projects/RESEARCH/InternalCoordination/ahullVectors.rds")
#randomSample = ahullVectors %>% group_by(folder, Frame) %>% nest(.key = "vectorField") %>% sample_n(20)
# ahullVectors = read_rds("D:/Mircea/Projects/RESEARCH/InternalCoordination/ahullVectors.rds") %>% filter(folder == "2013-07-01-run_1")
# ahullVectors = ahullVectors %>% group_by(folder, Frame) %>% nest()
# testField = ahullVectors$data[[1]]
#
#
# #write_rds(sampleField, "D:/Mircea/ScriptLibraries/R/Collective/real.rds")
# ggplot(sampleField, aes(X, Y)) + geom_point() + coord_fixed()
# sampleField = select(sampleField, -Frame)
# quiverPlot(sampleField, 1, colormapped = T)
# fluctuations = calculateFluctuationField(sampleField)
# quiverPlot(fluctuations, 1, colormapped = T)
