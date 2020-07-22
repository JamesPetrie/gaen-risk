
require(plyr)
require(data.table)


probAsymptomatic = 0.15
asymptomaticSheddingDurationCdf = c(0.054054054,0.094594595,0.12162162,0.148648649,0.189189189,0.21621621,0.256756757,0.310810811,0.378378378,0.432432432,0.486486486,0.540540541,0.554054054,0.621621622,0.662162162,0.662162162,0.689189189,0.702702703,0.756756757,0.783783784,0.824324324,0.851351351,0.864864865,0.891891892,0.905405405,0.945945946,0.959459459,0.959459459,0.959459459,0.972972973,0.986486,0.986486486,0.986486486,0.986486486,0.986486486,0.986486486,0.986486486,0.986486486,0.986486486,1)
incubationPdf = c(0,4E-05,0.011842,0.088541,0.181965,0.207344,0.174797,0.123761,0.081488,0.051057,0.031469,0.018734,0.011235,0.006786,0.00422,0.002518,0.001626,0.000978,0.000592,0.000364,0.000231,0.00014,0.000093,0.000062,0.00004,0.000025,0.000017,0.000011,0.000008)
incubationCdf = cumsum(incubationPdf)
#discountSchedule = probAsymptomatic*(1-asymptomaticDoneSheddingCdf) + (1- probAsymptomatic)*(1 - symptomaticSymptomStartCdf)


discountSchedule = c(1,0.99998,0.994059,0.9497885,0.858806,0.755134,0.660103392,0.586894919,0.533407703,0.494373128,0.463039432,0.438587189,0.416241392,
0.393207216,0.367287169,0.340932595,0.313997176,0.286927378,0.265554932,0.240765331,0.217746365,0.201059905,0.185435372,0.172969757,
0.156689676,0.141405162,0.124388311,0.108319101,0.094752304,0.081300662,0.070016527,0.056302622,0.044703284,0.036214683,0.030309399,
0.024554527,0.018833743,0.014769669)

testSensitivity = c(0.000161329, 0.015159562,0.101394973, 0.260243002, 0.429042926, 0.565863684, 0.659973336, 0.715982881, 0.742706949, 
                    0.748928452, 0.739946231, 0.720327199, 0.693023891, 0.66030651, 0.624367052, 0.586499145, 0.548150589, 0.510691298,
                    0.474791693, 0.44082903, 0.409543316)

doseResponseLambda = 1.71E-05 # transmission probability of single particle
attenuationDurationExpectedParticles = c(2.0182978,1.1507629,0.6651614) # with units of [particles per minute] for a normalized transmission risk

attenuationDurationWeights = doseResponseLambda*attenuationDurationExpectedParticles

riskLevels = c(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5)/100

riskConfig = list(
  attenuationDurationThresholds = c(50, 70),
  
  attenuationDurationWeights = attenuationDurationWeights,
  
  transmissionRiskValuesForLevels = c(0, 10, 10^(1+(2:7)/6)),
  
  discountSchedule = discountSchedule
)

#' Get the default risk configuration parameters
#' @keywords 
#' @export
#' @examples
#' getRiskConfig()
getRiskConfig = function(){
  return(riskConfig)
}

#' Randomly generate exposure data
#' @keywords 
#' @export
#' @examples
#' generateRandomExposures()
generateRandomExposures = function(count, currentDay, version){
  return(llply(1:count, function(x){
    list(attenuationDurations = c(sample(5*(1:6), 1),sample(5*(1:6), 1),sample(5*(1:6), 1)), transmissionRiskLevel = 6, day = sample((currentDay-14):currentDay, 1), ENVersion = version)
  }))
}

computeAttenuationScore = function(attenuation, riskConfig){
  roundedAttenuation = round(attenuation)
  return(riskConfig$attenuationScores[roundedAttenuation])
}

#' Compute risk
#' @keywords 
#' @export
#' @examples
#' computeExposureInfoRisk()
computeExposureInfoRisk = function(exposureInfo, riskConfig){
  score = as.numeric(exposureInfo$attenuationDurations %*% riskConfig$attenuationDurationWeights) # dot product
  transmissionFactor = riskConfig$transmissionRiskValuesForLevels[exposureInfo$transmissionRiskLevel + 1]
  score = score*transmissionFactor
  return((1 - exp(-score)))
}

# EN version 1.5
computeExposureWindowRisk = function(exposureWindow, riskConfig){
  score = 0
  for(scanResponse in exposureWindow$scans){
    score = score + scanResponse$duration*computeAttenuationScore(scanResponse$attenuation, riskConfig)
  }
  return((1 - exp(-score)))
}

combineRisks = function(risks){
  inverseRisks = 1 - risks
  return(1 - prod(inverseRisks))
}



getDayRisks = function(exposures){
  riskDt = data.table(Day = 0:60, DayTransmissionRisk = 0, InfectedRisk = 0)
  
  for(exposure in exposures){
    if(exposure$ENVersion == "1.5"){
      newRisk = computeExposureWindowRisk(exposure, riskConfig)
    }else{ # EN Version 1.1
      newRisk = computeExposureInfoRisk(exposure, riskConfig)
    }
    riskDt[Day == exposure$day, DayTransmissionRisk := combineRisks(c(DayTransmissionRisk, newRisk))]
  }
  
  
  for(day in riskDt$Day){
    if(riskDt[day == Day]$DayTransmissionRisk > 0){
      newInfectedRisk = riskDt[day == Day]$DayTransmissionRisk*riskConfig$discountSchedule
      riskDt[Day >= day & Day < (day + length(riskConfig$discountSchedule)), InfectedRisk := 1 - (1-InfectedRisk)*(1-newInfectedRisk)]
    }
  }
  return(riskDt)
}

getRiskLevel = function(risk){
  level = 0
  for(riskLevel in riskLevels){
    if(risk < riskLevel){
      return(level)
    }
    level = level + 1 
  }
  return(length(riskLevels))
}

