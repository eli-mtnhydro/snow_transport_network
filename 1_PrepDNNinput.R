# Copyright Eli Boardman 2024
# All Rights Reserved

library(terra)

setwd("C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/")

################################################################################
# SWE map

# SWE map from 3 m ASO depth + Mountain Hydrology density
# Aggregated to 100 m
# DNN imputation outside lidar survey extent - R^2 = 0.85
# WindRiverRange_SWE_100m_2024-05-31_Updated_Jul-27-2024_Extrapolated_Oct-20-2024.tif
# Manually clipped to desired extent around Torrey Creek watershed to cover snow area

SWEmap = rast("SWEreference.tif")

plot(SWEmap)

################################################################################
# Prepare other data

# Pre-merged SRTM data from other project
DEMfull = rast("F:/MountainHydrologyLLC/BoR_SnowForecasting/Glaciers/GlacierPaper/Data/DEM_30m_FullExtent.tif")
DEM = project(DEMfull,SWEmap,method="bilinear")
writeRaster(DEM,"DEM.tif")

# Imputed broad-band albedo
AlbedoFull = rast("F:/MountainHydrologyLLC/BoR_SnowForecasting/ASOflights/SnowOn_24-05-31/WindRiverRange_AlbedoBB_100m_2024-05-31_Extrapolated_Oct-21-2024.tif")
Albedo = project(AlbedoFull,SWEmap,method="bilinear") / 100
writeRaster(Albedo,"AlbedoReference.tif",overwrite=TRUE)

################################################################################
# Setup and calibrate DHSVM, then run counterfactual without pattern
# (Handled by standalone DHSVM pipeline)

##########
# Need to re-create snow redistribution map for entire gridMET domain, not just model area

SnowMapLocal = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/NormalSnowPatternMap.tif")
SnowPattern = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/NormalSnowPatternMapScaled.tif")
DEMext = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/NEDdownload_NED_1.tif")
DEMext = project(DEMext,SnowMapLocal)

SnowLowerElev = 3000 # m
SnowUpperElev = 3300 # m
UniformExtraPctMean = 0.01 # Fraction of mean value added to all cells

# Only consider area above upper snowline
UpperElevPtrs = which(DEMext[,] > SnowUpperElev)

# Rescale relative to mean of snowy area
MeanSnowVal = mean(unlist(SnowMapLocal[UpperElevPtrs]))
SnowMapScaled = SnowMapLocal / MeanSnowVal

SnowMapWeighted = SnowMapScaled

# Use constant mean value below snowline
SnowMapWeighted[which(DEMext[,] < SnowLowerElev)] = MeanSnowVal

# Use linear ramp between snow lines
IntermediatePtrs = which(DEMext[,] > SnowLowerElev & DEMext[,] < SnowUpperElev)
IntermediateWeights = (unlist(DEMext[IntermediatePtrs]) - SnowLowerElev) / (SnowUpperElev - SnowLowerElev)
SnowMapWeighted[IntermediatePtrs] = SnowMapScaled[IntermediatePtrs] * IntermediateWeights + MeanSnowVal * (1 - IntermediateWeights)

# Preserve outlying snow that may exist at lower elevations
SnowMapWeighted[,] = pmax(SnowMapWeighted[,], SnowMapScaled[,])

SnowMapWeighted = SnowMapWeighted + MeanSnowVal * UniformExtraPctMean

SnowMapWeightedClipped = project(SnowMapWeighted, SWEmap)

plot(SnowMapWeightedClipped - SnowPattern) # Same

writeRaster(SnowMapWeighted,paste0("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/NormalSnowPatternMapScaled",
                                     "_FullStationArea.tif"),overwrite=TRUE)

##########
# Smooth measured SWE map to capture large-scale orographic patterns
# Based on snow pattern redistribution map used for DHSVM calibration
WindowSize = 21
SnowPatternSmooth = focal(focal(SnowMapWeighted,
                                w=WindowSize,fun="mean",expand=TRUE),
                          w=WindowSize,fun="mean",expand=TRUE)
SnowPatternSmoothClipped = project(SnowPatternSmooth,SWEmap,method="near")

plot(SnowPatternSmoothClipped)

writeRaster(SnowPatternSmoothClipped,paste0("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/NormalSnowPatternMapScaled",
                                     "_Smoothed2x",WindowSize,".tif"),overwrite=TRUE)
writeRaster(SnowPatternSmoothClipped,paste0("E:/DHSVM_MountainHydrology/2_Calibration/Watersheds/TorreySnowRegion/ScenarioRuns/",
                                     "SnowPattern_Smoothed2x",WindowSize,".asc"),overwrite=TRUE)

##########
# Create station-normal files for each counterfactual

GridMetStationMap = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/GridMetStationMap.tif")
GridMetStationStats = read.csv("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/GridMetStationStats.csv")[,-1]

for (i in 1:length(GridMetStationStats$StationID)){
  
  StationID = GridMetStationStats$StationID[i]
  StationPtrs = which(GridMetStationMap[,]==StationID)
  StationMean = mean(unlist(SnowPatternSmooth[StationPtrs]))
  
  OutFilename = paste0("E:/DHSVM_MountainHydrology/2_Calibration/Watersheds/TorreySnowRegion/CalibrationRuns/gridforcing",
                       "_Smoothed2x",WindowSize,"/",
                       GridMetStationStats$StationName[i],".snowpattern")
  writeLines(paste(round(StationMean,3), collapse=" "), OutFilename)
}

####################################################################################################################################
# Run DHSVM scenarios

##########
# Account for difference between accumulation and measured SWE

SWEmapDHSVM = rast(paste0("E:/DHSVM_MountainHydrology/2_Calibration/Watersheds/TorreySnowRegion/ScenarioRuns/",
                            "SnowPattern","/Map.Snow.Swq.asc"))
crs(SWEmapDHSVM) = crs(SWEmap)

AccumMapDHSVM = rast(paste0("E:/DHSVM_MountainHydrology/2_Calibration/Watersheds/TorreySnowRegion/ScenarioRuns/",
                            "SnowPattern","/Map.Snow.SumAccum.asc"))
crs(AccumMapDHSVM) = crs(SWEmap)

DHSVMdiffAccumSWE = AccumMapDHSVM - SWEmapDHSVM

plot(DHSVMdiffAccumSWE,main="DHSVM with snow pattern: Accum - SWE")

AccumImplicit = SWEmap + DHSVMdiffAccumSWE

# Sometimes modeled SWE can be higher than accumulation due to refreezing rain
# So we need to threshold at 0
AccumImplicit[,] = pmax(0, AccumImplicit[,])

TimeEffect = AccumImplicit - SWEmap

plot(AccumImplicit)
plot(TimeEffect)
plot(TimeEffect / SWEmap,range=c(0,2))

writeRaster(TimeEffect,"AccumMinusSWE.tif",overwrite=TRUE)
writeRaster(AccumImplicit,"AccumulationReference.tif",overwrite=TRUE)

AccumMapData = as.data.frame(as.matrix(AccumImplicit,wide=TRUE))

write.table(AccumMapData,"MassDistribution_Reference.csv",
            sep=",",row.names=FALSE,col.names=FALSE)

##########
# Use counterfactual accumulation map to match implied measurement-based accumulation map

AccumImplicit = rast("AccumulationReference.tif")

CounterfactualNames = c("Smoothed2x21","Smoothed2x41","GridmetLanczos","Uniform")

for (CounterfactualName in CounterfactualNames){
  
  AccumMapNoDrift = rast(paste0("E:/DHSVM_MountainHydrology/2_Calibration/Watersheds/TorreySnowRegion/ScenarioRuns/",
                              CounterfactualName,"/Map.Snow.SumAccum.asc"))
  crs(AccumMapNoDrift) = crs(AccumImplicit)
  
  ScaleFactor = mean(AccumImplicit[,]) / mean(AccumMapNoDrift[,])
  AccumMapNoDrift = AccumMapNoDrift * ScaleFactor
  
  plot(AccumMapNoDrift,main=paste0(CounterfactualName," Scale: ",signif(ScaleFactor,3)))
  
  dir.create(CounterfactualName)
  dir.create(paste0(CounterfactualName,"/EpochProgressMaps"))
  
  writeRaster(AccumMapNoDrift,paste0("AccumCounterfactual_",CounterfactualName,".tif"),overwrite=TRUE)
  writeRaster(AccumMapNoDrift,paste0(CounterfactualName,"/AccumCounterfactual.tif"),overwrite=TRUE)
  AccumMapNoDriftData = as.data.frame(as.matrix(AccumMapNoDrift,wide=TRUE))
  write.table(AccumMapNoDriftData,paste0(CounterfactualName,"/MassDistribution_Counterfactual.csv"),
              sep=",",row.names=FALSE,col.names=FALSE)
}












