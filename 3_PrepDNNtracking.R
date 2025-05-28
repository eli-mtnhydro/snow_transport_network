# Copyright Eli Boardman 2024
# All Rights Reserved

library(terra)

setwd("C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/")

CounterfactualName = "Smoothed2x21_22.5deg_0.1loss"

SWEmap = rast("SWEreference.tif")
AccumMap = rast("AccumulationReference.tif")
DEM = rast("DEM.tif")
AccumMapNoDrift = rast(paste0("AccumCounterfactual_",sub("_.*","",CounterfactualName),".tif"))

setwd(CounterfactualName)

AccumMapMod = rast("DNN_AccumFinal.tif")
ModelPred = as.matrix(read.table("FinalStorage.csv",sep=","))
AccumMapNoDriftData = as.matrix(read.table("MassDistribution_Counterfactual.csv",sep=","))
NetLoss = as.matrix(read.table("NetLoss.csv",sep=","))
NetLossMap = rast("DNN_NetLoss.tif")

EfficiencyVal = as.numeric(readLines("Efficiency.txt"))
ScaleVal = as.numeric(readLines("InputScale.txt"))

AccumMapNoDrift = AccumMapNoDrift * ScaleVal
AccumMapNoDriftData = AccumMapNoDriftData * ScaleVal

################################################################################
# Create necessary inputs

# Add back loss so that particle-tracking simulation is mass-conserving
# NetLossShift = NetLoss * 0
# NetLossShift[,2:ncol(NetLoss)] = NetLoss[,1:(ncol(NetLoss)-1)]
# ModelPredAdj = ModelPred + NetLossShift

NetExport = AccumMapNoDriftData - ModelPred
NetExport[,] = pmax(0, NetExport[,])

NetGain = ModelPred - AccumMapNoDriftData
NetGain[,] = pmax(0, NetGain[,])

image(t(apply(NetGain, 2, rev)))
image(t(apply(NetExport, 2, rev)))

write.table(NetExport,"NetExport.csv",
            sep=",",row.names=FALSE,col.names=FALSE)

##########
# Same thing, but with rasters for easier visualization

# Add back loss so that particle-tracking simulation is mass-conserving
# NetLossMapShift = NetLossMap * 0
# NetLossMapShift[,] = as.numeric(t(NetLossShift))
# AccumMapModAdj = AccumMapMod + NetLossMapShift

SWEexport = AccumMapNoDrift - AccumMapMod
SWEexport[,] = pmax(0, SWEexport[,])

SWEimport = AccumMapMod - AccumMapNoDrift
SWEimport[,] = pmax(0, SWEimport[,])

SWEkeep = AccumMapNoDrift
SWEkeep[,] = pmin(AccumMapNoDrift[,], AccumMapMod[,])

plot(SWEexport)
plot(SWEimport)
plot(SWEkeep)

writeRaster(SWEexport,"DNN_SWEexport.tif",overwrite=TRUE)
writeRaster(SWEimport,"DNN_SWEimport.tif",overwrite=TRUE)
writeRaster(SWEkeep,"DNN_SWEkeep.tif",overwrite=TRUE)

# Two other equivalent ways of defining SWEkeep:
plot((AccumMapNoDrift - SWEexport) - SWEkeep)
plot((AccumMapMod - SWEimport) - SWEkeep)

##########
# Re-scale transfer fractions so that they represent the flux of net export only

OutgoingFracUp = as.matrix(read.table("TransferFracs_Up.csv",sep=","))
OutgoingFracAcross = as.matrix(read.table("TransferFracs_Across.csv",sep=","))
OutgoingFracDown = as.matrix(read.table("TransferFracs_Down.csv",sep=","))

OutgoingFracUpExport = OutgoingFracUp
OutgoingFracAcrossExport = OutgoingFracAcross
OutgoingFracDownExport = OutgoingFracDown

OutgoingFracTotal = OutgoingFracUp + OutgoingFracAcross + OutgoingFracDown

NetExportPtrs = (NetExport[,] > 0 & OutgoingFracTotal > 0)
NetGainPtrs = (NetGain[,] > 0 & OutgoingFracTotal > 0)

# Cells with net export from initial storage must first export ANY imported snow
OutgoingFracUpExport[NetExportPtrs] = OutgoingFracUp[NetExportPtrs] / OutgoingFracTotal[NetExportPtrs]
OutgoingFracAcrossExport[NetExportPtrs] = OutgoingFracAcross[NetExportPtrs] / OutgoingFracTotal[NetExportPtrs]
OutgoingFracDownExport[NetExportPtrs] = OutgoingFracDown[NetExportPtrs] / OutgoingFracTotal[NetExportPtrs]

# Cells with net gain ONLY export a fraction of the imported flux necessary to match learned flux
ExportFlux = pmin(1,OutgoingFracTotal[,]) * ModelPred / (1 - pmin(1,OutgoingFracTotal[,]))
OutgoingFracTotalExport = ExportFlux / (ExportFlux + NetGain)
OutgoingFracRatio = OutgoingFracTotalExport / OutgoingFracTotal
OutgoingFracUpExport[NetGainPtrs] = OutgoingFracUp[NetGainPtrs] * OutgoingFracRatio[NetGainPtrs]
OutgoingFracAcrossExport[NetGainPtrs] = OutgoingFracAcross[NetGainPtrs] * OutgoingFracRatio[NetGainPtrs]
OutgoingFracDownExport[NetGainPtrs] = OutgoingFracDown[NetGainPtrs] * OutgoingFracRatio[NetGainPtrs]

OutgoingFracTotalExport = OutgoingFracUpExport + OutgoingFracAcrossExport + OutgoingFracDownExport

which(!is.finite(OutgoingFracTotalExport))
max(OutgoingFracTotalExport) - 1
min(OutgoingFracTotalExport)

# OversizePtrs = (OutgoingFracTotal > 1)
# OutgoingFracUpExport[OversizePtrs] = OutgoingFracUpExport[OversizePtrs] / OutgoingFracTotalExport[OversizePtrs]
# OutgoingFracAcrossExport[OversizePtrs] = OutgoingFracAcrossExport[OversizePtrs] / OutgoingFracTotalExport[OversizePtrs]
# OutgoingFracDownExport[OversizePtrs] = OutgoingFracDownExport[OversizePtrs] / OutgoingFracTotalExport[OversizePtrs]
# 
# OutgoingFracTotalExport = OutgoingFracUpExport + OutgoingFracAcrossExport + OutgoingFracDownExport

image(t(apply(OutgoingFracTotal, 2, rev)))
image(t(apply(OutgoingFracTotalExport, 2, rev)))

write.table(OutgoingFracUpExport,"TransferFracs_ExportOnly_Up.csv",
            sep=",",row.names=FALSE,col.names=FALSE)
write.table(OutgoingFracAcrossExport,"TransferFracs_ExportOnly_Across.csv",
            sep=",",row.names=FALSE,col.names=FALSE)
write.table(OutgoingFracDownExport,"TransferFracs_ExportOnly_Down.csv",
            sep=",",row.names=FALSE,col.names=FALSE)

OutgoingExportFracMap = 0 * SWEmap
OutgoingExportFracMap_Up = 0 * SWEmap
OutgoingExportFracMap_Across = 0 * SWEmap
OutgoingExportFracMap_Down = 0 * SWEmap
OutgoingExportFracMap[,] = as.numeric(t(OutgoingFracTotalExport))
OutgoingExportFracMap_Up[,] = as.numeric(t(OutgoingFracUpExport))
OutgoingExportFracMap_Across[,] = as.numeric(t(OutgoingFracAcrossExport))
OutgoingExportFracMap_Down[,] = as.numeric(t(OutgoingFracDownExport))

writeRaster(OutgoingExportFracMap,"DNN_TransferFrac_ExportOnly.tif",overwrite=TRUE)
writeRaster(OutgoingExportFracMap_Across,"DNN_TransferFrac_ExportOnly_X.tif",overwrite=TRUE)
writeRaster(OutgoingExportFracMap_Up - OutgoingExportFracMap_Down,"DNN_TransferFrac_ExportOnly_Y.tif",overwrite=TRUE)

################################################################################
# Save masks for fate tracking

if (FALSE){
  
  ##########
  # Torrey Creek stream gauge
  
  Mask = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreyCreek_30m/Watershed_TorreyCreekGage.tif")
  plot(Mask)
  
  # Bilinear interpolation, then clamp to 0-1
  Mask = project(Mask,SWEmap)
  Mask[,] = round(Mask[,])
  Mask[which(!is.finite(Mask[,]))] = 0
  plot(Mask)
  unique(Mask[,])
  
  MaskName = "TorreyWatershed"
  
  writeRaster(Mask,paste0("../Mask_",MaskName,".tif"),overwrite=TRUE)
  MaskData = as.data.frame(as.matrix(Mask,wide=TRUE))
  write.table(MaskData,paste0("../Mask_",MaskName,".csv"),
              sep=",",row.names=FALSE,col.names=FALSE)
  
  ##########
  # Continental Glacier (CG1 stream gauge of V&V)
  
  Mask = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/ContinentalGlacierStream/Watershed.tif")
  plot(Mask)
  
  # Bilinear interpolation, then clamp to 0-1
  Mask = project(Mask,SWEmap)
  Mask[,] = round(Mask[,])
  Mask[which(!is.finite(Mask[,]))] = 0
  plot(Mask)
  unique(Mask[,])
  
  MaskName = "ContinentalGlacier"
  
  writeRaster(Mask,paste0("../Mask_",MaskName,".tif"),overwrite=TRUE)
  MaskData = as.data.frame(as.matrix(Mask,wide=TRUE))
  write.table(MaskData,paste0("../Mask_",MaskName,".csv"),
              sep=",",row.names=FALSE,col.names=FALSE)
  
}

################################################################################
# Test one-hot precipitation pattern and distance

if (FALSE){
  # ModelPredOneHot = as.matrix(read.table("FinalStorage_OneHotExample.csv",sep=","))
  # image(t(apply(log10(ModelPredOneHot), 2, rev)))
  # OneHotMap = 0 * SWEmap
  # OneHotMap[,] = as.numeric(t(ModelPredOneHot))
  # writeRaster(OneHotMap,"DNN_SWEfinal_OneHotExample.tif",overwrite=TRUE)
  # 
  # DistData = as.matrix(read.table("DistX145Y80.csv",sep=","))
  # DistMap = 0 * SWEmap
  # DistMap[,] = as.numeric(t(DistData)) * res(DistMap)[1]
  # writeRaster(DistMap,"DistanceExampleX145Y80.tif",overwrite=TRUE)
  
  ModelPredHot = as.matrix(read.table("TotalImport_OneHotSum.csv",sep=","))
  image(t(apply(ModelPredHot, 2, rev)))
  OneHotMap = 0 * SWEmap
  OneHotMap[,] = as.numeric(t(ModelPredHot))
  OneHotMap = OneHotMap + SWEkeep
  OneHotErr = OneHotMap - AccumMapMod
  
  plot(OneHotErr)
  mean(abs(OneHotErr[,])) / mean(AccumMapMod[,]) # 10^-8
  
  writeRaster(OneHotErr,"DNN_OneHotErr.tif",overwrite=TRUE)
}



