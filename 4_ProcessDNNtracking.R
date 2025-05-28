# Copyright Eli Boardman 2024
# All Rights Reserved

library(sf)
library(terra)
library(smoothr)

setwd("C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/")

CounterfactualName = "Smoothed2x21_22.5deg_0.01loss"

SWEmap = rast("SWEreference.tif")
AccumMap = rast("AccumulationReference.tif")
DEM = rast("DEM.tif")
AccumMapNoDrift = rast(paste0("AccumCounterfactual_",sub("_.*","",CounterfactualName),".tif"))

setwd(CounterfactualName)

EfficiencyVal = as.numeric(readLines("Efficiency.txt"))
ScaleVal = as.numeric(readLines("InputScale.txt"))

SWEmapMod = rast("DNN_SWEfinal.tif")
AccumMapMod = rast("DNN_AccumFinal.tif")
ModelPred = as.matrix(read.table("FinalStorage.csv",sep=","))
AccumMapNoDriftData = as.matrix(read.table("MassDistribution_Counterfactual.csv",sep=","))

AccumMapNoDrift = AccumMapNoDrift * ScaleVal
AccumMapNoDriftData = AccumMapNoDriftData * ScaleVal

################################################################################
# Transport distance results

TransportDist = as.matrix(read.table("AvgTransportDist_ForNetImport.csv",sep=","))

image(t(apply(TransportDist, 2, rev)))

TransportDistMapImport = 0 * SWEmap
TransportDistMapImport[,] = as.numeric(t(TransportDist))
TransportDistMapImport = TransportDistMapImport * res(TransportDistMapImport)[1]

TransportDistMapImport[which(TransportDistMapImport[,] < 0)] = 0
TransportDistMapImport[which(!is.finite(TransportDistMapImport[,]))] = 0

plot(TransportDistMapImport)

writeRaster(TransportDistMapImport,"AvgTransportDistance_ForNetImport.tif",overwrite=TRUE)

##########
# Weighted average transport distance of final SWE

TransportDistMap = 0 * SWEmap

EffectiveImportFrac = (AccumMapMod - AccumMapNoDrift) / AccumMapMod

max(EffectiveImportFrac[,])
plot(EffectiveImportFrac,range=c(0,1))

ImportPtrs = which(EffectiveImportFrac[,] > 0)
TransportDistMap[ImportPtrs] = TransportDistMapImport[ImportPtrs] * EffectiveImportFrac[ImportPtrs]

TransportDistMap[which(!is.finite(TransportDistMap[,]))] = 0

plot(TransportDistMap)

writeRaster(EffectiveImportFrac,"EffectiveImportFrac.tif",overwrite=TRUE)
writeRaster(TransportDistMap,"AvgTransportDistance.tif",overwrite=TRUE)

################################################################################
# Transport distance histograms

DistHistData = read.csv("TransportDistanceHistogram.csv")
head(DistHistData)

DistHistData$AvgDist[!is.finite(DistHistData$AvgDist)] = 0

# Add columns with final accumulation and final SWE
FinalAccumTotalData = as.data.frame(AccumMapMod,xy=TRUE)
FinalSWEdata = as.data.frame(SWEmapMod,xy=TRUE)
RasterRes = res(SWEmapMod)[1]
FinalAccumTotalData$xy = paste0((FinalAccumTotalData$x - ext(SWEmapMod)[1] - RasterRes/2) / RasterRes,
                                (ext(SWEmapMod)[4] - FinalAccumTotalData$y - RasterRes/2) / RasterRes)
FinalSWEdata$xy = paste0((FinalSWEdata$x - ext(SWEmapMod)[1] - RasterRes/2) / RasterRes,
                         (ext(SWEmapMod)[4] - FinalSWEdata$y - RasterRes/2) / RasterRes)

DistHistData$FinalAccumTotal = FinalAccumTotalData[match(paste0(DistHistData$x,DistHistData$y),FinalAccumTotalData$xy),3]
DistHistData$FinalSWE = FinalSWEdata[match(paste0(DistHistData$x,DistHistData$y),FinalSWEdata$xy),3]

DistHistData$RetainedInitial = DistHistData$FinalAccumTotal - DistHistData$ImportTotal

write.csv(DistHistData,"TransportDistanceHistogram.csv")

plot(unlist(DistHistData[which.max(DistHistData$FinalSWE),grepl("ImportMinDist",names(DistHistData))]),type="l")

ImportDistCols = which(grepl("ImportMinDist",names(DistHistData)))
DistMins = as.numeric(sub("ImportMinDist","",names(DistHistData)[ImportDistCols]))
DistMins = DistMins * res(SWEmap)[1]
DistMaxs = c(DistMins[-1],Inf)

DistHistDataByDepth = data.frame()
DepthBinLimitMax = 6
DepthBinSpacing = 1
DepthBinLimits = seq(DepthBinSpacing,DepthBinLimitMax,DepthBinSpacing)
for (DepthBinLimit in DepthBinLimits){
  
  Ptrs = which(DistHistData$FinalSWE >= (DepthBinLimit - DepthBinSpacing) &
                 DistHistData$FinalSWE < DepthBinLimit)
  nCells = length(Ptrs)
  
  DataChunk = data.frame(MinDepth=(DepthBinLimit - DepthBinSpacing),
                         MaxDepth=DepthBinLimit,
                         nCells=nCells,
                         AvgSWE=sum(DistHistData$FinalAccumTotal[Ptrs]) / nCells,
                         AvgAccum=sum(DistHistData$FinalAccumTotal[Ptrs]) / nCells,
                         AvgImport=sum(DistHistData$ImportTotal[Ptrs]) / nCells,
                         AvgRetained=sum(DistHistData$RetainedInitial[Ptrs]) / nCells,
                         DistMin=DistMins,
                         DistMax=DistMaxs,
                         DistContrib=unname(unlist(colSums(DistHistData[Ptrs,ImportDistCols]))) / nCells)
  
  DataChunk$CumulativeAccum = DataChunk$AvgRetained[1] + cumsum(DataChunk$DistContrib)
  DataChunk$CumulativeAccumPct = DataChunk$CumulativeAccum / DataChunk$AvgAccum
  
  DataChunk$DistContribPct = DataChunk$DistContrib / DataChunk$AvgAccum
  DataChunk$DistContribPctCumulative = cumsum(DataChunk$DistContribPct)
  
  DataChunk$DCPctImport = DataChunk$DistContrib / DataChunk$AvgImport
  DataChunk$DCPctImportCumulative = cumsum(DataChunk$DCPctImport)
  
  DistHistDataByDepth = rbind(DistHistDataByDepth, DataChunk)
}

DistHistDataPlot = DistHistDataByDepth
DistHistDataPlot$TransportDist = DistHistDataPlot$DistMin + res(SWEmap)[1]/2
DistHistDataPlot$TransportDist[DistHistDataPlot$DistMin==0] = 0
DistHistDataPlot$TransportDist[DistHistDataPlot$DistMin==max(DistHistDataPlot$DistMin)] = max(DistHistDataPlot$DistMin) + res(SWEmap)[1]
DistHistDataPlot$DepthRange = paste0(DistHistDataPlot$MinDepth," - ",DistHistDataPlot$MaxDepth," m")

write.csv(DistHistDataPlot,"DistanceHistogramPlotData.csv")

################################################################################
# Accumulation area for each pixel

ImportAreaData = read.csv("ImportAreaMaps.csv")
ImportAreaData$x = ImportAreaData$x + 1 # Python to R
ImportAreaData$y = ImportAreaData$y + 1 # Python to R

hist(ImportAreaData$ImportAreaMin1mm)
max(ImportAreaData$ImportAreaMin1mm) * 100^2 / 1000^2
max(ImportAreaData$ImportAreaMin10mm) * 100^2 / 1000^2
max(ImportAreaData$ImportAreaMin100mm) * 100^2 / 1000^2

ImportAreaData1mm = 0 * ModelPred
ImportAreaData10mm = 0 * ModelPred
ImportAreaData100mm = 0 * ModelPred

for (i in 1:nrow(ImportAreaData)){
  x = ImportAreaData$x[i]
  y = ImportAreaData$y[i]
  ImportAreaData1mm[y,x] = ImportAreaData$ImportAreaMin1mm[i]
  ImportAreaData10mm[y,x] = ImportAreaData$ImportAreaMin10mm[i]
  ImportAreaData100mm[y,x] = ImportAreaData$ImportAreaMin100mm[i]
}

ConversionFact = res(SWEmap)[1]^2 / 1000^2 # Grid cells to km^2

ImportAreaMap1mm = 0 * SWEmap
ImportAreaMap10mm = 0 * SWEmap
ImportAreaMap100mm = 0 * SWEmap
ImportAreaMap1mm[,] = as.numeric(t(ImportAreaData1mm)) * ConversionFact
ImportAreaMap10mm[,] = as.numeric(t(ImportAreaData10mm)) * ConversionFact
ImportAreaMap100mm[,] = as.numeric(t(ImportAreaData100mm)) * ConversionFact

plot(ImportAreaMap1mm)

writeRaster(ImportAreaMap1mm,"ImportAreaMin1mm.tif",overwrite=TRUE)
writeRaster(ImportAreaMap10mm,"ImportAreaMin10mm.tif",overwrite=TRUE)
writeRaster(ImportAreaMap100mm,"ImportAreaMin100mm.tif",overwrite=TRUE)

plot(TransportDistMap[,], ImportAreaMap1mm[,])
cor(as.numeric(TransportDistMap[,]), as.numeric(ImportAreaMap1mm[,])) # 0.899

plot(AccumMapMod[,], ImportAreaMap1mm[,])
cor(as.numeric(AccumMapMod[,]), as.numeric(ImportAreaMap1mm[,])) # 0.657

PlotData = data.frame(SWEdepth=unname(unlist(SWEmapMod[,])),
                      TotalAccum=unname(unlist(AccumMapMod[,])),
                      LocalSnowfall=unname(unlist(AccumMapNoDrift[,])),
                      NetImport=unlist(pmax(0,(SWEmapMod - AccumMapNoDrift)[,])),
                      TransportDist=pmax(0,unlist(TransportDistMap[,])),
                      TransportDistImport=pmax(0,unlist(TransportDistMapImport[,])),
                      ImportAreaMin1mm=unname(unlist(ImportAreaMap1mm[,])),
                      ImportAreaMin10mm=unname(unlist(ImportAreaMap10mm[,])),
                      ImportAreaMin100mm=unname(unlist(ImportAreaMap100mm[,])),
                      Elevation=unname(unlist(DEM[,])))

write.csv(PlotData,"AvgDistancePlotData.csv")

################################################################################
# Watershed contribution results

MaskNames = c("TorreyWatershed","ContinentalGlacier")

for (MaskName in MaskNames){
  
  Mask = rast(paste0("../Mask_",MaskName,".tif"))
  
  MaskContrib = as.matrix(read.table(paste0("ExportContribToMask_",MaskName,".csv"),sep=","))
  MaskExit = as.matrix(read.table(paste0("ExportExitFromMask_",MaskName,".csv"),sep=","))
  
  image(t(apply(MaskContrib, 2, rev)))
  image(t(apply(MaskExit, 2, rev)))
  
  MaskContribFrac = 0 * SWEmap + 1
  MaskContribMapRaw = 0 * SWEmap
  MaskExitMapRaw = 0 * SWEmap
  MaskContribMapRaw[,] = as.numeric(t(MaskContrib))
  MaskExitMapRaw[,] = as.numeric(t(MaskExit))
  
  plot(MaskContribMapRaw)
  plot(MaskExitMapRaw)
  
  ##########
  # Import and export depth
  
  # First, only consider positive inputs from outside the mask
  MaskContribMap = MaskContribMapRaw * (1 - Mask)
  
  # Now compute negative contributions (export) from inside the mask
  SWEexport = rast("DNN_SWEexport.tif")
  NetExportInsidePtrs = which(MaskExitMapRaw[,] > 0 & Mask[,] == 1)
  MaskContribMap[NetExportInsidePtrs] = -1 * MaskExitMapRaw[NetExportInsidePtrs]
  
  ##########
  # Import and export fraction
  
  # Snow originating in cells with a net gain just stays put
  NetGainPtrs = which(AccumMapMod[,] >= AccumMapNoDrift[,])
  MaskContribFrac[NetGainPtrs] = Mask[NetGainPtrs]
  
  # Snow exported from outside the mask is divided by the total accumulation at the origin cell
  NetExportOutsidePtrs = which(AccumMapMod[,] < AccumMapNoDrift[,] & Mask[,] == 0)
  MaskContribFrac[NetExportOutsidePtrs] = MaskContribMapRaw[NetExportOutsidePtrs] / AccumMapNoDrift[NetExportOutsidePtrs]
  
  # Remaining snow must first be added to the export from cells inside the mask
  NetExportInsidePtrs = which(MaskExitMapRaw[,] > 0 & Mask[,] == 1)
  MaskContribFrac[NetExportInsidePtrs] = 1 - (MaskExitMapRaw[NetExportInsidePtrs] / AccumMapNoDrift[NetExportInsidePtrs])
  
  plot(MaskContribMap,range=c(-1.5,1.5))
  plot(MaskContribFrac)
  
  writeRaster(MaskContribMap,paste0("MaskContribTotalDepth_",MaskName,".tif"),overwrite=TRUE)
  writeRaster(MaskContribFrac,paste0("MaskContribFrac_",MaskName,".tif"),overwrite=TRUE)
  
  ##########
  # Polygons at various thresholds
  
  dir.create(paste0("MaskPolygons_",MaskName))
  
  RawPolygons = as.polygons(Mask,values=TRUE)
  RawPolygons = RawPolygons[RawPolygons[[1]]==1]
  RawPolygons = st_geometry(st_as_sf(RawPolygons))
  SmoothPolygonsBase = smooth(RawPolygons,method="ksmooth",smoothness=2)
  
  plot(SmoothPolygonsBase)
  
  st_write(SmoothPolygonsBase,paste0("MaskPolygons_",MaskName,"/",MaskName,".shp"),
           append=FALSE)
  
  for (Pcontrib in c(0.01,0.1,0.5,0.9,0.99)){
    
    MaskThresh = 0 * Mask
    MaskThresh[which(MaskContribFrac[,] >= Pcontrib)] = 1
    
    RawPolygons = as.polygons(MaskThresh,values=TRUE)
    RawPolygons = RawPolygons[RawPolygons[[1]]==1]
    RawPolygons = st_geometry(st_as_sf(RawPolygons))
    SmoothPolygons = smooth(RawPolygons,method="ksmooth",smoothness=2)
    
    plot(SmoothPolygons,main=Pcontrib)
    
    st_write(SmoothPolygons,paste0("MaskPolygons_",MaskName,"/",MaskName,"_Min",Pcontrib*100,"pct.shp"),
             append=FALSE)
    
    MaskThreshExport = MaskThresh * 0
    MaskThreshExport[which(MaskThresh[,]==0 & Mask[,]==1)] = 1
    MaskThreshImport = MaskThresh * 0
    MaskThreshImport[which(MaskThresh[,]==1 & Mask[,]==0)] = 1
    
    RawPolygons = as.polygons(MaskThreshExport,values=TRUE)
    SmoothPolygons = smooth(st_geometry(st_as_sf(RawPolygons[RawPolygons[[1]]==1])),
                            method="ksmooth",smoothness=2)
    st_write(SmoothPolygons,paste0("MaskPolygons_",MaskName,"/",MaskName,"_Min",Pcontrib*100,"pct_Export.shp"),
             append=FALSE)
    
    RawPolygons = as.polygons(MaskThreshImport,values=TRUE)
    SmoothPolygons = smooth(st_geometry(st_as_sf(RawPolygons[RawPolygons[[1]]==1])),
                            method="ksmooth",smoothness=2)
    st_write(SmoothPolygons,paste0("MaskPolygons_",MaskName,"/",MaskName,"_Min",Pcontrib*100,"pct_Import.shp"),
             append=FALSE)
  }
  
  if (FALSE){
    # Create smoothed stream network
    
    StreamsRaw = st_read("E:/DHSVM_MountainHydrology/1_Setup/ChannelData/ChannelNetwork_NHD_HU8-10080001/ChannelNetwork_NHD_HU8-10080001.shp")
    StreamsRaw = st_transform(StreamsRaw,"EPSG:32612")
    StreamsRaw = st_intersection(st_zm(StreamsRaw),st_read("MaskPolygons_TorreyWatershed/TorreyWatershed.shp"))
    StreamsRaw = st_geometry(StreamsRaw)
    
    st_write(StreamsRaw,"../Extra/TorreyCreekStreamsNHD.shp",append=FALSE)
    
  }
  
  ##########
  # Area within a given threshold
  
  MaskContribThreshData = data.frame(Pcontrib=c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99),
                                     Area_km2=NA)
  RastRes = res(Mask)[1]
  TotalArea_km2 = length(which(Mask[,]==1)) * RastRes^2 / (1000^2)
  for (i in 1:nrow(MaskContribThreshData)){
    MaskContribThreshData$Area_km2[i] = length(which(MaskContribFrac[,] > MaskContribThreshData$Pcontrib[i])) * RastRes^2 / (1000^2)
  }
  plot(MaskContribThreshData$Pcontrib,MaskContribThreshData$Area_km2 / TotalArea_km2,type="o")
  lines(c(-1,2),c(1,1))
  lines(c(0.5,0.5),c(0,2))
  
  write.csv(MaskContribThreshData,paste0("MaskContribThreshData_",MaskName,".csv"))
  
  ImportTot = sum(pmax(0,MaskContribMap[,]))
  ExportTot = sum(pmin(0,MaskContribMap[,]))
  FinalTot = sum((SWEmapMod * Mask)[,])
  ImportTot / FinalTot
  ExportTot / FinalTot
  
}

################################################################################
# Visualizations

if (FALSE){
  ValidMask = 0 * SWEmap + 1
  
  CanopyMap = rast("E:/DHSVM_MountainHydrology/1_Setup/Watersheds/TorreySnowRegion/CanopyCover.tif")
  
  ValidMask[which(DEM[,] < 3000)] = 0
  ValidMask[which(CanopyMap[,] > 0)] = 0
  
  writeRaster(ValidMask,"../ValidMask.tif",overwrite=TRUE)
}


if (FALSE){
  
  RawPolygons = st_read("Extra/PhotoLocationsRaw.shp")
  RawPolygons = st_geometry(RawPolygons)
  SmoothPolygons = smooth(RawPolygons,method="ksmooth",smoothness=1)
  
  plot(RawPolygons)
  plot(SmoothPolygons,add=TRUE)
  
  st_write(SmoothPolygons,"Extra/PhotoLocationsSmooth.shp",
           append=FALSE)
  
}


