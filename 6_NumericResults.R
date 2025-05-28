# Copyright Eli Boardman 2024
# All Rights Reserved

library(terra)
library(sf)

setwd("C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/")

CounterfactualName = "Smoothed2x21_22.5deg_0.01loss"

setwd(CounterfactualName)

TrainingProgressData = read.csv("TrainingProgress.csv")
TransportData = read.csv("AvgDistancePlotData.csv")[,-1]
DistHistData = read.csv("DistanceHistogramPlotData.csv")[,-1]
DistHistDataFull = read.csv("TransportDistanceHistogram.csv")

AccumMap = rast("../AccumulationReference.tif")
AccumMapNoDrift = rast(paste0("../AccumCounterfactual_",sub("_.*","",CounterfactualName),".tif"))
SWEmap = rast("../SWEreference.tif")
AccumMapMod = rast("DNN_AccumFinal.tif")
SWEmapMod = rast("DNN_SWEfinal.tif")
SWEerr = rast("DNN_SWEerror.tif")
OutgoingFlux = rast("DNN_TransferFlux.tif")
OutgoingFrac = rast("DNN_TransferFracTotal.tif")
NetExport = rast("DNN_SWEexport.tif")
NetImport = rast("DNN_SWEimport.tif")
TransportDist = rast("AvgTransportDistance.tif")
# TransportDistImport = rast("AvgTransportDistance_ForNetImport.tif")
ImportFrac = rast("EffectiveImportFrac.tif")
AccumArea1mm = rast("ImportAreaMin1mm.tif")
AccumArea10mm = rast("ImportAreaMin10mm.tif")

MaskTorrey = rast("../Mask_TorreyWatershed.tif")
MaskFracTorrey = rast("MaskContribFrac_TorreyWatershed.tif")
MaskAbsTorrey = rast("MaskContribTotalDepth_TorreyWatershed.tif")
AreaDataTorrey = read.csv("MaskContribThreshData_TorreyWatershed.csv")[,-1]

MaskContinental = rast("../Mask_ContinentalGlacier.tif")
MaskFracContinental = rast("MaskContribFrac_ContinentalGlacier.tif")
MaskAbsContinental = rast("MaskContribTotalDepth_ContinentalGlacier.tif")
AreaDataContinental = read.csv("MaskContribThreshData_ContinentalGlacier.csv")[,-1]

DEM = rast("../DEM.tif")
WindXY = rast("../WindNinja/AbsWindXY.tif")
WindShelter = rast("../DEM_SX_270deg.tif")

MassVals = read.table("../SnowPattern_Mass.Balance",header=TRUE)

CanopyMap = rast("../CanopyCover.tif")

##########
# Make Winstral Sx map

# library(whitebox)
# wbt_horizon_angle(dem="../DEM.tif",
#                   output="../DEM_SX_270deg.tif",
#                   azimuth=270,
#                   max_dist=NULL,
#                   wd=getwd())

################################################################################
# Sublimation

plot(MassVals$Swq,type="l")
plot(cumsum(MassVals$SnowVaporFlux),type="l")

sum(MassVals$SnowVaporFlux) / 0.4197
max(MassVals$Swq)

sum(MassVals$Precip.m.)
sum(MassVals$Snow.m.) # Note that total accumulation is higher because of re-freezing rain
sum(MassVals$SnowVaporFlux)

# Amount above treeline

plot(SWEmap)
plot(CanopyMap)

sum(unlist(SWEmap[which(CanopyMap[,] == 0)]))
sum(unlist(SWEmap[which(CanopyMap[,] > 0)]))
12035 / (12035 + 2189)

sum(unlist(SWEmap[which(CanopyMap[,] == 0 & SWEmap[,] > 1)]))
sum(unlist(SWEmap[which(CanopyMap[,] > 0 & SWEmap[,] > 1)]))
3905.7 / (3905.7 + 3.3)

sum(unlist(SWEmap[which(CanopyMap[,] > 0 & SWEmap[,] > 2)]))

################################################################################
# Overall redistribution flux

mean(OutgoingFrac[,]) # 0.54

TrainingProgressData$FluxCost[nrow(TrainingProgressData)] # 0.92
mean(OutgoingFlux[,]) # 0.80
quantile(OutgoingFlux[,],c(0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999,0.9999))

mean(AccumMapMod[,]) # 0.42
plot(AccumMap / AccumMapNoDrift)
quantile(unlist((AccumMap / AccumMapNoDrift)[which(SWEmap[,] > 3)]))

# Relationship between topographic shelter and import/export
median(unlist(WindShelter[which(NetExport[,] > 0)]))
median(unlist(WindShelter[which(NetImport[,] > 0)]))

median(unlist(WindShelter[which(NetExport[,] > 0.5)]))
median(unlist(WindShelter[which(NetImport[,] > 0.5)]))
median(unlist(WindShelter[which(NetImport[,] > 1)]))
median(unlist(WindShelter[which(NetImport[,] > 2)]))

# Relationship between wind speed and redistribution flux
plot(WindXY[,], OutgoingFlux[,])
cor.test(WindXY[,], OutgoingFlux[,])

quantile(WindXY[,],c(0.1,0.25,0.5,0.75,0.9))
median(unlist(OutgoingFlux[which(WindXY[,] < quantile(WindXY[,],c(0.25)))]))
median(unlist(OutgoingFlux[which(WindXY[,] > quantile(WindXY[,],c(0.75)))]))

# mean(unlist(NetExport[which(WindXY[,] < quantile(WindXY[,],c(0.25)))]))
# mean(unlist(NetExport[which(WindXY[,] > quantile(WindXY[,],c(0.75)))]))

# Snowfall vs. net accumulation patterns
sd(unlist(AccumMapNoDrift[,]))
median(unlist(AccumMapNoDrift[,]))
max(unlist(AccumMapNoDrift[,]))

sd(unlist(AccumMapMod[,]))
median(unlist(AccumMapMod[,]))
max(unlist(AccumMapMod[,]))

# SWE error
sqrt(mean(unlist(SWEerr[,])^2))
1 - sum((unlist(SWEerr[,]))^2) / sum((SWEmap[,] - mean(SWEmap[,]))^2)

################################################################################
# Transport distance and source areas

# Each cell is already mass-weighted, now compute mass-weighted average across cells
mean(unlist(TransportDist[,] * AccumMapMod[,])) / mean(unlist(AccumMapMod[,]))
Ptrs = which(SWEmapMod[,] > 1)
mean(unlist(TransportDist[Ptrs] * AccumMapMod[Ptrs])) / mean(unlist(AccumMapMod[Ptrs]))
Ptrs = which(SWEmapMod[,] > 3)
mean(unlist(TransportDist[Ptrs] * AccumMapMod[Ptrs])) / mean(unlist(AccumMapMod[Ptrs]))

DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==1 & DistHistData$TransportDist==0]
1 - DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==2 & DistHistData$TransportDist==0]

1 - DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==2 & DistHistData$DistMax==1000]
1 - DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==3 & DistHistData$DistMax==1000]
1 - DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==4 & DistHistData$DistMax==1000]
1 - DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==5 & DistHistData$DistMax==1000]
1 - DistHistData$CumulativeAccumPct[DistHistData$MaxDepth==6 & DistHistData$DistMax==1000]

# Fraction of all snow (by mass) originating > 1 km or > 2 km
DistCols = which(as.numeric(sub("ImportMinDist","",names(DistHistDataFull))) >= 10) # or >=20
DepthRows = which(DistHistDataFull$FinalSWE > 4)
sum(DistHistDataFull[DepthRows,DistCols]) / sum(DistHistDataFull$FinalAccumTotal[DepthRows])
length(DepthRows)

cor.test(TransportData$SWEdepth,TransportData$TransportDist)

sort(TransportData$TransportDist[TransportData$SWEdepth > 3])
range(TransportData$TransportDist[TransportData$SWEdepth > 1])

Ptrs = which(TransportData$SWEdepth > 1.9 & TransportData$SWEdepth < 2.1)
cor.test(TransportData$Elevation[Ptrs],TransportData$TransportDist[Ptrs])

TransportData[TransportData$TransportDist > 1800,]

cor(TransportData$TransportDist, TransportData$ImportAreaMin1mm)

median(TransportData$ImportAreaMin1mm[TransportData$SWEdepth > 2])
median(TransportData$ImportAreaMin10mm[TransportData$SWEdepth > 2])
median(TransportData$ImportAreaMin100mm[TransportData$SWEdepth > 2])

median(TransportData$ImportAreaMin1mm[TransportData$SWEdepth > 3])
median(TransportData$ImportAreaMin10mm[TransportData$SWEdepth > 3])
median(TransportData$ImportAreaMin100mm[TransportData$SWEdepth > 3])

median(TransportData$ImportAreaMin1mm[TransportData$SWEdepth > 5])
median(TransportData$ImportAreaMin10mm[TransportData$SWEdepth > 5])
median(TransportData$ImportAreaMin100mm[TransportData$SWEdepth > 5])

TransportData[TransportData$ImportAreaMin1mm > 2.5,]

sort(TransportData$ImportAreaMin1mm[TransportData$SWEdepth > 3])

cor.test(TransportData$ImportAreaMin1mm[TransportData$SWEdepth > 3], TransportData$ImportAreaMin10mm[TransportData$SWEdepth > 3])
cor.test(TransportData$ImportAreaMin1mm[TransportData$SWEdepth > 3], TransportData$ImportAreaMin100mm[TransportData$SWEdepth > 3])
length(which(TransportData$SWEdepth > 3))

TransportData[TransportData$ImportAreaMin100mm > 0.1,]
mean(TransportData$TransportDist[TransportData$ImportAreaMin100mm > 0.1])
range(TransportData$TotalAccum[TransportData$ImportAreaMin100mm > 0.1])
mean(TransportData$TransportDist[TransportData$SWEdepth > 2.1 & TransportData$SWEdepth < 4.5])

################################################################################
# Snowshed boundaries

plot(MaskFracTorrey)
plot(MaskFracContinental)

print(AreaDataTorrey)
print(AreaDataContinental)

# Torrey

sum(MaskTorrey[,]) * 100^2 / 1000^2 # 124, but really 124.5-->125 based on the polygon before coarsening to 100 m resolution

length(which(MaskFracTorrey[,] > 0.5 & MaskTorrey[,] == 0)) * 100^2 / 1000^2
length(which(MaskFracTorrey[,] > 0.1 & MaskTorrey[,] == 0)) * 100^2 / 1000^2
length(which(MaskFracTorrey[,] > 0.01 & MaskTorrey[,] == 0)) * 100^2 / 1000^2

length(which(MaskFracTorrey[,] < 0.5 & MaskTorrey[,] == 1)) * 100^2 / 1000^2
length(which(MaskFracTorrey[,] < 0.9 & MaskTorrey[,] == 1)) * 100^2 / 1000^2
length(which(MaskFracTorrey[,] < 0.99 & MaskTorrey[,] == 1)) * 100^2 / 1000^2

sum(pmax(0,MaskAbsTorrey[,])) / sum((MaskTorrey * AccumMapMod)[,])
sum(pmin(0,MaskAbsTorrey[,])) / sum((MaskTorrey * AccumMapMod)[,])

sum(MaskAbsTorrey[,]) * 100^2 / 1000^3
1000 * sum(MaskAbsTorrey[,]) / sum(MaskTorrey[,])
sum(MaskAbsTorrey[,]) / sum((MaskTorrey * AccumMapMod)[,])

# Continental

sum(MaskContinental[,]) * 100^2 / 1000^2 # 1.98, consistent with 1.95-->2 from polygon

sum(pmax(0,MaskAbsContinental[,])) / sum((MaskContinental * AccumMapMod)[,])

length(which(MaskFracContinental[,] > 0.01 & MaskContinental[,] == 0)) * 100^2 / 1000^2

# This only considers the trans-divide flux that ends up in the catchment
# (sum(pmax(0,MaskAbsContinental[,])) * 100^2 * 1000) / 2700 # Approx. mass (kg) flux per length of Continental Divide along import area

# True trans-divide flux, only picking 1 cell per row (to avoid double-counting the W-to-E flux)
ContinentalDivide = vect("../Extra/ContinentalDivideFluxBoundary.shp")
plot(OutgoingFlux)
plot(ContinentalDivide,add=TRUE)
CDfluxVals = extract(OutgoingFlux,ContinentalDivide)[,2]
mean(CDfluxVals)
mean(CDfluxVals) * 100^2 * 1000 / 100 # m SWE --> m^3 water --> kg --> kg per m along Divide
length(CDfluxVals) * 100 / 1000
sum(CDfluxVals) * 100^2 / 1000^3

sum(CDfluxVals) * 100^2 / as.numeric(difftime(as.POSIXct("6/1/2024",format="%m/%d/%Y"),as.POSIXct("10/1/2023",format="%m/%d/%Y"),units="s"))
# 0.10 m^3/s

TotalImportVolm3 = sum(pmax(0,MaskAbsTorrey[,])) * 100^2
TotalImportVolm3 / 1000^3 # 3.40x10^-3

TorreyStreamDataFull = read.csv("../TorreyCreek_DailyGageRecord_2022-2024.csv")[,-1]
TorreyStreamDataFull$Date = as.Date(TorreyStreamDataFull$Date,format="%Y-%m-%d")
# Subset to WY 2024 only for consistency with 2024 snow accumulation season
TorreyStreamData = TorreyStreamDataFull[TorreyStreamDataFull$Date > as.Date("2023-09-30",format="%Y-%m-%d"),]
plot(TorreyStreamData$Date, TorreyStreamData$Streamflow_cfs,type="l")
mean(TorreyStreamData$Streamflow_cfs) / (3.28084^3) # 1.51 m^3/s
OctToMayPtrs = which(as.numeric(format(TorreyStreamData$Date,"%m")) %in% c(10:12,1:5))
mean(TorreyStreamData$Streamflow_cfs[OctToMayPtrs]) / (3.28084^3) # 0.52 m^3/s

as.numeric(difftime(as.POSIXct("6/1/2024",format="%m/%d/%Y"),as.POSIXct("10/1/2023",format="%m/%d/%Y"),units="d"))
length(OctToMayPtrs) # Same number of days

TotalImportVolm3 / as.numeric(difftime(as.POSIXct("6/1/2024",format="%m/%d/%Y"),as.POSIXct("10/1/2023",format="%m/%d/%Y"),units="s"))
# 0.161 m^3/s (5.7 cfs)
0.1614 / 1.5051 # vs. annual streamflow
0.1614 / 0.5216 # vs. streamflow on equivalent time period (Oct-May)

TotalImportVolm3 / (60*60*24 * sum(TorreyStreamData$Streamflow_cfs) / (3.28084^3))

# Super interesting:
# Imported Snow / (Imported Snow + Net Glacier Mass Loss) = 45%
# In zero-net-mass-loss counterfactual, Torrey Creek overproduces August-September streamflow in Green/Pine/Bull by 45%
# Seems like the extra late-summer streamflow is closely related to the imported snow, which makes sense!

################################################################################
# Sensitivity analysis - sublimation

head(TransportData)

TransportData0pct = read.csv("../Smoothed2x21/AvgDistancePlotData.csv")[,-1]
TransportData1pct = read.csv("../Smoothed2x21_22.5deg_0.01loss/AvgDistancePlotData.csv")[,-1]
TransportData10pct = read.csv("../Smoothed2x21_22.5deg_0.1loss/AvgDistancePlotData.csv")[,-1]

Ptrs = which(TransportData$SWEdepth > 3)

Mean0pct = mean(TransportData0pct$TransportDist[Ptrs] * TransportData0pct$TotalAccum[Ptrs]) / mean(TransportData0pct$TotalAccum[Ptrs])
Mean1pct = mean(TransportData1pct$TransportDist[Ptrs] * TransportData1pct$TotalAccum[Ptrs]) / mean(TransportData1pct$TotalAccum[Ptrs])
Mean10pct = mean(TransportData10pct$TransportDist[Ptrs] * TransportData10pct$TotalAccum[Ptrs]) / mean(TransportData10pct$TotalAccum[Ptrs])

Mean0pct / Mean1pct
Mean10pct / Mean1pct

DistHistData10pct = read.csv("../Smoothed2x21_22.5deg_0.1loss/DistanceHistogramPlotData.csv")[,-1]

1 - DistHistData10pct$CumulativeAccumPct[DistHistData10pct$MaxDepth==6 & DistHistData10pct$DistMax==1000]

MaskAbsTorrey0pct = rast("../Smoothed2x21/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm3_0pct = sum(pmax(0,MaskAbsTorrey0pct[,])) * 100^2 / 1000^3

MaskAbsTorrey1pct = rast("../Smoothed2x21_22.5deg_0.01loss/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm3_1pct = sum(pmax(0,MaskAbsTorrey1pct[,])) * 100^2 / 1000^3

MaskAbsTorrey10pct = rast("../Smoothed2x21_22.5deg_0.1loss/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm3_10pct = sum(pmax(0,MaskAbsTorrey10pct[,])) * 100^2 / 1000^3

################################################################################
# Sensitivity analysis - pattern

# Distance sensitivity

head(TransportData)

TransportData2kmKernel = read.csv("../Smoothed2x21/AvgDistancePlotData.csv")[,-1]
TransportData4kmKernel = read.csv("../Smoothed2x41/AvgDistancePlotData.csv")[,-1]
TransportDataLanczos = read.csv("../GridmetLanczos/AvgDistancePlotData.csv")[,-1]
TransportDataUniform = read.csv("../Uniform/AvgDistancePlotData.csv")[,-1]

Ptrs = which(TransportData$SWEdepth > 3)
MeanRef = mean(TransportData2kmKernel$TransportDist[Ptrs] * TransportData2kmKernel$TotalAccum[Ptrs]) / mean(TransportData2kmKernel$TotalAccum[Ptrs])
Mean4kmKernel = mean(TransportData4kmKernel$TransportDist[Ptrs] * TransportData4kmKernel$TotalAccum[Ptrs]) / mean(TransportData4kmKernel$TotalAccum[Ptrs])
MeanLanczos = mean(TransportDataLanczos$TransportDist[Ptrs] * TransportDataLanczos$TotalAccum[Ptrs]) / mean(TransportDataLanczos$TotalAccum[Ptrs])
MeanUniform = mean(TransportDataUniform$TransportDist[Ptrs] * TransportDataUniform$TotalAccum[Ptrs]) / mean(TransportDataUniform$TotalAccum[Ptrs])

Mean4kmKernel / MeanRef
MeanLanczos / MeanRef
MeanUniform / MeanRef

# Watershed sensitivity

MaskAbsTorrey2kmKernel = rast("../Smoothed2x21/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm3 = sum(pmax(0,MaskAbsTorrey2kmKernel[,])) * 100^2 / 1000^3

MaskAbsTorrey4kmKernel = rast("../Smoothed2x41/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm34kmKernel = sum(pmax(0,MaskAbsTorrey4kmKernel[,])) * 100^2 / 1000^3

MaskAbsTorreyLanczos = rast("../GridmetLanczos/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm3Lanczos = sum(pmax(0,MaskAbsTorreyLanczos[,])) * 100^2 / 1000^3

MaskAbsTorreyUniform = rast("../Uniform/MaskContribTotalDepth_TorreyWatershed.tif")
TotalImportVolkm3Uniform = sum(pmax(0,MaskAbsTorreyUniform[,])) * 100^2 / 1000^3

TotalImportVolkm34kmKernel / TotalImportVolkm3
TotalImportVolkm3Lanczos / TotalImportVolkm3
TotalImportVolkm3Uniform / TotalImportVolkm3

# TransportData$TransportDistDiffLanczos = TransportDataLanczos$TransportDist / TransportData$TransportDist
# 
# plot(TransportData$TransportDist,TransportData$TransportDistDiffLanczos,ylim=c(0,1))
# 
# sd(TransportData$TransportDistDiffLanczos[TransportData$TransportDist > 1000],na.rm=TRUE)





