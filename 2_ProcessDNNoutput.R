# Copyright Eli Boardman 2024
# All Rights Reserved

library(terra)

setwd("C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/")

CounterfactualName = "Smoothed2x21_22.5deg_0.1loss"

AccumMapData = as.matrix(read.table("MassDistribution_Reference.csv",sep=","))
AccumMap = rast("AccumulationReference.tif")
AccumMinusSWE = rast("AccumMinusSWE.tif")
SWEmap = rast("SWEreference.tif")
AccumMapNoDrift = rast(paste0("AccumCounterfactual_",sub("_.*","",CounterfactualName),".tif"))

setwd(CounterfactualName)
AccumMapNoDriftData = as.matrix(read.table("MassDistribution_Counterfactual.csv",sep=","))

EfficiencyVal = as.numeric(readLines("Efficiency.txt"))
ScaleVal = as.numeric(readLines("InputScale.txt"))

AccumMapNoDrift = AccumMapNoDrift * ScaleVal
AccumMapNoDriftData = AccumMapNoDriftData * ScaleVal

################################################################################
# Get results from PyTorch

ModelPred = as.matrix(read.table("FinalStorage.csv",sep=","))
OutgoingFlux = as.matrix(read.table("TotalFlux.csv",sep=","))
OutgoingFrac = as.matrix(read.table("TransferFracs.csv",sep=","))
OutgoingFracUp = as.matrix(read.table("TransferFracs_Up.csv",sep=","))
OutgoingFracAcross = as.matrix(read.table("TransferFracs_Across.csv",sep=","))
OutgoingFracDown = as.matrix(read.table("TransferFracs_Down.csv",sep=","))
NetLoss = as.matrix(read.table("NetLoss.csv",sep=","))

max(rbind(OutgoingFracUp/OutgoingFracAcross,OutgoingFracDown/OutgoingFracAcross),na.rm=TRUE)
min(rbind(OutgoingFracUp/OutgoingFracAcross,OutgoingFracDown/OutgoingFracAcross),na.rm=TRUE)

image(t(apply(log10(OutgoingFlux + 1), 2, rev)))

sqrt(mean(unlist(ModelPred[,] - AccumMapData[,])^2))
1 - sum((ModelPred[,] - AccumMapData[,])^2) / sum((AccumMapData[,] - mean(AccumMapData[,]))^2)
Ptrs = which(AccumMapData[,] > 1)
1 - sum((ModelPred[Ptrs] - AccumMapData[Ptrs])^2) / sum((AccumMapData[Ptrs] - mean(AccumMapData[Ptrs]))^2)

AccumMapMod = 0 * AccumMap
OutgoingFluxMap = 0 * AccumMap
OutgoingFracMap = 0 * AccumMap
OutgoingFracMap_Up = 0 * AccumMap
OutgoingFracMap_Across = 0 * AccumMap
OutgoingFracMap_Down = 0 * AccumMap
NetLossMap = 0 * AccumMap
AccumMapMod[,] = as.numeric(t(ModelPred))
OutgoingFluxMap[,] = as.numeric(t(OutgoingFlux))
OutgoingFracMap[,] = as.numeric(t(OutgoingFrac))
OutgoingFracMap_Up[,] = as.numeric(t(OutgoingFracUp))
OutgoingFracMap_Across[,] = as.numeric(t(OutgoingFracAcross))
OutgoingFracMap_Down[,] = as.numeric(t(OutgoingFracDown))
NetLossMap[,] = as.numeric(t(NetLoss))

ErrMap = AccumMapMod - AccumMap
plot(ErrMap,range=c(-0.16,0.16))
plot(abs(ErrMap),range=c(0,0.16))
cor(ErrMap[,],AccumMapMod[,])

plot(AccumMap,range=c(0,5.5))
plot(AccumMapMod,range=c(0,5.5))
plot((AccumMapMod - AccumMap),range=c(-1,1))
hist(abs(AccumMapMod - AccumMap))
plot(ModelPred[,], AccumMapData[,])
plot(AccumMapData[,], ModelPred[,] - AccumMapData[,])

# Bias vs. depth
SWEerrs = ModelPred[,] - AccumMapData[,]
SWEdepths = AccumMapData[,]
BiasSpacing = 0.1
SWEerrBiasData = data.frame(DepthMin=seq(0,6-BiasSpacing,BiasSpacing),
                            DepthMax=seq(BiasSpacing,6,BiasSpacing),
                            MeanErr=NA,
                            RMSE=NA)
for (i in 1:nrow(SWEerrBiasData)){
  SWEerrsSubset = SWEerrs[ModelPred > SWEerrBiasData$DepthMin[i] & 
                            ModelPred < SWEerrBiasData$DepthMax[i]]
  SWEerrBiasData$MeanErr[i] = mean(SWEerrsSubset)
  SWEerrBiasData$RMSE[i] = sqrt(mean(SWEerrsSubset^2))
}
mean(abs(SWEerrBiasData$MeanErr),na.rm=TRUE)
mean(SWEerrBiasData$RMSE,na.rm=TRUE)
plot(SWEerrBiasData$DepthMin,SWEerrBiasData$MeanErr,type="o")
lines(c(-10,10),c(0,0))
plot(SWEerrBiasData$DepthMin,SWEerrBiasData$RMSE,type="o")
lines(c(-10,10),c(0.12,0.12))

par(mfrow=c(2,2))
plot(OutgoingFracMap,range=c(0,1))
plot(OutgoingFracMap_Up,range=c(0,1))
plot(OutgoingFracMap_Across,range=c(0,1))
plot(OutgoingFracMap_Down,range=c(0,1))
par(mfrow=c(1,1))

AcrossRatio = abs(OutgoingFracMap_Up - OutgoingFracMap_Down) / OutgoingFracMap_Across
OutgoingAngleMap = atan(AcrossRatio / (1 + AcrossRatio)) * 180 / pi
OutgoingAngleMap = OutgoingAngleMap * sign(OutgoingFracMap_Up - OutgoingFracMap_Down)
plot(OutgoingAngleMap)
plot(OutgoingFluxMap)

SWEmapMod = AccumMapMod - AccumMinusSWE
SWEmapMod[,] = pmax(0, SWEmapMod[,])

writeRaster(AccumMapMod,"DNN_AccumFinal.tif",overwrite=TRUE)
writeRaster(SWEmapMod,"DNN_SWEfinal.tif",overwrite=TRUE)
writeRaster(SWEmapMod - SWEmap,"DNN_SWEerror.tif",overwrite=TRUE)
writeRaster(OutgoingFluxMap,"DNN_TransferFlux.tif",overwrite=TRUE)
writeRaster(OutgoingFracMap,"DNN_TransferFracTotal.tif",overwrite=TRUE)
writeRaster(OutgoingFracMap_Across,"DNN_TransferFrac_X.tif",overwrite=TRUE)
writeRaster(OutgoingFracMap_Up,"DNN_TransferFrac_Yup.tif",overwrite=TRUE)
writeRaster(OutgoingFracMap_Down,"DNN_TransferFrac_Ydown.tif",overwrite=TRUE)
writeRaster(OutgoingAngleMap,"DNN_TransferAngle.tif",overwrite=TRUE)
writeRaster(NetLossMap,"DNN_NetLoss.tif",overwrite=TRUE)

################################################################################
# Validation with slow but easy to understand version of the network structure

if (FALSE){
  
  nX = ncol(AccumMap)
  nY = nrow(AccumMap)
  
  InputVals = AccumMapNoDriftData
  FinalVals = InputVals * 0
  TransferVals = InputVals * 0
  FluxVals = InputVals * 0
  
  for (x in 1:nX){
    for (y in 1:nY){
      
      if (x > 1){
        TotalAvail = InputVals[y,x] + TransferVals[y,x-1]
      } else {
        TotalAvail = InputVals[y,x]
      }
      
      OutgoingFracVal = OutgoingFracUp[y,x] + OutgoingFracAcross[y,x] + OutgoingFracDown[y,x]
      FinalVals[y,x] = TotalAvail * (1 - OutgoingFracVal)
      
      if (y > 1){
        TransferVals[(y-1),x] = TransferVals[(y-1),x] + TotalAvail * OutgoingFracUp[y,x] * EfficiencyVal
      }
      TransferVals[y,x] = TransferVals[y,x] + TotalAvail * OutgoingFracAcross[y,x] * EfficiencyVal
      if (y < nY){
        TransferVals[(y+1),x] = TransferVals[(y+1),x] + TotalAvail * OutgoingFracDown[y,x] * EfficiencyVal
      }
    }
  }
  
  # Mass balance
  TotalInput = sum(InputVals[,])
  FinalStorage = sum(FinalVals[,])
  LeftDomain = sum(TransferVals[,nX])
  TotalFlux = sum(TransferVals[,])
  print(paste("In:",round(TotalInput),"Storage:",round(FinalStorage),
              "Left:",round(LeftDomain),
              "Loss:",round(TotalInput - FinalStorage - LeftDomain),
              "Flux:",round(TotalFlux)))
  
  TotalFlux / TotalInput
  
  image(t(apply(log10(TransferVals + 1), 2, rev)))
  image(t(apply(FinalVals, 2, rev)))
  
  Diff = ModelPred - FinalVals
  image(t(apply(Diff, 2, rev)))
  hist(unlist(Diff))
  
  DiffTransfer = OutgoingFlux - TransferVals
  image(t(apply(DiffTransfer, 2, rev)))
  hist(unlist(DiffTransfer))
}

################################################################################
# Plot progress

if (FALSE){
  
  AvailEpochs = list.files("EpochProgressMaps/")
  AvailEpochs = AvailEpochs[grepl("FinalStorage",AvailEpochs)]
  AvailEpochs = as.numeric(sub("Epoch","",sub("_FinalStorage.csv","",AvailEpochs)))
  AvailEpochs = sort(AvailEpochs)
  
  AccumMapMod = 0 * AccumMap
  OutgoingFluxMap = 0 * AccumMap
  
  par(mfrow=c(2,1))
  for (nEpoch in 30000){
    ModelPred = as.matrix(read.table(paste0("EpochProgressMaps/Epoch",nEpoch,"_FinalStorage.csv"),sep=","))
    OutgoingFlux = as.matrix(read.table(paste0("EpochProgressMaps/Epoch",nEpoch,"_TotalFlux.csv"),sep=","))
    
    AccumMapMod[,] = as.numeric(t(ModelPred))
    OutgoingFluxMap[,] = pmin(14,as.numeric(t(OutgoingFlux)))
    
    plot(AccumMapMod,main=paste0("Storage Map Epoch ",nEpoch),range=c(0,6))
    plot(OutgoingFluxMap,main=paste0("Transport Flux Epoch ",nEpoch),range=c(0,14))
  }
  par(mfrow=c(1,1))
  
  TrainingProgressData = read.csv("TrainingProgress.csv")
  NudgePtrs = which(TrainingProgressData$RMSE < 0.1 & TrainingProgressData$RMSEbinned == 0 & TrainingProgressData$Bias < 0.01)
  plot(TrainingProgressData$Epoch,TrainingProgressData$FluxCost,type="l",ylim=c(0,2))
  lines(c(-1e9,1e9),c(0.9,0.9))
  abline(v=20000)
  abline(v=TrainingProgressData[NudgePtrs,"Epoch"],
         col="yellow")
  plot(TrainingProgressData$Epoch[10001:nrow(TrainingProgressData)],
       TrainingProgressData$FluxCost[10001:nrow(TrainingProgressData)] /
         TrainingProgressData$FluxCost[1:(nrow(TrainingProgressData)-10000)],
       type="l",ylim=c(0.95,1.01))
  lines(c(-1e9,1e9),c(1,1))
  
  plot(TrainingProgressData$Epoch,TrainingProgressData$RMSE,type="l")
  abline(v=TrainingProgressData[NudgePtrs,"Epoch"],
         col="yellow")
  plot(TrainingProgressData$RMSEbinned,type="l",ylim=c(0,0.165))
  plot(TrainingProgressData$RMSEbinned / TrainingProgressData$RMSE,type="l",ylim=c(2,2.5))
  plot(TrainingProgressData$Epoch,TrainingProgressData$Bias,type="l",ylim=c(0,0.1))
  lines(c(-1e9,1e9),c(0.01,0.01))
  
  # Pareto frontier
  FluxOrdered = sort(TrainingProgressData$FluxCost)
  RMSEordered = TrainingProgressData$RMSE[order(TrainingProgressData$FluxCost)]
  plot(FluxOrdered,RMSEordered,type="l")
}








