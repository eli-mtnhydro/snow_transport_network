# Copyright Eli Boardman 2024
# All Rights Reserved

library(ggplot2)
library(ggh4x)

setwd("C:/Users/board/Desktop/DriftNetwork/Data/TorreySnowRegion/")

CounterfactualName = "Smoothed2x21_22.5deg_0.01loss"

setwd(CounterfactualName)

TransportDataPlot = read.csv("AvgDistancePlotData.csv")[,-1]
DistHistDataPlot = read.csv("DistanceHistogramPlotData.csv")[,-1]

TransportDataPlot$DepthRange = paste0(floor(TransportDataPlot$SWEdepth)," - ",floor(TransportDataPlot$SWEdepth)+1," m")

DistHistDataPlot$DepthRange = paste0(" ",DistHistDataPlot$DepthRange)

################################################################################

g1 = ggplot() + 
  geom_point(data=TransportDataPlot,
             aes(x=SWEdepth, y=TransportDist, fill=Elevation, size=SWEdepth),
             shape=21,
             color="black",alpha=0.75,stroke=0.5) +
  scale_size(breaks=seq(1,5),limits=c(0,6),range=c(0,6)) +
  scale_fill_gradientn(colors=c("black","#0a0054","#570f6e","#e6315b","yellow","white"),
                       limits=c(3000,4000),breaks=seq(3000,4000,250),expand=c(0,0)) +
  scale_x_continuous(limits=c(0,5.5),breaks=seq(0,5),
                     # labels=c("0",paste0(breaks=seq(1,5)," m")),
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(0,2200),breaks=seq(0,2000,500),
                     # labels=c("0",paste0(breaks=seq(1,2,1)," km")),
                     labels=seq(0,2,0.5),
                     expand=c(0,0),
                     minor_breaks=seq(0,2500,100), guide="axis_minor") +
  guides(y.sec=guide_axis(title=NULL)) +
  labs(x="Final SWE, m",
       y="Mean Contributing Dist., km",
       fill="Elevation, m") +
  theme_bw() +
  theme(axis.text.x=element_text(color="black",size=32),
        axis.text.y=element_text(color="black",size=32),
        axis.text.y.right=element_blank(),
        axis.title.x=element_text(size=40,face="bold",margin=margin(0.25,0,0,0,"cm")),
        axis.title.y=element_text(size=40,face="bold",margin=margin(0,0.25,0,0,"cm")),
        axis.ticks.length=unit(0.5,"cm"),
        axis.ticks=element_line(color="black",linewidth=1),
        plot.margin=margin(1,0.25,0.5,0.5,"cm"),
        panel.border=element_rect(fill=NA,color="black",linewidth=2),
        panel.grid=element_blank(),
        panel.spacing=unit(1,"cm"),
        legend.title=element_text(size=52,face="bold",margin=margin(0,0,0.25,0,"cm")),
        legend.title.align=0.5,
        legend.text=element_text(size=40,hjust=0.5),
        legend.key.width=unit(2,"cm"),
        legend.spacing.y=unit(0.5,"cm"),
        legend.text.align=0,
        legend.margin=margin(0,0,0,0.25,"cm"),
        plot.title=element_blank()) +
  # guides(fill=guide_colorbar(barheight=unit(4,"cm"),ticks.colour=NA,byrow=TRUE,order=1),
  #        size="none")
  guides(fill="none",
         size="none")

print(g1)

ggsave("DepthVsTravelDistance.png", plot=g1,
       width=unit(4000/300,"in"), height=unit(2750/300,"in"), dpi=300)

g2 = ggplot() + 
  geom_violin(data=data.frame(DepthRange=sub(" m","",TransportDataPlot$DepthRange),
                              ImportArea=c(TransportDataPlot$ImportAreaMin1mm,
                                           TransportDataPlot$ImportAreaMin10mm,
                                           TransportDataPlot$ImportAreaMin100mm),
                              MinContrib=c(rep(" >1 mm",nrow(TransportDataPlot)),
                                           rep(" >10 mm",nrow(TransportDataPlot)),
                                           rep(" >100 mm SWE",nrow(TransportDataPlot)))),
              aes(x=DepthRange, y=ImportArea, fill=MinContrib),
              position="dodge", scale="width", width=0.9, linewidth=1, color="black") +
  scale_fill_manual(values=c("gray90","gray50","black")) +
  scale_y_continuous(limits=c(0,4),breaks=seq(0,4,1),
                     # labels=c("0",expression("1 km"^2),expression("2 km"^2),expression("3 km"^2)),
                     expand=c(0,0),
                     minor_breaks=seq(0,4,0.25), guide="axis_minor") +
  guides(y.sec=guide_axis(title=NULL)) +
  labs(x="Final SWE, m",
       y=expression(bold("Contributing Source Area, km"^2)),
       fill="Min Contribution") +
  theme_bw() +
  theme(axis.text.x=element_text(color="black",size=32),
        axis.text.y=element_text(color="black",size=32),
        axis.text.y.right=element_blank(),
        axis.title.x=element_text(size=40,face="bold",margin=margin(0.25,0,0,0,"cm")),
        axis.title.y=element_text(size=40,face="bold",margin=margin(0,0.5,0,0,"cm")),
        axis.ticks.length=unit(0.5,"cm"),
        axis.ticks=element_line(color="black",linewidth=1),
        plot.margin=margin(1,0.25,0.5,0.5,"cm"),
        plot.background=element_rect(fill="transparent",color=NA_character_),
        panel.border=element_rect(fill=NA,color="black",linewidth=2),
        panel.grid=element_blank(),
        panel.spacing=unit(1,"cm"),
        legend.title=element_blank(),
        legend.title.align=0.5,
        legend.text=element_text(size=32,hjust=0.5),
        legend.key.width=unit(1.5,"cm"),
        legend.key.spacing.x=unit(1.5,"cm"),
        legend.text.align=0,
        legend.margin=margin(0,0,0,0.25,"cm"),
        legend.position=c(0.47,0.9),
        legend.background=element_blank(),
        plot.title=element_blank()) +
  guides(fill=guide_legend(byrow=TRUE,
                            nrow=1,title.position="top",order=1,reverse=FALSE))

print(g2)

ggsave("AccumulationAreaDistributions.png", plot=g2,
       width=unit(4000/300,"in"), height=unit(3000/300,"in"), dpi=300)

g3 = ggplot() + 
  geom_line(data=DistHistDataPlot,
            aes(x=TransportDist, y=CumulativeAccumPct, color=DepthRange),
            linewidth=3,lineend="round") +
  scale_color_manual(values=c("lightsteelblue2","skyblue2","steelblue1","dodgerblue3","navy","black")) +
  scale_x_continuous(limits=c(0,3001),breaks=seq(0,3000,500),
                     labels=c(seq(0,2.5,0.5),">3  "),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1.001),breaks=seq(0,1,0.1),
                     labels=paste0(seq(0,1,0.1)*100,"%"),expand=c(0,0),
                     sec.axis=dup_axis(name=NULL,labels=NULL)) +
  labs(x="Contributing Distance, km",
       y="Fraction of Net Accumulation",
       color="Final SWE") +
  theme_bw() +
  theme(axis.text.x=element_text(color="black",size=32),
        axis.text.y=element_text(color="black",size=32),
        axis.title.x=element_text(size=40,face="bold",margin=margin(0.25,0,0,0,"cm")),
        axis.title.y=element_text(size=40,face="bold",margin=margin(0,0,0,0,"cm")),
        axis.ticks.length=unit(0.5,"cm"),
        axis.ticks=element_line(color="black",linewidth=1),
        plot.margin=margin(1,1,0.5,0.5,"cm"),
        plot.background=element_rect(fill="transparent",color=NA_character_),
        panel.border=element_rect(fill=NA,color="black",linewidth=2),
        panel.grid=element_blank(),
        panel.spacing=unit(4,"cm"),
        legend.title=element_text(size=40,face="bold",margin=margin(0,0,0.5,0,"cm"),hjust=0.5,vjust=0.5),
        legend.title.align=0.5,
        legend.text=element_text(size=32,hjust=0.5,vjust=0.5),
        legend.key.width=unit(2.5,"cm"),
        legend.key.spacing.y=unit(0.5,"cm"),
        legend.text.align=0,
        legend.margin=margin(1,0,0,0,"cm"),
        legend.position=c(0.7,0.4),
        legend.background=element_blank(),
        plot.title=element_blank()) +
  guides(color=guide_legend(override.aes=c(linewidth=4),byrow=TRUE,
                            nrow=6,title.position="top",order=1,reverse=FALSE))

print(g3)

ggsave("CumulativeAccumulationFrac.png", plot=g3,
       width=unit(4000/300,"in"), height=unit(3000/300,"in"), dpi=300)

################################################################################
# Bonus plots

g2 = ggplot() + 
  geom_line(data=DistHistDataPlot[DistHistDataPlot$DistMin > 0,],
            aes(x=TransportDist, y=DistContrib, color=DepthRange),
            linewidth=3,lineend="round") +
  scale_color_manual(values=c("lightsteelblue2","skyblue2","steelblue1","dodgerblue3","navy","black")) +
  scale_x_continuous(limits=c(0,3000),breaks=seq(0,3000,500),
                     labels=c(seq(0,2500,500),">3000  "),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.42),breaks=seq(0,0.5,0.1),expand=c(0,0),
                     sec.axis=dup_axis(name=NULL,labels=NULL)) +
  labs(x="Transport Distance, m",
       y="Imported SWE, m",
       color="  Final SWE:  ") +
  theme_bw() +
  theme(axis.text.x=element_text(color="black",size=32),
        axis.text.y=element_text(color="black",size=32),
        axis.title.x=element_text(size=40,face="bold",margin=margin(0.5,0,0,0,"cm")),
        axis.title.y=element_text(size=40,face="bold",margin=margin(0,0,0,0,"cm")),
        axis.ticks.length=unit(0.5,"cm"),
        axis.ticks=element_line(color="black",linewidth=1),
        plot.margin=margin(0.25,1,0.25,0.5,"cm"),
        panel.border=element_rect(fill=NA,color="black",linewidth=2),
        panel.grid=element_blank(),
        panel.spacing=unit(4,"cm"),
        legend.title=element_text(size=52,face="bold",margin=margin(0,0,0,0,"cm"),hjust=0.5,vjust=0.5),
        legend.title.align=0.5,
        legend.text=element_text(size=40,hjust=0.5,vjust=0.5),
        legend.key.width=unit(2,"cm"),
        legend.spacing.x=unit(0.25,"cm"),
        legend.text.align=0,
        legend.margin=margin(1,0,0,0,"cm"),
        legend.position="bottom",
        plot.title=element_blank()) +
  guides(color=guide_legend(override.aes=c(linewidth=4),byrow=TRUE,
                            nrow=1,title.position="left",order=1,reverse=FALSE))

print(g2)

g3 = ggplot() + 
  geom_line(data=DistHistDataPlot,
            aes(x=TransportDist, y=CumulativeAccumPct, color=DepthRange),
            linewidth=3,lineend="round") +
  scale_color_manual(values=c("lightsteelblue2","skyblue2","steelblue1","dodgerblue3","navy","black")) +
  scale_x_continuous(limits=c(0,3001),breaks=seq(0,3000,500),
                     labels=c(seq(0,2500,500),">3000  "),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,1.001),breaks=seq(0,1,0.1),
                     labels=paste0(seq(0,1,0.1)*100,"%"),expand=c(0,0),
                     sec.axis=dup_axis(name=NULL,labels=NULL)) +
  labs(x="Transport Distance, m",
       y="Fraction of Accumulation",
       color="  Final SWE:  ") +
  theme_bw() +
  theme(axis.text.x=element_text(color="black",size=32),
        axis.text.y=element_text(color="black",size=32),
        axis.title.x=element_text(size=40,face="bold",margin=margin(0.5,0,0,0,"cm")),
        axis.title.y=element_text(size=40,face="bold",margin=margin(0,0,0,0,"cm")),
        axis.ticks.length=unit(0.5,"cm"),
        axis.ticks=element_line(color="black",linewidth=1),
        plot.margin=margin(0.25,1,0.25,0,"cm"),
        panel.border=element_rect(fill=NA,color="black",linewidth=2),
        panel.grid=element_blank(),
        panel.spacing=unit(4,"cm"),
        legend.title=element_text(size=52,face="bold",margin=margin(0,0,0,0,"cm"),hjust=0.5,vjust=0.5),
        legend.title.align=0.5,
        legend.text=element_text(size=40,hjust=0.5,vjust=0.5),
        legend.key.width=unit(2,"cm"),
        legend.spacing.x=unit(0.25,"cm"),
        legend.text.align=0,
        legend.margin=margin(1,0,0,0,"cm"),
        legend.position="bottom",
        plot.title=element_blank()) +
  guides(color=guide_legend(override.aes=c(linewidth=4),byrow=TRUE,
                            nrow=1,title.position="left",order=1,reverse=FALSE))

print(g3)

g = patchwork::wrap_plots(list(g2,g3), nrow=1, guides="collect") &
  theme(legend.position="bottom")
print(g)

ggsave("SnowTransportComboPlot.png", plot=g,
       width=unit(26.667,"in"), height=unit(10,"in"), dpi=300)

g4 = ggplot() + 
  geom_line(data=DistHistDataPlot[DistHistDataPlot$CumulativeAccumPct <= 0.99 &
                                    DistHistDataPlot$CumulativeAccumPct > 0.90,],
            aes(x=TransportDist, y=CumulativeAccum, color=DepthRange, linetype="90% - 99%"),
            linewidth=2,lineend="round") +
  geom_line(data=DistHistDataPlot[DistHistDataPlot$CumulativeAccumPct <= 0.90 &
                                    DistHistDataPlot$CumulativeAccumPct > 0.75,],
            aes(x=TransportDist, y=CumulativeAccum, color=DepthRange, linetype="75% - 90%"),
            linewidth=2,lineend="round") +
  geom_line(data=DistHistDataPlot[DistHistDataPlot$CumulativeAccumPct <= 0.75,],
            aes(x=TransportDist, y=CumulativeAccum, color=DepthRange, linetype="0% - 75%"),
            linewidth=2,lineend="round") +
  scale_color_manual(values=c("lightsteelblue2","skyblue2","steelblue1","dodgerblue3","navy","black")) +
  scale_linetype_manual(values=c("0% - 75%"="solid",
                                 "75% - 90%"="longdash",
                                 "90% - 99%"="12")) +
  scale_x_continuous(limits=c(0,3000),breaks=seq(0,3000,500),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,5.5),breaks=seq(0,6,1),expand=c(0,0)) +
  labs(x="Transport Distance, m",
       y="Cumulative Accumulation, m",
       color="  Final SWE",
       linetype="\nFraction of Total") +
  theme_bw() +
  theme(axis.text.x=element_text(color="black",size=24),
        axis.text.y=element_text(color="black",size=24),
        axis.title.x=element_text(size=32,margin=margin(0.5,0,0,0,"cm")),
        axis.title.y=element_text(size=32,margin=margin(0,0.5,0,0,"cm")),
        axis.ticks.length=unit(0.5,"cm"),
        axis.ticks=element_line(color="black",linewidth=1),
        plot.margin=margin(1,0.5,0.5,0.5,"cm"),
        panel.background=element_rect("white","black", linewidth=2),
        panel.grid=element_blank(),
        panel.spacing=unit(1,"cm"),
        legend.title=element_text(size=32,margin=margin(0,0,0.25,0,"cm")),
        legend.text=element_text(size=24,hjust=0.5),
        legend.key.width=unit(2,"cm"),
        legend.spacing.y=unit(0.5,"cm"),
        legend.text.align=0,
        legend.margin=margin(0,0,0,0.25,"cm"),
        plot.title=element_blank()) +
  guides(color=guide_legend(override.aes=c(linewidth=4),byrow=TRUE,order=1,reverse=TRUE),
         linetype=guide_legend(override.aes=c(linewidth=2),byrow=TRUE,order=2))

print(g4)

ggsave("Histogram_CumulativeTotal.png", plot=g4,
       width=unit(12,"in"), height=unit(8,"in"), dpi=300)











