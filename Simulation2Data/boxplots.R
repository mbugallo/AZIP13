###############################################################################
###############################################################################
###
###       Created on Thursday May 19 2022
###       Directory  /Users/mariabugalloporto/Desktop/UMH_Elche/MariaBugallo/Datos EPF2014-16/simu2                                              
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(RColorBrewer)

B50=read.csv("simu2_B50.csv")
B100=read.csv("simu2_B100.csv")
B200=read.csv("simu2_B200.csv")
B400=read.csv("simu2_B400.csv")
B500=read.csv("simu2_B500.csv")
B600=read.csv("simu2_B600.csv")


# RB
boxplot(B100$RBijk_100, B200$RBijk_200, B400$RBijk_400, B500$RBijk_500, B600$RBijk_600, xlab=' ', ylab=' ', 
	xaxt='n', yaxt='n', pch=20, col=rep('white',5), ylim=c(-60,145))
axis(1, cex.axis=1.6, at=c(0.4,seq(1,5, by=1),5.6), labels=c('','B=100', 'B=200','B=400', 'B=500', 'B=600',''))
axis(2, cex.axis=1.6)

# RRE
boxplot(B100$RREijk_100, B200$RREijk_200, B400$RREijk_400, B500$RREijk_500, B600$RREijk_600, xlab=' ', ylab=' ', 
	xaxt='n', yaxt='n', pch=20, col=rep('white',5), ylim=c(0,500))
axis(1, cex.axis=1.6, at=c(0.4,seq(1,5, by=1),5.6), labels=c('','B=100', 'B=200','B=400', 'B=500', 'B=600',''))
axis(2, cex.axis=1.6)

# 

MSE_PI=read.csv("MSE_PI.csv")$rmse_muPIr^2

sum(abs(B600$RBijk_600*MSE_PI/100)) # AB

sum(abs(B600$RREijk_600*MSE_PI/100)) # RE