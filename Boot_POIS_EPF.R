###############################################################################
###############################################################################
###
###       Created on Monday April 11 2022
###       Directory  /Users/mariabugalloporto/Desktop/UMH_Elche/MariaBugallo/Datos_EPF2006_2016                                                 
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(glmmTMB)
options(warn=-1)

set.seed(2509)

# Bootstrap Confidence intervals for model parameters
# Relative root-mean square error (RRMSE)

########################################
######## PART 1 
########################################

## Datos EPF2016 + 4 EPA2016 trimestrales
EPF_2016aux=read.csv("EPF_EPA_2016.csv", sep=',', dec='.', header=T) 

EPF_2016aux$prov=factor(EPF_2016aux$prov)
EPF_2016aux$sex=factor(EPF_2016aux$sex)
EPF_2016aux$age4=factor(EPF_2016aux$age4)
attach(EPF_2016aux)

fit_zipoisson0<-glmmTMB(floor(h1total)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=poisson)

x1=as.matrix(data.frame(rep(1,416)))
beta1=as.matrix((data.frame(fit_zipoisson0$sdr$par.fixed[5])))
phi1=sqrt(VarCorr(fit_zipoisson0)[[2]][[1]][[1]])

x2=as.matrix(data.frame(rep(1,416), EPF_2016aux[c(13, 19, 20)]))
beta2=as.matrix((data.frame(fit_zipoisson0$sdr$par.fixed[1:4])))
phi2=sqrt(VarCorr(fit_zipoisson0)[[1]][[1]][[1]])


########################################
######## PART 2: SIMULATION 0 
########################################

t <- proc.time()

R=1000 # R=10^3 SIMULACIONES
u1_2r=matrix(ncol=416, nrow=R)
u2_2r=matrix(ncol=416, nrow=R)

p_ijk=matrix(ncol=416, nrow=R); lambda_ijk=matrix(ncol=416, nrow=R); mu_ijk=matrix(ncol=416, nrow=R)
z_ijk=matrix(ncol=416, nrow=R); y_ijk=matrix(ncol=416, nrow=R)
# y_ijkPOIS=matrix(ncol=416, nrow=R)

for (r in 1:R){  	
	u1_2r[r,]=rep(rnorm(4,mean=0,sd=1), each=104)
	u2_2r[r,]=rnorm(416,mean=0,sd=1) 	
	}

# Generamos los datos de nuestro modelo
for (l in 1:416){
	for (r in 1:R){  	
		# 1.2 Generate p_ijk^r and lambda_ijk^r
		p_ijk[r,l]=exp(beta1+phi1*u1_2r[r,l]) / (1+exp(beta1+phi1*u1_2r[r,l]))
		lambda_ijk[r,l]=exp(x2[l,]%*%beta2+phi2*u2_2r[r,l])
		mu_ijk[r,l]=Nd[l]*(1-p_ijk[r,l])*lambda_ijk[r,l]
		
		z_ijk[r,l]<-rbinom(1,size=1, prob=p_ijk[r,l])   # z_ijk* Be(p_ijk[r,l])
		
		# y_ijkPOIS[r,l]<-rpois(1,lambda=Nd[l]*lambda_ijk[r,l])
		
		if(z_ijk[r,l]==1){
			y_ijk[r,l]=0
		} else{
			y_ijk[r,l]=rpois(1,lambda=Nd[l]*lambda_ijk[r,l]) }		}	
			}	

# count_rAZIP=rep(NA, R); count_rPOIS=rep(NA, R)
# for (r in 1:R){ 
	# count_rPOIS[R] = length(y_ijkPOIS[r, y_ijkPOIS[r,]==0 ])
	# count_rAZIP[r] = length(y_ijk[r, y_ijk[r,]==0 ])
	# }

beta1r=matrix(ncol=1, nrow=R)
beta2r=matrix(ncol=4, nrow=R)
phi1r=matrix(ncol=1, nrow=R)
phi2r=matrix(ncol=1, nrow=R)

mu_PIr=matrix(ncol=416, nrow=R)

# Realizamos el estudio de simulacion en si

for (r in 1:R){  
	print(r)
	EPF_2016aux$y_ijk_sim=y_ijk[r,]
	fit_zipoissonr<-glmmTMB(floor(y_ijk_sim)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=poisson)
	
	beta1r[r,1]=as.matrix((data.frame(fit_zipoissonr$sdr$par.fixed[5])))
	beta2r[r,]=as.matrix((data.frame(fit_zipoissonr$sdr$par.fixed[1:4])))

	phi1r[r]=sqrt(VarCorr(fit_zipoissonr)[[2]][[1]][[1]])
	phi2r[r]=sqrt(VarCorr(fit_zipoissonr)[[1]][[1]][[1]])

	mu_PIr[r,]=fitted(fit_zipoissonr)
}

proc.time()-t

# CI Bootstrap
round(quantile(beta1r, prob=c(0.025, 0.975)),3)
round(quantile(beta2r[,1], prob=c(0.025, 0.975)),3)
round(quantile(beta2r[,2], prob=c(0.025, 0.975)),3)
round(quantile(beta2r[,3], prob=c(0.025, 0.975)),3)
round(quantile(beta2r[,4], prob=c(0.025, 0.975)),3)

round(quantile(phi1r, prob=c(0.025, 0.975)),3)
round(quantile(phi2r, prob=c(0.025, 0.975)),3)

# RMSE and RRMSE PREDICCIONES
rmse_muPIr=sqrt(colMeans( t(t(mu_PIr-mu_ijk)/Nd)^2 ,na.rm=T))

rrmse_muPIr=100*rmse_muPIr/(fitted(fit_zipoisson0)/Nd)

# Hajek
rmse_Hajek=sqrt(colMeans(t(t(mu_ijk)/Nd-t(y_ijk)/Nd)^2,na.rm=T))

rrmse_Hajek=100*rmse_Hajek/(fitted(fit_zipoisson0)/Nd)

# Save information
df=data.frame(rmse_muPIr, rrmse_muPIr, rmse_Hajek, rrmse_Hajek)
colnames(df)=c('RMSE_PlugIN', 'RRMSE_PlugIN', 'RMSE_Hajek', 'RRMSE_Hajek')
write.csv(df,"rmse_PlugIN.csv", row.names = FALSE)


########################################
######## PART 3: Plots
########################################

EPF_2016aux=read.csv("EPF_EPA_2016.csv", sep=',', dec='.', header=T) 

PlugIN=read.csv("rmse_PlugIN.csv", header=T, sep=',')
rmse_muPIr=PlugIN$RMSE_PlugIN
rrmse_muPIr=PlugIN$RRMSE_PlugIN

rmse_muHajek=PlugIN$RMSE_Hajek
rrmse_muHajek=PlugIN$RRMSE_Hajek

i=sort(nd, index.return=T)$ix

# 1. RMSE - RVAR
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(sqrt(h1var)[EPF_2016aux$sex==1], rmse_muPIr[EPF_2016aux$sex==1], main=' ', xlab='RVAR Hajek', 
	ylab='RMSE IN', col='darkblue', pch=19, cex.main=1.7, cex.lab=1.7, cex.axis=1.7, xlim=c(0, 0.3), ylim=c(0,0.3))
points(sqrt(h1var)[EPF_2016aux$sex==6], rmse_muPIr[EPF_2016aux$sex==6], pch=3, col='darkblue')	
new <- data.frame(x = seq(0,0.3,0.01))
lines(new$x, new$x, col='red', lwd=3, lty=2)	
legend('topleft', col='darkblue', c('Men', 'Women'), cex=1.5, pt.cex = 1.4, pch=c(19,3))


# 2. RMSE - RMSE
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(rmse_muHajek[EPF_2016aux$sex==1], rmse_muPIr[EPF_2016aux$sex==1], main=' ', xlab='RMSE Hajek', 
	ylab='RMSE IN', col='darkblue', pch=19, cex.main=1.7, cex.lab=1.7, cex.axis=1.7, xlim=c(0, 0.3), ylim=c(0,0.3))
points(rmse_muHajek[EPF_2016aux$sex==6], rmse_muPIr[EPF_2016aux$sex==6], pch=3, col='darkblue')	
new <- data.frame(x = seq(0,0.3,0.01))
lines(new$x, new$x, col='red', lwd=3, lty=2)	


