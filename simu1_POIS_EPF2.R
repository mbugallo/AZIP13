###############################################################################
###############################################################################
###
###       Created on Monday April 11 2022
###       Directory  /Users/mariabugalloporto/Desktop/UMH_Elche/MariaBugallo/Datos_EPF2006_2016                                                 
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(glmmTMB)
library(sae)
options(warn=-1)
source('EBP_ESP.R')

set.seed(1605)

# Simulations under the EPF2016 scenario

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
# beta1=log(1/0.8-1) # p=0.2
# beta1=0 # p =0.5
phi1=sqrt(VarCorr(fit_zipoisson0)[[2]][[1]][[1]])

x2=as.matrix(data.frame(rep(1,416), EPF_2016aux[c(13, 19, 20)]))
beta2=as.matrix((data.frame(fit_zipoisson0$sdr$par.fixed[1:4])))
phi2=sqrt(VarCorr(fit_zipoisson0)[[1]][[1]][[1]])

S=2000

########################################
######## PART 2: SIMULATION 1 
########################################

# Check the behaviour of the fitting algoritms

t <- proc.time()

R=1000
u1_2r=matrix(ncol=416, nrow=R)
u2_2r=matrix(ncol=416, nrow=R)

for (r in 1:R){  	
	u1_2r[r,]=rep(rnorm(4,mean=0,sd=1), each=104)
	u2_2r[r,]=rnorm(416,mean=0,sd=1) 	
	}

# 1.2 Generate p_ijk^r and lambda_ijk^r
p_ijk=matrix(ncol=416, nrow=R); lambda_ijk=matrix(ncol=416, nrow=R); mu_ijk=matrix(ncol=416, nrow=R)
z_ijk=matrix(ncol=416, nrow=R); y_ijk=matrix(ncol=416, nrow=R)

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
    
    if(z_ijk[r,l]==1){
      y_ijk[r,l]=0
    } else{
      y_ijk[r,l]=rpois(1,lambda=Nd[l]*lambda_ijk[r,l]) }		}	
}

beta1r=matrix(ncol=1, nrow=R)
beta2r=matrix(ncol=4, nrow=R)
phi1r=matrix(ncol=1, nrow=R)
phi2r=matrix(ncol=1, nrow=R)

mu_PIr <- mu_PIr_POISSON <- mu_PIr_POISSON1 <- mu_SP <- 
	mu_ESP <- mu_EBLUP <- mu_NB <- matrix(ncol=416, nrow=R)

# Realizamos el estudio de simulacion en si

for (r in 1:R){  
	print(r)
	EPF_2016aux$y_ijk_sim=y_ijk[r,]
	fit_zipoissonr<-glmmTMB(floor(y_ijk_sim)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=poisson)
	
	beta1r[r]=as.matrix((data.frame(fit_zipoissonr$sdr$par.fixed[5])))
	beta2r[r,]=as.matrix((data.frame(fit_zipoissonr$sdr$par.fixed[1:4])))

	phi1r[r]=sqrt(VarCorr(fit_zipoissonr)[[2]][[1]][[1]])
	phi2r[r]=sqrt(VarCorr(fit_zipoissonr)[[1]][[1]][[1]])

	mu_PIr[r,]=fitted(fit_zipoissonr)
	mu_PIr_POISSON[r, ] = exp( x2 %*% beta2r[r,] + phi2r[r] * ranef(fit_zipoissonr)[[1]][[1]][[1]] )*Nd
	aux = 1/(1+exp(beta1r[r,1]))
	mu_PIr_POISSON1[r, ] = rep(aux, 416)*mu_PIr_POISSON[r, ]
	
	result=EBP_ESP(x1,x2, beta1, beta2, phi1, phi2, y_ijk[r,], age4, S) 
	mu_SP[r, ]=result[,2]

	empirical_result=EBP_ESP(x1,x2, beta1r[r], beta2r[r,], phi1r[r], phi2r[r], y_ijk[r,], age4, S) 
	mu_ESP[r, ]=empirical_result[,2]
	
	EBLUP <- eblupFH(y_ijk[r,] ~ 1+edu3+ec2+ec3+ offset(log(floor(Nd))), h1var)
	mu_EBLUP[r,]=EBLUP$eblup
	
	fit_ziNBr <- glmmTMB(floor(y_ijk_sim)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=nbinom2)
	mu_NB[r,]=fitted(fit_ziNBr)
	
	}

proc.time()-t


###### MEDIDAS ABSOLUTAS ###### 
# BIAS PARAMETROS
sesgo_beta1r=colMeans(beta1r-rep(beta1,R)) 
sesgo_beta2r=colMeans(beta2r-t(matrix(beta2, length(beta2), R)))
sesgo_phi1r=colMeans(phi1r-rep(phi1,R)) 
sesgo_phi2r=colMeans(phi2r-rep(phi2,R)) 

print('###### MEDIDAS ABSOLUTAS ###### ')
print('BIAS PARAMETROS')
c(sesgo_beta1r, sesgo_beta2r, sesgo_phi1r, sesgo_phi2r)
print(' ')

# RMSE PARAMETROS
rmse_beta1r=sqrt(colMeans((beta1r-rep(beta1,R))^2))
rmse_beta2r=sqrt(colMeans((beta2r-rep(beta2,R))^2))
rmse_phi1r=sqrt(colMeans((phi1r-rep(phi1,R))^2))
rmse_phi2r=sqrt(colMeans((phi2r-rep(phi2,R))^2))

print('RMSE PARAMETROS')
c(rmse_beta1r, rmse_beta2r, rmse_phi1r, rmse_phi2r)
print(' ')

# BIAS PREDICCIONES
sesgo_muPIr=colMeans(mu_PIr-mu_ijk, na.rm=T)
sesgo_muESP=colMeans(mu_ESP-mu_ijk, na.rm=T)
sesgo_muPIr_POISSON=colMeans(mu_PIr_POISSON-mu_ijk, na.rm=T)
sesgo_muPIr_POISSON1=colMeans(mu_PIr_POISSON1-mu_ijk, na.rm=T)
sesgo_muSP=colMeans(mu_SP-mu_ijk, na.rm=T)
sesgo_muEBLUP=colMeans(mu_EBLUP-mu_ijk, na.rm=T)
sesgo_muNB=colMeans(mu_NB-mu_ijk, na.rm=T)

# RMSE PREDICCIONES
rmse_muPIr=sqrt(colMeans((mu_PIr-mu_ijk)^2))
rmse_muESPr=sqrt(colMeans((mu_ESP-mu_ijk)^2,na.rm=T))
rmse_muPIr_POISSON=sqrt(colMeans((mu_PIr_POISSON-mu_ijk)^2))
rmse_muPIr_POISSON1=sqrt(colMeans((mu_PIr_POISSON1-mu_ijk)^2))
rmse_muSPr=sqrt(colMeans((mu_SP-mu_ijk)^2,na.rm=T))
rmse_muEBLUPr=sqrt(colMeans((mu_EBLUP-mu_ijk)^2,na.rm=T))
rmse_muNBr=sqrt(colMeans((mu_NB-mu_ijk)^2,na.rm=T))

# ABIAS PREDICCIONES
abias_muPIr=mean(abs(sesgo_muPIr))
abias_muESPr=mean(abs(sesgo_muESP), na.rm=T)
abias_muPIr_POISSON=mean(abs(sesgo_muPIr_POISSON))
abias_muPIr_POISSON1=mean(abs(sesgo_muPIr_POISSON1))
abias_muSPr=mean(abs(sesgo_muSP), na.rm=T)
abias_muEBLUPr=mean(abs(sesgo_muEBLUP), na.rm=T)
abias_muNBr=mean(abs(sesgo_muNB), na.rm=T)

print('ABIAS PREDICCIONES')
c(abias_muPIr, abias_muESPr, abias_muPIr_POISSON, abias_muPIr_POISSON1, 
	abias_muSPr, abias_muEBLUPr, abias_muNBr)
print(' ')

# ARMSE PREDICCIONES
armse_muPIr=mean(rmse_muPIr)
armse_muESPr=mean(rmse_muESPr, na.rm=T)
armse_muPIr_POISSON=mean(rmse_muPIr_POISSON)
armse_muPIr_POISSON1=mean(rmse_muPIr_POISSON1)
armse_muSPr=mean(rmse_muSPr, na.rm=T)
armse_muEBLUPr=mean(rmse_muEBLUPr, na.rm=T)
armse_muNBr=mean(rmse_muNBr, na.rm=T)

print('ARMSE PREDICCIONES')
c(armse_muPIr,  armse_muESPr, armse_muPIr_POISSON, armse_muPIr_POISSON1, 
	armse_muSPr, armse_muEBLUPr, armse_muNBr)
print(' ')

###### MEDIDAS RELATIVAS ###### 
# RBIAS PARAMETROS
rsesgo_beta1r=100*sesgo_beta1r/abs(beta1)
rsesgo_beta2r=100*sesgo_beta2r/abs(beta2)
rsesgo_phi1r=100*sesgo_phi1r/abs(phi1)
rsesgo_phi2r=100*sesgo_phi2r/abs(phi2)

print('###### MEDIDAS RELATIVAS ###### ')
print('RBIAS PARAMETROS')
c(rsesgo_beta1r, rsesgo_beta2r, rsesgo_phi1r, rsesgo_phi2r)
print(' ')

# RRMSE PARAMETROS
rrmse_beta1r=100*rmse_beta1r/abs(beta1)
rrmse_beta2r=100*rmse_beta2r/abs(beta2)
rrmse_phi1r=100*rmse_phi1r/abs(phi1)
rrmse_phi2r=100*rmse_phi2r/abs(phi2)

print('RRMSE PARAMETROS')
c(rrmse_beta1r, rrmse_beta2r, rrmse_phi1r, rrmse_phi2r)
print(' ')

# RBIAS PREDICCIONES
rsesgo_muPIr=100*sesgo_muPIr/colMeans(mu_ijk)
rsesgo_muESPr=100*sesgo_muESP/colMeans(y_ijk)
rsesgo_muPIr_POISSON=100*sesgo_muPIr_POISSON/colMeans(mu_ijk)
rsesgo_muPIr_POISSON1=100*sesgo_muPIr_POISSON1/colMeans(mu_ijk)
rsesgo_muSPr=100*sesgo_muSP/colMeans(mu_ijk)
rsesgo_muEBLUPr=100*sesgo_muEBLUP/colMeans(mu_ijk)
rsesgo_muNBr=100*sesgo_muNB/colMeans(mu_ijk)

# RRMSE PREDICCIONES
rrmse_muPIr=100*rmse_muPIr/colMeans(mu_ijk)
rrmse_muESPr=100*rmse_muESPr/colMeans(mu_ijk)
rrmse_muPIr_POISSON=100*rmse_muPIr_POISSON/colMeans(mu_ijk)
rrmse_muPIr_POISSON1=100*rmse_muPIr_POISSON1/colMeans(mu_ijk)
rrmse_muSPr=100*rmse_muSPr/colMeans(mu_ijk)
rrmse_muEBLUPr=100*rmse_muEBLUPr/colMeans(mu_ijk)
rrmse_muNBr=100*rmse_muNBr/colMeans(mu_ijk)

# RABIAS PREDICCIONES
rabias_muPIr=mean(abs(rsesgo_muPIr))
rabias_muESPr=mean(abs(rsesgo_muESPr)[!is.infinite(rsesgo_muESPr)], na.rm=T)
rabias_muPIr_POISSON=mean(abs(rsesgo_muPIr_POISSON))
rabias_muPIr_POISSON1=mean(abs(rsesgo_muPIr_POISSON1))
rabias_muSPr=mean(abs(rsesgo_muSPr)[!is.infinite(rsesgo_muSPr)], na.rm=T)
rabias_muEBLUPr=mean(abs(rsesgo_muEBLUPr))
rabias_muNBr=mean(abs(rsesgo_muNBr))

print('RABIAS PREDICCIONES')
c(rabias_muPIr, rabias_muESPr, rabias_muPIr_POISSON, rabias_muPIr_POISSON1, 
	rabias_muSPr, rabias_muEBLUPr, rabias_muNBr)
print(' ')

# RARMSE PREDICCIONES
rarmse_muPIr=mean(rrmse_muPIr, na.rm=T)
rarmse_muESPr=mean(rrmse_muESPr[!is.infinite(rrmse_muESPr)], na.rm=T)
rarmse_muPIr_POISSON=mean(rrmse_muPIr_POISSON, na.rm=T)
rarmse_muPIr_POISSON1=mean(rrmse_muPIr_POISSON1, na.rm=T)
rarmse_muSPr=mean(rrmse_muSPr[!is.infinite(rrmse_muSPr)], na.rm=T)
rarmse_muEBLUPr=mean(rrmse_muEBLUPr, na.rm=T)
rarmse_muNBr=mean(rrmse_muNBr, na.rm=T)

print('RARMSE PREDICCIONES')
c(rarmse_muPIr, rarmse_muESPr, rarmse_muPIr_POISSON, rarmse_muPIr_POISSON1, 
	rarmse_muSPr, rarmse_muEBLUPr, rarmse_muNBr)
print(' ')

# Save information
df=data.frame(rmse_muPIr, rmse_muESPr)
colnames(df)=c('PlugIN', 'ESP')
write.csv(df,"rmse_predictors.csv", row.names = FALSE)


dff=data.frame(rrmse_muPIr, rrmse_muESPr)
colnames(dff)=c('RRPlugIN', 'RRESP')
write.csv(dff,"rrmse_predictors.csv", row.names = FALSE)

