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
phi1=sqrt(VarCorr(fit_zipoisson0)[[2]][[1]][[1]])

x2=as.matrix(data.frame(rep(1,416), EPF_2016aux[c(13, 19, 20)]))
beta2=as.matrix((data.frame(fit_zipoisson0$sdr$par.fixed[1:4])))
phi2=sqrt(VarCorr(fit_zipoisson0)[[1]][[1]][[1]])


########################################
######## PART 3: SIMULATION 2 
########################################

# Studies the behavior of the parametric bootstrap estimator of the MSE of a predictor mu_yijk^hat of mu_yijk

rmse_predictors=read.csv("rmse_predictors.csv",header=T, sep=',')
rmse_muPIr=rmse_predictors$PlugIN
# rmse_muPIr=rmse_predictors$ESP

t <- proc.time()

MSE_ijk=rmse_muPIr^2 # empirical MSE (Simulation 1)
R=500 # R = 500
B=50 # B = 50, 100, 200, 400, 500

u1_2r=matrix(ncol=416, nrow=R)
u2_2r=matrix(ncol=416, nrow=R)

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

# Realizamos el estudio de simulacion en si (para cada r=1,...,R, estos datos los reescribimos)
pijk_rb=matrix(nrow=B, ncol=416)
lambdaijk_rb=matrix(nrow=B, ncol=416) 
muijk_rb=matrix(nrow=B, ncol=416) 
zijk_rb=matrix(nrow=B, ncol=416) 
yijk_rb=matrix(nrow=B, ncol=416)
muijk_rb_est=matrix(nrow=B, ncol=416) 

mse=matrix(nrow=R, ncol=416)

S=2000
	
for (r in 1:R){  
  print(r)
	EPF_2016aux$y_ijk_sim=y_ijk[r,]
	fit_zipoissonr<-glmmTMB(floor(y_ijk_sim)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=poisson)

	# Reescribense en cada r (en 1:R):
	u1_2b=matrix(ncol=416, nrow=B); u2_2b=matrix(ncol=416, nrow=B)
	for(b in 1:B){
		u1_2b[b,]=rep(rnorm(4,mean=0,sd=1), each=104)
		u2_2b[b,]=rnorm(416,mean=0,sd=1)}
	
	for(b in 1:B){ 
		pijk_rb[b,]=exp(rep(beta1, 416)+phi1*u1_2b[b,]) / (1 + exp(rep(beta1, 416)+phi1*u1_2b[b,]))
		lambdaijk_rb[b,]=exp(x2%*%beta2+phi2*u2_2b[b, ])
		muijk_rb[b,]=Nd*(1-pijk_rb[b,])*lambdaijk_rb[b,]
		for (l in 1:416){
			zijk_rb[b,l]=rbinom(1,size=1, prob=pijk_rb[b,l]) 
				if(zijk_rb[b,l]==1){
					yijk_rb[b,l]=0
				} else{
					yijk_rb[b,l]=rpois(1,lambda=Nd[l]*lambdaijk_rb[b,l])
				}	}
		EPF_2016aux$yijk_simb=yijk_rb[b,]
		fit_zipoissonb<-glmmTMB(floor(yijk_simb)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=poisson)	
	
    # beta1rb=as.matrix((data.frame(fit_zipoissonb$sdr$par.fixed[5])))
    # phi1rb=sqrt(VarCorr(fit_zipoissonb)[[2]][[1]][[1]])

    # beta2rb=as.matrix((data.frame(fit_zipoissonb$sdr$par.fixed[1:4])))
    # phi2rb=sqrt(VarCorr(fit_zipoissonb)[[1]][[1]][[1]])

		muijk_rb_est[b,]=fitted(fit_zipoissonb)  }
    #EBP_ESP(x1, x2, beta1rb, beta2rb, phi1rb, phi2rb, yijk_rb[b,], age4, S)[,2] 	}
				
		mse[r,]=colMeans((muijk_rb_est-muijk_rb)^2, na.rm=T)	# calculo del mse*_ijk^r
}

proc.time()-t

MSE_ijk_replicate=t(replicate(expr=MSE_ijk,n=R))

###### MEDIDAS ABSOLUTAS ###### 

print('Bijk: colMeans(mse*-MSE_ijk)')
Bijk=colMeans(mse-MSE_ijk_replicate, na.rm=T); # Bijk
print(' ')

print('REijk: sqrt(colMeans((mse*-MSE)^2))')
REijk=sqrt(colMeans((mse-MSE_ijk_replicate)^2,na.rm=T)); # REijk
print(' ')

print('AB: sum(abs(Bijk))')
AB=sum(abs(Bijk)); AB
print(' ')

print('RE: sum(abs(RE))')
RE=sum(REijk); RE
print(' ')

###### MEDIDAS RELATIVAS ###### 

RBijk=100*Bijk/MSE_ijk
RREijk=100*REijk/MSE_ijk

print('AARB: sum(abs(RBijk)')
AARB=mean(abs(RBijk)); AARB
print(' ')

print('ARRE: sum(abs(RREijk))')
ARRE=mean(RREijk); ARRE
print('')

##### INFORMACION #####
df=data.frame(RREijk, RBijk)
colnames(df)=c('RREijk_50', 'RBijk_50')
write.csv(df,"simu2_B50.csv", row.names = FALSE)




