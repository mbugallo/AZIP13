###############################################################################
###############################################################################
###
###       Created on Monday October 17 2022
###       Directory  /Users/mariabugalloporto/Desktop/UMH_Elche/MariaBugallo/1Paper_Singleperson/Write_Results                                               
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(glmmTMB)
library(KernSmooth)

########################################
######## PART 1 
########################################

## Datos EPF2016 + 4 EPA2016 trimestrales
EPF_2016aux=read.csv("EPF_EPA_2016.csv", sep=',', dec='.', header=T) 

EPF_2016aux$prov=factor(EPF_2016aux$prov)
EPF_2016aux$sex=factor(EPF_2016aux$sex)
EPF_2016aux$age4=factor(EPF_2016aux$age4)
attach(EPF_2016aux)


########################################
######## PART 2 & 3
########################################

fit_zipoisson0<-glmmTMB(floor(h1total)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	ziformula=~1+(1|age4), data=EPF_2016aux, family=poisson)
summary(fit_zipoisson0)
round(confint(fit_zipoisson0),3) 

# Modelo log-lineal de Poisson no inflado en el 0
fit_poisson0<-glmmTMB(floor(h1total)~1+edu3+ec2+ec3+offset(log(floor(Nd)))+(1|prov:sex:age4), 
	 data=EPF_2016aux, family=poisson)
summary(fit_poisson0)

########################################
######## MODELO FINAL: fit_zipoisson0
########################################

sigma_f=sqrt(VarCorr(fit_zipoisson0)[[1]][[1]][[1]])
ud_f<-ranef(fit_zipoisson0)[[1]][[1]][[1]]/sigma_f  # modes of the random effects conditional model

sigma_z=sqrt(VarCorr(fit_zipoisson0)[[2]][[1]][[1]])
ud_z<-ranef(fit_zipoisson0)[[2]][[1]][[1]]/sigma_z  # modes of the random effects z-i model


########################################
######## PART 4
########################################

predic_values=fitted(fit_zipoisson0)
predic_values_PO=fitted(fit_poisson0)

########### Residuos brutos
res.bt= h1total- predic_values
res.bt_PO= h1total-predic_values_PO

########### Residuos estandarizados
res.st=(res.bt-mean(res.bt))/sd(res.bt)
res.st_PO=(res.bt_PO-mean(res.bt_PO))/sd(res.bt_PO)

########################################

# 1. Standarized residuals vs domain index
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(1:416, rep(0, 416), type='l', xlab='Domain index', 
	ylab='Standarized residuals', main=' ', cex.main=1.7, ylim=c(-6,6), cex.lab=1.7, cex.axis=1.7)
points((1:416)[EPF_2016aux$sex==1], res.st[EPF_2016aux$sex==1], pch=19, col='darkblue')	
points((1:416)[EPF_2016aux$sex==6], res.st[EPF_2016aux$sex==6], pch=3, col='darkblue')
new <- data.frame(x = seq(1, 416))
lines(new$x, rep(0, 416), col='red', lwd=3, lty=2)	
legend('topright', col='darkblue', c('Men', 'Women'), cex=1.5, pt.cex = 1.4, pch=c(19,3))


# 2. Standarized residuals vs predicted proportions
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot((predic_values/Nd)[EPF_2016aux$sex==1], res.st[EPF_2016aux$sex==1], 
	xlab='Predicted proportions', ylab='Standarized residuals', col='darkblue',  main=' ',
	pch=19, cex.main=1.7, ylim=c(-6,6), xlim=c(0,1), cex.lab=1.7, cex.axis=1.7)
points((predic_values/Nd)[EPF_2016aux$sex==6], res.st[EPF_2016aux$sex==6], pch=3, col='darkblue')
new <- data.frame(x = seq(0,1,0.01))
lines(new$x, rep(0, 101), col='red', lwd=3, lty=2)	


# 3. Standarized residuals vs log predicted proportions
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(log(predic_values/Nd)[EPF_2016aux$sex==1 & res.st > -6], res.st[EPF_2016aux$sex==1 & res.st > -6], 
	xlab='Log predicted proportions', ylab='Standarized residuals', col='darkblue',  main=' ',
	pch=19, cex.main=1.7, ylim=c(-6,6), xlim=c(-3.5, 0), cex.lab=1.7, cex.axis=1.7)
points(log(predic_values/Nd)[EPF_2016aux$sex==6 & res.st > -6], res.st[EPF_2016aux$sex==6 & res.st > -6], pch=3, col='darkblue')
new <- data.frame(x = seq(-3.5,0,0.1))
lines(new$x, rep(0, 36), col='red', lwd=3, lty=2)

# 4. Boxplot standarized residuals
par(mar=c(4, 4.3, 4, 4), xpd=T)
boxplot(res.st~EPF_2016aux$prov, main='Province', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
boxplot(res.st~EPF_2016aux$sex, main='Sex', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)
boxplot(res.st~EPF_2016aux$age4, main='Age group', ylab='', xlab='', cex.main=1.8, cex.lab=1.8, cex.axis=1.8)

########################################

########### Null counts
# 1
j=sort(h1n[h1n==0], index.return=T)$ix

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(h1mean[j], type='o', pch=19, col="black", lty=1, ylim=c(0,1),lwd=1.2, xaxt='n', yaxt='n',
	xlab='Count', ylab=' Proportion', axes = T, cex.lab=1.8, cex=0.7, main='', cex.main=1.8)
axis(1, cex.axis=1.8, at=c(0, 28), labels=c('0', '1'))
axis(2, cex.axis=1.8) 
points(predic_values[j]/Nd[j], type='o', pch=4, col='red', lty=2, cex=0.7) 
points(predic_values_PO[j]/Nd[j], type='p', pch=20, col='blue', lty=2, cex=0.7) 
legend('topright', c('Hajek', 'IN', 'IN0'), col = c('black', 'red', 'blue'), 
lty = c(1, 2), pch=c(20, 4, 20), bty = "n", cex=1.8, pt.cex = 1.8, lwd=2)

# 2
j=sort(h1n[h1n<5], index.return=T)$ix

par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(h1mean[j], type='o', pch=19, col="black", lty=1, ylim=c(0,1),lwd=1.2,
	xlab='Count', ylab=' Proportion', xaxt='n', yaxt='n', cex.lab=1.8, cex=0.7, main='', cex.main=1.8)
axis(1, cex.axis=1.8, at=c(0, 28, 69, 109, 137, 173), labels=c('0', '1', '2', '3', '4', '5'))
axis(2, cex.axis=1.8) 
points(predic_values[j]/Nd[j], type='o', pch=4, col='red', lty=2, cex=0.7) 
points(predic_values_PO[j]/Nd[j], type='p', pch=20, col='blue', lty=2, cex=0.7) 
legend('topright', c('Hajek', 'IN', 'IN0'), col = c('black', 'red', 'blue'), 
lty = c(1, 2), pch=c(20, 4, 20), bty = "n", cex=1.8, pt.cex = 1.8, lwd=2)


# 3
par(mar=c(4, 4.3, 4, 4), xpd=T)
pred.df=data.frame('prop'=c(res.st[h1n<5], res.st_PO[h1n<5]), 
	'occ'=rep(h1n[h1n<5], 2), 'id'= c(rep(1, 173), rep(2, 173)))
boxplot(pred.df$prop ~ pred.df$id:pred.df$occ, col=c('red', 'blue'), xaxt='n', yaxt='n', xlab='Count',
	ylab='Standarized residuals', cex.lab=1.8, cex.axis=1.8, pch=19, ylim=c(-6,2))
axis(1, at=c(0, 1.5, 3.5, 5.5, 7.5, 9.5, 11), labels=c('', '0', '1', '2', '3', '4', ''),
	cex.lab=1.85, cex.axis=1.8)	; axis(2, cex.axis=1.8) 
legend(-0.8,3.3, legend=c("aZIP13","GLMM Poisson"), lty=c(1,1), col=c('red', 'blue'),
	horiz=T, cex=1.7, lwd=5, bty = "n")	


########################################
######## PART 5
########################################	

# Proporciones de hogares unipersonales
# 1. Domain index
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(h1mean, type='l', col='gray40', main=' ', ylab='Proportion', xaxt = "n", xlab='Age group', lwd=1.5, 
	cex.main=1.8, yaxt = "n", cex.lab=1.7, ylim=c(0,1)); 
axis(1, at=c(0, 52, 104, 156, 208, 260, 312, 364, 416), labels=c('', 1, '', 2, '', 3, '', 4, ''), cex.axis=1.7); 
axis(2, cex.axis=1.7); lines(predic_values/Nd, type='l', col='red')	
legend('topleft', col=c('gray40','red'), c('Hajek', 'IN'), bty = "n", cex=2, pt.cex = 2, lwd=2)


# 2. Sample size
par(mar=c(4, 4.3, 4, 4), xpd=T)
i=sort(nd, index.return=T)$ix
plot(h1mean[i], type='l', col='gray40', main=' ', ylab='Proportion', xaxt = "n", xlab='Sample size', lwd=1.5, 
	cex.main=1.8, yaxt = "n", cex.lab=1.7, ylim=c(0,1)); 
points(predic_values[i]/Nd[i], type='l', col='red') 	
axis(1, cex.axis=1.7, at=round(seq(1,416,length=5),0), labels=floor(sort(nd, index.return=T)$x)[round(seq(1,416,length=5),0)])
axis(2, cex.axis=1.7) 


# 3. PI vs Hajek
par(mar=c(4, 4.3, 4, 4), xpd=T)
plot(h1mean[EPF_2016aux$sex==1], predic_values[EPF_2016aux$sex==1]/Nd[EPF_2016aux$sex==1], pch=19, xlab='Hajek', 
	ylab='PI', main=' ', cex.main=1.7, xlim=c(0,1), ylim=c(0,1), cex.lab=1.7, cex.axis=1.7, col='darkblue')
points(h1mean[EPF_2016aux$sex==6], (predic_values/Nd)[EPF_2016aux$sex==6], pch=3, col='darkblue')
new <- data.frame(x = seq(0,1, by=0.01))
lines(new$x, new$x, col='red', lwd=3, lty=2)	
fit0 <- locpoly(h1mean, predic_values/Nd, bandwidth = 1.3, degree=1) 	
lines(fit0, lwd=3)
legend('topleft', col='darkblue', c('Men', 'Women'), cex=1.5, pt.cex = 1.4, pch=c(19,3))



########################################
######## PART 6: Se crean los csv para hacer los mapas en otro script
########################################

	# SEXOS
EPF_2016auxp=cbind(EPF_2016aux[,c(1,2,3)], predic_values, predic_values/Nd)

datos_m=EPF_2016auxp[EPF_2016aux[,2]==6,] # mujeres
datos_m1=datos_m[datos_m[,3]==1,] # mujeres del primer grupo de edad

datos_h=EPF_2016auxp[EPF_2016aux[,2]==1,] # hombres
datos_h1=datos_h[datos_h[,3]==1,] # hombres del primer grupo de edad

used=cbind(1:52,datos_m1['predic_values/Nd'], datos_h1['predic_values/Nd']) 
names(used)=c('province', 'women_age1', 'men_age1')

write.table(used,"provincias_age1.txt", sep=";",col.names = TRUE, row.names = FALSE)

	# EDADES
a1=aggregate(predic_values, by=list(prov, age4), sum)
a2=aggregate(Nd, by=list(prov, age4), sum)

used2=cbind(1:52, a1[a1[,2]==1,3]/a2[a2[,2]==1,3], a1[a1[,2]==2,3]/a2[a2[,2]==2,3], 
	a1[a1[,2]==3,3]/a2[a2[,2]==3,3], a1[a1[,2]==4,3]/a2[a2[,2]==4,3])
names(used2)=c('province', 'age1', 'age2', 'age3', 'age4')

write.table(used2,"provincias_ages.txt", sep=";",col.names = TRUE, row.names = FALSE)


########################################
######## PART 7 (NO SE USA)
########################################
# Binomial negativa

# variance = phi mu
fit_zinbinom1 <- update(fit_zipoisson0,family=nbinom1)
summary(fit_zinbinom1) 

fit_zinbinom2 <- update(fit_zipoisson0,family=nbinom2)
summary(fit_zinbinom2)

# Hurdle models treat zero-count and non-zero outcomes as two completely separate categories
fit_hnbinom1 <- update(fit_zipoisson0, ziformula=~1.,data=EPF_2016aux,family=truncated_nbinom1)
summary(fit_hnbinom1)


########################################
######## PART 8
########################################

##### Otras predicciones

x1=as.matrix(data.frame(rep(1,416)))
beta1=as.matrix((data.frame(fit_zipoisson0$sdr$par.fixed[5])))
phi1=sqrt(VarCorr(fit_zipoisson0)[[2]][[1]][[1]])

x2=as.matrix(data.frame(rep(1,416), EPF_2016aux[c(13, 19, 20)]))
beta2=as.matrix((data.frame(fit_zipoisson0$sdr$par.fixed[1:4])))
phi2=sqrt(VarCorr(fit_zipoisson0)[[1]][[1]][[1]])


####  The empirical best predictor EBP and the empirical simplify predictor ESP
source('EBP_ESP.R')
t <- proc.time()
resultado=EBP_ESP(x1,x2, beta1, beta2, phi1, phi2, h1total, age4, S=2000) 
proc.time()-t

predict_ebp=resultado[,1]
predict_esp=resultado[,2]


# Comparemos su comportamiento de los predictores (EBP, ESP, Hajek y Plug-In)

i=sort(nd, index.return=T)$ix
# 1.
plot(h1mean[i]*Nd[i], type='o', pch=20, axes=F, xlab='', ylab='Count variable', main='', 
	cex.main=1.6, cex=0.65, lty=1, col='black', cex.lab=1.3)
axis(1, cex.axis=1.55, at=round(seq(1,416,length=10),0), labels=floor(sort(nd, index.return=T)$x)[round(seq(1,416,length=10),0)])
axis(2, cex.axis=1.55) 
points(predic_values[i], type='o', pch=3, col="blue", lty=2, cex=0.5) 
points(predict_esp[i], type='o', pch=4, col="orange", lty=2, cex=0.5) 

legend( "top", c('Hajek', 'Plug-In', 'ESP'), col = c("black", 'blue', "orange"), 
lty = c(1, 2, 2), pch=c(20, 25, 3), bty = "n", cex=1.2, pt.cex = 1.3)

# 2.
plot(h1mean[i], type='o', pch=20, axes=F, xlab='', ylab='Proportion', main='', cex.main=1.6, cex=0.65, lty=1, 
	col='black', cex.lab=1.3, ylim=c(0,1))
axis(1, cex.axis=1.55, at=round(seq(1,416,length=10),0), labels=floor(sort(nd, index.return=T)$x)[round(seq(1,416,length=10),0)])
axis(2, cex.axis=1.55) 
points(predic_values[i]/Nd[i], type='o', pch=3, col="blue", lty=2, cex=0.5) 
points(predict_esp[i]/Nd[i]*0.9, type='o', pch=4, col="orange", lty=2, cex=0.5) 

legend( "topright", c('Hajek', 'Plug-In', 'ESP'), col = c("black", 'blue', "orange"), 
lty = c(1, 2, 2), pch=c(20, 25, 3), bty = "n", cex=1.2, pt.cex = 1.3)





