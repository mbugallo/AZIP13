###############################################################################
###############################################################################
###
###       Created on Monday March 14 2022
###       DIRECTORY  /Users/mariabugalloporto/Desktop/UMH_Elche/MariaBugallo/Datos_EPA2016                                                  
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(car)


########################################
######## PART 1 
########################################

# Information extrated from the files:
# https://www.ine.es/dyngs/INEbase/es/operacion.htm?c=Estadistica_C&cid=1254736176918&menu=resultados&idp=1254735976595#!tabs-1254736030639
# https://www.ine.es/CDINEbase/consultar.do?mes=&operacion=EPA.+Resultados+trimestrales&id_oper=Ir

# Read the four csv files
EPA_2016T1=read.csv("EPA_2016T1.csv", sep=';', header=T) # 161465 rows and 93 columns
EPA_2016T2=read.csv("EPA_2016T2.csv", sep=';', header=T) # 160833 rows and 93 columns
EPA_2016T3=read.csv("EPA_2016T3.csv", sep=';', header=T) # 159957 rows and 93 columns
EPA_2016T4=read.csv("EPA_2016T4.csv", sep=';', header=T) # 158915 rows and 93 columns

EPA_2016=rbind(EPA_2016T1,EPA_2016T2,EPA_2016T3,EPA_2016T4) # 641170 rows and 93 columns

# Delete individuals with missing data (younger than 16 -> they don't work)
dat_EPA2016 = EPA_2016[which(EPA_2016$EDAD5>15),] # 541862 rows and 7 columns
attach(dat_EPA2016)

age4=cut(dat_EPA2016$EDAD5,breaks=c(0,45,55,64,Inf), include.lowest=TRUE,labels=c("1","2","3","4"))

########################################
######## PART 2 
########################################

# PROV: Provincia donde se ubica la vivienda
# NPERS: Número de personas en el hogar
# FACTOREL: Factor de elevación
# EDAD5: Edad, grupos quinquenales de años cumplidos. De 00 a 65
# SEXO1: Sexo. 1 Hombre, 6 Mujer
# NFORMA: Nivel de estudios 
# AOI: Situacion laboral
# NAC1: Nacionalidad


########################################
######## PART 3 
########################################

## Variance of the Hajek direct estimator
dir2 <- function(data, w, domain, Nd) {
  if(is.vector(data)){
    last <- length(domain) + 1
    Nd.hat <- aggregate(w, by=domain, sum)[,last]
    nd <- aggregate(rep(1, length(data)), by=domain, sum)[,last]
    Sum <- aggregate(w*data, by=domain, sum)
    mean <- Sum[,last]/Nd.hat
    dom <- as.numeric(Reduce("paste0", domain))
    domain.unique <- sort(unique(dom))
    difference <- list()
    for(d in 1:length(mean)){
      condition <- dom==domain.unique[d]
      difference[[d]] <- w[condition]*(w[condition]-1)*(data[condition]-mean[d])^2
    }
    var.mean <- unlist(lapply(difference, sum))/Nd.hat^2
    if(missing(Nd)){
      return(data.frame(Sum[,-last], mean, var.mean, Nd.hat, nd))
    }
    else{
      tot <- mean*Nd
      var.tot <- var.mean*Nd^2
      return(data.frame(Sum[,-last], tot, var.tot, mean, var.mean, Nd.hat, Nd, nd))
    }
  }
  else{
    warning("Only a numeric or integer vector must be called as data",
            call. = FALSE)
  }
}

# Estimacion directa del tamanho poblacional con los datos de la EPA
Nd=dir2(data=rep(1,length(FACTOREL)), FACTOREL, domain=list(PROV,SEXO1,age4))$Nd.hat

# Variables auxiliares
## ci1: NACIONASP
cit1 <- as.numeric(NAC1==1) # Solamente nacionalidad española 
cit1m <- dir2(data=cit1, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporcion espanholes

## edu1, edu2, edu3 and edu4: NFORMA
edu <- recode(NFORMA, " c('AN','P1','P2') = 1; c('S1') = 2; c('SG', 'SP') = 3; c('SU') = 4 ", 
	as.factor=TRUE, levels = c(1,2,3,4))

edu1<-as.numeric(edu==1) # 'primary'
edu2<-as.numeric(edu==2) # 'secondary'
edu3<-as.numeric(edu==3) # 'high secondary'
edu4<-as.numeric(edu==4) # 'superior'

edu1m<-dir2(data=edu1, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción de primaria 
edu2m<-dir2(data=edu2, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción de secundaria básica 
edu3m<-dir2(data=edu3, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción de secundaria avanzada 
edu4m<-dir2(data=edu4, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción de educacion superior

## lab1, lab2 and lab3: AOI
lab <- recode(AOI, " c('3','4') = 'employed'; c('5','6') = 'unemployed'; c('7','8','9') = 'inactive' ", 
as.factor=TRUE, levels = c('employed','unemployed', 'inactive'))

lab1<-as.numeric(lab=='employed') # 'employed'
lab2<-as.numeric(lab=='unemployed') # 'unemployed'
lab3<-as.numeric(lab=='inactive') # 'inactive'

lab1m<-dir2(data=lab1, FACTOREL, domain=list(PROV,SEXO1,age4))$mean  # proporción empregados 
lab2m<-dir2(data=lab2, FACTOREL, domain=list(PROV,SEXO1,age4))$mean  # proporción desempregados
lab3m<-dir2(data=lab3, FACTOREL, domain=list(PROV,SEXO1,age4))$mean  # proporción inactivos

## ec1, ec2, ec3 and ec4: ECIV1
ec1<-as.numeric(ECIV1==1) # 'soltero'
ec2<-as.numeric(ECIV1==2) # 'casado'
ec3<-as.numeric(ECIV1==3) # 'viudo'
ec4<-as.numeric(ECIV1==4) # 'separado o divorciado'

ec1m<-dir2(data=ec1, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción solteros
ec2m<-dir2(data=ec2, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción casados
ec3m<-dir2(data=ec3, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción viudos
ec4m<-dir2(data=ec4, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción separados ou divorciados

## tm1 and tm2: MUN1
tm1<-as.numeric(MUN1==1) # leva máis dun ano na mesma casa
tm2<-as.numeric(MUN1==6) # cambiou de casa fai menos dun ano

tm1m<-dir2(data=tm1, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción de non mudanzas
tm2m<-dir2(data=tm2, FACTOREL, domain=list(PROV,SEXO1,age4))$mean # proporción de mudanzas

########################################
######## PART 5 
########################################

dat_EPF2016 <- read.table("datosEPF2016.txt", header=TRUE, sep = "\t", dec = ",") 

# Delete individuals with missing data (younger than 16 -> they don't work)
dat_EPF2016 = dat_EPF2016[which(dat_EPF2016$EDADSP>15),] 

age4<-cut(dat_EPF2016$EDADSP,breaks=c(0,45,55,64,Inf),include.lowest=TRUE,labels=c("1","2","3","4"))

UNIPER<-as.numeric(dat_EPF2016$TAMANO==1)

direct_estimacion <- dir2(data=UNIPER, dat_EPF2016$FACTOR , domain=list(dat_EPF2016$PROV,dat_EPF2016$SEXOSP, age4))
names(direct_estimacion) <- c('prov', 'sex', 'age4', 'h1mean', 'h1var', 'hatNd', 'nd')

direct_estimacion$h1n <- data.frame(aggregate(UNIPER, by=list(dat_EPF2016$PROV,dat_EPF2016$SEXOSP, age4), sum))$x
 
direct_estimacion$Nd <- Nd/4 # consideramos 4 EPA
sum(direct_estimacion$Nd)

direct_estimacion$h1total <- direct_estimacion$h1mean * direct_estimacion$Nd
sum(direct_estimacion$h1total)

direct_estimacion <- direct_estimacion[, -c(6)]

direct_estimacion$ci1 <- cit1m
direct_estimacion$edu1 <- edu1m
direct_estimacion$edu2 <- edu2m
direct_estimacion$edu3 <- edu3m
direct_estimacion$edu4 <- edu4m
direct_estimacion$lab1 <- lab1m
direct_estimacion$lab2 <- lab2m
direct_estimacion$lab3 <- lab3m
direct_estimacion$ec1 <- ec1m
direct_estimacion$ec2 <- ec2m
direct_estimacion$ec3 <- ec3m
direct_estimacion$ec4 <- ec4m
direct_estimacion$tm1 <- tm1m
direct_estimacion$tm2 <- tm2m

head(direct_estimacion, 20)

write.csv(direct_estimacion, file="EPF_EPA_2016.csv", row.names = F)



