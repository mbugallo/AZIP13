

###########  Crear Mapas ####################

rm(list=ls(all.names=TRUE))              

library(maptools)
library(RColorBrewer)

##########################################


GroupClassification <- function(data,datacompare,intervals)
{
   n = length(data)
   group = matrix(0,nrow=n) 
   ninterv = length(intervals)
   for (i in 1:n)
   {
      for (j in 1:ninterv)
         if (datacompare[i]<intervals[j])
         {
            group[i]= intervals[j]
            break
         }
   }
   result = split(data,group)
   return (result)
}
 
PrintSpainMap <- function(pathmap,datos,colors,titlemap,textlegend,eliminarprov)
{

   m <- matrix(c(1,1,1,2),2,2)
   layout(m, widths=c(1.5, 1), heights=c(1.5, 1), respect=F)

   xName <- readShapePoly(pathmap, IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))     

   xName$datos <- NA
   for (i in 1:length(colors))
      xName$datos[datos[[i]]] <- colors[i]
 
   xSC <- xName[xName$ESP_PROV_I < 35 | xName$ESP_PROV_I >38 | xName$ESP_PROV_I==36 | xName$ESP_PROV_I ==37,]
   plot(xSC,  xlab="",  col=xSC$datos, axes=F)

  title(titlemap, line=-1.3, cex.main=1.5)
  legend( "topright", textlegend, col = 1, pt.bg = colors, pch=21, bty = "n", cex=1.2, pt.cex = 2.5)     

   xC <- xName[xName$ESP_PROV_I == 35 | xName$ESP_PROV_I == 38,]
   plot(xC,  xlab="",  col=xC$datos)
   box()
}


##############################################################

## Define file names
file_est   <- "provincias_age4.txt" 
pathmap    <- "spainmap/esp_prov.shp"

c_dom     <- 1     # domain column number
c_estML   <- 2   # column number

datos <- read.table(file=file_est, header=TRUE, sep=';', dec='.')

nprov <- nrow(datos)
datos <- datos[1:(nprov-2),]  # elimino Ceuta y Melilla

dom     <- datos[,c_dom]
estML   <- datos[,c_estML]

# Intervalos  
intervals_prop <- c(0.35, 0.5, 0.65, 0.75, 1)  

# Colores 
colorsprop <- c(brewer.pal(7,"GnBu")[c(4,5,6,7)] )
 
# Leyendas 
legend_prop <- expression("35 - 50", "50 - 65", "65 - 75", 'over 75') 


result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"Young women: Age group 1",legend_prop,eliminarprov)

##############################################################
pathmap    <- "spainmap/esp_prov.shp"
PlugIN=read.csv("rmse_PlugIN.csv", header=T, sep=',')
rmse_muPIr=PlugIN$RMSE_PlugIN
rrmse_muPIr=PlugIN$RRMSE_PlugIN/100
dom=seq(1,50, by=1)

estML1=rrmse_muPIr[1:50] # hombres jovenes
estML2=rrmse_muPIr[53:102] # mujeres jovenes
#estML1=rrmse_muPIr[105:154] #hombres mid
#estML2=rrmse_muPIr[157:206] # mujeres mid

#estML1=rrmse_muPIr[209:258] #hombres adult
#estML2=rrmse_muPIr[261:310] # mujeres adult
#estML1=rrmse_muPIr[313:362] #hombres old
#estML2=rrmse_muPIr[365:414] # mujeres old

estML=estML1

# Intervalos  
intervals_prop <- c(0.10, 0.15, 0.20, 0.30, 0.4, 1)

# Colores 
colorsprop <- c(brewer.pal(9,"YlOrRd")[c(2, 4, 5, 6, 7, 9)] )
# Leyendas 
legend_prop <- expression('under 10%',  "10 - 15 %", "15 - 20 %",  "20 - 30 %", "30 - 40 %", "over 40 %") 
result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"RRMSE Young men: Age group 1",legend_prop,eliminarprov)



