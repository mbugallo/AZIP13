
# Productorio en r=1,..., I y t=1,...,J. Solo depende de k=1,...,K (condition=(age4==k))
Prodk_EBP = function(condition, den1, num1, u1_2s, u2_2s, h1total){
	suma=rep(0, 104)
	for (m in 1:104){
		if(h1total[condition][m]!=0){ suma[m]=sum(log(1: (h1total[condition][m])))
			}	
		}
	return( prod( 1/den1[, condition] * (((h1total[condition]==0)*1) * (den1[, condition]-1 + exp(-Nd[condition]*num1[,condition])) 
			+ (1-(h1total[condition]==0)*1) * (exp(-suma+h1total[condition]* log(num1[,condition])
			- Nd[condition]*num1[,condition] + h1total[condition]*log(Nd[condition])	)))))
}

Prodk_ESP = function(l, den1, num1, u1_2s, u2_2s, h1total){
	if(h1total[l]==0){  suma=0  } else {  suma=sum(log(1:h1total[l])) 
			}	
	return( 1/den1[, l] * (((h1total[l]==0)*1) * (den1[, l]-1 + exp(-Nd[l]*num1[,l])) 
			+ (1-(h1total[l]==0)*1) * (exp(-suma+h1total[l]* log(num1[,l])
			- Nd[l]*num1[,l] + h1total[l]*log(Nd[l])))))	
}


# Empirical best predictor EBP and empirical simplify predictor ESP
EBP_ESP <- function(x1, x2, beta1, beta2, phi1, phi2, h1total, age4, S){
	u1_2s=matrix(0,ncol=416, nrow=2*S)
	u2_2s=matrix(ncol=416, nrow=2*S)

	for (s in 1:S){  	
		u1_2s[s,]=rep(rnorm(4,mean=0,sd=1), each=104)
		u1_2s[S+s,]=-u1_2s[s,]
		u2_2s[s,]=rnorm(416,mean=0,sd=1)
		u2_2s[S+s,]=-u2_2s[s,]
	}
	
	num1=exp( matrix(rep(x2%*%beta2, each=2*S),ncol=416) + phi2*u2_2s)
	den1=1+exp( matrix(rep(x1%*%beta1, each=2*S), ncol=416) + phi1*u1_2s)
	
	# Calculamos todos los exp{x2,ijkβ2 +φ2u2,ijk} (num1) y 1 + exp{x1,ijkβ1 + φ1u1,k} (den1)
	# Calculamos tambien los terminos del productorio en r=1,..., I y t=1,...,J.	    
	
	prod_ebp=cbind(Prodk_EBP(age4==1, den1, num1, u1_2s, u2_2s, h1total), Prodk_EBP(age4==2,den1, num1, u1_2s, u2_2s, h1total),
				Prodk_EBP(age4==3,den1, num1, u1_2s, u2_2s, h1total), Prodk_EBP(age4==4,den1, num1, u1_2s, u2_2s, h1total)	)
	predict_ebp = Nd * colMeans(num1 / den1 * prod_ebp[, age4]) / colMeans( 1 / den1 * prod_ebp[,age4] )
	
	prod_esp=matrix(nrow=2*S, ncol=416)
	for (l in 1:416){
		prod_esp[, l]=Prodk_ESP(l,den1, num1, u1_2s, u2_2s, h1total)
	}
	predict_esp = Nd * colMeans(num1 / den1^2 * prod_esp) / colMeans( 1 / den1 * prod_esp )
	
	return(cbind(predict_ebp, predict_esp))
}

