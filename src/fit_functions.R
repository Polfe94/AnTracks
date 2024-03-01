library(nls.multstart)

## Sum of squares
ssq <- function(obs, est) sum((obs-est)^2)

## Negative exponential trend
nexp <- function(x, a, k){
	fx <- a*exp(-k*x)
	fx
}

## Negative double exponential trend
doub_nexp <- function(x, a, k, a2, k2){
	fx <- a*exp(-k*x) + a2*exp(-k2*x)
	fx
}

## Negative double exponential trend with a baseline
doub_nexpC <- function(x, a, k, a2, k2, C){
	fx <- a*exp(-k*x) + a2*exp(-k2*x) + C
	fx
}

## Negative exponential trend with a baseline
nexpC <- function(x, a, k, C){
	fx <- a*exp(-k*x) + C
	#if(any(fx < ylim[1L] | ylim[2L] < fx)) fx <- rep(NA_real_, length(fx))
	fx
}

## Decaying exponential
fit_nexp <- function(data, a_0, k_0, alims, klims, plot=F){
	
	## Set lower starting parameters
	start_lower <- c(a=a_0*.1, k=k_0*.8)
	
	## Set upper starting parameters
	start_upper <- c(a=a_0*1.2, k=k_0*1.2)
	
	## Fit an exponential model with additive noise
	mod <- nls_multstart(y ~ nexp(x=x, a=a, k=k), data=data, supp_errors = 'Y', start_lower=start_lower,
			     start_upper=start_upper, iter=3e3, control=nls.control(minFactor=1/8192, maxiter=1e4),
			     convergence_count=100L, lower=c(a=alims[1L], k=klims[1L]),
			     upper=c(a=alims[2L], k=klims[2L]),
			     na.action=na.omit)
	
	## In case the model is not NULL
	if(!is.null(mod)){
		
		## Store model parameters in res
		a <- coef(mod)[['a']]
		k <- coef(mod)[['k']]
		
		## Compute the AIC
		aic <- AIC(mod)
		
		## Compute the sum of squares
		err <- ssq(data[['y']], nexp(x=data[['x']], a=a, k=k))
		
		## Useful plot
		if(plot){
			fx <- 10^seq(floor(log10(min(data[['x']]))), ceiling(log10(max(data[['x']]))), length.out=300)
			fy <- a*exp(-k*fx)
			plot(log10(data[['x']]), data[['y']], pch=20, cex=data[['modelweights']], xlim=range(log10(data[['x']])),
			     ylim=range(c(data[['y']], fy)), #ylim=range(c(data[['y']], a*exp(-k*fx) + C)),
			     col=adjustcolor('black', alpha.f=0.25), xlab='t', ylab='rate')
			lines(log10(fx), fy, col='red', lwd=2)
		}
		
		## Return the model data
		data.table(a=a, k=k, aic=aic, ssq=err)
	}
	
}

## Decaying exponential with a constant baseline
fit_nexpC <- function(data, a_0, k_0, C_0, alims, klims, Clims, plot=F){
	
	## Set lower starting parameters
	start_lower <- c(a=a_0*.1, k=k_0*.8, C=min(C_0*c(.8, 1.2)))
	
	## Set upper starting parameters
	start_upper <- c(a=a_0*1.2, k=k_0*1.2, C=max(C_0*c(.8, 1.2)))
	
	## Fit an exponential model with additive noise
	mod <- nls_multstart(y ~ nexpC(x=x, a=a, k=k, C=C), data=data, supp_errors='Y', start_lower=start_lower,
			     start_upper=start_upper, iter=3e3, control=nls.control(minFactor=1/8192, maxiter=1e4),
			     convergence_count=100L, lower=c(a=alims[1L], k=klims[1L], C=Clims[1L]),
			     upper=c(a=alims[2L], k=klims[2L], C=Clims[2L]),
			     na.action=na.omit)
	
	## In case the model is not NULL
	if(!is.null(mod)){
		
		## Store model parameters in res
		a <- coef(mod)[['a']]
		k <- coef(mod)[['k']]
		C <- coef(mod)[['C']]
		
		## Compute the AIC
		aic <- AIC(mod)
		
		## Compute the sum of squares
		err <- ssq(data[['y']], nexpC(x=data[['x']], a=a, k=k, C=C))
		
		## Useful plot
		if(plot){
			fx <- 10^seq(floor(log10(min(data[['x']]))), ceiling(log10(max(data[['x']]))), length.out=300)
			fy <- a*exp(-k*fx) + C
			plot(log10(data[['x']]), data[['y']], pch=20, cex=data[['modelweights']], xlim=range(log10(data[['x']])),
			     ylim=range(c(data[['y']], fy)), #ylim=range(c(data[['y']], a*exp(-k*fx) + C)),
			     col=adjustcolor('black', alpha.f=0.25), xlab='t', ylab='rate')
			lines(log10(fx), fy, col='red', lwd=2)
		}
		
		## Return the model data
		data.table(a=a, k=k, C=C, aic=aic, ssq=err)
	}
	
}


## Decaying double exponential
## Do not set convergence_count=FALSE. Otherwise the code will crash very easily.
fit_doub_nexp <- function(data, a_0, k_0, a2_0, k2_0, alims, klims, a2lims, k2lims, ylim, plot=F){
	
	## Get the starting values
	start_lower <- c(a=a_0*.8, k=k_0*.8, a2=a2_0*.8, k2=k2_0*.8)
	start_upper <- c(a=a_0*1.2, k=k_0*1.2, a2=a2_0*1.2, k2=k2_0*1.2)
	
	## Get the lower and upper limits of the parameter space
	lower <- c(a=alims[1L], k=klims[1L], a2lims=a2lims[1L], k2lims=k2lims[1L])
	upper <- c(a=alims[2L], k=klims[2L], a2lims=a2lims[2L], k2lims=k2lims[2L])
	
	## Fit an exponential model with additive noise
	mod <- nls_multstart(y ~ doub_nexp(x=x, a=a, k=k, a2=a2, k2=k2), data=data, iter=3e3, start_lower=start_lower, 
			     start_upper=start_upper, supp_errors='Y', control=nls.control(minFactor=1/8192, maxiter=1e4), 
			     convergence_count=100L, lower=lower, upper=upper, na.action=na.omit)
	
	## In case the model is not NULL
	if(!is.null(mod)){
		
		## Store model parameters in res
		a <- coef(mod)[['a']]
		k <- coef(mod)[['k']]
		a2 <- coef(mod)[['a2']]
		k2 <- coef(mod)[['k2']]
		
		## Compute the AIC
		aic <- AIC(mod)
		
		## Compute the sum of squares
		err <- ssq(data[['y']], doub_nexp(x=data[['x']], a=a, k=k, a2=a2, k2=k2))
		
		## Useful plot
		if(plot){
			fx <- 10^seq(floor(log10(min(data[['x']]))), ceiling(log10(max(data[['x']]))), length.out=300)
			fy <- a*exp(-k*fx) + a2*exp(-k2*fx)
			plot(log10(data[['x']]), data[['y']], pch=20, cex=data[['modelweights']], xlim=range(log10(data[['x']])),
			     ylim=range(c(data[['y']], fy)), #ylim=range(c(data[['y']], a*exp(-k*fx) + a2*exp(-k2*fx) + C)), 
			     col=adjustcolor('black', alpha.f=0.25), xlab='log10 t', ylab='rate')
			lines(log10(fx), fy, col='red', lwd=2)
		}
		
		## Return the model data
		data.table(a=a, k=k, a2=a2, k2=k2, aic=aic, ssq=err)
		
	}
	
}


## Decaying double exponential with a constant baseline.
## Do not set convergence_count=FALSE. Otherwise the code will crashes very easily.
fit_doub_nexpC <- function(data, a_0, k_0, a2_0, k2_0, C_0, alims, klims, a2lims, k2lims, 
			   Clims, ylim, plot=F){
	
	## Get the starting values
	start_lower <- c(a=a_0*.8, k=k_0*.8, a2=a2_0*.8, k2=k2_0*.8, C=min(C_0*c(.8, 1.2)))
	start_upper <- c(a=a_0*1.2, k=k_0*1.2, a2=a2_0*1.2, k2=k2_0*1.2, C=max(C_0*c(.8, 1.2)))
	
	## Get the lower and upper limits of the parameter space
	lower <- c(a=alims[1L], k=klims[1L], 
		   a2lims=a2lims[1L], k2lims=k2lims[1L], C=Clims[1L])
	upper <- c(a=alims[2L], k=klims[2L], 
		   a2lims=a2lims[2L], k2lims=k2lims[2L], C=Clims[2L])
	
	## Fit an exponential model with additive noise
	mod <- nls_multstart(y ~ doub_nexpC(x=x, a=a, k=k, a2=a2, k2=k2, C=C), data=data, iter=3e3, start_lower=start_lower, 
			     start_upper=start_upper, supp_errors='Y', control=nls.control(minFactor=1/8192, maxiter=1e4), 
			     convergence_count=100L, lower=lower, upper=upper, na.action=na.omit)
	
	## In case the model is not NULL
	if(!is.null(mod)){
		
		## Store model parameters in res
		a <- coef(mod)[['a']]
		k <- coef(mod)[['k']]
		a2 <- coef(mod)[['a2']]
		k2 <- coef(mod)[['k2']]
		C <- coef(mod)[['C']]
		
		## Compute the AIC
		aic <- AIC(mod)
		
		## Compute the sum of squares
		err <- ssq(data[['y']], doub_nexpC(x=data[['x']], a=a, k=k, a2=a2, k2=k2, C=C))
		
		## Useful plot
		if(plot){
			fx <- 10^seq(floor(log10(min(data[['x']]))), ceiling(log10(max(data[['x']]))), length.out=300)
			fy <- a*exp(-k*fx) + a2*exp(-k2*fx) + C
			plot(log10(data[['x']]), data[['y']], pch=20, cex=data[['modelweights']], xlim=range(log10(data[['x']])),
			     ylim=range(c(data[['y']], fy)), #ylim=range(c(data[['y']], a*exp(-k*fx) + a2*exp(-k2*fx) + C)), 
			     col=adjustcolor('black', alpha.f=0.25), xlab='log10 t', ylab='rate')
			lines(log10(fx), fy, col='red', lwd=2)
		}
		
		## Return the model data
		data.table(a=a, k=k, a2=a2, k2=k2, C=C, aic=aic, ssq=err)
		
	}
	
}
