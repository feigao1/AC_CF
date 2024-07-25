### Testing with Counterfactual Placebo ###
Rej_const <- function(lambdaP_est,lambdaA_est,Hist_PA_RR,alpha=0.025){
	if (lambdaA_est$Est==0) {
		Test = (log(Hist_PA_RR$RR)) / sqrt(Hist_PA_RR$var_logRR)
	} else{
		Test = (log(lambdaP_est$Est) - log(lambdaA_est$Est) - log(Hist_PA_RR$RR)) / sqrt(lambdaA_est$var_log + lambdaP_est$var_log+ Hist_PA_RR$var_logRR)
	}
	return(list(Rej = Test<qnorm(alpha),Test = Test))
}

Rej_PA <- function(lambdaP_est, lambdaA_est, alpha=0.025){
	if (lambdaP_est$Est==0) return(list(Rej=F,Test=0))
	if (lambdaA_est$Est==0) return(list(Rej=T,Test=0))
	Test = (log(lambdaP_est$Est) - log(lambdaA_est$Est)) / sqrt(lambdaP_est$var_log+lambdaA_est$var_log)
	return(list(Rej = Test>-qnorm(alpha),Test = Test))
}

Rej_RAE <- function(lambdaP_est, lambdaE_est, lambdaA_est, gamma,alpha=0.025){
	if (Rej_PA(lambdaP_est, lambdaA_est, alpha=0.025)$Rej == F) return(list(Rej=F,Test = 0, RejAP = F))
	if (lambdaA_est$Est==0) return(list(Rej=F,Test=0,RejAP = T))
	if (lambdaE_est$Est==0) return(list(Rej=T,Test=0,RejAP = T))
	Test = (log(lambdaP_est$Est)*(1-gamma) - log(lambdaE_est$Est) + gamma*log(lambdaA_est$Est)) / sqrt((1-gamma)^2*lambdaP_est$var_log+lambdaE_est$var_log+gamma^2*lambdaA_est$var_log)
	return(list(Rej = Test>-qnorm(alpha),Test = Test, RejAP =F))
}

Rej_PE <- function(lambdaP_est, lambdaE_est, zeta, alpha=0.025){
	if (lambdaP_est$Est==0) return(list(Rej=F,Test=0))
	if (lambdaE_est$Est==0) return(list(Rej=T,Test=0))
	Test = (log(lambdaP_est$Est) - log(lambdaE_est$Est) - zeta) / sqrt(lambdaP_est$var_log+lambdaE_est$var_log)
	return(list(Rej = Test>-qnorm(alpha), Test = Test))
}

Rej_RAE95 <- function(lambdaP_est, lambdaE_est, lambdaA_est, gamma,alpha=0.025){
	lambdaP95 = lambdaP_est
	lambdaP95$Est = exp(log(lambdaP_est$Est) + qnorm(0.025) * sqrt(lambdaP_est$var_log))
	lambdaP95$var_log = 0
	return(Rej_RAE(lambdaP95, lambdaE_est, lambdaA_est, gamma,alpha))
}

Rej_PAdouble <- function(lambdaP_est, lambdaA_est, alpha=0.025){
	if (lambdaA_est$Est==0) return(list(Rej=F,Test=0))
	Test = (log(lambdaP_est$Est) - log(lambdaA_est$Est)) / sqrt(2*lambdaP_est$var_log+lambdaA_est$var_log)
	return(list(Rej = Test>-qnorm(alpha),Test = Test))
}

Rej_RAE_double <- function(lambdaP_est, lambdaE_est, lambdaA_est, gamma,alpha=0.025){
	if (Rej_PA(lambdaP_est, lambdaA_est, alpha=0.025)$Rej ==F) return(list(Rej=F,Test = 0))
	if (lambdaE_est$Est==0) return(list(Rej=T,Test=0))
	Test = (log(lambdaP_est$Est)*(1-gamma) - log(lambdaE_est$Est) + gamma*log(lambdaA_est$Est)) / sqrt((1-gamma)^2*2*lambdaP_est$var_log+lambdaE_est$var_log+gamma^2*lambdaA_est$var_log)
	return(list(Rej = Test>-qnorm(alpha),Test = Test))
}

### Recency Testing Functions ###
para_set <- function(lambda0,p,u,v,MDRI,FRR,Time,sd_MDRI,sd_FRR,tau){
  return(list(lambda0=lambda0, p=p, u=u, v=v, OmegaT=MDRI, sigmaOmegaT=sd_MDRI, betaT=FRR, sigmabetaT=sd_FRR, Time=Time, tau=tau))
}
generate_rec <- function(Npos, para){
	p = para$p; lambda0 = para$lambda0; 
	OmegaT = para$OmegaT; sigmaOmegaT = para$sigmaOmegaT; betaT = para$betaT; sigmabetaT = para$sigmabetaT; Time = para$Time
	PR = betaT + lambda0*(1-p)/p*(OmegaT-betaT*Time)
	NR = rbinom(1,Npos,PR)
	para$betaT = rnorm(1,betaT,sigmabetaT); para$OmegaT = rnorm(1,OmegaT,sigmaOmegaT)
	return(list(Npos = Npos, NR = NR, para=para))
}
rec_est <- function(N,Npos,NR,para){
	pest = Npos/N; uest = 1; Nneg = N - Npos
	betaT = para$betaT; OmegaT = para$OmegaT; Time = para$Time
	lambdaest = (NR/uest - betaT*Npos) / (Nneg*(OmegaT-betaT*Time))
	para_est = list(lambda0=lambdaest, p=pest, u=uest, OmegaT=OmegaT, sigmaOmegaT = para$sigmaOmegaT, betaT=betaT, sigmabetaT = para$sigmabetaT, Time=Time)
	varlog = rec_gamma(para_est,N)$varlog
	Est=lambdaest; SE = sqrt(varlog)*lambdaest;
	CI = c(lambdaest-qnorm(0.975)*SE, lambdaest+qnorm(0.975)*SE)
	CI_log = c(lambdaest * exp(-qnorm(0.975)*sqrt(varlog)), lambdaest * exp(qnorm(0.975)*sqrt(varlog)))
	return(list(Est=Est,SE = SE, var_log = varlog, CI = CI, CI_log = CI_log))
}
rec_gamma<-function(para, N=NULL){
	lambda0 = para$lambda0; p = para$p; u = para$u; OmegaT = para$OmegaT; betaT = para$betaT; Time = para$Time
	sigma2Omega = para$sigmaOmegaT^2;sigma2beta = para$sigmabetaT^2
	PR = betaT + lambda0 *(1-p)/p*(OmegaT-betaT*Time)
	gamma00 = (PR*(1-PR)/(PR-betaT)^2/u + 1/(1-p) + (1-p*u)*sigma2beta/(PR-betaT)^2/u)/p
	gamma01 = sigma2Omega/(OmegaT-betaT*Time)^2 + sigma2beta*((OmegaT-PR*Time)/(PR-betaT)/(OmegaT-betaT*Time))^2
	return(list(gamma00 = gamma00, gamma01 = gamma01, varlog = gamma00/N+gamma01))
}

NI_CF_samp <- function(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha = 0.025, beta = 0.8){
	samplesize <-function(N){
		x = (gammastar-gamma)*(log(lambdaP)-log(lambdaA))/sqrt(((1-gamma)^2*cP0 + cE + gamma^2*cA)/N+(1-gamma)^2*cP1)
		y = (log(lambdaP)-log(lambdaA))/sqrt((cP0+cA)/N+cP1)
		pnorm(qnorm(alpha) + x) + pnorm(qnorm(alpha)+y) - 1 - beta
	}
	
	N = try(uniroot(samplesize,interval=Nrange)$root, TRUE)
	return(N)
}


Onearm_CF_samp <- function(zeta, zetastar, cP0, cP1, cE, alpha = 0.025, beta = 0.8){
	N = (cP0 + cE)/(((zetastar - zeta)/(qnorm(beta)-qnorm(alpha)))^2-cP1)
	return(N)
}


NI_CF_samp_9595 <- function(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha = 0.025, beta = 0.8){
	samplesize <-function(N){
	  Vgamma = (cE + gamma^2*cA)/N
	  sp2 = cP0/N+cP1
	  x1 = (sqrt(Vgamma) +(1-gamma)*sqrt(sp2))/sqrt(Vgamma+(1-gamma)^2*sp2)
	  x2 = (gammastar-gamma)*(log(lambdaP)-log(lambdaA))/sqrt(Vgamma+(1-gamma)^2*sp2)
		y1 = (sqrt(sp2)+sqrt(cA/N))/sqrt(sp2+cA/N)
		y2 = (log(lambdaP)-log(lambdaA))/sqrt(sp2+cA/N)
		pnorm(qnorm(alpha)*x1 + x2) + pnorm(qnorm(alpha)*y1+y2) - 1 - beta
	}
	N = try(uniroot(samplesize,interval=Nrange)$root, TRUE)
	return(N)
}

power_func_9595 <-function(N, gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, alpha = 0.025){
    Vgamma = (cE + gamma^2*cA)/N
    sp2 = cP0/N+cP1
    x1 = (sqrt(Vgamma) +(1-gamma)*sqrt(sp2))/sqrt(Vgamma+(1-gamma)^2*sp2)
    x2 = (gammastar-gamma)*(log(lambdaP)-log(lambdaA))/sqrt(Vgamma+(1-gamma)^2*sp2)
    y1 = (sqrt(sp2)+sqrt(cA/N))/sqrt(sp2+cA/N)
    y2 = (log(lambdaP)-log(lambdaA))/sqrt(sp2+cA/N)
    power_AS = pnorm(qnorm(alpha)*y1+y2)
    power_RAE = pnorm(qnorm(alpha)*x1 + x2)
  return(list(power_AS = power_AS, power_RAE = power_RAE))
}

Analytical_9595 <- function(N, gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, alpha = 0.025){
  rAP = sqrt(cA/N / (cP0/N+cP1))
  rEA = sqrt(cE/cA)
  alpha1 = pnorm(qnorm(alpha)*(sqrt(rEA^2+gamma^2*rAP^2)+(1-gamma))/(sqrt(rEA^2+gamma^2*rAP^2+(1-gamma))))
  alpha2 = pnorm(qnorm(alpha)*(rAP+1)/(sqrt(rAP^2+1)))
  power = power_func_9595(N, gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, alpha = 0.025)
  return(list(type1=min(alpha1,alpha2),power = power$power_AS + power$power_RAE-1))
}

NI_CF_samp_2sd <- function(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha = 0.025, beta = 0.8){
	cP0 = cP0*2; cP1 = cP1*2;
	samplesize <-function(N){
		x = (gammastar-gamma)*(log(lambdaP)-log(lambdaA))/sqrt(((1-gamma)^2*cP0 + cE + gamma^2*cA)/N+(1-gamma)^2*cP1)
		y = (log(lambdaP)-log(lambdaA))/sqrt((cP0+cA)/N+cP1)
		pnorm(qnorm(alpha) + x) + pnorm(qnorm(alpha)+y) - 1 - beta
	}
	
	N = try(uniroot(samplesize,interval=Nrange)$root, TRUE)
	return(N)
}

NI_CF_samp_fixP <- function(gamma,gammastar,lambdaA,lambdaE, lambdaP, sigmaP2, alpha = 0.025, beta = beta){ # PY enrolled
	cA = 1/(lambdaA/2); cE = 1/(lambdaE/2);
	effectsize = (gamma-gammastar)*(log(lambdaA)-log(lambdaP))
	numer = cE + gamma^2*cA
	vZ = 1 # varZ(para, lambdaE,lambdaA,gamma,gammastar)
	denom = (effectsize/(qnorm(beta,0,sqrt(vZ))-qnorm(alpha)))^2 - (1-gamma)^2*sigmaP2
	return(numer/denom)
}
varZ <-function(para, lambdaE,lambdaA,gamma,gammastar){
	p=para$p; u= para$u; OmegaT = para$OmegaT;betaT=para$betaT; v=para$v;tau =para$tau; Time = para$Time
	lambdaP = para$lambda0
	PR = betaT + lambdaP *(1-p)/p*(OmegaT-betaT*Time)
	varAB = var_AB(p,u,PR,OmegaT,betaT,v,lambdaE,lambdaA,tau,gamma)
	tA = (gamma-gammastar)*(log(lambdaA)-log(lambdaP))
	tB = (1-gamma)^2*(PR*(1-PR)/(p*u*(PR-betaT)^2) + 1/p) + 1/(1-p) + 1/((1-p)*v*lambdaE*tau/2) + gamma^2/((1-p)*v*lambdaA*tau/2)
	c = rep(0,2); c[1] = 1/sqrt(tB); c[2] = -tA/(2*tB^(3/2))
	return(t(c)%*%varAB%*%c)
}


var_AB <- function(p,u,PR,OmegaT,betaT,v,lambdaE,lambdaA,tau,gamma){
	a = matrix(0,8,2)
	a[1,1] = -1/(p*u*(PR-betaT))*(1-gamma)
	a[2,1] = 1/(p*u)*(1-gamma)
	a[3,1] = -1/(p*(1-p))*(1-gamma)
	a[4,1] = 1/OmegaT*(1-gamma)
	a[5,1] = 1/((1-p)*v*lambdaE*tau)*(-1)
	a[6,1] = -1/((1-p)*v)*(1-gamma)
	a[8,1] = 1/((1-p)*v*lambdaA*tau)*gamma
	a[1,2] = -2*PR*(1-PR)/(p*u)^2/(PR-betaT)^3*(1-gamma)^2
	a[2,2] = PR^2/(p*u*(PR-betaT))^2*(1-gamma)^2
	a[3,2] = -(p^2+(1-p)^2)/(p*(1-p))^2*(1-gamma)^2
	a[4,2] = 0
	a[5,2] = -1/((1-p)*v*lambdaE*tau/2)^2*(-1)^2
	a[6,2] = 0
	a[7,2] = (1-2*PR)/(p*u*(PR-betaT))^2*(1-gamma)^2
	a[8,2] = -1/((1-p)*v*lambdaA*tau/2)^2*(-1)^2
	cov_mat = varW(p,u,PR,OmegaT, betaT, v,lambdaE,lambdaA,tau)
	return(t(a)%*%cov_mat%*%a)
}
#### NR - N+test*wbetaT, N+test, N+, wOmega -wbeta*T, NeventE, N-enroll, NR, NeventA
varW <-function(p,u,PR,OmegaT,betaT,v,lambdaE,lambdaA,tau){
	variance = matrix(0,8,8)
	variance[1,1] = p*u*PR*(1-PR) + p*u*(1-p*u)*(PR-betaT)^2
	variance[2,2] = p*u*(1-p*u)
	variance[3,3] = p*(1-p)
	variance[5,5] = (1-p)*v*lambdaE*tau*(1+lambdaE*tau*(1-v)+lambdaE*tau*v*p)/4
	variance[6,6] = (1-p)*v*(1-v+p*v)
	variance[7,7] = p*u*PR*(1-PR*u*p)
	variance[8,8] = (1-p)*v*lambdaA*tau*(1+lambdaA*tau*(1-v)+lambdaA*tau*v*p)/4
	variance[1,2] = variance[2,1] = p*u*(1-p*u)*(PR-betaT)
	variance[1,3] = variance[3,1] = p*(1-p)*u*(PR-betaT)
	variance[1,5] = variance[5,1] = -p*(1-p)*(PR-betaT)*u*v*lambdaE*tau/2
	variance[1,6] = variance[6,1] = -p*(1-p)*(PR-betaT)*u*v
	variance[1,7] = variance[7,1] = p*u*PR*(1-PR) + p*u*(1-p*u)*(PR-betaT)*PR
	variance[1,8] = variance[8,1] = -p*(1-p)*(PR-betaT)*u*v*lambdaA*tau/2
	variance[2,3] = variance[3,2] = p*(1-p)*u
	variance[2,5] = variance[5,2] = -p*(1-p)*u*v*lambdaE*tau/2
	variance[2,6] = variance[6,2] = -p*(1-p)*u*v
	variance[2,7] = variance[7,2] = p*u*(1-p*u)*PR
	variance[2,8] = variance[8,2] = -p*(1-p)*u*v*lambdaA*tau/2
	variance[3,5] = variance[5,3] = -p*(1-p)*v*lambdaE*tau/2
	variance[3,6] = variance[6,3] = -p*(1-p)*v
	variance[3,7] = variance[7,3] = p*(1-p)*u*PR
	variance[3,8] = variance[8,3] = -p*(1-p)*v*lambdaA*tau/2
	variance[5,6] = variance[6,5] = (1-p)*v*(1-v+p*v)*lambdaE*tau/2
	variance[5,7] = variance[7,5] = -p*(1-p)*PR*u*v*lambdaE*tau/2
	variance[8,6] = variance[6,8] = (1-p)*v*(1-v+p*v)*lambdaA*tau/2
	variance[8,7] = variance[7,8] = -p*(1-p)*PR*u*v*lambdaA*tau/2
	variance[6,7] = variance[7,6] = -p*(1-p)*PR*u*v
	variance[5,8] = variance[8,5] = (1-p)*v*(1-v+p*v)*lambdaA*lambdaE*tau^2/4
	return(variance)
}

### Find alpha for correlated N(0,1)s ###
library("mvtnorm")
findcv = function(sigma, alpha=.025, npoint = 5){ # find bivariate critical values under null (mean=c(0,0)), assume sigma has been standardized so variances are 1
	f = function(q2,q1,mean,sigma,alpha){pmvnorm(upper=c(q1,q2),mean=mean,sigma=sigma)-alpha}
	# focus on lower tail of two-sided test 
	b1 = qnorm(alpha)
	b2 = qmvnorm(alpha,tail="lower.tail",mean=c(0,0),sigma=sigma)$quantile #equicoordinate quantile
	q1 = seq(b1,b2,length.out = npoint+1)[-1]
	q2 = rep(NA,length(q1))
	for (i in 1:length(q1)){
		q2[i] = try(uniroot(f,interval=c(b1,5),mean=c(0,0),sigma=sigma,alpha=alpha,q1=q1[i])$root, TRUE)
	}
	x = c(q1,rev(q2)); y = c(q2,rev(q1))
	list(xy=matrix(c(x,y),ncol=2),alpha=matrix(c(pnorm(x),pnorm(y)),ncol=2))
}