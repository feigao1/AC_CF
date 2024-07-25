### Incidence Estimation ###
inc_est <- function(Event,PY){
	Est = Event / PY; varlog = 1/Event; SE = sqrt(varlog)* Est
	CI = c(Est-qnorm(0.975)*SE, Est +qnorm(0.975)*SE)
	CI_log = c(Est * exp(-qnorm(0.975)*sqrt(varlog)), Est * exp(qnorm(0.975)*sqrt(varlog)))
	return(list(Est=Est, SE = SE, var_log = varlog, CI = CI, CI_log = CI_log))
}

### Two arm Trial ###
Two_arm_res<-function(n,lambda1,lambda2){
	Event1 = rpois(1,lambda1*n/2); 	Event2 = rpois(1,lambda2*n/2)
	lambda1_est = inc_est(Event1, n/2); lambda2_est = inc_est(Event2, n/2)
	RR = lambda1_est$Est/lambda2_est$Est
	sd_logRR = sqrt(lambda1_est$var_log + lambda2_est$var_log)
	CI_RR = c(RR*exp(qnorm(0.025)* sd_logRR),RR*exp(qnorm(0.975)* sd_logRR))
	return(list(lambda1 = lambda1_est, lambda2 = lambda2_est, RR = RR, var_logRR = sd_logRR^2, CI_RR = CI_RR))
}

### One arm Trial ###
One_arm_res<-function(n,lambda1){
	Event1 = rpois(1,lambda1*n);	lambda1_est = inc_est(Event1, n)
	return(lambda1_est)
}

##### Non-inferiority Trial Design #####
NI_samp = function(delta,deltastar,lambdaA, alpha=0.025,beta = 0.9){
	return(round(2*(1+1/exp(deltastar))/lambdaA/((delta-deltastar)/(qnorm(beta)-qnorm(alpha)))^2))
}

Rej_NI <- function(lambdaE_est, lambdaA_est, delta,alpha=0.025){
	if (lambdaE_est$Est==0) return(list(Rej=T,Test=0))
	if (lambdaA_est$Est==0) return(list(Rej=F,Test=0))
	sigmaEA = sqrt(lambdaE_est$var_log + lambdaA_est$var_log)
	Test = (log(lambdaE_est$Est) - log(lambdaA_est$Est) - delta) / sigmaEA
	return(list(Rej = Test<qnorm(alpha), Test = Test))
}

Analytical_RAE <- function(delta, deltastar, gamma, gammastar, lambdaP, lambdaA, sigmaPA02, alpha, beta){
	sigmaEA = (delta-deltastar)/(qnorm(beta)-qnorm(alpha))
	x = sigmaEA / (1-gamma) / sqrt(sigmaPA02)
	type1 = pnorm(qnorm(alpha)*(x+1)/sqrt(x^2+1))
	x = ((1-gammastar)*(log(lambdaP)-log(lambdaA))-deltastar)/sigmaEA
	power = pnorm(qnorm(beta)-x)
	return(list(type1=type1,power=power))
}