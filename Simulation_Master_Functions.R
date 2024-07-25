source('NI_func.R'); source('CF_func.R')
#### Incidene Specification ####
Inc_spec_083 <- function(){
	R_PA = 2.2; lambdaP = 0.03; lambdaP0 = 0.05; 
	gamma = 0.5; deltastar = log(0.75)
	lambdaA = lambdaP/R_PA; lambdaA0 = lambdaP0/R_PA
	gammastar = 1-deltastar / (log(R_PA))
	n_hist = 3610; 
	return(list(lambdaP = lambdaP, lambdaP0 = lambdaP0, lambdaA = lambdaA, lambdaA0 = lambdaA0, gamma = gamma, gammastar = gammastar,n_hist = n_hist))
}
Inc_spec_CAB<- function(){
	R_PA = 1/0.07; lambdaP = 0.03; lambdaP0 = 0.05; 
	gamma = 0.5; deltastar = log(1)
	lambdaA = lambdaP/R_PA; lambdaA0 = lambdaP0/R_PA
	gammastar = 1-deltastar / (log(R_PA))
	n_hist = 820*2; 
	return(list(lambdaP = lambdaP, lambdaP0 = lambdaP0, lambdaA = lambdaA, lambdaA0 = lambdaA0, gamma = gamma, gammastar = gammastar,n_hist = n_hist))
}


#### Analytical Type-1 Error and Power ####
## NI Design ##
Anal_NI <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=123){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		res = NULL
		for (rep in 1:nrep){
			hist_PC = Two_arm_res(n_hist,lambdaP0,lambdaA0)
			delta = (1-gamma) * log(hist_PC$CI_RR[1])
			res_rep = Analytical_RAE(delta, deltastar, gamma, gammastar, lambdaP, lambdaA, sigmaPA02, alpha, beta)
			res = rbind(res, c(res_rep$type1,res_rep$power))
		}
		Result = rbind(Result,c(alpha,beta,apply(res,2,mean)))
	}
	colnames(Result) = c('Design_type1', 'Design_power','Anal_type1','Anal_power')
	return(Result)
}


Anal_cACCF_FU <- function(Inc_spec, beta_set = c(0.8,0.9)){
	alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		cP0 = 0; cP1=1/(n_hist/2*lambdaP); ## based on trial with size n_hist/2
		lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
		cE = 2/lambdaE; cA = 2/lambdaA
		Ntrial = ceiling(NI_CF_samp_9595(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,20000), alpha, beta)) # two arms
		res = Analytical_9595(Ntrial, gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, alpha)
		Result = rbind(Result,c(alpha,beta,res))
	}
	colnames(Result) = c('Design_type1', 'Design_power','Anal_type1','Anal_power')
	return(Result)
}

Anal_cACCF_RA <- function(Inc_spec, beta_set = c(0.8,0.9)){
	alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		for (tau in c(1,2)){
			para = para_set(lambdaP,p=0.15,u=1,v=1,MDRI=142/365.25,FRR=.01,Time=2,sd_MDRI=.07*142/365.25,sd_FRR=.25*.01,tau=tau)
			cP = rec_gamma(para);
			cP0 = cP$gamma00 *(1-para$p)*para$tau; cP1 = cP$gamma01
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			cE = 2/lambdaE; cA = 2/lambdaA
			Ntrial = ceiling(NI_CF_samp_9595(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,20000), alpha, beta)/2/tau)*2*tau # two arms
			res = Analytical_9595(Ntrial, gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, alpha)
			Result = rbind(Result,c(alpha,beta,tau, res))
		}
	}
	colnames(Result) = c('Design_type1', 'Design_power','Design_FU','Anal_type1','Anal_power')
	return(Result)
}

#### Empirical Type-1 Error and Power ####
## NI Design ##
Empir_NI <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		Nres = Rej_NI0 = Rej_NIa = Rej_RAE0 = Rej_RAEa = NULL
		for (rep in 1:nrep){
			hist_PC = Two_arm_res(n_hist,lambdaP0,lambdaA0)
			delta = (1-gamma) * log(hist_PC$CI_RR[1])
			if (is.na(delta)) {rep = rep-1; print('NA margin!'); break}
			N = ceiling(NI_samp(delta, deltastar,lambdaA, alpha,beta))
			Nres = c(Nres, N)
			
			### H0 ###
			lambdaE = lambdaA * exp(delta)
			Data_H0 = Two_arm_res(N,lambdaE,lambdaA)
			res = Rej_NI(Data_H0$lambda1, Data_H0$lambda2, delta)
			Rej_NI0 = c(Rej_NI0,res$Rej)
			
			### H1 ###
			lambdaE = lambdaA * exp(deltastar)
			Data_Ha = Two_arm_res(N,lambdaE,lambdaA)
			res = Rej_NI(Data_Ha$lambda1, Data_Ha$lambda2, delta)
			Rej_NIa = c(Rej_NIa,res$Rej)
			
			### K0 ###
			lambdaE = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
			Data_K0 = Two_arm_res(N,lambdaE,lambdaA)
			res = Rej_NI(Data_K0$lambda1, Data_K0$lambda2, delta)
			Rej_RAE0 = c(Rej_RAE0,res$Rej)
	
			### K1 ###
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			Data_Ka = Two_arm_res(N,lambdaE,lambdaA)
			res = Rej_NI(Data_Ka$lambda1, Data_Ka$lambda2, delta)
			Rej_RAEa = c(Rej_RAEa, res$Rej)
		}
		Result = rbind(Result,c(alpha, beta, mean(Nres), mean(Rej_NI0), mean(Rej_NIa), mean(Rej_RAE0), mean(Rej_RAEa)))
	}
	colnames(Result) = c('Design_type1','Design_power','Total_PY','Type1_NI','Power_NI','Type1_RAE','Power_RAE')
	return(Result)
}

## AC-CF Design with Follow-up CF ##
Empir_ACCF_FU <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		cP0 = 0; cP1=1/(n_hist/2*lambdaP); ## based on trial with size n_hist/2
		lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
		cE = 2/lambdaE; cA = 2/lambdaA
		Ntrial = ceiling(NI_CF_samp(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha, beta)) # two arms
		Rej_null = Rej_alt = NULL
		for (rep in 1:nrep){
			lambdaP_est = One_arm_res(n_hist/2,lambdaP)
		
			### Type-1 error on Absolute Efficacy ###
			lambdaE = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
			Data_K0 = Two_arm_res(Ntrial,lambdaE,lambdaA)
			res2 = Rej_RAE(lambdaP_est, Data_K0$lambda1, Data_K0$lambda2, gamma,alpha)
			Rej_null = rbind(Rej_null,c(res2$Rej))
	
			### Power on Absolute Efficacy ###
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			Data_Ka = Two_arm_res(Ntrial,lambdaE,lambdaA)
			res2 = Rej_RAE(lambdaP_est, Data_Ka$lambda1, Data_Ka$lambda2, gamma,alpha)
			Rej_alt = rbind(Rej_alt,c(res2$Rej))
		}
		Result = rbind(Result,c(alpha, beta,Ntrial,apply(Rej_null,2,mean),apply(Rej_alt,2,mean)))
	}
	colnames(Result) = c('Design_type1','Design_power','Total_PY','Type1_RAE','Power_RAE')
	return(Result)
}



## AC-CF Design with Recency-testing CF ##
Empir_ACCF_RA <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		for (tau in c(1,2)){
			para = para_set(lambdaP,p=0.15,u=1,v=1,MDRI=142/365.25,FRR=.01,Time=2,sd_MDRI=.07*142/365.25,sd_FRR=.25*.01,tau=tau)
			cP = rec_gamma(para);
			cP0 = cP$gamma00 *(1-para$p)*para$tau; cP1 = cP$gamma01
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			cE = 2/lambdaE; cA = 2/lambdaA
			Nneg = ceiling(NI_CF_samp(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha, beta)/2/tau)*2 # two arms
			Rej_null = Rej_alt = N_res = NULL
			for (rep in 1:nrep){
				N = round(Nneg/(1-para$p)); Npos = N - Nneg
				data_rec = generate_rec(Npos, para)
				lambdaP_est = rec_est(N,Npos,data_rec$NR,data_rec$para)
				N_res = rbind(N_res,c(N, Npos,data_rec$NR, Nneg,Nneg*tau))
			
				### Type-1 error on Absolute Efficacy ###
				lambdaE = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
				Data_K0 = Two_arm_res(Nneg*tau,lambdaE,lambdaA)
				res2 = Rej_RAE(lambdaP_est, Data_K0$lambda1, Data_K0$lambda2, gamma,alpha)
				Rej_null = rbind(Rej_null,c(res2$Rej))
	
				### Power on Absolute Efficacy ###
				lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
				Data_Ka = Two_arm_res(Nneg*tau,lambdaE,lambdaA)
				res2 = Rej_RAE(lambdaP_est, Data_Ka$lambda1, Data_Ka$lambda2, gamma,alpha)
				Rej_alt = rbind(Rej_alt,c(res2$Rej))
			}
			Result = rbind(Result,c(alpha, beta,tau,apply(N_res,2,mean),apply(Rej_null,2,mean),apply(Rej_alt,2,mean)))
		}
	}
	colnames(Result) = c('Design_type1','Design_power','tau','Nscreen','Npos','Nrec','Nneg','Total_PY','Type1_RAE','Power_RAE')
	return(Result)
}

## Conservative AC-CF Design with Follow-up CF ##
Empir_cACCF_FU <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		cP0 = 0; cP1=1/(n_hist/2*lambdaP); ## based on trial with size n_hist/2
		lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
		cE = 2/lambdaE; cA = 2/lambdaA
		Ntrial = ceiling(NI_CF_samp_9595(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,20000), alpha, beta)) # two arms
		Rej_null = Rej_alt = NULL
		for (rep in 1:nrep){
			lambdaP_est = One_arm_res(n_hist/2,lambdaP)
		
			### Type-1 error on Absolute Efficacy ###
			lambdaE = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
			Data_K0 = Two_arm_res(Ntrial,lambdaE,lambdaA)
			res2 = Rej_RAE95(lambdaP_est, Data_K0$lambda1, Data_K0$lambda2, gamma,alpha)
			Rej_null = rbind(Rej_null,c(res2$Rej))
	
			### Power on Absolute Efficacy ###
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			Data_Ka = Two_arm_res(Ntrial,lambdaE,lambdaA)
			res2 = Rej_RAE95(lambdaP_est, Data_Ka$lambda1, Data_Ka$lambda2, gamma,alpha)
			Rej_alt = rbind(Rej_alt,c(res2$Rej))
		}
		Result = rbind(Result,c(alpha, beta,Ntrial,apply(Rej_null,2,mean),apply(Rej_alt,2,mean)))
	}
	colnames(Result) = c('Design_type1','Design_power','Total_PY','Type1_RAE','Power_RAE')
	return(Result)
}



## Conservative AC-CF Design with Recency-testing CF ##
Empir_cACCF_RA <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	Result = NULL
	for (beta in beta_set){
		for (tau in c(1,2)){
			para = para_set(lambdaP,p=0.15,u=1,v=1,MDRI=142/365.25,FRR=.01,Time=2,sd_MDRI=.07*142/365.25,sd_FRR=.25*.01,tau=tau)
			cP = rec_gamma(para);
			cP0 = cP$gamma00 *(1-para$p)*para$tau; cP1 = cP$gamma01
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			cE = 2/lambdaE; cA = 2/lambdaA
			Nneg = ceiling(NI_CF_samp_9595(gamma, gammastar, lambdaP, lambdaA, cP0, cP1, cE, cA, Nrange=c(1,20000), alpha, beta)/2/tau)*2 # two arms
			Rej_null = Rej_alt = N_res = NULL
			for (rep in 1:nrep){
				N = round(Nneg/(1-para$p)); Npos = N - Nneg
				data_rec = generate_rec(Npos, para)
				lambdaP_est = rec_est(N,Npos,data_rec$NR,data_rec$para)
				N_res = rbind(N_res,c(N, Npos,data_rec$NR, Nneg,Nneg*tau))
			
				### Type-1 error on Absolute Efficacy ###
				lambdaE = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
				Data_K0 = Two_arm_res(Nneg*tau,lambdaE,lambdaA)
				res2 = Rej_RAE95(lambdaP_est, Data_K0$lambda1, Data_K0$lambda2, gamma,alpha)
				Rej_null = rbind(Rej_null,c(res2$Rej))
	
				### Power on Absolute Efficacy ###
				lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
				Data_Ka = Two_arm_res(Nneg*tau,lambdaE,lambdaA)
				res2 = Rej_RAE95(lambdaP_est, Data_Ka$lambda1, Data_Ka$lambda2, gamma,alpha)
				Rej_alt = rbind(Rej_alt,c(res2$Rej))
			}
			Result = rbind(Result,c(alpha, beta,tau,apply(N_res,2,mean),apply(Rej_null,2,mean),apply(Rej_alt,2,mean)))
		}
	}
	colnames(Result) = c('Design_type1','Design_power','tau','Nscreen','Npos','Nrec','Nneg','Total_PY','Type1_RAE','Power_RAE')
	return(Result)
}

## OneArm-CF Design with Follow-up CF ##
Empir_OneArmCF_FU <- function(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP = Inc_spec$lambdaP; lambdaA = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP/lambdaA; deltastar = (1-gammastar)* (log(R_PA))
	zeta = gamma * log(R_PA); zetastar = gammastar * log(R_PA)
	Result = NULL
	for (beta in beta_set){
		cP0 = 0; cP1=1/(n_hist/2*lambdaP); ## based on trial with size n_hist/2
		lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
		cE = 1/lambdaE
	Ntrial = ceiling(Onearm_CF_samp(zeta, zetastar, cP0, cP1, cE, alpha, beta)) # one arm
		Rej_null = Rej_alt = NULL
		for (rep in 1:nrep){
			lambdaP_est = One_arm_res(n_hist/2,lambdaP)
		
			### Type-1 error on Absolute Efficacy ###
			lambdaE = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
			lambdaE_est = One_arm_res(Ntrial,lambdaE)
			res2 = Rej_PE(lambdaP_est, lambdaE_est, zeta,alpha)
			Rej_null = rbind(Rej_null,c(res2$Rej))
	
			### Power on Absolute Efficacy ###
			lambdaE = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
			lambdaE_est = One_arm_res(Ntrial,lambdaE)
			res2 = Rej_PE(lambdaP_est, lambdaE_est, zeta,alpha)
			Rej_alt = rbind(Rej_alt,c(res2$Rej))
		}
		Result = rbind(Result,c(alpha, beta, Ntrial,apply(Rej_null,2,mean),apply(Rej_alt,2,mean)))
	}
	colnames(Result) = c('Design_type1','Design_power','Total_PY','Type1_AE','Power_AE')
	return(Result)
}


Violation <- function(Inc_spec, lambdaP_grid, lambdaA_grid, Method = c('NI','AC-CF','cAC-CF','OneArm-CF'), beta = 0.8, nrep=10000, seed=1234){
	set.seed(seed); alpha = 0.025
	lambdaP0 = Inc_spec$lambdaP0; lambdaA0 = Inc_spec$lambdaA0; n_hist = Inc_spec$n_hist
	sigmaPA02 = (1/(lambdaP0*n_hist/2) + 1/(lambdaA0*n_hist/2))
	lambdaP_design = Inc_spec$lambdaP; lambdaA_design = Inc_spec$lambdaA
	gamma = Inc_spec$gamma; gammastar = Inc_spec$gammastar
	R_PA = lambdaP_design/lambdaA_design; deltastar = (1-gammastar)* (log(R_PA))
	zeta = gamma * log(R_PA); zetastar = gammastar * log(R_PA)
	Result = NULL
	if (Method %in% c('AC-CF','cAC-CF')){
		cP0 = 0; cP1=1/(n_hist/2* lambdaP_design); ## based on trial with size n_hist/2
		lambdaE_design = exp(log(lambdaP_design) - gammastar * (log(lambdaP_design)-log(lambdaA_design)))
		cE = 2/lambdaE_design; cA = 2/lambdaA_design
		if (Method =='AC-CF'){
			PY = ceiling(NI_CF_samp(gamma, gammastar, lambdaP_design, lambdaA_design, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha, beta)) # two arms
		} else {
			PY = ceiling(NI_CF_samp_9595(gamma, gammastar, lambdaP_design, lambdaA_design, cP0, cP1, cE, cA, Nrange=c(1,10000), alpha, beta)) # two arms
		} 
	} else {
		if (Method =='OneArm-CF') {
			cP0 = 0; cP1=1/(n_hist/2* lambdaP_design); ## based on trial with size n_hist/2
			lambdaE_design = exp(log(lambdaP_design) - gammastar * (log(lambdaP_design)-log(lambdaA_design)))
			cE = 1/lambdaE_design
			PY = ceiling(Onearm_CF_samp(zeta, zetastar, cP0, cP1, cE, alpha, beta)) # one arm
		}
	}
	for (lambdaP in lambdaP_grid){
		for (lambdaA in lambdaA_grid){
			Rej_null = Rej_alt = PY_NI = NULL
			for (rep in 1:nrep){
				if (Method=='NI'){
					hist_PC = Two_arm_res(n_hist,lambdaP0,lambdaA0)
					delta = (1-gamma) * log(hist_PC$CI_RR[1])
					if (is.na(delta)) {rep = rep-1; print('NA margin!'); break}
					PY = ceiling(NI_samp(delta, deltastar,lambdaA_design, alpha,beta))
					PY_NI = c(PY_NI,PY)
				} else {
					lambdaP_est = One_arm_res(n_hist/2,lambdaP_design)
				}
				### Type-1 error on Absolute Efficacy ###
				lambdaE_null = exp(log(lambdaP) - gamma * (log(lambdaP)-log(lambdaA)))
				lambdaE_alt = exp(log(lambdaP) - gammastar * (log(lambdaP)-log(lambdaA)))
				if (Method=='OneArm-CF'){
					lambdaE_null = exp(log(lambdaP) - zeta)
					lambdaE_est = One_arm_res(PY,lambdaE_null)
					res1 = Rej_PE(lambdaP_est, lambdaE_est, zeta, alpha)
					lambdaE_alt = exp(log(lambdaP) - zetastar)
					lambdaE_est = One_arm_res(PY,lambdaE_alt)
					res2 = Rej_PE(lambdaP_est, lambdaE_est, zeta, alpha)
				} else {
					Data_K0 = Two_arm_res(PY,lambdaE_null,lambdaA)
					Data_Ka = Two_arm_res(PY,lambdaE_alt,lambdaA)
					if (Method=='AC-CF'){
						res1 = Rej_RAE(lambdaP_est, Data_K0$lambda1, Data_K0$lambda2, gamma,alpha)
						res2 = Rej_RAE(lambdaP_est, Data_Ka$lambda1, Data_Ka$lambda2, gamma,alpha)
					} else {
						if (Method=='cAC-CF'){
							res1 = Rej_RAE95(lambdaP_est, Data_K0$lambda1, Data_K0$lambda2, gamma,alpha)
							res2 = Rej_RAE95(lambdaP_est, Data_Ka$lambda1, Data_Ka$lambda2, gamma,alpha)
						} else{ # NI
							res1 = Rej_NI(Data_K0$lambda1, Data_K0$lambda2, delta)
							res2 = Rej_NI(Data_Ka$lambda1, Data_Ka$lambda2, delta)
						}
					}
				}
				Rej_null = rbind(Rej_null,c(res1$Rej))
				Rej_alt = rbind(Rej_alt,c(res2$Rej))
			}
			if (Method=='NI') PY = mean(PY_NI)
			Result = rbind(Result,c(alpha, beta, PY,lambdaP,lambdaA,apply(Rej_null,2,mean),apply(Rej_alt,2,mean)))
		}
	}
	colnames(Result) = c('alpha','beta','PY','lambdaP','lambdaA','Type_1','Power')
	return(Result)
}
