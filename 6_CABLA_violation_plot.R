library('ggplot2')
library('ggpubr')
library('dplyr')

Violatin_plot_prepare <- function(data){
	data = data[data$lambdaP > data$lambdaA,] ## Viable Designs
	data$Type_1_cut = cut(data$Type_1,breaks = c(0,0.01, 0.025,0.05,0.2,0.5,1),right = F)
	data$Power_cut = cut(data$Power,breaks = c(0,0.3,0.5,0.8,0.9,1),right = T)
	return(data)
}

Result_NI = Violatin_plot_prepare(data.frame(readRDS(file = "Violation_Result_CAB_NI.rds")))
Result_ACCF = Violatin_plot_prepare(data.frame(readRDS(file = "Violation_Result_CAB_ACCF.rds")))
Result_cACCF = Violatin_plot_prepare(data.frame(readRDS(file = "Violation_Result_CAB_cACCF.rds")))

#### Conservativeness Cutoff Approximation ####
Result_NI_cons = Result_NI[Result_NI$Type_1<=0.025,]; 
Result_NI_cons = Result_NI_cons[order(Result_NI_cons$lambdaP),]
NI_cons <- Result_NI_cons %>%
	group_by(lambdaA) %>%
	summarize(lambdaP = min(lambdaP))
NI_cons_cut = lm(lambdaA ~ 0 + lambdaP, data = NI_cons)$coef; 1- NI_cons_cut # 0.8320723

Result_ACCF_cons = Result_ACCF[Result_ACCF$Type_1<=0.025,]; 
ACCF_cons <- Result_ACCF_cons %>%
	group_by(lambdaA) %>%
	summarize(lambdaP = min(lambdaP))
ACCF_cons_cut = mean(ACCF_cons$lambdaP); ACCF_cons_cut # 0.02909756

Result_cACCF_cons = Result_cACCF[Result_cACCF$Type_1<=0.025,]; 
cACCF_cons <- Result_cACCF_cons %>%
	group_by(lambdaA) %>%
	summarize(lambdaP = min(lambdaP))
cACCF_cons_cut = mean(cACCF_cons$lambdaP); cACCF_cons_cut # 0.02329268

R_PA = 10;
pdf('Simu_Violation_CAB.pdf',height=4,width=12,onefile=TRUE)
p1 = ggplot(data = Result_NI, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Type_1_cut), colour = "white") + scale_fill_brewer() + geom_abline(intercept = 0, slope = 1/R_PA, color="red", linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + theme(legend.position = "none") + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3)  + geom_abline(intercept = 0, slope = NI_cons_cut, color="yellow", linetype="dashed", size=1.5)
p2 = ggplot(data = Result_ACCF, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Type_1_cut), colour = "white") + scale_fill_brewer() + geom_vline(xintercept=0.03,color='red', linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + theme(legend.position = "none") + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3) + geom_vline(xintercept= ACCF_cons_cut,color='yellow', linetype="dashed", size=1.5)
p3 = ggplot(data = Result_cACCF, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Type_1_cut), colour = "white") + scale_fill_brewer() + geom_vline(xintercept=0.03,color='red', linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3) + geom_vline(xintercept= cACCF_cons_cut,color='yellow', linetype="dashed", size=1.5)

ggarrange(p1,p2,p3,labels = c("(a)","(b)","(c)"), ncol = 3, nrow = 1,widths=c(1,1,1.35))
dev.off()


R_PA = 10;
pdf('Simu_Violation_CAB_slide.pdf',height=4,width=9,onefile=TRUE)
p1 = ggplot(data = Result_NI, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Type_1_cut), colour = "white") + scale_fill_brewer() + geom_abline(intercept = 0, slope = 1/R_PA, color="red", linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + theme(legend.position = "none") + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3)  + geom_abline(intercept = 0, slope = NI_cons_cut, color="yellow", linetype="dashed", size=1.5)
p2 = ggplot(data = Result_ACCF, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Type_1_cut), colour = "white") + scale_fill_brewer() + geom_vline(xintercept=0.03,color='red', linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + theme(legend.position = "none") + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3) + geom_vline(xintercept= ACCF_cons_cut,color='yellow', linetype="dashed", size=1.5)
p3 = ggplot(data = Result_cACCF, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Type_1_cut), colour = "white") + scale_fill_brewer() + geom_vline(xintercept=0.03,color='red', linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3) + geom_vline(xintercept= cACCF_cons_cut,color='yellow', linetype="dashed", size=1.5)

ggarrange(p2,p3, ncol = 2, nrow = 1,widths=c(1,1.35))
dev.off()

#### Power Cutoff Approximation ####
Result_NI_cons = Result_NI[Result_NI$Power>=0.5,]; range(Result_NI_cons$lambdaP) # [1] 0.014 0.050
Result_cACCF_cons = Result_cACCF[Result_cACCF$Power>=0.5,]; range(Result_cACCF_cons$lambdaP) # [1] 0.01 0.05

pdf('Simu_Violation_power_CAB.pdf',height=4,width=12,onefile=TRUE)
p1 = ggplot(data = Result_NI, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Power_cut), colour = "white") + scale_fill_brewer() + geom_abline(intercept = 0, slope = 1/R_PA, color="red", linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY)))+ theme(legend.position = "none") + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3)
p2 = ggplot(data = Result_ACCF, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Power_cut), colour = "white") + scale_fill_brewer() + geom_vline(xintercept=0.03,color='red', linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + theme(legend.position = "none") + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3)
p3 = ggplot(data = Result_cACCF, aes(x = lambdaP, y = lambdaA)) + geom_tile(aes(fill = Power_cut), colour = "white") + scale_fill_brewer() + geom_vline(xintercept=0.03,color='red', linetype="dashed", size=1.5) + labs(fill='',x = expression(lambda[P] (cases/PY)), y= expression(lambda[A] (cases/PY))) + geom_point(aes(x=0.03,y=0.03/R_PA),colour = 'red',shape=17,size = 3)
ggarrange(p1,p2,p3,labels = c("(a)","(b)","(c)"), ncol = 3, nrow = 1,widths=c(1,1,1.35))
dev.off()