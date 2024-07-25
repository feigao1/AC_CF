source('Simulation_Master_Functions.R')
Inc_spec = Inc_spec_083()
lambdaP_grid = seq(0.01,0.05,0.001); lambdaA_grid = lambdaP_grid/Inc_spec$lambdaP*Inc_spec$lambdaA

res = Violation(Inc_spec, lambdaP_grid, lambdaA_grid, Method = c('NI'), beta = 0.8, nrep=10000, seed=1234)
saveRDS(res, file = "Violation_Result_NI.rds")
res = Violation(Inc_spec, lambdaP_grid, lambdaA_grid, Method = c('cAC-CF'), beta = 0.8, nrep=10000, seed=1234)
saveRDS(res, file = "Violation_Result_cACCF.rds")
res = Violation(Inc_spec, lambdaP_grid, lambdaA_grid, Method = c('AC-CF'), beta = 0.8, nrep=10000, seed=1234)
saveRDS(res, file = "Violation_Result_ACCF.rds")
res = Violation(Inc_spec, lambdaP_grid, lambdaA_grid, Method = c('OneArm-CF'), beta = 0.8, nrep=10000, seed=1234)
saveRDS(res, file = "Violation_Result_Onearm_CF.rds")