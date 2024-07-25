source('Simulation_Master_Functions.R')
Inc_spec = Inc_spec_CAB()

#### NI Design ####
Anal_NI(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=123)
    # Design_type1 Design_power  Anal_type1 Anal_power
# [1,]        0.025          0.8 0.003658476        0.8
# [2,]        0.025          0.9         NaN        NaN
Empir_NI(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY   Type1_NI  Power_NI    Type1_RAE Power_RAE
# [1,]        0.025          0.8 16738.04 0.02410849 0.7940733 0.0005022602 0.7845304
# [2,]        0.025          0.9 22355.79 0.02330000 0.8871000 0.0014000000 0.8907000


#### AC-CF Design with FU CF####
Empir_ACCF_FU(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY Type1_RAE Power_RAE
# [1,]        0.025          0.8     5074    0.0209     0.824
# [2,]        0.025          0.9     6858    0.0218     0.923

#### Conservative AC-CF Design with FU CF####
Anal_cACCF_FU(Inc_spec, beta_set = c(0.8,0.9))
     # Design_type1 Design_power Anal_type1  Anal_power
# [1,] 0.025        0.8          0.005139408 0.8000547 
# [2,] 0.025        0.9          0.004381345 0.9000183 
Empir_cACCF_FU(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY Type1_RAE Power_RAE
# [1,]        0.025          0.8     6378    0.0070    0.8266
# [2,]        0.025          0.9     8606    0.0058    0.9251