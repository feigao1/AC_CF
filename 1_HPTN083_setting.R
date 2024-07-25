source('Simulation_Master_Functions.R')
Inc_spec = Inc_spec_083()

#### NI Design ####
Anal_NI(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=123)
     # Design_type1 Design_power  Anal_type1 Anal_power
# [1,]        0.025          0.8 0.004096519        0.8
# [2,]        0.025          0.9 0.003583260        0.9
Empir_NI(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY Type1_NI Power_NI Type1_RAE Power_RAE
# [1,]        0.025          0.8 12016.30   0.0238   0.8053    0.0034    0.8006
# [2,]        0.025          0.9 16190.12   0.0245   0.9035    0.0025    0.9036


#### AC-CF Design with FU CF####
Empir_ACCF_FU(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY Type1_RAE Power_RAE
# [1,]        0.025          0.8     4942    0.0207    0.8439
# [2,]        0.025          0.9     6554    0.0220    0.9210

#### AC-CF Design with RA CF####
Empir_ACCF_RA(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power tau Nscreen Npos    Nrec Nneg Total_PY Type1_RAE Power_RAE
# [1,]        0.025          0.8   1    6391  959 69.7346 5432     5432    0.0216    0.8354
# [2,]        0.025          0.8   2    3922  588 42.7996 3334     6668    0.0206    0.8175
# [3,]        0.025          0.9   1    8080 1212 88.0998 6868     6868    0.0208    0.9208
# [4,]        0.025          0.9   2    4939  741 53.8185 4198     8396    0.0214    0.9022

#### Conservative AC-CF Design with FU CF####
Anal_cACCF_FU(Inc_spec, beta_set = c(0.8,0.9))
     # Design_type1 Design_power Anal_type1  Anal_power
# [1,] 0.025        0.8          0.002788086 0.8000419 
# [2,] 0.025        0.9          0.002863521 0.9000162 
Empir_cACCF_FU(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY Type1_RAE Power_RAE
# [1,]        0.025          0.8     8205    0.0038    0.8217
# [2,]        0.025          0.9    10938    0.0033    0.8991

#### Conservative AC-CF Design with RA CF####
Anal_cACCF_RA(Inc_spec, beta_set = c(0.8,0.9))
     # Design_type1 Design_power Design_FU Anal_type1  Anal_power
# [1,] 0.025        0.8          1         0.002789359 0.8000477 
# [2,] 0.025        0.8          2         0.003073366 0.8001337 
# [3,] 0.025        0.9          1         0.002798741 0.9000139 
# [4,] 0.025        0.9          2         0.003123023 0.9000361 

Empir_cACCF_RA(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power tau Nscreen Npos     Nrec  Nneg Total_PY Type1_RAE Power_RAE
# [1,]        0.025          0.8   1    9725 1459 106.1020  8266     8266    0.0031    0.8344
# [2,]        0.025          0.8   2    6158  924  67.2117  5234    10468    0.0031    0.8247
# [3,]        0.025          0.9   1   11920 1788 129.8169 10132    10132    0.0043    0.9178
# [4,]        0.025          0.9   2    7518 1128  82.0802  6390    12780    0.0029    0.9029

#### OneArm-CF Design with RA CF####
Empir_OneArmCF_FU(Inc_spec, beta_set = c(0.8,0.9), nrep=10000, seed=1234)
     # Design_type1 Design_power Total_PY Type1_AE Power_AE
# [1,]        0.025          0.8     2398   0.0218   0.8251
# [2,]        0.025          0.9     3792   0.0258   0.9017
