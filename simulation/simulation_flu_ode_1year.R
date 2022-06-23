source("../R/simulate_ode.R")

simulation_flu_ode_1year_base <- simulate_ode(npi=0, iend=21)
simulation_flu_ode_1year_20 <- simulate_ode(iend=21)
simulation_flu_ode_1year_40 <- simulate_ode(npi=0.4, iend=21)
simulation_flu_ode_1year_60 <- simulate_ode(npi=0.6, iend=21)
simulation_flu_ode_1year_80 <- simulate_ode(npi=0.8, iend=21)

save("simulation_flu_ode_1year_base",
     "simulation_flu_ode_1year_20",
     "simulation_flu_ode_1year_40",
     "simulation_flu_ode_1year_60",
     "simulation_flu_ode_1year_80",
     file="simulation_flu_ode_1year.rda")
