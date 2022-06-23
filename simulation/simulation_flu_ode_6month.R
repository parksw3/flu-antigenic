source("../R/simulate_ode.R")

simulation_flu_ode_6month_base <- simulate_ode(npi=0)
simulation_flu_ode_6month_20 <- simulate_ode()
simulation_flu_ode_6month_40 <- simulate_ode(npi=0.4)
simulation_flu_ode_6month_60 <- simulate_ode(npi=0.6)
simulation_flu_ode_6month_80 <- simulate_ode(npi=0.8)

save("simulation_flu_ode_6month_base",
     "simulation_flu_ode_6month_20",
     "simulation_flu_ode_6month_40",
     "simulation_flu_ode_6month_60",
     "simulation_flu_ode_6month_80",
     file="simulation_flu_ode_6month.rda")
