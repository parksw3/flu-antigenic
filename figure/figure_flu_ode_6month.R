library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())

load("../simulation/simulation_flu_ode_6month.rda")

n <- 300

simulation_all <- list(
  simulation_flu_ode_6month_base,
  simulation_flu_ode_6month_20,
  simulation_flu_ode_6month_40,
  simulation_flu_ode_6month_60,
  simulation_flu_ode_6month_80
)

simulation_all2 <- simulation_all %>%
  lapply(function(x) x[,c(1, (n+2):(2*n+1))]) %>%
  bind_rows(.id="sim") %>%
  gather(key, value, -sim, -time) %>%
  mutate(
    strain=as.numeric(gsub("I", "", key))
  )

simulation_all3 <- simulation_all2 %>%
  filter(time >= 18) %>%
  mutate(
    sim=factor(sim, levels=1:5,
               labels=c("No reduction", paste0(1:4*20, "% reduction")))
  )

simulation_all4 <- simulation_all3 %>%
  filter(strain >= 150, strain <= 280) %>%
  mutate(
    value2 = ifelse(value <= 1e-8, 1e-8, value)
  )

g1 <- ggplot(simulation_all4) +
  geom_tile(aes(time, strain, fill=log10(value2))) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("Strains", limits=c(150, 280), expand=c(0, 0)) +
  scale_fill_viridis_c("log(I)", breaks=c(-2, -4, -6, -8),
                       labels=c("-2", "-4", "-6", "< -8")) +
  facet_wrap(~sim, nrow=1)

ggsave("figure_flu_ode_6month_1.pdf", g1, width=12, height=3)
ggsave("figure_flu_ode_6month_1.png", g1, width=12, height=3)

simulation_all_r <- simulation_all3 %>%
  group_by(time, sim) %>%
  mutate(
    total=sum(value),
    rvalue=value/total
  ) %>%
  filter(strain >= 150, strain <= 280)

g2 <- ggplot(simulation_all_r) +
  geom_tile(aes(time, strain, fill=rvalue)) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("Strains", limits=c(150, 280), expand=c(0, 0)) +
  scale_fill_viridis_c("Relative frequency") +
  facet_wrap(~sim, nrow=1)

ggsave("figure_flu_ode_6month_2.pdf", g2, width=12, height=3)
ggsave("figure_flu_ode_6month_2.png", g2, width=12, height=3)

g3 <- ggplot(simulation_all_r) +
  geom_hline(yintercept=max(filter(simulation_all_r, sim=="No reduction")$total), lty=2) +
  geom_bar(aes(time, value, fill=strain), stat="identity", width=1/300) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("I", limits=c(0, 0.14), expand=c(0, 0)) +
  scale_fill_viridis_c() +
  facet_wrap(~sim, ncol=1) +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave("figure_flu_ode_6month_3.pdf", g3, width=12, height=6)
ggsave("figure_flu_ode_6month_3.png", g3, width=12, height=6)

simulation_all_diversity <- simulation_all3 %>%
  group_by(time, sim) %>%
  mutate(
    total=sum(value),
    rvalue=value/total
  ) %>%
  summarize(
    diversity=exp(-sum(rvalue * log(rvalue))),
    total=sum(value)
  )

g4 <- ggplot(simulation_all_diversity) +
  geom_line(aes(time, total*80+4), lty=2, col=2) +
  geom_line(aes(time, diversity)) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("First-order Hill diversity", expand=c(0, 0),
                     sec.axis = sec_axis(trans=~(.-4)/80,
                                         name="I")) +
  scale_fill_viridis_c() +
  facet_wrap(~sim, ncol=1) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line.y.right = element_line(color="red"),
    axis.ticks.y.right = element_line(color="red"),
    axis.text.y.right = element_text(color="red"),
    axis.title.y.right = element_text(color="red")
  )

ggsave("figure_flu_ode_6month_4.pdf", g4, width=12, height=6)

simulation_erate <- simulation_all3 %>%
  group_by(sim, time) %>%
  summarize(
    mean=weighted.mean(strain, value),
    total=sum(value)
  ) %>%
  arrange(sim, time) %>%
  mutate(
    erate=c(diff(mean), NA)*300
  )

g5 <- ggplot(simulation_erate) +
  geom_line(aes(time, total*500), lty=2, col=2) +
  geom_line(aes(time, erate)) +
  scale_x_continuous("Time (years)", expand=c(0, 0)) +
  scale_y_continuous("Evolution rate (per year)", expand=c(0, 0),
                     sec.axis = sec_axis(trans=~./500,
                                         name="I")) +
  scale_fill_viridis_c() +
  facet_wrap(~sim, ncol=1) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line.y.right = element_line(color="red"),
    axis.ticks.y.right = element_line(color="red"),
    axis.text.y.right = element_text(color="red"),
    axis.title.y.right = element_text(color="red")
  )

ggsave("figure_flu_ode_6month_5.pdf", g5, width=12, height=6)
