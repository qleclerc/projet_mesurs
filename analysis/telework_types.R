
# quick plot to show different "shapes" of telework impact on chronic disease risk

library(ggplot2)
library(cowplot)

res = data.frame(alpha = seq(0,1,0.01))

res$lin_inc = res$alpha*2
res$lin_dec = (1-res$alpha)*2
res$U = (res$alpha-0.5)^2*2
res$rev_U = (1-(res$alpha-0.5)^2)*2

custom_theme = theme(axis.title = element_text(size = 12),
                     axis.ticks = element_blank(),
                     axis.text = element_blank(),
                     plot.title = element_text(hjust = 0.5, face = "bold"))

pa = ggplot(res) +
  geom_line(aes(alpha, lin_inc), linewidth = 1) +
  theme_classic() +
  custom_theme +
  labs(x = "Telework frequency", y = "Risk of chronic disease", title = "Linear increasing")

pb = ggplot(res) +
  geom_line(aes(alpha, lin_dec), linewidth = 1) +
  theme_classic() +
  custom_theme +
  labs(x = "Telework frequency", y = "", title = "Linear decreasing")

pc = ggplot(res) +
  geom_line(aes(alpha, U), linewidth = 1) +
  theme_classic() +
  custom_theme +
  labs(x = "Telework frequency", y = "", title = "U-shaped")

pd = ggplot(res) +
  geom_line(aes(alpha, rev_U), linewidth = 1) +
  theme_classic() +
  custom_theme +
  labs(x = "Telework frequency", y = "", title = "Reverse U-shaped")

plot_grid(pa,pb,pc,pd, nrow=1,
          labels = c("a)"))

ggsave(here::here("figures", "figx.png"), height = 3)
