

setwd("/Users/Hause/Dropbox/Working Projects/Akina preregistered report/simulations")

dt1 <- fread()


dt2 <- dt1[, .(power_b1 = mean(sig_x1, na.rm = T), power_b2 = mean(sig_x2, na.rm = T), r_x1 = mean(r_x1), r_x2 = mean(r_x2)),
           keyby = .(N.test, b1.test, b2.test, varInt.test, varSlope_b1.test, varSlope_b2.test, varResid.test)]

# dt2 %>% fwrite("/Users/hause/desktop/power_simulations_avg.csv")

dt2[, .(r_x1 = min(r_x1), r_x2 = min(r_x2)), keyby = .(b1.test, b2.test)]

dt2[, mean(power_b1)]
dt2[, unique(b2.test)]
dt2[, mean(power_b2)]
dt2[, min(power_b2)]

dt2[N.test == 90, hist(power_b1)]

dt2[N.test == 60 & power_b1 < .60, ]

dt2[N.test == 60 & power_b1 > .90 & power_b2 > .90, .(min_r_x1 = min(r_x1), min_r_x2 = min(r_x2))]

dt2[power_b1 > .90 & power_b2 > .90, .(min_r_x1 = min(r_x1), min_r_x2 = min(r_x2)), by = N.test]

dt2[, .(min_r_x1 = min(r_x1),
        min_r_x2 = min(r_x2)),
    by = N.test]

avg <- dt2[, .(min_r_x1 = min(r_x1),
        min_r_x2 = min(r_x2),
        pow_b1 = mean(power_b1),
        pow_b2 = mean(power_b2)),
    keyby = .(N.test, b1.test)]
avg[, f := es(r = min_r_x1)$f, by = .(N.test, b1.test)]


ggplot(avg, aes(f, pow_b1, col = as.factor(N.test))) +
  geom_point() +
  # geom_line() +
  geom_hline(yintercept = 0.90, linetype = 'dashed') +
  labs(x = "Mininum simulated effect size (Cohen's f)", y = "Statistical power (alpha = 0.02)", col = "Sample size")
ggsave("simulations.jpg", dpi = 200, width = 7, height = 5)
