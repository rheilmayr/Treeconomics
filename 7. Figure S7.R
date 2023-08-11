














#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Transect plots ------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# pred_dat <- plot_dat %>% 
#   select(cwd.q, pet.q, 
#          rwi_change = rwi_change, 
#          rwi_change_lb = rwi_change_lb, 
#          rwi_change_ub = rwi_change_ub) %>% 
#   mutate(scenario = "rwi_change")
# 
# pclim_dat <- plot_dat %>% 
#   select(cwd.q, pet.q, 
#          rwi_change = rwi_change_pclim, 
#          rwi_change_lb = rwi_change_pclim_lb, 
#          rwi_change_ub = rwi_change_pclim_ub) %>% 
#   mutate(scenario = "rwi_change_pclim")
# 
# 
# transect_dat <- pred_dat %>% 
#   rbind(pclim_dat) %>% 
#   mutate(scenario = fct_relevel(scenario, "rwi_change", "rwi_change_pclim"))
# 
# transect_1 <- transect_dat %>% 
#   filter(pet.q == 0) %>%
#   group_by(cwd.q, scenario) %>% 
#   summarise(rwi_change = mean(rwi_change),
#             rwi_change_lb = mean(rwi_change_lb),
#             rwi_change_ub = mean(rwi_change_ub)) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_ribbon(aes(ymin = rwi_change_lb,
#                   ymax = rwi_change_ub,
#                   fill = scenario),
#               alpha = 0.2) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-1, 0.2)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted")) +
#   scale_fill_manual(name = "Scenario",
#                     labels = c("Full model", 
#                                "Variable shift in climate,\nconstant sensitivity"), 
#                     values = c("dark blue", "dark red", "dark green")) +
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std above mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.75),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# transect_1
# 
# 
# transect_2 <- transect_dat %>% 
#   filter(pet.q == 2) %>% 
#   ggplot(aes(x = cwd.q, y = rwi_change, group = scenario, color = scenario)) +
#   geom_ribbon(aes(ymin = rwi_change_lb,
#                   ymax = rwi_change_ub,
#                   fill = scenario),
#               alpha = 0.2) +
#   geom_line(size = 2) +
#   theme_bw(base_size = 20)+
#   ylim(c(-1.1, 0.2)) +
#   xlim(c(-2, 2)) +
#   scale_linetype_manual(values=c("solid", "dotted", "dotted"))+
#   scale_fill_manual(name = "Scenario",
#                     labels = c("Full model", 
#                                "Variable shift in climate,\nconstant sensitivity"), 
#                     values = c("dark blue", "dark red", "dark green")) +
#   scale_color_manual(name = "Scenario",
#                      labels = c("Full model", 
#                                 "Variable shift in climate,\nconstant sensitivity"), 
#                      values = c("dark blue", "dark red", "dark green")) +
#   ggtitle("Historic PET = 1 std below mean") +
#   ylab("Predicted change in RWI") +
#   xlab("Historic CWD (Deviation from species mean)") +
#   theme(legend.position = c(.18,.25),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank()) +
#   geom_hline(yintercept = 0, linetype = "dashed", size = 1)
# transect_2
# 
#   locator <- rwi_bin + 
#   theme_bw(base_size = 20)+
#   theme(legend.position = c(.18,.83),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=18),
#         legend.background = element_blank())+
#   geom_hline(yintercept = 1, size = 1) + 
#   geom_hline(yintercept = -1, size = 1)
# 
# locator | transect_1 / transect_2




transect_0 <- plot_dat %>% 
  filter(pet.q == 0) %>%
  group_by(cwd.q) %>% 
  summarise(rwi_change = mean(rwi_dif),
            rwi_change_lb = mean(rwi_dif_lb),
            rwi_change_ub = mean(rwi_dif_ub)) %>% 
  ggplot(aes(x = cwd.q, y = rwi_change)) +
  geom_ribbon(aes(ymin = rwi_change_lb,
                  ymax = rwi_change_ub),
              alpha = 0.2) +
  geom_line(size = 2) +
  theme_bw(base_size = 20)+
  ylim(c(-0.3, 0.2)) +
  xlim(c(-2, 2)) +
  # ggtitle("Historic PET = historic species mean") +
  # ylab("Predicted difference in RWI change - neutral model vs ourse") +
  xlab("Historic CWD\n(Deviation from species mean)") +
  ylab("Difference in predicted RWI changes by 2100\n(Heterogeneous sensitivity model - Neutral model)") +
  theme(legend.position = c(.18,.75),
        legend.text = element_text(size=13),
        legend.title = element_text(size=18),
        legend.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)
transect_0
ggsave(paste0(wdir, "figures\\", "a5_dif_pred.svg"), transect_0, width = 8, height = 8)
