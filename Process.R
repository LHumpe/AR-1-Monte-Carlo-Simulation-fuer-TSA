source("Functions.R")

if (!require('ggplot2')) install.packages('ggplot2'); 
if (!require('doParallel')) install.packages('doParallel'); 
if (!require('gridExtra')) install.packages('gridExtra'); 
if (!require('kableExtra')) install.packages('kableExtra');
if (!require('matlib')) install.packages('matlib');

library('ggplot2')
library('doParallel')
library('gridExtra')
library('kableExtra')

no_cores <- 6
registerDoParallel(cores=no_cores)  

set.seed(123)

mu = 0
sigma = 1 
theorethical_phi = seq(from = -0.9, to = 0.9, by = 0.1) 
N = c(25, 50, 100, 500, 1000) 


result_frame = data.frame(phi = c(0),
                          N = c(0),
                          phi_est = c(0),
                          bias = c(0),
                          se_phi = c(0))

i = 1
for (ac in theorethical_phi){
  sim_result <- foreach(i = N, .packages = 'matlib') %dopar% mc.sim_pacf(i, 
                                                                         ac, 
                                                                         mu, 
                                                                         sigma)
  for (t in 1:length(N)){
    result_frame[i,] <- c(ac, N = N[t], 
                          phi_est = sim_result[[t]][[1]] , 
                          bias = sim_result[[t]][[2]], 
                          se_phi = sim_result[[t]][[3]])
    i <- i + 1
  }
}


ggplot(result_frame, aes(x = phi, y = se_phi, group = factor(N))) +
  geom_line() +
  geom_point(aes(shape = factor(N))) +
  scale_y_continuous(expression(SD ~ group("(",hat(phi),")"))) + 
  scale_x_continuous(expression(phi), 
                     breaks = theorethical_phi, 
                     labels = as.character(theorethical_phi)) + 
  scale_shape_discrete(labels = lapply(N, 
                                       function(i) paste("T=", i, sep=" "))) +
  labs(shape="") + 
  theme_minimal()


par(mfrow=c(1,2))

bias_plot <- ggplot(data = result_frame, aes(x = phi, 
                                             y = bias, 
                                             group = factor(N))) + 
  geom_line() +
  geom_point(aes(shape = factor(N))) +
  scale_y_continuous(expression(Bias~of~hat(phi))) + 
  scale_x_continuous(expression(phi), 
                     breaks=theorethical_phi, 
                     labels = as.character(theorethical_phi)) +
  scale_shape_discrete(labels = lapply(N,function(i) paste("T=", i, sep=" "))) +
  labs(shape = "") + 
  theme_minimal()

est_plot <- ggplot(result_frame, aes(x = phi, 
                                     y = phi_est, 
                                     group = factor(N))) +
  geom_line() +
  geom_point(aes(shape = factor(N))) + 
  scale_y_continuous(expression(hat(phi)), 
                     breaks = theorethical_phi, 
                     labels = as.character(theorethical_phi)) +
  scale_x_continuous(expression(phi), 
                     breaks = theorethical_phi, 
                     labels = as.character(theorethical_phi)) +
  scale_shape_discrete(labels = lapply(N,function(i) paste("T=", i, sep=" "))) +
  labs(shape = "") + 
  theme_minimal()

grid.arrange(est_plot, bias_plot, nrow = 2)

tbl_reshape <- reshape(result_frame, 
                       idvar = "phi", 
                       timevar = "N", 
                       direction = "wide")

kable_column_names = c("Phi", rep(c("Phi_Est", "Bias", "Phi_Sd"), 5))

kable(tbl_reshape, digits = 3,
      "latex",
      row.names = F, 
      col.names = kable_column_names, 
      booktabs = T, 
      escape = T, 
      align = 'c') %>% 
  add_header_above(c("","T=25"=3,"T=50"=3,"T=100"=3,"T=500"=3,"T=1000"=3)) %>%
  kable_styling(latex_options = "hold_position") %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = T) %>%
  collapse_rows(columns = 1)

