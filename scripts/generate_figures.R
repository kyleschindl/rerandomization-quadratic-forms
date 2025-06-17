# scripts/generate_figures.R
#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Generate all figures
# -----------------------------------------------------------------------------

library(here)
source(here("R", "replication_functions.R"))
library(reshape2)
library(ggplot2)
library(haven)
library(fastDummies)
library(dplyr)

### Figure 1
set.seed(42)
simulated_nu = lapply(1:10000, function(i) simulate_nu_df(n=500,
                                                          d=50,
                                                          tau=1,
                                                          conc=1,
                                                          beta='uniform'))

combined_array <- abind(simulated_nu, along = 3)
nu.df <- data.frame(apply(combined_array, c(1, 2), mean))
plt.df = melt(nu.df, id.vars='idx')

p1 <- ggplot(plt.df, aes(x=idx, y=value, color=variable)) +
  geom_point(size = 0.5) +
  geom_line() + 
  theme_minimal() +
  xlab('Ordered Eigenvalue') +
  ylab('Variance reduction factor') +
  ylim(0, 1.01) +
  labs(color='Method') + 
  scale_color_discrete(labels = c("Mahalanobis",
                                  "Euclidean",
                                  "Sq. Euclidean",
                                  'Kaiser',
                                  'Wtd. Eigenvalue'))

ggsave(filename=here('output', 'figures', 'figure1.png'),
       plot=p1, width=6, height=4, units='in')

write.csv(nu.df,
          file=here('output', 'tables', 'averaged_nus.csv'),
          row.names = FALSE)

### Figure 2
set.seed(42)

d_list = c(5, 25, 50, 75, 100, 150, 200, 250)
gamma_vals = c(0.05, 0.5, 1)

combinations <- expand.grid(d=d_list,gamma=gamma_vals)
n = 500
rerands = 2500
tau = 1
p = 0.5
alpha = 0.01

set.seed(42)
results = list()
for(i in 1:length(combinations$d)) {
  
  results[[i]] = simulate_var_study(n=n,
                                    d=combinations[i,]$d,
                                    tau=tau,
                                    conc=combinations[i,]$gamma,
                                    p=p,
                                    alpha=alpha,
                                    rerands=rerands,
                                    a_sims=250)
}

csv_list <- lapply(seq_len(nrow(combinations)), function(i) {
  d_i     <- combinations$d[i]
  gamma_i <- combinations$gamma[i]
  var.df  <- read.csv(
    here('output', 'calculations', sprintf('var_%s_%s_alpha_01_a2500.csv', d_i, gamma_i))
  )
  var.df[, 1:7] <- sqrt(n) * var.df[, 1:7]
  var.df[, 1:7] <- var.df[, 1:7] / var.df$Complete.Randomization
  subset(var.df, select = -Complete.Randomization)
})





results.combined <- do.call(rbind, csv_list)

write.csv(results.combined,
          file=here('output', 'tables', 'simulation_study.csv'),
          row.names = FALSE)

plt.results <- reshape2::melt(results.combined,
                              id.vars = c("d", "Concentration", "beta"))

p2 <- ggplot(plt.results, aes(x = as.numeric(d),
                              y = as.numeric(value),
                              color = variable,
                              group = variable)) +
  geom_point(size = 0.85) +
  geom_line() +
  facet_grid(beta ~ Concentration, scales = "free_y") +
  theme_minimal() +
  xlab("Number of Covariates")


ggsave(filename=here('output', 'figures', 'figure2.png'),
       plot=p2, width=6, height=4, units='in')


### Figure 3

data = read_dta(here('data', 'DataMentoring.dta'))
data = data.frame(data)

#the paper only keeps the following data:
data = subset(data, Pairs_NoTwin_w1 == 0)

X = subset(data, select = c(male_w1,
                            age_w1,
                            migrant_w1,
                            books_w1,
                            grade_math_gr_w1,
                            grade_german_gr_w1,
                            grade_english_gr_w1,
                            privatteach_w1,
                            parentssupport_w1,
                            bfi_con_w1,
                            bfi_neur_w1,
                            bfi_open_w1,
                            bfi_extra_w1,
                            bfi_agree_w1,
                            HigherSES))

X = dummy_cols(X, select_columns = c("grade_math_gr_w1",
                                     "grade_german_gr_w1",
                                     "grade_english_gr_w1"),
               remove_selected_columns = TRUE,
               remove_first_dummy = TRUE)

# Technically the indicators for 2 and 3 are still missing.
# For simplicity, we can impute 0s:
X$grade_math_gr_w1_2 = ifelse(is.na(X$grade_math_gr_w1_2), 0,
                              X$grade_math_gr_w1_2)

X$grade_math_gr_w1_3 = ifelse(is.na(X$grade_math_gr_w1_3), 0,
                              X$grade_math_gr_w1_3)

X$grade_german_gr_w1_2 = ifelse(is.na(X$grade_german_gr_w1_2), 0,
                                X$grade_german_gr_w1_2)

X$grade_german_gr_w1_3 = ifelse(is.na(X$grade_german_gr_w1_3), 0,
                                X$grade_german_gr_w1_3)

X$grade_english_gr_w1_2 = ifelse(is.na(X$grade_english_gr_w1_2), 0,
                                 X$grade_english_gr_w1_2)

X$grade_english_gr_w1_3 = ifelse(is.na(X$grade_english_gr_w1_3), 0,
                                 X$grade_english_gr_w1_3)

treatment = data$treatment_w1
y = data$I_lmprospect_w2

xyData <- data.frame(cbind(y, treatment, X))
xyData <- na.omit(xyData)


plot.pca = prcomp(as.matrix(xyData %>% dplyr::select(-c(y, treatment))), scale = TRUE)
plot.pca$sdev^2 / sum(plot.pca$sdev^2)

d = ncol(xyData %>% dplyr::select(-c(y, treatment)))

p3 <- ggplot() +
  geom_point(aes(x=1:d, y=plot.pca$sdev^2 / sum(plot.pca$sdev^2))) +
  geom_line(aes(x=1:d, y=plot.pca$sdev^2 / sum(plot.pca$sdev^2))) +
  ylim(0, 1) +
  theme_minimal()

ggsave(filename = here('output', 'figures', 'figure3.png'),
       plot=p3, width=6, height=4, units='in')

scree.df = data.frame(cbind(d=1:d, var.explained=plot.pca$sdev^2 / sum(plot.pca$sdev^2)))

write.csv(scree.df,
          file=here('output', 'tables', 'scree_dataframe.csv'),
          row.names = FALSE)
