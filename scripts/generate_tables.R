# scripts/generate_tables.R
#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Generate all tables
# -----------------------------------------------------------------------------

library(here)
source(here("R", "replication_functions.R"))
library(dplyr)
library(haven)
library(fastDummies)
library(SuperLearner)
library(ranger)
library(gbm)
library(glmnet)

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

d = ncol(X)
p = sum(data$treatment) / nrow(xyData)
n = nrow(xyData)
alpha = 0.01

# Constructing potential outcomes
sl.lib <- c("SL.ranger", "SL.gam", "SL.gbm", "SL.glm", "SL.glmnet")

treated = xyData[xyData$treatment == 1,]
untreated = xyData[xyData$treatment == 0,]

treated = treated %>% dplyr::select(-treatment)
untreated = untreated %>% dplyr::select(-treatment)

mu_treated <- SuperLearner(Y=treated$y, X=treated %>% dplyr::select(-y), SL.library=sl.lib)
mu_untreated <- SuperLearner(Y=untreated$y, X=untreated %>% dplyr::select(-y), SL.library=sl.lib)

psihat_treated1 <- mu_treated$SL.predict
psihat_treated0 <- predict(mu_untreated, treated[,2:(d+1)])$pred

psihat_untreated1 <- predict(mu_treated, untreated[,2:(d+1)])$pred
psihat_untreated0 <- mu_untreated$SL.predict

treated_sim <- cbind(treated, y1=psihat_treated1, y0=psihat_treated0, treatment=1)
untreated_sim <- cbind(untreated, y1=psihat_untreated1, y0=psihat_untreated0, treatment=0)

sim.df <- rbind(treated_sim, untreated_sim)
sim.df = sim.df %>% dplyr::select(-y)
sim.df <- sim.df[order(as.numeric(rownames(sim.df))), ]
sim.df = sim.df %>% dplyr::select(-treatment)


# Scaling covariates and preparing PCA based methods
covariates <- data.frame(scale(sim.df[,1:d]))
xyData <- cbind(covariates, sim.df[,(d+1):(d+2)]) #Adding back in y0, y1

inv_cov_matrix = solve(cov(as.matrix(covariates)))
eigs = n * p * (1 - p) * eig(cov(as.matrix(covariates)))

pca = prcomp(as.matrix(covariates), scale = TRUE)
var_explained = pca$sdev^2 / sum(pca$sdev^2)

avg_explained = mean(var_explained)
kaiser_k = min(max(which(var_explained > avg_explained)), d)
weighted_mk = get_weighted_k_mahalanobis(eigs=eigs, d=d, alpha=alpha)

pcs = prcomp(covariates, scale=TRUE)$x

kaiser_df = data.frame(cbind(pcs[,1:kaiser_k], y0=xyData$y0, y1=xyData$y1))
weighted_df = data.frame(cbind(pcs[,1:weighted_mk], y0=xyData$y0, y1=xyData$y1))

kaiser_inv = solve(cov(as.matrix(kaiser_df[, 1:kaiser_k])))
weighted_inv = solve(cov(as.matrix(weighted_df[, 1:weighted_mk])))


p=0.5
rerands = 10000
set.seed(42)

md_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_MD, sims=10000, alpha=alpha, p=p)
ed_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_ED, sims=10000, alpha=alpha, p=p)
sed_a <- get_threshold(data=xyData, inv_cov=inv_cov_matrix, metric=get_SED, sims=10000, alpha=alpha, p=p)
kaiser_a <- get_threshold(data=kaiser_df, inv_cov=kaiser_inv, metric=get_MD, sims=10000, alpha=alpha, p=p)
weighted_a <- get_threshold(data=weighted_df, inv_cov=weighted_inv, metric=get_MD, sims=10000, alpha=alpha, p=p)

complete_var_rerands <- lapply(1:rerands, function(z) randomize(xyData, p=p))
md_rerands <- lapply(1:rerands, function(z) get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=md_a, metric=get_MD, p=p))
ed_rerands <- lapply(1:rerands, function(z) get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=ed_a, metric=get_ED, p=p))
sed_rerands <- lapply(1:rerands, function(z) get_rerand(data=xyData, inv_cov=inv_cov_matrix, threshold=sed_a, metric=get_SED, p=p))
kaiser_rerands <- lapply(1:rerands, function(z) get_rerand(data=kaiser_df, inv_cov=kaiser_inv, threshold=kaiser_a, metric=get_MD, p=p))
weighted_rerands <- lapply(1:rerands, function(z) get_rerand(data=weighted_df, inv_cov=weighted_inv, threshold=weighted_a, metric=get_MD, p=p))

complete.var <- sd(unlist(lapply(complete_var_rerands, FUN = getMeanDiff)))
md.var = sd(unlist(lapply(md_rerands, FUN = getMeanDiff)))
ed.var = sd(unlist(lapply(ed_rerands, FUN = getMeanDiff)))
sed.var = sd(unlist(lapply(sed_rerands, FUN = getMeanDiff)))
kaiser.var = sd(unlist(lapply(kaiser_rerands, FUN = getMeanDiff)))
weighted.var = sd(unlist(lapply(weighted_rerands, FUN = getMeanDiff)))

var.df = data.frame(cbind('Complete Randomization'=complete.var,
                          'Mahalanobis'=md.var,
                          'Euclidean'=ed.var,
                          'Sq_Euclidean'=sed.var,
                          'Kaiser'=kaiser.var,
                          'Weighted'=weighted.var))

write.csv(var.df, file=here('output', 'tables', 'table2.csv'), row.names = FALSE)