



library(lavaan)
library(igraph)
library(Matrix)
library(dplyr)

source('sem ram func.R')

###############################################################################

load('..\\data\\velgerundersokelsen2013.rdata')

velger2013 %>% 
	group_by(Sex) %>% 
	mutate_all(funs(. - mean(.))) %>% 
	mutate_all(funs(. / sd(.))) %>% 
	ungroup() -> velger2013_scaled 

velger2013_scaled %>% 
	group_by(Sex) %>% 
	summarise(n())

###############################################################################
model_string <- '
intr =~ Interessert + Prate + Bryseg
participation =~ DiskusjonInternett + Demo + Opprop

participation ~ intr

Lettbestemt ~ intr + participation
'

summary(sem(model_string, data=velger2013_scaled, group='Sex'))


###############################################################################



# Assymetric relations.
graph_df_fb <- data.frame(from=c(rep('intr', 3), rep('participation', 3), 'intr', 'intr', 'participation'),
													 to= c('Interessert', 'Prate', 'Bryseg', 'DiskusjonInternett', 'Demo', 'Opprop', 'participation', 'Lettbestemt', 'Lettbestemt'))

model_igraph_fb <- graph_from_data_frame(graph_df_fb)

plot(model_igraph_fb)



ram_matrices_vg <- make_ram_matrices(model_igraph_fb, latent_variables = c('intr', 'participation'))

# reproduce results from lavaan.
est1 <- estimate_ram_group(ram_matrices_vg, velger2013_scaled, s_equal = FALSE,
													 group_var='Sex', difference_param=FALSE)


est1$res_params[[1]]
est1$res_params[[2]]

est1$ram_matrices_est[[1]][['A']]
est1$ram_matrices_est[[2]][['A']]




est_lasso1_1 <- estimate_ram_group_lasso(ram_matrices_vg, velger2013_scaled,
																				 group_var = 'Sex', n_lambdas = 35,
																				 verbose=TRUE,
																				 refit_var = TRUE,
																				 refit_cov=TRUE,
																				 refit_means = TRUE)


plot_lasso_trace(est_lasso1_1, plot_min_lambda=FALSE, plot_lambda = FALSE,
								 xlim=c(0.31,0))


########################################
######## CV
########################################

#velger2013_scaled %>% 
#	mutate(IDX=1:nrow(.)) -> velger2013_scaled

cv_seed <- 65468
set.seed(cv_seed)
est_lasso1_1_cv <- cv_ram_group_lasso(ram_matrices_vg, velger2013_scaled,
															 n_folds=2, n_cvs = 10,
															 group_var = 'Sex', n_lambdas = 15,
															 refit_var = TRUE,
															 refit_cov=TRUE,
															 refit_means = TRUE)

save(list=c('est_lasso1_1_cv', 'cv_seed'), file='cv_velgunder_20171120.rdata')

plot_cv_results(est_lasso1_1_cv)


plot(est_lasso1_1_cv$lambdas_cv,
		 est_lasso1_1_cv$mean_cv_fit, type='l',
		 ylab='Mean out of sample error', xlab='lambda')


best_lambda <- est_lasso1_1_cv$lambdas_cv[which.min(est_lasso1_1_cv$mean_cv_fit)]

# 
# est_lasso1_1_cv$mean_cv_fit + apply(est_lasso1_1_cv$lambda_fits_cv, MARGIN = 2, sd)
# 
# 
# est_lasso1_1_cv$lambdas_cv[which.min(est_lasso1_1_cv$mean_cv_fit)]




est_lasso1_2 <- estimate_ram_group_lasso(ram_matrices_vg, velger2013_scaled,
																				 group_var = 'Sex', lambda_values = best_lambda,
																				 verbose=TRUE,
																				 refit_var = TRUE,
																				 refit_cov=TRUE,
																				 refit_means = TRUE)

est_ram_best <- get_estimated_ram_matrices(est_lasso1_2)
est_ram_best[[1]]

est_lasso1_2$sem_res_init$ram_matrices_est[[2]][['A']]


est1$ram_matrices_est[[2]][['A']]
est1$ram_matrices_est[[1]][['A']]

est1$ram_matrices_est[[1]][['A']] + est1$ram_matrices_est[[2]][['A']] 

est_ram_best[[1]][[1]][['A']]
est_ram_best[[1]][[2]][['A']]






