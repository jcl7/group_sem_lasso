
###############################################################################
# This is the script that analyses the Affect Valuation Theory data
# in the paper. The resampling analysis of these data can be found
# in the AVT bootstrap folder.
###############################################################################



library(lavaan)
library(igraph)
library(Matrix)
library(dplyr)
library(reshape2)
library(ggplot2)
library(xtable)


source('sem ram func.R')


###################################################
### Affect Valuation Theory model
##################################################
# http://web.stanford.edu/class/psych253/section/section_8/section8.html#question-c-group-comparisons
# http://journals.sagepub.com/doi/pdf/10.1111/j.1745-6916.2007.00043.x

dta <- read.csv("..\\data\\jt-data1.csv")

cc_idx <- complete.cases(dta[,c('ideahap', 'cultatt',  'temperatt', 'actuhap', 'rigoract', 'depress')])

dta <- dta[cc_idx,]

grplab1 <- c("EA", "AA", "CH") # European American, Asian American, Hong Kong Chinese.
dta$group <- grplab1[dta$group]

#dta <- dta[dta$group != 'AA',] # Remove Asian Americans.

dta %>% 
	group_by(group) %>% 
	summarise(N=n())

# SCALE AND CENTERERR!!!!!!
dta %>%  group_by(group) %>% 
	mutate_all(funs(. - mean(.))) %>% 
	mutate_all(funs(. / sd(.))) %>% 
	ungroup() -> dta_scaled 


dta_scaled %>% 
	group_by(group) %>% 
	summarise(n())

avt.model1 <- '
# regressions 
ideahap ~ cultatt + temperatt
actuhap ~ cultatt + temperatt
rigoract ~ ideahap + actuhap + cultatt + temperatt
depress ~ ideahap + actuhap + cultatt + temperatt

# variances and covariances
# X ~~ Y
rigoract ~~ 0*depress
ideahap ~~ actuhap
cultatt ~~ temperatt
'

avt_res <- sem(avt.model1, fixed.x = FALSE, data = dta_scaled, group='group')
summary(avt_res)


# Assymetric relations.
graph_df_avt <- data.frame(from=c('cultatt', 'temperatt', 'cultatt', 'temperatt', 'ideahap', 'actuhap', 'cultatt', 'temperatt', 'ideahap', 'actuhap', 'cultatt', 'temperatt'),
															to=c('ideahap', 'ideahap', 'actuhap', 'actuhap', 'rigoract', 'rigoract', 'rigoract', 'rigoract', 'depress', 'depress', 'depress', 'depress'))

# covariances
graph_df_avt <- rbind(graph_df_avt,
												 data.frame(from=c('ideahap', 'actuhap', 'cultatt', 'temperatt'),
					 										to=c('actuhap', 'ideahap', 'temperatt', 'cultatt')))


model_igraph_avt <- graph_from_data_frame(graph_df_avt)

#write_graph(model_igraph_avt, file = 'model_igraph_avt1', format = 'dot')

model_plot_layout <- layout_with_sugiyama(model_igraph_avt)
plot(model_igraph_avt, layout=model_plot_layout$layout)

#save(list=c('model_plot_layout', 'model_igraph_avt'), file='avt_model_igraph.rdata')


ram_matrices_avt <- make_ram_matrices(model_igraph_avt)


est1 <- estimate_ram_group(ram_matrices_avt, dta_scaled,
													 group_var='group',
													 difference_param=FALSE)

# Group 1: EA
est1$ram_matrices_est[[1]][['A']]
est1$ram_matrices_est[[1]][['M']]
est1$ram_matrices_est[[1]][['S']]

# Group 2: AA
est1$ram_matrices_est[[2]][['A']]
est1$ram_matrices_est[[2]][['M']]
est1$ram_matrices_est[[2]][['S']]

# Group 3: CH
est1$res_params[[1]]
est1$ram_matrices_est[[3]][['A']]
est1$ram_matrices_est[[3]][['M']]
est1$ram_matrices_est[[3]][['S']]

lapply(est1$model_sigmas, FUN= function(x) eigen(x)$values)
lapply(est1$obs_cov, FUN= function(x) eigen(x)$values)



### LASSO!!

est_lasso_avt <- estimate_ram_group_lasso(ram_matrices_avt, dta_scaled,
																				 group_var = 'group', n_lambdas = 20,
																				 verbose=TRUE,
																				 refit_var = TRUE,
																				 refit_cov=TRUE,
																					refit_means = TRUE)


#save(est_lasso_avt, file='avt_3g_20180418.rdata')

trace_grp2 <- cbind(tidy_lasso_trace(est_lasso_avt, grp=2), Group='Asian Americans') 
trace_grp3 <- cbind(tidy_lasso_trace(est_lasso_avt, grp=3), Group='Hong Kong Chinese') 

trace_lasso_both_grps <- bind_rows(trace_grp2, trace_grp3)


ggplot(trace_lasso_both_grps, aes(Lambda, Delta_value, group=Coef)) +
	geom_line() +
	facet_wrap(~ Group) +
	theme_bw() +
	xlab('Lambda') +
	ylab('Delta parameter value') +
	geom_abline(intercept = 0, slope=0) + 
	scale_x_continuous(limits=c(0.48, 0), trans='reverse') +
	theme(panel.grid.major.x = element_blank(),
				panel.grid.major.y = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.minor.y = element_blank()) 
	

################################################################
# Cross validation
#################################################################


cv_seed <- 3803881
set.seed(cv_seed)

cv_res_avt <- cv_ram_group_lasso2(ram_matrices_avt, dta_scaled,
															group_var = 'group', n_cvs=5, n_lambdas_cv=20,
															refit_var = TRUE,
															refit_cov=TRUE,
															refit_means = TRUE)

#save(list=c('cv_res_avt', 'cv_seed'), file='cv_avt_3g_20180222.rdata')

#### plot out-of-sample errors
avg_grp_cv_error <- apply(cv_res_avt$lambda_fits_cv, MARGIN = c(1, 2), mean)

grp_labels <- c('GRP.2' = 'Asian Americans', 'GRP.3'='Hong Kong Chinese')

data.frame(Lambda=cv_res_avt$lambdas_cv, GRP=avg_grp_cv_error) %>% 
	melt(id.vars='Lambda') %>% 
	filter(variable != 'GRP.1') %>% 
	ggplot(aes(Lambda, value)) +
	geom_line(lwd=1.1) +
	facet_wrap(~variable, scales='free_y', labeller = as_labeller(grp_labels)) +
	ylab('Mean out-of-sample error') +
	theme_bw() +
	theme(panel.grid.major.x = element_blank(),
				panel.grid.major.y = element_blank(),
				panel.grid.minor.x = element_blank(),
				panel.grid.minor.y = element_blank()) 




##### fit the model on the complete data set, using the optimal lambdas.

load('cv_avt_3g_20180222.rdata')

avg_grp_cv_error <- apply(cv_res_avt$lambda_fits_cv, MARGIN = c(1, 2), mean)

best_lambda_grp2 <- cv_res_avt$lambdas_cv[which.min(avg_grp_cv_error[,2])]
best_lambda_grp3 <- cv_res_avt$lambdas_cv[which.min(avg_grp_cv_error[,3])]

print(c(best_lambda_grp2, best_lambda_grp3))
best_lambdas <- unique(c(best_lambda_grp2, best_lambda_grp3))

est_lasso_avt_best_lambda <- estimate_ram_group_lasso(ram_matrices_avt, dta_scaled,
																					group_var = 'group',
																					lambda_values = best_lambdas,
																					verbose=TRUE,
																					refit_var = TRUE,
																					refit_cov=TRUE,
																					refit_means = TRUE)


get_nice_par_table(est_lasso_avt_best_lambda, which_group = 1)

#save(est_lasso_avt_best_lambda, file='avt_3g_optimal_lambda_20180222.rdata')


grp2_lasso_bestl_idx <- which(est_lasso_avt_best_lambda$lambdas == best_lambdas)
get_estimated_ram_matrices(est_lasso_avt_best_lambda)[[grp2_lasso_bestl_idx]][[2]][['A']]


grp3_lasso_bestl_idx <- which(est_lasso_avt_best_lambda$lambdas == best_lambdas)
get_estimated_ram_matrices(est_lasso_avt_best_lambda)[[grp3_lasso_bestl_idx]][[3]][['A']]






dim(est_lasso_avt_best_lambda$sem_res_init$ram_matrices$A)

aii <- arrayInd(est_lasso_avt_best_lambda$sem_res_init$asym_idx,
				 .dim=dim(est_lasso_avt_best_lambda$sem_res_init$ram_matrices$A))

# from
colnames(est_lasso_avt_best_lambda$sem_res_init$ram_matrices$A)[aii[,1]]

# to 
rownames(est_lasso_avt_best_lambda$sem_res_init$ram_matrices$A)[aii[,2]]

# par est
get_estimated_ram_matrices(est_lasso_avt_best_lambda)[[1]][[2]][['A']][est_lasso_avt_best_lambda$sem_res_init$asym_idx]



tab_grp1 <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 1)

tab_grp2_lasso <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 2)
tab_grp2 <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 2, which_lambda = 0)

tab_grp3_lasso <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 3)
tab_grp3 <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 3, which_lambda = 0)


xtable(
	cbind(BETA1=tab_grp1, DELTA_grp2=tab_grp2[,3], DELTA_LASSO_grp2=tab_grp2_lasso[,3],
				DELTA_grp3=tab_grp3[,3], DELTA_LASSO_grp3=tab_grp3_lasso[,3])
)




