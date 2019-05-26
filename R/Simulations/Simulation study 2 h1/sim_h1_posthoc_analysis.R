

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(lavaan)
library(xtable)

source('..\\..\\sem ram func.R')


nns <- c(142, 300, 500, 1000)



###################################################
# Compare with the true model
###################################################

load('H1_TRUEMODEL_partable.rdata')

# ustart : the paramters used to simulate data.
mypar %>% 
	select(lhs, op, rhs, group, ustart) %>% 
	filter(op == '~') %>% 
	mutate(group = paste('group_', group, sep='')) %>% 
	dcast(lhs + op + rhs ~ group) %>% 
	mutate(true_diff = group_2 - group_1,
				 lhs = as.character(lhs),
				 rhs = as.character(rhs),
				 Parameter = paste(lhs, rhs, sep = ' ~ '))   -> true_parameters


xtable(true_parameters[,c(1, 3,4,5, 6)])


all_param_ests <- data.frame()


for (nn in 1:length(nns)){
	for (ii in 1:100){
		
		print(ii)
		fname <- sprintf('cv_avt_sim_20180625_n%d_%04d.rdata', nns[nn], ii)
		load(fname)

		get_nice_par_table(est_lasso_avt_best_lambda, which_group=2) %>% 
			mutate(N = nns[nn],
						 simno = ii) %>% 
			bind_rows(all_param_ests, .) -> all_param_ests
		
	}
}



head(all_param_ests)


all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ ')) %>% 
	left_join(true_parameters %>% select(Parameter, true_diff), by='Parameter') %>% 
	mutate(SqE = (true_diff - Estimate)^2,
				 est_error = true_diff - Estimate) -> all_param_ests_tmp

head(all_param_ests_tmp, 15)


all_param_ests_tmp %>% 
	group_by(Parameter, N) %>% 
	summarise(MSE = mean(SqE),
						Mean_bias = mean(est_error)) %>% 
	left_join(true_parameters %>% select(Parameter, true_diff))



all_param_ests_tmp %>% 
	group_by(Parameter, N) %>% 
	summarise(Prop0 = sum(Estimate == 0) / 100,
						Mean_est = mean(Estimate),
						MSE = mean(SqE),
						Mean_bias = mean(est_error),
						true_diff = mean(true_diff)) %>% 
	left_join(true_parameters %>% select(Parameter, true_diff), by='Parameter') -> ddd



varnames <- c('actuhap' = 'Actual Affect', 'depress' = 'Depression',
							'ideahap' = 'Ideal Affect', 'rigoract'='Rigorous Activties',
							'cultatt'='Attitudes', 'temperatt'='Temparement')

column_order <- c("Response", "Predictor", "group_1", #"group_2", 
									"true_diff", "142_Prop0", 
									"142_MSE", "142_Mean_bias", "300_Prop0","300_MSE", "300_Mean_bias",
									"500_Prop0","500_MSE" , "500_Mean_bias","1000_Prop0", "1000_MSE", "1000_Mean_bias")


ddd %>% 
	select(Parameter, N, Prop0, MSE, Mean_bias) %>% 
	melt(id.vars=c('Parameter', 'N')) %>% 
	dcast(Parameter ~ N + variable) %>% 
	separate(Parameter, into=c('Response', 'Predictor')) %>% 
	left_join(true_parameters, by=c('Response'='lhs', 'Predictor'='rhs')) %>% 
	dplyr::select(-op, -Parameter) %>% 
	mutate(Response = varnames[Response],
				 Predictor = varnames[Predictor]) %>% 
	arrange(Predictor) %>% 
	select(column_order) %>% 
	xtable()
	



all_param_ests_tmp %>%
	mutate(N_ = as.factor(N/2)) %>% 
	ggplot(aes(x=N_, y=Estimate)) +
	geom_hline(aes(yintercept=true_diff), color='red') +
	geom_violin(draw_quantiles = 0.5) +
	facet_wrap( ~ Parameter) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
				panel.grid.minor = element_blank(),
				panel.background = element_blank()) +
	xlab('Sample Size') +
	ylab(expression(delta ~ 'estimate')) +
	labs(title='Simulation study 2')

ggsave('simstudy2.png')
#


####################################################################################
####################################################################################
####################################################################################


for (nn in 1:length(nns)){
	for (ii in 1:100){
		
		print(ii)
		fname <- sprintf( 'cv_avt_sim_20180625_n%d_%04d.rdata', nns[nn], ii)
		load(fname)
		#print(best_cv_lamba)
		
		mypar0 %>% 
			left_join(get_nice_par_table(est_lasso_avt_best_lambda, which_group=2),
								by = c('lhs'='Response', 'rhs'='Predictor')) %>% 
			rename(True_param = est) %>% 
			filter(!is.na(Estimate)) %>% 
			mutate(Error = True_param - Estimate,
						 NN = nns[nn],
						 simno = ii) %>% 
			bind_rows(all_param_ests, .) -> all_param_ests
		
		
		
	}
}


all_param_ests %>% 
	filter(lhs == 'ideahap',
				 rhs == 'cultatt') %>% 
	ggplot(aes(y=Estimate, x=NN)) +
	geom_point() +
	geom_hline(aes(yintercept=True_param))



all_param_ests %>% 
	filter(lhs == 'rigoract',
				 rhs == 'cultatt') %>% 
	ggplot(aes(y=Estimate, x=NN)) +
	geom_point() +
	geom_hline(aes(yintercept=True_param))





#################################################
# Explore lambda values across sample sizes
#################################################


best_lambdasss <- array(dim=c(length(nns), 100, 2))



for (nn in 1:length(nns)){
	for (ii in 1:100){
		
		print(ii)
		fname <- sprintf( 'cv_avt_sim_20180625_n%d_%04d.rdata', nns[nn], ii)
		load(fname)
		
		best_lambdasss[nn, ii,] <- best_cv_lamba
		
		if (length(best_cv_lamba) == 1){
			best_lambda_grp2 <- best_cv_lamba
			best_lambda_grp3 <- best_cv_lamba
		} else {
			best_lambda_grp2 <- best_cv_lamba[1]
			best_lambda_grp3 <- best_cv_lamba[2]
		}
		
		bl_idx_grp2 <- which(est_lasso_avt_best_lambda$lambdas == best_lambda_grp2)
		bl_idx_grp3 <- which(est_lasso_avt_best_lambda$lambdas == best_lambda_grp3)
		
		
	}
}



boxplot(best_lambdasss[1,,1])
boxplot(best_lambdasss[2,,1])
boxplot(best_lambdasss[3,,1])
boxplot(best_lambdasss[4,,1])


melt(best_lambdasss) %>% 
	rename(N = Var1, 
				 Iter = Var2,
				 Group = Var3) %>% 
	mutate(N = sprintf('%04d', nns[N])) -> dta_best_lambdas



ggplot(dta_best_lambdas, aes(y=value, x=as.character(N))) +
	geom_boxplot() +
	facet_wrap(~Group) +
	theme_bw()



for (ii in 1:100){
	
	
	
	print(ii)
	fname <- sprintf( 'cv_avt3g_boot_20180524_%04d.rdata', ii)
	load(fname)

	best_lambdasss[ii,] <- best_lambda_csv
	
	if (length(best_lambda_csv) == 1){
		best_lambda_grp2 <- best_lambda_csv
		best_lambda_grp3 <- best_lambda_csv
	} else {
		best_lambda_grp2 <- best_lambda_csv[1]
		best_lambda_grp3 <- best_lambda_csv[2]
	}
	
	
	bl_idx_grp2 <- which(est_lasso_avt_best_lambda$lambdas == best_lambda_grp2)
	bl_idx_grp3 <- which(est_lasso_avt_best_lambda$lambdas == best_lambda_grp3)

	if (length(bl_idx_grp2) != 0){
		pt2 <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 2, which_lambda = bl_idx_grp2)
	} else {
		# implies thatlambda = 0 was the best. All parameters not 0.
		pt2_non_zero$Counts <- pt2_non_zero$Counts + 1
	}
	
	
	if (length(bl_idx_grp3) != 0){
		pt3 <- get_nice_par_table(est_lasso_avt_best_lambda, which_group = 3, which_lambda = bl_idx_grp3)
	} else {
		# implies thatlambda = 0 was the best. All parameters not 0.
		pt3_non_zero$Counts <- pt3_non_zero$Counts + 1
	}
	
	if (ii == 1){
		pt2_non_zero <- pt2
		pt2_non_zero$Counts <- 0
		
		pt3_non_zero <- pt3
		pt3_non_zero$Counts <- 0
	}
	
	if (length(bl_idx_grp2) != 0){
		pt2_non_zero$Counts <- pt2_non_zero$Counts + as.numeric(pt2$Estimate != 0)
	}
	
	if (length(bl_idx_grp3) != 0){
		pt3_non_zero$Counts <- pt3_non_zero$Counts + as.numeric(pt3$Estimate != 0)
	}

}




mean(best_lambdasss[,1])
median(best_lambdasss[,1])

mean(best_lambdasss[,2])
median(best_lambdasss[,2])

quantile(best_lambdasss[,1], probs = c(0.25, 0.75))
quantile(best_lambdasss[,2], probs = c(0.25, 0.75))


hist(best_lambdasss[,1])
hist(best_lambdasss[,2])



pt2_non_zero$Prop <- pt2_non_zero$Counts / 100
pt3_non_zero$Prop <- pt3_non_zero$Counts / 100


pt2_non_zero # AA
pt3_non_zero # HC



