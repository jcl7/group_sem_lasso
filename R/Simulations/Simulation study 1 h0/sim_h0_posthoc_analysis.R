

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(reshape2)
library(lavaan)
library(xtable)

source('..\\..\\sem ram func.R')



#### Load true parameters
load('sem sim avt h0 TRUE MODEL.rdata')

parTable(avt_res) %>% 
	select(lhs, op, rhs, est) %>% 
	filter(op == '~') %>% 
	mutate(Parameter = paste(lhs, rhs, sep = ' ~ ')) -> true_parameters




##### Load simulation results.

nns <- c(142, 300, 500, 1000)

all_param_ests <- data.frame()


for (nn in 1:length(nns)){
	for (ii in 1:100){
		
		print(ii)
		fname <- sprintf('cv_avt_sim_20180625_n%d_%04d.rdata', nns[nn], ii)
		load(fname)
		#print(best_cv_lamba)
		
		get_nice_par_table(est_lasso_avt_best_lambda, which_group=2) %>% 
			mutate(N = nns[nn],
						 simno = ii) %>% 
			bind_rows(all_param_ests, .) -> all_param_ests
		
	}
}




all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ ')) %>% 
	left_join(true_parameters, by='Parameter') %>% 
	rename(group1 = est) %>% 
	group_by(Parameter, N) %>% 
	summarise(group1=mean(group1),
						Prop0 = sum(Estimate == 0) / 100,
						Mean_bias = mean(Estimate),
						MSE = mean(Estimate^2)) %>% 
	ungroup() %>% 
	select(Parameter, group1, N, Prop0, Mean_bias, MSE) -> params_summary


varnames <- c('actuhap' = 'Actual Affect', 'depress' = 'Depression',
							'ideahap' = 'Ideal Affect', 'rigoract'='Rigorous Activties',
							'cultatt'='Attitudes', 'temperatt'='Temparement')

column_order <- c("Response", "Predictor", "group1", 
									"true_diff", "142_Prop0", 
									"142_MSE", "142_Mean_bias", "300_Prop0","300_MSE", "300_Mean_bias",
									"500_Prop0","500_MSE" , "500_Mean_bias","1000_Prop0", "1000_MSE", "1000_Mean_bias")


params_summary %>% 
	melt(id.vars=c('Parameter', 'group1', 'N')) %>% 
	dcast(Parameter + group1 ~ N + variable) %>% 
	separate(Parameter, into=c('Response', 'Predictor')) %>% 
	mutate(true_diff = 0) %>% 
	select(column_order) %>% 
	mutate(Response = varnames[Response],
				 Predictor = varnames[Predictor]) %>% 
	arrange(Predictor) %>% 
	xtable()




all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ ')) %>% 
	group_by(Parameter, N) %>% 
	summarise(Prop = sum(Estimate == 0) / 100) %>% 
	arrange(N) %>% 
	dcast(Parameter ~ N) -> params_not0


all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ ')) %>% 
	group_by(Parameter, N) %>% 
	summarise(Prop = sum(Estimate == 0) / 100) %>% 
	summarise(Min = min(Prop),
						Max = max(Prop))


all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ ')) %>% 
	group_by(Parameter, N) %>% 
	summarise(Prop0 = sum(Estimate == 0) / 100,
						Mean = mean(Estimate),
						MSE = mean(Estimate^2),
						Median = median(Estimate),
						SD = sd(Estimate),
						lwr_95 = quantile(Estimate, 0.025),
						upr_95 = quantile(Estimate, 0.975)) -> params_summary


params_summary %>% 
	mutate(N2 = paste('N = ', N, sep='')) %>% 
	dcast(Parameter ~ N2, value.var = 'MSE')


params_summary %>% 
	mutate(N2 = paste('N = ', N, sep='')) %>% 
	dcast(Parameter ~ N2, value.var = 'Prop0')






	
all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ ')) %>% 
	ggplot(aes(x=N, y=Estimate)) +
	geom_point() +
	#geom_smooth(method = 'lm') +
	facet_wrap( ~ Parameter) +
	theme_bw()
	




all_param_ests %>% 
	mutate(Parameter = paste(Response, Predictor, sep = ' ~ '),
				 N_ = as.factor(N/2)) %>% 
	ggplot(aes(x=N_, y=Estimate)) +
	geom_violin() +
	facet_wrap( ~ Parameter) +
	theme_bw() +
	xlab('Sample size') +
	ylab(expression(delta ~ 'estimate'))
	













#######################################################################
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
















