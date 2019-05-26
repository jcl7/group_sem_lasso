
library(lavaan)
library(igraph)
library(Matrix)
library(dplyr)
library(reshape2)



source('..\\..\\sem ram func.R')

############################
# Simulation parameters
###########################

mydatestamp <- '20180625'

sim_nobs <- 150*2 # total number of observations
nsim <- 100


seed_file <- 'sim_seeds_02.txt' # todo
used_seed_file <- paste(seed_file, '_used.txt', sep='')



###################################################
### Affect Valuation Theory model
##################################################
# http://web.stanford.edu/class/psych253/section/section_8/section8.html#question-c-group-comparisons
# http://journals.sagepub.com/doi/pdf/10.1111/j.1745-6916.2007.00043.x

dta <- read.csv("..\\..\\..\\data\\jt-data1.csv")

cc_idx <- complete.cases(dta[,c('ideahap', 'cultatt',  'temperatt', 'actuhap', 'rigoract', 'depress')])

dta <- dta[cc_idx,]

grplab1 <- c("EA", "AA", "CH") # European American, Asian American, Hong Kong Chinese.
dta$group <- grplab1[dta$group]

table(dta$group)

dta %>% 
	group_by(group) %>% 
	summarise(N=n())

# SCALE AND CENTERERR!!!!!!
dta %>%  group_by(group) %>% 
	mutate_all(funs(. - mean(.))) %>% 
	mutate_all(funs(. / sd(.))) %>% 
	ungroup() -> dta_scaled 



###########################################
##### FIT model in lavaan,, without any grouping.
############################################

avt.model1 <- '
# regressions 
ideahap ~ cultatt + temperatt
actuhap ~ cultatt + temperatt
rigoract ~ ideahap + actuhap + cultatt + temperatt
depress ~ ideahap + actuhap + cultatt + temperatt

# variances and covariances
rigoract ~~ 0*depress
ideahap ~~ actuhap
cultatt ~~ temperatt
'


avt_res <- sem(avt.model1, fixed.x = FALSE, data = dta_scaled)
summary(avt_res)


# We can simulate data using the simulateData function from lavaan. 
# This can tak model syntax or a parameter table as input. 
# Model syntax with fixed parameters have their values stored in
# the lavaanified model data.frame variable ustart (and with free=0).
# We therefore need to set the estimates from the parameter table from the 
# fitted model to the ustart column. 

mypar <- parTable(avt_res)
mypar$ustart <- mypar$est

# # Estimates from the model fitted on real data and
# # the model fitted on simulated data are in agreement.
# set.seed(23942)
# sdd <- simulateData(mypar, fixed.x = FALSE)
# sim_res <- sem(avt.model1, fixed.x = FALSE, data=sdd)
# 
# plot(parTable(sim_res)$est, parTable(avt_res)$est, asp=1)
# abline(coef=c(0,1))

##########################################
# Specify model
########################################3

# Assymetric relations.
graph_df_avt <- data.frame(from=c('cultatt', 'temperatt', 'cultatt', 'temperatt', 'ideahap', 'actuhap', 'cultatt', 'temperatt', 'ideahap', 'actuhap', 'cultatt', 'temperatt'),
													 to=c('ideahap', 'ideahap', 'actuhap', 'actuhap', 'rigoract', 'rigoract', 'rigoract', 'rigoract', 'depress', 'depress', 'depress', 'depress'))

# covariances
graph_df_avt <- rbind(graph_df_avt,
											data.frame(from=c('ideahap', 'actuhap', 'cultatt', 'temperatt'),
																 to=c('actuhap', 'ideahap', 'temperatt', 'cultatt')))


model_igraph_avt <- graph_from_data_frame(graph_df_avt)
ram_matrices_avt <- make_ram_matrices(model_igraph_avt)




random_seeds <- as.numeric(readLines(seed_file))


#################################
# Commence simulation
################################


sim_count <- 1
seed_count <- 0

while (sim_count <= nsim){

	cat(sprintf('#### Simulation round %d of %d ####\n', sim_count, nsim))
	cat(sprintf('#### N = %d ####\n', sim_nobs))
	
	# set seed.
	seed_count <- seed_count + 1
	current_seed <- random_seeds[seed_count]
	set.seed(current_seed) 
	write(c(sim_count, current_seed), file=used_seed_file, append=TRUE, sep=',')
	
	
	## simulate two balanced groups that are identical
	simdta <- simulateData(mypar, fixed.x = FALSE, sample.nobs = sim_nobs)
	simdta$group <- sample(c(rep('grp1', sim_nobs/2), rep('grp2', sim_nobs/2)), replace=FALSE)
	
	# Scale and center the simulated data
	simdta %>%  group_by(group) %>% 
		mutate_all(funs(. - mean(.))) %>% 
		mutate_all(funs(. / sd(.))) %>% 
		ungroup() -> simdta_scaled 
	
	
	# Cross validate
	sim_cv_res_avt <- tryCatch({cv_ram_group_lasso(ram_matrices_avt, simdta_scaled, 
																						 group_var = 'group', n_cvs=5,
																						 n_lambdas_cv = 20)},
												 error = function(e) {print(e)
												 	return(0)}					 
	)
	
	
	if (is.numeric(sim_cv_res_avt)){
		print('Error in cv_ram_group_lasso. Simulating a new data set...')
		next
	}
	
	# fit optimal
	# average, then find the one with least average CV error.
	cv_mean_error <- apply(sim_cv_res_avt$lambda_fits_cv[,2,], MARGIN = 1, FUN=mean)
	best_cv_lamba <- sim_cv_res_avt$lambdas_cv[which.min(cv_mean_error)]
	
	
	print(sprintf('best_lambda_csv: %0.3f', best_cv_lamba))
	
	est_lasso_avt_best_lambda <- estimate_ram_group_lasso(ram_matrices_avt,
																												inndata = simdta_scaled,
																												group_var = 'group',
																												lambda_values = best_cv_lamba,
																												verbose=FALSE)
	
	fname <- sprintf('cv_avt_sim_%s_n%d_%0.4d.rdata', mydatestamp, sim_nobs, sim_count)
	
	save(list=c('sim_cv_res_avt', 'best_cv_lamba', 'est_lasso_avt_best_lambda', 'simdta_scaled', 'current_seed'), file=fname)
	
	
	cat(sprintf('-- Results saved as %s \n', fname))
	
	sim_count <- sim_count + 1
	
}






