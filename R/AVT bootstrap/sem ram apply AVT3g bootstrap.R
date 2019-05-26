


library(lavaan)
library(igraph)
library(Matrix)
library(dplyr)
library(reshape2)
library(ggplot2)


source("..\\sem ram func.R")

list.files()


random_seeds <- as.numeric(readLines('random_seed_bootstrap.txt'))

###################################################
# Bootstrap run parameters
####################################################


nboot <- 100

mydatestamp <- '20180524'


###################################################
### Affect Valuation Theory model
##################################################
# http://web.stanford.edu/class/psych253/section/section_8/section8.html#question-c-group-comparisons
# http://journals.sagepub.com/doi/pdf/10.1111/j.1745-6916.2007.00043.x

dta <- read.csv("..\\..\\data\\jt-data1.csv")

cc_idx <- complete.cases(dta[,c('ideahap', 'cultatt',  'temperatt', 'actuhap', 'rigoract', 'depress')])

dta <- dta[cc_idx,]

grplab1 <- c("EA", "AA", "CH") # European American, Asian American, Hong Kong Chinese.
dta$group <- grplab1[dta$group]

dta %>% 
	group_by(group) %>% 
	summarise(N=n())

# SCALE AND CENTERERR!!!!!!
dta %>%  group_by(group) %>% 
	mutate_all(funs(. - mean(.))) %>% 
	mutate_all(funs(. / sd(.))) %>% 
	ungroup() -> dta_scaled 


## Specify model

# Assymetric relations.
graph_df_avt <- data.frame(from=c('cultatt', 'temperatt', 'cultatt', 'temperatt', 'ideahap', 'actuhap', 'cultatt', 'temperatt', 'ideahap', 'actuhap', 'cultatt', 'temperatt'),
															to=c('ideahap', 'ideahap', 'actuhap', 'actuhap', 'rigoract', 'rigoract', 'rigoract', 'rigoract', 'depress', 'depress', 'depress', 'depress'))

# covariances
graph_df_avt <- rbind(graph_df_avt,
												 data.frame(from=c('ideahap', 'actuhap', 'cultatt', 'temperatt'),
					 										to=c('actuhap', 'ideahap', 'temperatt', 'cultatt')))


model_igraph_avt <- graph_from_data_frame(graph_df_avt)



ram_matrices_avt <- make_ram_matrices(model_igraph_avt)


################## BOOTSTRAP ####################
boot_count <- 1

seed_count <- 0

while (boot_count <= nboot){
  
  
	cat(sprintf('#### Bootsrap round %d of %d ####\n', boot_count, nboot))
	
	
  seed_count <- seed_count + 1
  current_seed <- random_seeds[seed_count]
  set.seed(current_seed)
  write(current_seed, file='used_seeds.txt', append=TRUE)
  
  cat(sprintf('#### Using seed number %d: %d ####\n', seed_count, current_seed))
  
  
	# resample the data set.
	sample_idx <- c(base::sample(which(dta_scaled$group == 'EA'), size=sum(dta_scaled$group == 'EA'), replace=TRUE),
									base::sample(which(dta_scaled$group == 'AA'), size=sum(dta_scaled$group == 'AA'), replace=TRUE),
									base::sample(which(dta_scaled$group == 'CH'), size=sum(dta_scaled$group == 'CH'), replace=TRUE))
	
	dta_scaled_boot <- dta_scaled[sample_idx,]
	stopifnot(all(dim(dta_scaled_boot) == dim(dta_scaled)))
	
	
	# Cross validate
	cv_res_avt <- tryCatch({cv_ram_group_lasso(ram_matrices_avt, dta_scaled_boot,
																	 group_var = 'group', n_cvs=5, n_lambdas_cv = 20)},
												error = function(e) {print(e)
												  return(0)}					 
	)
	
	if (is.numeric(cv_res_avt)){
	  print('Error in cv_ram_group_lasso2. Trying a new bootstrap sampling...')
	  next
	}
	
	# fit optimal
	# average, then find the one with least average error.
	best_lambda_idxx <- apply(apply(cv_res_avt$lambda_fits_cv, MARGIN = c(1,2), mean)[,-1], MARGIN=2, which.min)
	best_lambda_csv <- unique(cv_res_avt$lambdas_cv[best_lambda_idxx])
	print(sprintf('best_lambda_csv: %0.3f', best_lambda_csv))
	
	est_lasso_avt_best_lambda <- estimate_ram_group_lasso(ram_matrices_avt,
																												inndata = dta_scaled_boot,
																												group_var = 'group',
																												lambda_values = best_lambda_csv,
																												verbose=FALSE)
	
	
	# save 
	fname <- sprintf('cv_avt3g_boot_%s_%0.4d.rdata', mydatestamp, boot_count)
	
	save(list=c('cv_res_avt', 'best_lambda_csv', 'est_lasso_avt_best_lambda', 'dta_scaled_boot', 'current_seed'), file=fname)
	
	cat(sprintf('-- Results saved as %s \n', fname))
	
	boot_count <- boot_count + 1
	
}












