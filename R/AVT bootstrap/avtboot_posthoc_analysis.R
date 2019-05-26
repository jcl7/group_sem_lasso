


source("..\\sem ram func.R")


best_lambdasss <- matrix(ncol=2, nrow=100)


for (ii in 1:100){
	
	
	
	print(ii)
	fname <- sprintf( 'Result\\cv_avt3g_boot_20180524_%04d.rdata', ii)
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



