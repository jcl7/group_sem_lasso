



library(igraph)
library(Matrix)
library(dplyr)
library(numDeriv)

message('Version 0.7 (20190524)\n')


make_ram_matrices <- function(igraph_obj, latent_variables=NULL){
	
	
	## Make two graph objects, one with the asymmetric relations, one with the symmetric.
	# Mutual edges corresponds to covariances between edges.
	# Make undirected graph, keep only mutual edges.
	# Used later to make the symmetric matrix.
	model_igraph_sym <- as.undirected(igraph_obj, mode='mutual')
	
	# Remove mutual edges. Used later to make the asymmetric matrix.
	model_igraph_asym <- delete.edges(igraph_obj, which(is.mutual(igraph_obj)))
	
	# Define variable sets
	observed_variables <- V(model_igraph_asym)$name
	observed_variables <- observed_variables[!observed_variables %in% latent_variables]
	
	all_variables <- c(observed_variables, latent_variables)
	
	#exogenous_variables <- names(which(degree(igraph_obj, mode = 'in') == 0))
	exogenous_variables <- names(which(degree(model_igraph_asym, mode = 'in') == 0))
	observed_exogenous_variables <- exogenous_variables[exogenous_variables %in% observed_variables]
	
	
	# Identify indicator variables for the latent variables.
	indicator_variables <- list()
	if (!is.null(latent_variables)){
		for (ll in 1:length(latent_variables)){

			indvar <- ego(model_igraph_asym, nodes=latent_variables[ll], order=1, mode = "out", mindist = 1)
			indvar1 <- names(indvar[[1]])
			indicator_variables[[latent_variables[ll]]] <- indvar1[indvar1 %in% observed_variables][1]
			
		}
	}
	
	## Asymetric matrix. (directed relationships, regression coefs)
	Ma <- get.adjacency(model_igraph_asym)
	Ma <- t(Ma[all_variables,all_variables])
	
	stopifnot(all(colnames(Ma) == all_variables))
	
	## symmetric matrix (two headed arrows, ie variances and covariances)
	Ms_tmp <- get.adjacency(model_igraph_sym)
	Ms_tmp <- Ms_tmp[all_variables,all_variables] # Sorter rad og kollone i riktig rekkef?lge.

	
	Ms <- Diagonal(n=nrow(Ms_tmp))
	Ms <- Ms + Ms_tmp

	stopifnot(isSymmetric(Ms))
	stopifnot(all(colnames(Ms) == colnames(Ma)))
	stopifnot(all(colnames(Ms) == rownames(Ms)))
	

	## Filter matrix.
	Mf <- Diagonal(n=length(observed_variables))
	rownames(Mf) <- observed_variables
	
	Mf_lat <- sparseMatrix(i=NULL, j=NULL, dims=c(nrow(Mf), length(latent_variables)))
	Mf <- cbind(Mf, Mf_lat)
	colnames(Mf) <- all_variables
	
	## Mean vector (basicly a filter type vector indicating which means to estimate.)
	Mm <- rep(1, length(all_variables))
	names(Mm) <- all_variables
	Mm[names(Mm) %in% latent_variables] <- 0
	
	
	# Output
	out <- list(A=Ma, S=Ms, Filter=Mf, M=Mm,
							all_variables=all_variables, latent_variables=latent_variables,
							observed_variables=observed_variables,
							indicator_variables=indicator_variables,
							exogenous_variables=exogenous_variables,
							observed_exogenous_variables=observed_exogenous_variables)

	return(out)
	
}




obj_func_group_ccd <- function(pars, obs_cov, obs_mean, Amat, Smat, Fmat, Mvec,
													 sym_idx, asym_idx, mean_idx,
													 pars_skeleton, difference_param=FALSE,
													 s_equal = s_equal,
													 all_pars=NULL, which_param=0,
													 debug.=FALSE, groupwise=FALSE){
	
	stopifnot(length(obs_cov) == length(obs_mean))
	
	# if not all parameters should be estimated, update the all_pars vector with 
	# the values supplied in the pars argument (at the positions given in which_param). 
	if (length(which_param) == 1){
		ccd <- which_param > 0 & !is.null(all_pars)
		if (ccd){
			stopifnot(length(which_param) == length(pars))
			all_pars[which_param] <- pars
			pars <- all_pars
		}
	} else if (length(which_param) > 1 & !is.null(all_pars)) {
		# If there's more than 1, but not all, parameter to estimate.
		stopifnot(length(which_param) == length(pars))
		ccd <- FALSE
		all_pars[which_param] <- pars
		pars <- all_pars
	}
	
  if (length(which_param) == 1){
  	if (is.null(all_pars) & which_param != 0){
  		stop('Error: something wrong in the arguments. all_pars should be given')
  	}
  }
	
	# relist params
	pars <- relist(pars, skeleton=pars_skeleton)
	
	# Calculate the objective function (ie the likelihood).
	objectif <- numeric(length(obs_cov)) #0
	for (gg in 1:length(obs_cov)){
		
		# If difference parameterization is used, the parameter values
		# for the non-reference groups are the differences. To calculate 
		# the likelihood the differences need to be added to the parameters in the
		# reference group.
		if (gg > 1 & difference_param){
			Amat[asym_idx] <- pars[[1]][['asym']] + pars[[gg]][['asym']]
			
			if (s_equal) {
				Smat[sym_idx] <- pars[[1]][['sym']] 
			} else {
				Smat[sym_idx] <- pars[[gg]][['sym']] # difference param not applicable for S.
			}
			
			Mvec[mean_idx] <- pars[[1]][['mean']] + pars[[gg]][['mean']]
			
		} else {
			Amat[asym_idx] <- pars[[gg]][['asym']]

			if (s_equal){
				Smat[sym_idx] <- pars[[1]][['sym']]
			} else {
				Smat[sym_idx] <- pars[[gg]][['sym']]
			}
			
			Mvec[mean_idx] <- pars[[gg]][['mean']]
			
		}
		
		# exponentiate the variances to make sure they are postive.
		diag(Smat) <- exp(diag(Smat)) 
		Smat <- forceSymmetric(Smat, uplo = 'L')
		
		ia <- Diagonal(n=ncol(Amat)) - Amat
		ia_inv <- solve(ia, sparse=TRUE)
		fia_inv <- Fmat %*% ia_inv
		
		model_sigma <- fia_inv %*% Smat %*% t(fia_inv)
		model_mu <- fia_inv %*% Mvec 
		
		model_sigma_inv <- solve(model_sigma)
		
		mean_diff <- obs_mean[[gg]] - model_mu
		mean_obj <- t(mean_diff) %*% model_sigma_inv %*% mean_diff
		
		det_ms <- det(model_sigma)
		
		sigma_obj <- log(det_ms) + sum(diag(obs_cov[[gg]] %*% model_sigma_inv)) - log(det(obs_cov[[gg]])) - nrow(Fmat)
		
		#objectif <- objectif + as.numeric(sigma_obj + mean_obj)
		objectif[gg] <- as.numeric(sigma_obj + mean_obj)
	}
	
	if (groupwise){
		return(objectif)
	} else {
		return(sum(objectif))
	}
	
}



# wrapper around obj_grad() to make it work for grouped SEM.

obj_grad_grp_ccd <- function(pars, obs_cov, obs_mean, Amat, Smat, Fmat, Mvec,
												 sym_idx, asym_idx, mean_idx,
												 pars_skeleton, difference_param=FALSE,
												 s_equal = s_equal,
												 all_pars=NULL, which_param=0,
												 debug. = FALSE){
	
	
	if (length(which_param) == 1){
		ccd <- which_param > 0 & !is.null(all_pars)
		if (ccd){
			all_pars[which_param] <- pars
			pars <- all_pars
			
			n_grp <- length(obs_cov)
			n_par <- length(all_pars) / n_grp
			
			# convert the absolute paramter index to the within-group index.
			which_grp <- (floor((which_param-1) / n_par)+1)
			if (which_grp >= 2){
				which_param_new <- which_param - n_par*(which_grp-1)
			} else {
				which_param_new <- which_param # no conversion needed when it's the first group.
			}

		}
	} else if (length(which_param) > 1){
		stopifnot(length(which_param) == length(pars))
		ccd <- FALSE
		all_pars[which_param] <- pars
		pars <- all_pars
		
		n_grp <- length(obs_cov)
		n_par_tot <- length(all_pars)
		n_par <- n_par_tot / n_grp
		
		which_param <- sort(which_param)
		
		which_param_list <- vector('list', n_grp)
		for (ww in 1:length(which_param)){
			which_grp <- floor((which_param[ww]-1) / n_par)+1
			which_param_new <- which_param[ww] - (n_par * (which_grp-1))
	
			which_param_list[[which_grp]] <- c(which_param_list[[which_grp]], which_param_new)
		}

	}

	# relist params
	pars <- relist(pars, pars_skeleton)
	
	#gradient_vector <- numeric(0)
	gradient_list <- vector('list', length(obs_cov)) 
	
	for (gg in 1:length(obs_cov)){ # iterate over groups
		
		if (gg > 1 & difference_param){
			pars[[gg]][['asym']] <- pars[[1]][['asym']] + pars[[gg]][['asym']]
			
			pars[[gg]][['mean']] <- pars[[1]][['mean']] + pars[[gg]][['mean']]
		}
		
		if (s_equal){
			pars[[gg]][['sym']] <- pars[[1]][['sym']]
		}
		
		parsvec <- unlist(pars[[gg]])
		
		if (ccd){
			
			# Calculate the gradient for a single parameter.
			if (gg != which_grp){next}
			
			objg <- obj_grad(pars=parsvec, obs_cov[[gg]], obs_mean[[gg]],
											 Amat, Smat, Fmat, Mvec, sym_idx, asym_idx, mean_idx,
											 pars_skeleton=pars[[gg]], which_param=round(which_param_new))

			return(objg)
			
		} else if (length(which_param) > 1){
			# Calculate the gradient for just a subset (within a group) of the parameters
			# Used when updating auxillary parameters (more than 1).
			
			if (is.null(which_param_list[[gg]])){next}
			objg <- obj_grad(pars=parsvec, obs_cov[[gg]], obs_mean[[gg]],
											 Amat, Smat, Fmat, Mvec, sym_idx, asym_idx, mean_idx,
											 pars_skeleton=pars[[gg]], which_param=which_param_list[[gg]])

			
		}	else {
			# Calculate the gradient for all parameters
			objg <- obj_grad(pars=parsvec, obs_cov[[gg]], obs_mean[[gg]],
											 Amat, Smat, Fmat, Mvec, sym_idx, asym_idx, mean_idx,
											 pars_skeleton=pars[[gg]])
		}
		
		# Add the gradients to the gradient_list
		grad_lst_tmp <- relist(objg, pars_skeleton[[1]])
		
		if (gg > 1 & s_equal){
			gradient_list[[1]][['sym']] <- gradient_list[[1]][['sym']] + grad_lst_tmp[['sym']]
			grad_lst_tmp[['sym']] <- numeric(0)
		}
		
		
		gradient_list[[gg]] <- grad_lst_tmp
		
	}

	# make gradient_list a vector before returning.
	gradient_vector <- unlist(gradient_list)

	# When doing CCD and only a subset of parameters, there will be NA's in 
	# the gradient vector. Remove those.
	gradient_vector <- gradient_vector[!is.na(gradient_vector)]

	
	return(gradient_vector)
	
}


# The gradient for a single group.
obj_grad <- function(pars, obs_cov, obs_mean, Amat, Smat, Fmat, Mvec,
										 sym_idx, asym_idx, mean_idx, pars_skeleton,
										 browse_at=0, which_param=0){
	
	gradient_vector <- numeric(length(pars))
	
	# relist params
	pars_list <- relist(pars, pars_skeleton)
	
	Amat[asym_idx] <- pars_list$asym
	Smat[sym_idx] <- pars_list$sym
	Mvec[mean_idx] <- pars_list$mean
	
	asym_idxm <- arrayInd(asym_idx, .dim=dim(Amat))
	sym_idxm <- arrayInd(sym_idx, .dim=dim(Smat))

	diag(Smat) <- exp(diag(Smat)) # exponentiate the variances to always make them postive.
	Smat <- forceSymmetric(Smat, uplo = 'L')
	
	# Precompute B, C and E, model sigma etc.
	ia <- Diagonal(n=nrow(Amat)) - Amat
	ia_inv <- solve(ia, sparse=TRUE) ## B matrix (eq 7)
	fia_inv <- Fmat %*% ia_inv
	
	model_sigma <- fia_inv %*% Smat %*% t(fia_inv)
	model_sigma_inv <- solve(model_sigma)
	model_mu <- fia_inv %*% Mvec 
	
	CC <- Diagonal(n=nrow(model_sigma)) - model_sigma_inv %*% obs_cov # eq 8
	
	EE <- ia_inv %*% Smat %*% t(ia_inv) # eq 9
	
	b_ <- obs_mean - model_mu # eq 10.
	
	# Iterate over all parameters.
	pars_new_matrix_idx <- cumsum(sapply(pars_list, length))
	
	for (ii in 1:length(pars)){
		
		# Check if gradient is supposed to be calculated for the current parameter.
		if (length(which_param) == 1){
			if (which_param > 0 & ii != which_param){
				next
			}
		} else if (length(which_param) > 1){
			if (!(ii %in% which_param)){
				next
			}
		}
		

		# Keep track of which matrix we are in.
		which_matrix <- which.min(pars_new_matrix_idx < ii) # 1 (A), 2 (S) or 3 (M)
		

		if (ii == browse_at){browser()}
		
		# Set two of the gradient matrices for A, B and M to 0. 
		if (which_matrix == 1){
			# We are in A. Set the gradient of S and M to 0.
			Sgrad <- sparseMatrix(i = 1, j = 1,
														dims=dim(Smat), x=0)
			Mgrad <- rep(0, length(Mvec))
			
			idxm_idx <- ii
			Agrad <- sparseMatrix(i = asym_idxm[idxm_idx,1], j = asym_idxm[idxm_idx,2],
														dims=dim(Amat), x=1)
			
		} else if (which_matrix == 2) {
			# We are in S. Set the gradient of A and M to 0.
			Agrad <- sparseMatrix(i = 1, j = 1,
														dims=dim(Amat), x=0)
			Mgrad <- rep(0, length(Mvec))
			
			idxm_idx <- ii - pars_new_matrix_idx[which_matrix-1]
			Sgrad <- sparseMatrix(i = sym_idxm[idxm_idx,1], j = sym_idxm[idxm_idx,2],
														dims=dim(Smat), x=1)
			
		} else if (which_matrix == 3){
			# We are in M. Set the gradient of A and S to 0.
			Agrad <- sparseMatrix(i = 1, j = 1,
														dims=dim(Amat), x=0)
			Sgrad <- sparseMatrix(i = 1, j = 1,
														dims=dim(Smat), x=0)
			
			idxm_idx <- ii - pars_new_matrix_idx[which_matrix-1]
			Mgrad <- rep(0, length(Mvec))
			Mgrad[idxm_idx] <- 1
		}
		
		# sigma gradient (eq 15)
		sigma_grad <- symop(fia_inv %*% Agrad %*% EE %*% t(Fmat)) + (fia_inv %*% Sgrad %*% t(fia_inv))
		#mu_grad <- tryCatch({(fia_inv %*% Agrad %*% ia_inv %*% Mvec) + fia_inv %*% Mgrad}, error=browser())
		mu_grad <- (fia_inv %*% Agrad %*% ia_inv %*% Mvec)
		#if (ii == 28){browser()}
		mu_grad <- mu_grad + (fia_inv %*% Mgrad)
		
		term1 <- sum(diag(model_sigma_inv %*% sigma_grad %*% CC)) 
		term2 <- ((t(b_) %*% model_sigma_inv %*% sigma_grad) + 2*t(mu_grad)) %*% model_sigma_inv %*% b_
		gradient_vector[ii] <- term1 - as.numeric(term2)
		
	}
	
	if (length(which_param) == 1){
	  if (which_param != 0){
	    return(gradient_vector[which_param])
	  } else {
	    return(gradient_vector)
	  }
	} else {
	  return(gradient_vector[which_param])
	}
	
}



tr <- function(MM){
	sum(diag(MM))
}

symop <- function(MM){
	MM + t(MM)
}

obj_grad_numeric <- function(par, ...){
	grad(func = obj_func_group, x=par, ..., method='simple')
}


soft_treshold <- function(xx, lambda=0){
	stopifnot(length(xx)==1)
	sign(xx) * max(0, (abs(xx) - lambda))
}


# Small function that gives the index of the diagonal elements of a sqaure matrix
diag_idx <- function(matr){
	stopifnot(nrow(matr) == ncol(matr))
	which(diag(nrow=nrow(matr), ncol=ncol(matr)) == 1)
}



estimate_ram_group <- function(ram_matrices, inndata, group_var, s_equal=FALSE,
															 difference_param=FALSE, robust_covar=FALSE){
	

	stopifnot(ram_matrices$observed_variables %in% colnames(inndata))
	stopifnot(!ram_matrices$latent_variables %in% colnames(inndata))
	stopifnot(group_var %in% colnames(inndata))
		
	group_labels <- as.character(unique(inndata[[group_var]]))
	n_groups <- length(group_labels)
	
	# Calculate the observed mean vector and covariance matrix .
	obs_cov <- list()
	obs_mean <- list()
	for (ii in 1:n_groups){
		grp_idx <- inndata[[group_var]] == group_labels[ii]
		
		if (robust_covar){
			obs_cov[[group_labels[ii]]] <- robust_cov(inndata[grp_idx,ram_matrices$observed_variables])
		} else {
			obs_cov[[group_labels[ii]]] <- cov(inndata[grp_idx,ram_matrices$observed_variables])
			
		}
		
		obs_mean[[group_labels[ii]]] <- colMeans(inndata[grp_idx,ram_matrices$observed_variables])
	}


	# first indicators, fixed to 1. row and column in the A matrix.
	if (!is.null(ram_matrices$latent_variables)){
		first_inidc_idx <- matrix(ncol=2, nrow=length(ram_matrices$latent_variables))
		for (ll in 1:length(ram_matrices$latent_variables)){
			first_indic <- ram_matrices$indicator_variables[[ram_matrices$latent_variables[ll]]][1]
			first_inidc_idx[ll,1] <- which(first_indic == ram_matrices$all_variables)
			first_inidc_idx[ll,2] <- which(ram_matrices$latent_variables[ll] == ram_matrices$all_variables)
		}
		
		# set the latent indicators in A to 0, (to be set to 1 later.)
		for (ii in 1:nrow(first_inidc_idx)){
			rr <- first_inidc_idx[ii,1]
			cc <- first_inidc_idx[ii,2]
			ram_matrices$A[rr,cc] <- 0
		}
		
		
	}

	# Get indices for the parameters.
	asym_idx <- which(ram_matrices$A != 0)
	mean_idx <- which(ram_matrices$M != 0)
	sym_idx <- which(lower.tri(ram_matrices$S, diag = TRUE) & ram_matrices$S != 0)
	
	# set the latent indicators in A to 1, again
	if (!is.null(ram_matrices$latent_variables)){
		for (ii in 1:nrow(first_inidc_idx)){
			rr <- first_inidc_idx[ii,1]
			cc <- first_inidc_idx[ii,2]
			ram_matrices$A[rr,cc] <- 1
		}
	}
	
	param_list <- list()
	# estimates_Amat <- list() #??
	
	for (gg in 1:n_groups){
		
		if (gg > 1 & s_equal){
			# The sym-parameters are only specified for the first group.
			param_list[[gg]] <- list(asym=rep(0.01, length(asym_idx)),
															 sym=numeric(0),
															 mean = obs_mean[[gg]])
		} else {
			param_list[[gg]] <- list(asym=rep(0.01, length(asym_idx)),
															 sym=rep(0.01, length(sym_idx)),
															 mean = obs_mean[[gg]])
		}
	}
	
	params <- unlist(param_list)
	
	optim_ctrl <- list(maxit=500)
	
	optim_res <- optim(par=params, fn=obj_func_group_ccd, gr=obj_grad_grp_ccd,
										 obs_cov = obs_cov,
										 obs_mean = obs_mean,
										 Amat = ram_matrices$A,
										 Smat= ram_matrices$S,
										 Fmat = ram_matrices$Filter,
										 Mvec = ram_matrices$M,
										 sym_idx = sym_idx, asym_idx = asym_idx, mean_idx = mean_idx,
										 difference_param = difference_param,
										 s_equal = s_equal,
										 pars_skeleton=param_list, method='BFGS',
										 control=optim_ctrl)
	
	if (optim_res$convergence != 0){
		warning('Did not converge !')
	}

	
	res_params <- relist(optim_res$par, param_list)
	
	# Make RAM matrices filled with paramter estimates.
	ram_matrices_est <- vector('list', n_groups)
	for (gg in 1:n_groups){
		ram_matrices_est[[gg]] <- ram_matrices[c('A', 'S', 'M')]
		ram_matrices_est[[gg]][['A']][asym_idx] <- res_params[[gg]][['asym']]
		ram_matrices_est[[gg]][['S']][sym_idx] <- res_params[[gg]][['sym']]
		ram_matrices_est[[gg]][['M']][mean_idx] <- res_params[[gg]][['mean']]
		
		diag(ram_matrices_est[[gg]][['S']]) <- exp(diag(ram_matrices_est[[gg]][['S']]))
		ram_matrices_est[[gg]][['S']] <- forceSymmetric(ram_matrices_est[[gg]][['S']], uplo = 'L')
	}
	
	# calculate the model implied covariance matrices.
	model_sigmas <- vector('list', n_groups)
	for (gg in 1:n_groups){
		
		Smat <- ram_matrices_est[[gg]][['S']]
		if (difference_param & gg > 1){
			Amat <- ram_matrices_est[[1]][['A']] + ram_matrices_est[[gg]][['A']]
		} else {
			Amat <- ram_matrices_est[[gg]][['A']]
		}

		ia <- Diagonal(n=ncol(Amat)) - Amat
		ia_inv <- solve(ia, sparse=TRUE)
		fia_inv <- ram_matrices$Filter %*% ia_inv
		
		model_sigmas[[gg]] <- fia_inv %*% Smat %*% t(fia_inv)
	}
	
	
	out <- list(res_params=res_params, optim_res=optim_res,
							asym_idx=asym_idx, sym_idx=sym_idx, mean_idx=mean_idx,
							ram_matrices_est=ram_matrices_est, ram_matrices=ram_matrices,
							obs_cov=obs_cov, obs_mean=obs_mean, Filter=ram_matrices$Filter,
							model_sigmas=model_sigmas,
							difference_param=difference_param)
	
	class(out) <- 'sem_jcl'
	
	return(out)
	
}



estimate_ram_group_lasso <- function(ram_matrices, inndata, group_var,
																		 n_lambdas = 10, verbose=TRUE,
																		 lambda_values = NULL, refit_var=TRUE,
																		 refit_cov=TRUE, refit_means=TRUE,
																		 s_equal=FALSE, robust_covar = FALSE, 
																		 check_model_diagonal=FALSE){
	
	stopifnot(n_lambdas >= 0)
	n_lambdas <- ceiling(n_lambdas) # make integer.
	 
	
	# lambda_values can be used to manually specify which values of lambda to use.
	# By default this is NULL. In that case the lambda values are inferred automatically.
	if (!is.null(lambda_values)){
		stopifnot(all(lambda_values >= 0))
		
		
		if (any(lambda_values == 0.0)){
			# drop any lambda values == 0. This are estimated anyway. 
			lambda_values <- lambda_values[!lambda_values == 0.0]
		}
		
		n_lambdas <- length(lambda_values)
		
	} 
	
	
	
	
	# Give a helpfull message if s_equal = TRUE.
	if (s_equal & (refit_cov | refit_var)){
		refit_cov <- FALSE
		refit_var <- FALSE
		change_refit_cov_message <- 'Message: Since s_equal = TRUE, refit_cov and refit_var are set to be FALSE.\n'
		message(change_refit_cov_message)
		cat(change_refit_cov_message)
	}
	
	cat('Commence estimation of Structural Equation Model (SEM)\n')
	cat('Initial estimation\n')
	
	# initial estimation
	init_start_time <- Sys.time()
	sem_res_init <- estimate_ram_group(ram_matrices, inndata,
																		 group_var = group_var,
																		 difference_param = TRUE,
																		 s_equal=s_equal,
																		 robust_covar = robust_covar)
	
	init_end_time <- Sys.time()
	init_time <- as.numeric(init_end_time - init_start_time, units='secs')
	
	if (sem_res_init$optim_res$convergence != 0){
		cat('Intitial model estimation did not converge!!!')
	} else {
		cat(sprintf('Initial estimation successfully completed in %.2f seconds. \n', init_time))
	}
	
	## A quick test if model fit is not completely crazy.
	# Since the data is standardized, the covariance matric is actually a 
	# correlation matrix, implying diagonal elements = 1. This should be the case for 
	# the model implied covariance as well.
	if (check_model_diagonal){
		model_diagonal_ok <- all(abs(sapply(sem_res_init$model_sigmas, diag) - 1) < 0.01)
		if (!model_diagonal_ok){
			cat('Model sigma diagonals not OK!!!!!!!')
			stop()
		}
	}
	
	
	# if all user sumbitted lambda values are 0, do not fit.
	if (length(lambda_values) == 0){
		
		out <- list(lambdas=0.0, sem_res_init=sem_res_init, ram_matrices=ram_matrices,
								refit_var=refit_var, refit_cov=refit_cov, refit_means=refit_means)
		
		class(out) <- 'sem_lasso_jcl'
		
		warning('User supplied lamda = 0. Retunr only initial fit (equal to lambda=0).')
		return(out)
		
	}
		
	
	
	
	# Infer hte number of groups based on the initial estimation.
	ngroups <- length(sem_res_init$obs_cov)
	
	if (s_equal){
		for (gg in 1:ngroups){
			sem_res_init$res_params[[gg]][['sym']] <- sem_res_init$res_params[[1]][['sym']]
		}
	}
	
	# Number of parameters in each group
	n_pars <- length(unlist(sem_res_init$res_params[[1]]))
	
	params_init <- sem_res_init$res_params
	
	# set inital parameters to 0
	for (gg in 2:ngroups){
		params_init[[gg]]$asym <- rep(0, length(params_init[[gg]]$asym)) 
	}
	
	# Indexes for which parameters to estimate.
	
	# Index of which parameters to lasso (ie the asymetric parameters).
	which_asym <- grep('asym', names(unlist(params_init[[2]])))
	to_do <- rep(which_asym, times=ngroups-1) + rep((1:(ngroups-1))*n_pars, each=length(which_asym))
	#to_do <- n_pars + (which_asym * ((2:ngroups)-1))
	#browser()
	
	
	# Index of auxillary parameters, ie means and (co)variances.
	which_sym <- grep('^sym', names(unlist(params_init[[2]])))
	to_do_sym_idx <- rep(which_sym, times=ngroups-1) + rep((1:(ngroups-1))*n_pars, each=length(which_sym))
	# to_do_sym_idx <- n_pars + grep('(^sym)', names(unlist(params_init[[2]])))
	
	which_mean <-  grep('(^mean)', names(unlist(params_init[[2]]))) 
	to_do_mean <- rep(which_mean, times=ngroups-1) + rep((1:(ngroups-1))*n_pars, each=length(which_mean))


	to_do_aux <- sort(c(to_do_sym_idx, to_do_mean))
	
	
	## Lasso lambda-parameters.
	
	if (is.null(lambda_values)){
		# Lambda values not given by user. A set of useful lambda values is generated.
		max_lambda <- max(abs(unlist(sem_res_init$res_params[[2]][c('asym', 'mean')]))) * 1.01
		lambdas <- seq(max_lambda, 0+(max_lambda/n_lambdas), length.out = n_lambdas)
		cat('Generating useful lambda values\n')
		cat(sprintf('Number of lambda parameters: %d (range: %.3f - %.3f) \n', n_lambdas, min(lambdas), max(lambdas)))
	} else {
		max_lambda <- max(lambda_values)
		lambdas <- sort(lambda_values, decreasing = TRUE)
		cat('Lambda values specified by user.\n')
		cat(sprintf('Number of lambda parameters: %d (range: %.3f - %.3f) \n', n_lambdas, min(lambdas), max(lambdas)))
	}
	
	# Lists, vectors, to store results.
	param_est_lasso <- list()
	iteration_counts <- numeric(n_lambdas)
	
	params <- unlist(params_init) # make into vector to be used with optim / objective function.
	
	for (ll in 1:n_lambdas){
		
		if (ll > 1){
			params <- param_est_lasso[[ll-1]]
		}
		
		## Update SYM and MEAN parameters.
		if (length(to_do_aux) > 0){
			cat('Updating auxillary paramters (SYM and MEAN) \n')
			
			ores_sym <- optim(par=params[to_do_aux], fn=obj_func_group_ccd, gr=obj_grad_grp_ccd,
												obs_cov = sem_res_init$obs_cov,
												obs_mean = sem_res_init$obs_mean,
												Amat = ram_matrices$A,
												Smat = ram_matrices$S,
												Fmat = ram_matrices$Filter,
												Mvec = ram_matrices$M,
												sym_idx = sem_res_init$sym_idx, asym_idx = sem_res_init$asym_idx, mean_idx = sem_res_init$mean_idx,
												difference_param = TRUE, 
												which_param=to_do_aux, all_pars = params,
												s_equal=s_equal,
												pars_skeleton=params_init, method='BFGS',
												debug.=FALSE)
			
			params[to_do_aux] <- ores_sym$par
			
			if (ores_sym$convergence != 0){
				warning('Sym did not converge')
			}
		}

		
		cat(sprintf('=== LAMBDA: %.3f [%d of %d] === \n', lambdas[ll], ll, n_lambdas))
		
		for (jj in 1:15){ #iterate until convergence. Max iteration = 15.
			cat(sprintf('Iteration: %d \n', jj))
			cur_max_par_change <- 0
			
			for (pp in to_do){  # loop over parameters
				
				if (verbose){
					cat(sprintf('[%d] ', pp))
				}
				
				ores <- optim(par=params[pp], fn=obj_func_group_ccd, gr=obj_grad_grp_ccd,
											obs_cov = sem_res_init$obs_cov,
											obs_mean = sem_res_init$obs_mean,
											Amat = ram_matrices$A,
											Smat = ram_matrices$S,
											Fmat = ram_matrices$Filter,
											Mvec = ram_matrices$M,
											sym_idx = sem_res_init$sym_idx, asym_idx = sem_res_init$asym_idx, mean_idx = sem_res_init$mean_idx,
											difference_param = TRUE, which_param=pp, all_pars = params,
											s_equal=s_equal,
											pars_skeleton=params_init, method='BFGS',
											debug.=TRUE)
				
				
				# Apply penalty / soft tresholding.
				new_par <- ores$par
				new_par <- soft_treshold(new_par, lambda = lambdas[ll])
				
				par_change <- new_par - params[pp]
				cur_max_par_change <- max(cur_max_par_change, abs(par_change))
				
				if (verbose){
					cat(sprintf('New param: %f; prev=%f [change=%f] \n', new_par, params[pp], par_change))
				}
					
				params[pp] <- new_par
			}
			
			# Check for convergence.
			if (cur_max_par_change < 0.001){
				print('converged!')
				iteration_counts[ll] <- jj
				break
			}
			
		}
		
		param_est_lasso[[ll]] <- params
		
	}
	
	out <- list(param_est_lasso=param_est_lasso, lambdas=lambdas,
							sem_res_init=sem_res_init, iteration_counts=iteration_counts, 
							params_init=params_init, to_do=to_do, ram_matrices=ram_matrices,
							refit_var=refit_var, refit_cov=refit_cov, refit_means=refit_means)
	
	class(out) <- 'sem_lasso_jcl'
	
	return(out)
	
}
	

tidy_lasso_trace <- function(lasso_res, grp=2){
	n_lambdas <- length(lasso_res$lambdas)
	param_est_lasso <- lasso_res$param_est_lasso
	lambdas <- lasso_res$lambdas
	to_do <- lasso_res$to_do
	params_init <- lasso_res$params_init
	lambda0_est <- lasso_res$sem_res_init$res_params[[grp]][['asym']]
	
	trace_mat <- matrix(ncol=n_lambdas, nrow=length(to_do))
	#trace_mat <- matrix(ncol=n_lambdas+1, nrow=length(to_do))
	#trace_mat[,1] <- params_init[[grp]][['asym']]
	
	#colnames(trace_mat) <- sprintf('lambda_%.3f', lambdas)
	colnames(trace_mat) <- lambdas
	
	for (ii in 1:n_lambdas){
		trace_mat[,ii] <- relist(param_est_lasso[[ii]], skeleton = params_init)[[grp]][['asym']]
	}

	trace_mat <- cbind(trace_mat, '0'=lambda0_est)
	
	trace_mat <- melt(trace_mat)
	colnames(trace_mat) <- c('Coef', 'Lambda', 'Delta_value')
	trace_mat$Lambda <- as.numeric(trace_mat$Lambda)
	
	return(trace_mat)
	
}

# Make a trace plot of the parameter estimates from a lasso fit.
# Essentially plots parameter values as a function of lambda.
plot_lasso_trace <- function(lasso_res, plot_min_lambda=TRUE, plot_lambda=TRUE, grp=2,
														 xlab='Lambda', ylab='Delta parameter value', xlim=NULL,
														 main=''){
	
	n_lambdas <- length(lasso_res$lambdas)
	param_est_lasso <- lasso_res$param_est_lasso
	lambdas <- lasso_res$lambdas
	to_do <- lasso_res$to_do
	params_init <- lasso_res$params_init
	lambda0_est <- lasso_res$sem_res_init$res_params[[grp]][['asym']]
	
	trace_mat <- matrix(ncol=n_lambdas+1, nrow=length(to_do))
	trace_mat[,1] <- params_init[[grp]][['asym']]
	
	for (ii in 1:n_lambdas){
		trace_mat[,ii+1] <- relist(param_est_lasso[[ii]], skeleton = params_init)[[grp]][['asym']]
	}
	
	trace_mat <- cbind(trace_mat, lambda0_est)
	

	
	lambdas_plot <- c(lambdas[1]*1.05, lambdas, 0)
	
	if (is.null(xlim)){
		xlim <- rev(range(lambdas_plot))
	}
	
	plot(x=lambdas_plot, y=trace_mat[1,],
			 ylim=range(trace_mat), xlim=xlim, type='l',
			 ylab=ylab, xlab=xlab, main=main)
	for (ii in 2:nrow(trace_mat)){
		lines(x=lambdas_plot, y=trace_mat[ii,], type='l')
	}
	
	points(x=rep(0, length(lambda0_est)), y=lambda0_est)
	
	if (plot_lambda){
		par(xpd=FALSE) # https://stat.ethz.ch/pipermail/r-help/2011-April/276737.html
		
		for (ii in 1:length(lasso_res$lambdas)){
			abline(v=lasso_res$lambdas[ii], lty='dotdash', col='grey76')
		}
	}

	
	if (plot_min_lambda){
		abline(v=lambdas[length(lambdas)])
	}
	
	
}




# For a estimated model, how well does this fit a new data set?
# This function can take both a regular model fit, and a lasso fit. 
# lambda_idx: An integer specifying which lambda value to use, in the case of a lasso fit. 
ram_data_fit <- function(fitted_model, new_data, group_var=NULL, lambda_idx=0,
												 robust_covar=FALSE, s_equal=FALSE, groupwise=FALSE){
	
	# groupwise: Vil man ha fit for hver gruppe separat?
	
	if (is.null(group_var)){
		stop('')
	}

	group_labels <- as.character(unique(new_data[[group_var]]))
	n_groups <- length(group_labels)
	
	obs_cov <- list()
	obs_mean <- list()
	for (ii in 1:n_groups){
		grp_idx <- new_data[[group_var]] == group_labels[ii]
		
		if (robust_covar){
			obs_cov[[group_labels[ii]]] <- robust_cov(new_data[grp_idx, fitted_model$ram_matrices$observed_variables])
		} else {
			obs_cov[[group_labels[ii]]] <- cov(new_data[grp_idx, fitted_model$ram_matrices$observed_variables])
		}
		
		obs_mean[[group_labels[ii]]] <- colMeans(new_data[grp_idx, fitted_model$ram_matrices$observed_variables])
	}
	
	if (class(fitted_model) == 'sem_jcl'){
		
		if (lambda_idx != 0){
			warning('fitted_model not a lasso fit. Argument lambda_idx ignored.')
		}

		obj_func_group_ccd(pars = unlist(fitted_model$res_params),  
											 obs_cov = obs_cov,
											 obs_mean = obs_mean,
											 Amat = fitted_model$ram_matrices$A,
											 Smat= fitted_model$ram_matrices$S,
											 Fmat = fitted_model$ram_matrices$Filter,
											 Mvec = fitted_model$ram_matrices$M,
											 sym_idx = fitted_model$sym_idx,
											 asym_idx = fitted_model$asym_idx,
											 mean_idx = fitted_model$mean_idx,
											 pars_skeleton = fitted_model$res_params,
											 difference_param = fitted_model$difference_param,
											 groupwise = groupwise,
											 s_equal=s_equal)
		
	} else if (class(fitted_model) == 'sem_lasso_jcl'){
		
		stopifnot(lambda_idx > 0)
		stopifnot(lambda_idx <= length(fitted_model$lambdas))
		
		obj_func_group_ccd(pars = fitted_model$param_est_lasso[[lambda_idx]],  
											 obs_cov = obs_cov,
											 obs_mean = obs_mean,
											 Amat = fitted_model$sem_res_init$ram_matrices$A,
											 Smat= fitted_model$sem_res_init$ram_matrices$S,
											 Fmat = fitted_model$sem_res_init$ram_matrices$Filter,
											 Mvec = fitted_model$sem_res_init$ram_matrices$M,
											 sym_idx = fitted_model$sem_res_init$sym_idx,
											 asym_idx = fitted_model$sem_res_init$asym_idx,
											 mean_idx = fitted_model$sem_res_init$mean_idx,
											 pars_skeleton = fitted_model$sem_res_init$res_params,
											 difference_param = TRUE, groupwise=groupwise,
											 s_equal=s_equal) # difference_param is always TRUE for lasso.
		
	}
 	
	
}


# Skaff estimaerte RAM matriser fra lasso.
get_estimated_ram_matrices <- function(lasso_fit){
	
	stopifnot(class(lasso_fit) == 'sem_lasso_jcl')
	
	# first layer, lambda
	# second layer, group
	# third layer, matrix (A, S, M)
	
	idxx <- c('asym_idx', 'sym_idx', 'mean_idx')
	
	out <- list()
	
	# Iterate over lambdas.
	for (ll in 1:length(lasso_fit$lambdas)){
		ram_lst <- list()
		
		# par_lst <- relist(lasso_fit$param_est_lasso[[ll]], lasso_fit$params_init)
		
		if (length(lasso_fit$lambdas) == 1 & lasso_fit$lambdas == 0){
			par_lst <- lasso_fit$sem_res_init$res_params
		} else {
			par_lst <- relist(lasso_fit$param_est_lasso[[ll]], lasso_fit$sem_res_init$res_params)	
		}
		
		
		# iterate over groups.
		for (gg in 1:length(par_lst)){
			ram_lst[[gg]] <- list()
			mnames <- names(lasso_fit$sem_res_init$ram_matrices_est[[gg]])
			
			# Iterate over the three matrices.
			for (mm in 1:length(par_lst[[gg]])){
				
				# Get the ram matrix.
				ram_mat_tmp <- lasso_fit$sem_res_init$ram_matrices_est[[gg]][[mm]]
				
				# Fill the matrix with parameter estimates.
				ram_mat_tmp[lasso_fit$sem_res_init[[idxx[mm]]]] <- par_lst[[gg]][[mm]]
				
				# If the symetric matrix.
				if (mm == 2){
					diag(ram_mat_tmp) <- exp(diag(ram_mat_tmp))
					ram_mat_tmp <- forceSymmetric(ram_mat_tmp, uplo = 'L')
				}
				
				ram_lst[[gg]][[mnames[mm]]] <- ram_mat_tmp
				
			}
		}
		out[[ll]] <- ram_lst
	}
	return(out)
}


# Cross validation to tune the optimal lambda.
# n_folds: In how many parts should the data be split during CV? Recomended 2.
# n_cvs: How many times should the cross validation procedure be perfomed?
# n_lambdas_cv = How many values of lambda should be tried? 
# preproc_data: If TRUE, the data is scaled and centered groupwise. (recomended). 
# super_verbose: Should the detialed progress of all cross validation rounds be displayed? Not recomended.
# New version of cv_ram_group_lasso, to handle more than two groups.
cv_ram_group_lasso <- function(ram_matrices, inndata, group_var,
															 n_folds=2, n_cvs = 10, n_lambdas_cv = 15,
															 super_verbose=FALSE,
															 refit_var = TRUE, refit_cov = TRUE, refit_means = FALSE,
															 lambda_max=NULL, robust_covar=FALSE,
															 s_equal=FALSE){
	
	# Add index column to data. Usefull during CV.
	stopifnot(!('IDX' %in% colnames(inndata)))
	inndata %>% 
		mutate(IDX=1:nrow(.)) -> sem_data_cv
	

	# initial full data estimation, to determine reasonable lambda values
	if (is.null(lambda_max)){
		cat('Commence initial full data estimation to determine lambda_max. \n')
		cv_res_init <- estimate_ram_group(ram_matrices, sem_data_cv,
																			group_var = group_var, difference_param = TRUE, 
																			robust_covar=robust_covar,
																			s_equal=s_equal)
		
		if (cv_res_init$optim_res$convergence != 0){
			cat('Initial full data estimation did not converge!!!!! \n')
			cat(sprintf('Message: %s \n', cv_res_init$optim_res$message))
		} 
		
		# Use the full data estimation to determine which lambda values to use.
		max_lambda_cv <- max(abs(cv_res_init$res_params[[2]][['asym']])) * 1.05
		lambdas_cv <- seq(max_lambda_cv, 0+(max_lambda_cv/n_lambdas_cv), length.out = n_lambdas_cv)
		
	} else {
		stopifnot(lambda_max >= 0)
		lambdas_cv <- seq(lambda_max, 0+(lambda_max/n_lambdas_cv), length.out = n_lambdas_cv)
	}
	
	n_grp <- length(cv_res_init$obs_cov)
	
	# Matrices to save the results from CV in. 
	lambda_fits_cv <- array(dim=c(n_lambdas_cv+1, n_grp, n_cvs))
	#lambda_fits_cv <- matrix(ncol=n_lambdas_cv+1, nrow=n_cvs) # fitted values.
	
	cc <- 1 # Counter for CV rounds.
	iter_times <- c(n_cvs) # Timing of the CV rounds. Used to guess the remaining time. 
	
	lasso_fits <- vector(mode = 'list', length = n_cvs)
	cv_data_sets_train <- vector(mode = 'list', length = n_cvs)
	cv_data_sets_test <- vector(mode = 'list', length = n_cvs)
	
	
	cat(sprintf('COMMENCE CROSS VALIDATION', cc, n_cvs))
	
	while (cc <= n_cvs){
		
		cur_iter_start_time <- Sys.time()
		
		cat(sprintf('CV ROUND %d of %d \n', cc, n_cvs))
		
		if (cc > 1){
			exp_time_left <- mean(iter_times[1:(cc-1)]) * (n_cvs - cc + 1)
			cat(sprintf('ESTIMATED TIME REMAINING: %4.f minutes \n', exp_time_left))
			
		}
		
		# Split data in training and hold-out sets. Also center.
		sem_data_cv %>% 
			group_by_(group_var) %>% 
			sample_frac(0.5) %>% 
			mutate_all(funs(. - mean(.))) %>% 
			mutate_all(funs(. / sd(.))) %>% 
			ungroup() -> sem_data_cv_train
		
		sem_data_cv %>%
			filter(!IDX %in% sem_data_cv_train$IDX) %>% 
			group_by_(group_var) %>% 
			mutate_all(funs(. - mean(.))) %>% 
			mutate_all(funs(. / sd(.))) %>% 
			ungroup() -> sem_data_cv_test
		
		
		cv_data_sets_train[[cc]] <- sem_data_cv_train
		cv_data_sets_test[[cc]] <- sem_data_cv_test
		
		# fit model to training set.
		esttl_cv <- tryCatch({estimate_ram_group_lasso(ram_matrices, sem_data_cv_train,
																									 group_var = group_var,
																									 lambda_values=lambdas_cv,
																									 verbose=super_verbose,
																									 refit_var = refit_var,
																									 refit_cov=refit_cov,
																									 refit_means = refit_means, 
																									 robust_covar=robust_covar,
																									 s_equal=s_equal,
																									 check_model_diagonal=TRUE)
		}, error=function(e){
			print(e)
			return(0)}
		) 
		
		if (is.numeric(esttl_cv)){
			next
		}
		
		lasso_fits[[cc]] <- esttl_cv

		# Check fit on hold-out set.
		#lambda_fits_tmp <- numeric(length(esttl_cv$lambdas))
		lambda_fits_tmp <- matrix(nrow=length(esttl_cv$lambdas), ncol=n_grp)
		
		for (ii in 1:nrow(lambda_fits_tmp)){
			
			
			lambda_fits_tmp[ii,] <- ram_data_fit(esttl_cv, lambda_idx = ii,
																					new_data = sem_data_cv_test,
																					group_var = group_var, 
																					robust_covar=robust_covar,
																					groupwise = TRUE,
																					s_equal=s_equal)
		}
		# Fit when lambda = 0
		lambda_fits_tmp <- rbind(lambda_fits_tmp, ram_data_fit(esttl_cv$sem_res_init, 
																														new_data = sem_data_cv_test,
																														group_var = group_var,
																														robust_covar=robust_covar,
																														groupwise = TRUE,
																														s_equal=s_equal))
		
		# Save fits on hold-out set.
		lambda_fits_cv[, ,cc] <- lambda_fits_tmp
		
		cur_iter_end_time <- Sys.time()
		cur_iter_duration <- as.numeric(cur_iter_end_time - cur_iter_start_time, units='mins')
		iter_times[cc] <- cur_iter_duration
		
		cc <- cc + 1
	}
	
	# Add 0 to the lambda vector.
	lambdas_cv <- c(lambdas_cv, 0)
	
	
	# Find the lambda value giving the minimal CV error.
	# mean_cv_fit <- apply(lambda_fits_cv, MARGIN = 2, mean)
	# best_lambda <- lambdas_cv[which.min(mean_cv_fit)]
	
	out <- list(lambda_fits_cv=lambda_fits_cv, lambdas_cv=lambdas_cv,
							lasso_fits=lasso_fits, cv_data_sets_train=cv_data_sets_train,
							cv_data_sets_test=cv_data_sets_test,
							cv_res_init=cv_res_init)
							#mean_cv_fit=mean_cv_fit, best_lambda)
	return(out)
	
}



# make a data frame of the parameter estimates from a alsoo fit. Currently only works for A.
# which_lambda = 0 for the un-penalized estimates.
get_nice_par_table <- function(sem_lasso_res, which_matrix='A', which_lambda=1, which_group=1){
	
	stopifnot(which_matrix == 'A')
	
	if (which_lambda == 0){
		est_ram <- sem_lasso_res$sem_res_init$ram_matrices_est
		# Extract the desired matrix. 
		par_ram_mat <- est_ram[[which_group]][[which_matrix]]
	} else {
		est_ram <- get_estimated_ram_matrices(sem_lasso_res)
		# Extract the desired matrix. 
		par_ram_mat <- est_ram[[which_lambda]][[which_group]][[which_matrix]]
	}
	
	
	if (which_matrix == 'A'){
		idxx <- sem_lasso_res$sem_res_init$asym_idx # location of the parameters in the matrix.
	}
	
	aii <- arrayInd(idxx, .dim=dim(par_ram_mat))
	
	# to
	var_to <- colnames(sem_lasso_res$sem_res_init$ram_matrices[[which_matrix]])[aii[,1]]
	
	# from
	var_from <- rownames(sem_lasso_res$sem_res_init$ram_matrices[[which_matrix]])[aii[,2]]
	
	if (which_lambda == 0){
		params_est <- est_ram[[which_group]][[which_matrix]][idxx]
	} else {
		params_est <- est_ram[[which_lambda]][[which_group]][[which_matrix]][idxx]
	}
	
	out_df <- data.frame(Response=var_to, Predictor=var_from, Estimate=params_est)
	return(out_df)
}



# This function is a slight modification of the rsem.emusig function from the 
# rsem pacakge. This variant of the function does not deal with missing data
# and fixes an issue with a warning that say "Recycling array of length 1 in
# vector-array arithmetic is deprecated". 

robust_cov <- function(xx, varphi=0.1, max.it=1000){
	
	if (!is.matrix(xx)){
		xx <- as.matrix(xx)
	}
	
	ep <- 1e-06
	n <- dim(xx)[1]
	p <- dim(xx)[2]
	mu0 <- rep(0, p)
	sig0 <- diag(p)
	
	n_it <- 0
	dt <- 1
	if (varphi == 0) {
		ck <- 1e+11
		cbeta <- 1
	}
	else {
		prob <- 1 - varphi
		chip <- qchisq(prob, p)
		ck <- sqrt(chip)
		cbeta <- (p * pchisq(chip, p + 2) + chip * (1 - prob))/p
	}
	
	while (dt > ep && n_it <= max.it) {
		
		sumx <- rep(0, p)
		sumxx <- array(0, dim = c(p, p))
		sumw1 <- 0
		sumw2 <- 0
		
		sigin <- solve(sig0)
		for (i in 1:n) {
			xi <- xx[i, ]
			xi0 <- xi - mu0
			di2 <- as.numeric(xi0 %*% sigin %*% xi0) # fix the warning.
			di <- sqrt(di2)
			if (di <= ck) {
				wi1 <- 1
				wi2 <- 1/cbeta
			} else {
				wi1 <- ck/di
				wi2 <- wi1 * wi1/cbeta
			}
			sumw1 <- sumw1 + wi1
			xxi0 <- xi0 %*% t(xi0)
			sumx <- sumx + wi1 * xi
			sumxx <- sumxx + c(wi2) * xxi0
			sumw2 <- sumw2 + wi2
		}
		
		mu1 <- sumx/sumw1
		sig1 <- sumxx/n
		dt <- max(c(max(abs(mu1 - mu0)), max(abs(sig1 - sig0))))
		mu0 <- mu1
		sig0 <- sig1
		n_it <- n_it + 1
		
	}
	
	if (n_it >= max.it) 
		warning("The maximum number of iteration was exceeded. Please increase max.it in the input.")
	rownames(sig1) <- colnames(sig1)
	names(mu1) <- colnames(sig1)
	return(sig1)
	#list(mu = mu1, sigma = sig1, max.it = n_it)
	
}



