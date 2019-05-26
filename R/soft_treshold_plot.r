

################################################################################
# This script creates figure 2 in the paper, with the soft tresholding function.
################################################################################




soft_treshold <- function(xx, lambda=0){
	sign(xx) * pmax(0, (abs(xx) - lambda))
}



xx <- seq(-3, 3, length.out = 100)


soft_treshold(xx, lambda =7)


my_lambda <- 1

lambda_text <- expression(lambda)
minuslambda_text <- expression(-lambda)
st_text <- expression(paste(y, '=', S[lambda](x)))
yx_text <- expression(paste(y, '=', x))


plot(xx, soft_treshold(xx, lambda = my_lambda), type='l', lwd=3,
		 axes=FALSE, ylab='', xlab='')
lines(xx, xx, lty='dotdash', lwd=3)

text(x=my_lambda, y=-0.2, lambda_text, cex = 1.7)
text(x=-my_lambda-0.1, y=0.2, minuslambda_text, cex = 1.7)

text(x=my_lambda+1.25, y=0.7, st_text, cex = 1.3)
text(x=my_lambda+0.25, y=1.7, yx_text, cex = 1.3)

abline(h=0, v=0)


