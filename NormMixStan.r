# Claas Heuer, January 2016
#
# Mixture of K normals in Stan
#
# Adopted from here: https://gist.github.com/breakbee/1456db5dad0132c3f116

library(rstan)

modelstringHet <- "
data {
	int<lower=1> N;
	int<lower=1> K;
	real y[N];
}
parameters {
	simplex[K] theta;
	real mu[K];
	real<lower=0> tau[K];
}
model {
	real ps[K];
	for (i in 1:K){
		mu[i] ~ normal(0, 1.0e+2);
	}
	for(i in 1:N){
		for(j in 1:K){
			ps[j] <- log(theta[j]) + normal_log(y[i], mu[j], tau[j]);
		}
		increment_log_prob(log_sum_exp(ps));
	}

	tau ~ cauchy(0,5);

}
"

modelstringHom <- "
data {
	int<lower=1> N;
	int<lower=1> K;
	real y[N];
}
parameters {
	simplex[K] theta;
	real mu[K];
	real<lower=0> tau;
}
model {
	real ps[K];
	for (i in 1:K){
		mu[i] ~ normal(0, 1.0e+2);
	}
	for(i in 1:N){
		for(j in 1:K){
			ps[j] <- log(theta[j]) + normal_log(y[i], mu[j], tau);
		}
		increment_log_prob(log_sum_exp(ps));
	}

	tau ~ cauchy(0,5);

}
"


NormMixStan <- function(y, K = 2, model = "het", niter = 1000, burnin = 500, verbose = FALSE, ...) {

	if(anyNA(y)) stop("No NAs allowed in y")
	if(!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector")

	N = length(y)

	# make the dat list
	dat <- list(y = y, N = length(y), K = K)

	if(!model %in% c("het","hom")) stop("model must be 'het' or 'hom'")
	if(model == "het") modelstring = modelstringHet
	if(model == "hom") modelstring = modelstringHom

	model <- stan(model_name="Gauss mixture", model_code = modelstring, data=dat,
		      iter = niter, warmup = burnin, verbose = verbose, chains=1, ...)

	out <- model@sim$samples[[1]]
	out <- as.data.frame(out)

	return(model)

}

