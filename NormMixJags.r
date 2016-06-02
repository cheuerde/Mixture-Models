# Claas Heuer, January 2016
#
# Mixture of K normals in JAGS
#
# Adopted from here: http://doingbayesiandataanalysis.blogspot.com/2012/06/mixture-of-normal-distributions.html

	library(rjags)

	# The Jags model with unequal variances between X and Y (default)
	modelstringHet <- "
	model {
		# Likelihood:
		for( i in 1 : N ) {
			y[i] ~ dnorm( mu[i] , tau[i] ) 
			mu[i] <- muOfClust[ clust[i] ]
			tau[i] <- tauOfClust[clust[i]]
			clust[i] ~ dcat( pClust[1:K] )
		}
		# Prior:
		for ( clustIdx in 1: K ) {

			muOfClust[clustIdx] ~ dnorm( 0 , 1.0E-10 )
			tauOfClust[clustIdx] ~ dgamma( 0.01 , 0.01 )

		}

		pClust[1:K] ~ ddirch( onesRepK )
	}
	"

	# The Jags model with equal variances between X and Y (could make sense actually)
	modelstringHom <- "
	model {
		# Likelihood:
		for( i in 1 : N ) {

			y[i] ~ dnorm( mu[i] , tauOfClust ) 
			mu[i] <- muOfClust[ clust[i] ]
			clust[i] ~ dcat( pClust[1:K] )

		}
		# Prior:
		tauOfClust ~ dgamma( 0.01 , 0.01 )
		for ( clustIdx in 1: K ) {

			muOfClust[clustIdx] ~ dnorm( 0 , 1.0E-10 )

		}

		pClust[1:K] ~ ddirch( onesRepK )
	}
	"

	NormMixJags <- function(y, K = 2, model = "het", niter = 1000, burnin = 500, verbose = FALSE) {

		if(anyNA(y)) stop("No NAs allowed in y")
		if(!is.numeric(y) | !is.vector(y)) stop("y must be a numeric vector")

		N = length(y)

		# the cluster indicator
		clust <- rep(NA, N)

		# apparently we have to assign at least two datapoints to clusters
		clust[which.min(y)] <- 1
		clust[which.max(y)] <- 2

		onesRepK <- rep(1, K)

		# make the dat list
		dat <- list(y = y, N = length(y), K = K, clust = clust, onesRepK = onesRepK)

		if(!model %in% c("het","hom")) stop("model must be 'het' or 'hom'")
		if(model == "het") modelstring = modelstringHet
		if(model == "hom") modelstring = modelstringHom

		model <- jags.model(textConnection(modelstring), data = dat, quiet = verbose)

		# run the Gibbs Sampler
		# burnin
		update(model, n.iter = burnin, progress.bar = "none")

		# and store samples
		progress.bar = "text"
		if(!verbose) progress.bar = "none"
		output <- coda.samples(model=model,variable.names=c("muOfClust","tauOfClust","clust","pClust" ),
				       n.iter=niter - burnin,thin=1, progress.bar = progress.bar)

		return(output[[1]])

	}




