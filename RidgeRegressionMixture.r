# Claas Heuer, 2016
# 
# Inspired from: https://darrenjw.wordpress.com/2012/11/20/getting-started-with-bayesian-variable-selection-using-jags-and-rjags/

if (!require("pacman")) install.packages("pacman")
pacman::p_load(rjags)


RidgeRegressionMixture <- function(y, X = NULL, Z = NULL, k = 1) {

	# some checks
	if(!(is.numeric(y) | is.numeric(y) | is.vector(y))) stop("y must be a numeric/integer vector")
	N = length(Y)
	
	if(k < 1 | !is.integer(k)) stop("k must be an integer >= 2")

	if(is.null(X)) X <- array(1, dim = c(N, 1))
	Nx <- ncol(X)
	if(!is.matrix(X)) stop("X must be a numeric matrix")
	if(nrow(X) != length(y)) stop("X must have as many rows as elements in y")

	hasZ <- TRUE
	is(is.null(Z)) hasZ = FALSE
	if(hasZ & nrow(Z) != length(y)) stop("Z must have as many rows as elements in y")
	if(hasZ & !is.matrix(Z)) stop("Z must be a numeric matrix")
	Nz <- NULL
	if(hasZ) Nz <- ncol(Z)

	# Must have at least one data point with fixed assignment 
	# to each cluster, otherwise some clusters will end up empty:
	Nclust = k 
	clust = rep(NA,N) 

	for(i in 1:Nclust) y[
	clust[which.min(y)]=1 # smallest value assigned to cluster 1
	clust[which.max(y)]=2 # highest value assigned to cluster 2 

	# make data.frame
	dataList = list(
			y = y ,
			N = N ,
			Nclust = k ,
			clust = clust ,
			onesRepNclust = rep(1,Nclust),
			Nz = Nz,
			Nx = Nx,
			X = X,
			Z = Z
			)

	data = list(y = y, X=Z ,n = n,p = p, betaPriorAlpha = betaPriorAlpha, betaPriorBeta = betaPriorBeta)
	init = list(alpha = 0, betaT = rep(0,p), pind = 0, ind = rep(0,p), tauB = 1)

	modelstringOLS <- "

	model {

		# Likelihood:
		for( i in 1 : N ) {

			y[i] ~ dnorm(inprod(X[i,], beta[clust[i],]), tauE[clust[i]])

			clust[i] ~ dcat( pClust[1:Nclust] )

		}


		# Prior:
		for (clustIdx in 1:Nclust) {

			# fixed effects
			for(bIdx in 1:Nx) {

				beta[clustIdx,bIdx] ~ dnorm(0, 1.0E-10)

			}


			tauE[clustIdx] ~ dgamma( 0.01 , 0.01 )

		}

		pClust[1:Nclust] ~ ddirch( onesRepNclust )

	}
	"

	modelstringMME <- "

	model {

		# Likelihood:
		for( i in 1 : N ) {

			y[i] ~ dnorm(inprod(X[i,], beta[clust[i],]) + inprod(Z[i,], u[clust[i],]), tauE[clust[i]])

			clust[i] ~ dcat( pClust[1:Nclust] )

		}


		# Prior:
		for (clustIdx in 1:Nclust) {

			# fixed effects
			for(bIdx in 1:Nx) {

				beta[clustIdx,bIdx] ~ dnorm(0, 1.0E-10)

			}

			# random effects
			for(uIdx in 1:Nz) {

				u[clustIdx,uIdx] ~ dnorm(0, tauU[clustIdx])

			}

			tauE[clustIdx] ~ dgamma( 0.01 , 0.01 )
			tauU[clustIdx] ~ dgamma( 0.01 , 0.01 )

		}

		pClust[1:Nclust] ~ ddirch( onesRepNclust )

	}
	"



	# Generate random data from known parameter values:
	set.seed(47405)
	trueM1 = 100
	N1 = 200
	trueM2 = 145 # 145 for first example below; 130 for second example
	N2 = 200
	trueSD = 15
	effsz = abs( trueM2 - trueM1 ) / trueSD
	y1 = rnorm( N1 ) 
	y1 = (y1-mean(y1))/sd(y1) * trueSD + trueM1
	y2 = rnorm( N2 ) 
	y2 = (y2-mean(y2))/sd(y2) * trueSD + trueM2
	y = c( y1 , y2 ) 
	N = length(y)

	# Must have at least one data point with fixed assignment 
	# to each cluster, otherwise some clusters will end up empty:
	Nclust = 2
	clust = rep(NA,N) 
	clust[which.min(y)]=1 # smallest value assigned to cluster 1
	clust[which.max(y)]=2 # highest value assigned to cluster 2 

	X <- array(1, dim = c(N,1))
	Nx <- ncol(X)

	Nz <- 10
	Z <- array(rnorm(N * Nz), dim = c(N, Nz))

	dataList = list(
			y = y ,
			N = N ,
			Nclust = Nclust ,
			clust = clust ,
			onesRepNclust = rep(1,Nclust),
			Nz = Nz,
			Nx = Nx,
			X = X,
			Z = Z
			)




	model=jags.model(textConnection(modelstring), data=dataList)

	update(model,n.iter=1000)

	output=coda.samples(model=model,
			    variable.names=c("pClust", "tauE", "tauU", "beta", "u", "clust"),
			    n.iter=1000,thin=1)

	onames <- attributes(output[[1]])$dimnames[[2]]
	onames <- gsub("\\[.*\\]","",onames)

	keepvar <- c("pClust","tauE", "tauU")

	plot(output[, onames %in% keepvar])
