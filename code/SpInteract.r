library(spBayes);library(MASS);library(MBA);library(fields)


################################
##Spatial multivariate poisson
################################
##Generate some data
n <- 100 # number of locations
q <- 3   # number of 'species'

#make some locations
coords <- cbind(runif(n,1,100),runif(n,1,100))  

## generate the knots at which the spatial effects are captured
if(n<=100) knots=coords
if(n>100)  knots=as.matrix(makegrid(SpatialPoints(coords),n=100))

## spatial decay parameters (one for each species, here all the same)
theta <- rep(3/50,q) 

### Now make species interaction matrix
nltr <- q*(q+1)/2  #number of values in lower triangle
A <- matrix(0,q,q)
A[lower.tri(A,TRUE)] <- rep(c(1,-1),nltr/2) #rnorm(nltr,1,1)
K <- A%*%t(A)

## cor(K) is the interaction matrix as we think of it
## see http://blue.for.msu.edu/JBC_10/SC/slides/MultivariateSpatial.pdf
cor(K)

## Pure error  - here all zeros, can change this as needed to explore implications of noiser data...
Psi <- diag(0,q)

## Construct multivariate predictive process covariance matrices
## following example in spBayes spMvGLM function
## in this case the knots are the same as the points, so essentially this is just a 'normal' point-process model...
c1 <- mkSpCov(coords=coords, K=K, Psi=Psi,theta=theta,cov.model="exponential")

## generate the spatial-species effects for each species in each location
## note this is a n x q vector 
w <- mvrnorm(1,rep(0,nrow(c1)),c1)

## make some simple environmental data, here all ones to indicate constant environment.
## Given m univariate design matrices, the function mkMvX creates a multivariate design matrix
X <- mkMvX(list(matrix(1,n,1), matrix(1,n,1), matrix(1,n,1)))

## specify betas (species respond differently to the same environment)
beta <- c(-1,0,1)

## generate the response data - counts of each species at each point
y <- rpois(n*q, exp(X%*%beta+w))
## convert that to a location x species matrix of abundances
y=matrix(y,ncol=q,byrow=T)

## y is now our species abundance data for three species
head(y)

## build spatial dataframe for plotting
#d=data.frame(x=coords[,1],y=coords[,2],spid=spid,abun=y)
#coordinates(d)=c("x","y")

## look at species locations and abundances
scale=max(y)/8  #simple scale to make circles pleasing, adjust if too big or too small
plot(coords,cex=y[,1]/scale,col="red",main="Example species abundances at various locations generated with species interaction matrix",ylab="Y",xlab="X")
points(coords,cex=y.2/scale,col="green")
points(coords,cex=y.3/scale,col="blue")
points(knots,pch=3)
legend("bottomright",legend=c("Species 1","Species 2","Species 3"),text.col=c("red","green","blue"))
## look at the interaction matrix
cor(K)

## see how species abundance is correlated 
splom(cbind(y.1,y.2,y.3))

############################################################################
############################################################################
## now fit spMvGLM to recover species interactions and spatial effects from the data

##Specify starting values and collect samples. For
##a true analysis, several longer chains should be
##run.
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
beta.starting <- coefficients(glm(y~X-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~X-1, family="poisson"))))
n.samples <- 100000

m.1 <- spMvGLM(list(y.1~1,y.2~1,y.3~1), family="poisson", coords=coords,
               starting= list("beta"=beta.starting, "phi"=rep(0.06,q),
                 "A"=A.starting, "w"=0),
               tuning=list("beta"=beta.tuning, "phi"=rep(0.01,q),
                 "A"=rep(0.005,nltr), "w"=0.001),
               priors=list("beta.Flat", "phi.Unif"=rep(c(0.03, 0.3),q),
                 "K.IW"=list(q+1, diag(0.1,q))),
               cov.model="exponential",
               n.samples=n.samples, sub.sample=c(n.samples/2,n.samples,10),
               verbose=TRUE, n.report=500)

######################################
### Summarize posteriors

m.1$p.samples[,paste("phi_",1:q,sep="")] <- 3/m.1$p.samples[,paste("phi_",1:q,sep="")]

print(summary(mcmc(m.1$p.samples)))
beta.hat <- colMeans(m.1$p.samples[,1:q])
w.hat <- rowMeans(m.1$sp.effects)
y.hat <- exp(X%*%beta.hat+w.hat)
y.hat.1 <- y.hat[seq(1,length(y.hat),q)]
y.hat.2 <- y.hat[seq(2,length(y.hat),q)]
y.hat.3 <- y.hat[seq(3,length(y.hat),q)]

## get posterior species correlation matrix
Khat <- matrix(0,q,q)
Khat[lower.tri(Khat,TRUE)] <- colMeans(m.1$p.samples[,grep("K",colnames(m.1$p.samples))])
Khat <- Khat%*%t(Khat)
cor(Khat) # posterior matrix
cor(K)  # original matrix

### CI on interaction matrix
for(q in c(.025,.5,.975)) {
  Khat <- matrix(0,q,q)
  Khat[lower.tri(Khat,TRUE)]
  apply(m.1$p.samples[,grep("K",colnames(m.1$p.samples))],2,function(i) quantile(i,q))


summary(lm((as.vector(cor(K))~as.vector(cor(Khat)))))

##Take a look at spatial patterns
par(mfrow=c(3,2))
surf <- mba.surf(cbind(coords,y.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.1, cex=1, col="blue")
surf <- mba.surf(cbind(coords,y.hat.1),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, add=TRUE)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat.1,0), cex=1, col="blue")
surf <- mba.surf(cbind(coords,y.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.2, cex=1, col="blue")
surf <- mba.surf(cbind(coords,y.hat.2),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat.2,0), cex=1, col="blue")
surf <- mba.surf(cbind(coords,y.3),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=y.3, cex=1, col="blue")
surf <- mba.surf(cbind(coords,y.hat.3),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(y.hat.3,0), cex=1, col="blue")



