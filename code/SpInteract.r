library(spBayes);library(MASS);library(MBA);library(fields)
library(reshape);library(latticeExtra)

################################
##Spatial multivariate poisson
################################
##Generate some data
n <- 100 # number of locations
q <- 3   # number of 'species'

#make some locations
coords <- cbind.data.frame(x=runif(n,1,100),y=runif(n,1,100))  

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
c1 <- mkSpCov(coords=as.matrix(coords), K=K, Psi=Psi,theta=theta,cov.model="exponential")

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
spind=rep(1:q,each=n)

## convert that to ("wide") location x species matrix of abundances
yw=matrix(y,ncol=q,byrow=T)
## yw is now our species abundance data for three species
head(yw)

### convert to presence absence
yw2=ifelse(yw>0,1,0)
head(yw2)

## look at species locations and abundances
scale=max(yw)/8  #simple scale to make circles pleasing, adjust if too big or too small
plot(coords,cex=yw[,1]/scale,col="red",main="Example species abundances at various locations generated with species interaction matrix",ylab="Y",xlab="X")
points(coords,cex=yw[,2]/scale,col="green")
points(coords,cex=yw[,3]/scale,col="blue")
legend("bottomright",legend=c("Species 1","Species 2","Species 3"),text.col=c("red","green","blue"))
## look at the interaction matrix
cor(K)

## see how species abundance is correlated 
splom(yw)

############################################################################
############################################################################
## now fit spMvGLM to recover species interactions and spatial effects from the data

##Specify starting values and collect samples. For
##a true analysis, several longer chains should be
##run.
A.starting <- diag(1,q)[lower.tri(diag(1,q), TRUE)]
beta.starting <- coefficients(glm(y~X-1, family="poisson"))
beta.tuning <- t(chol(vcov(glm(y~X-1, family="poisson"))))
n.samples <- 10000

## Make formulas  - need to update this to match environmental data 
frmls=lapply(1:q,function(i) as.formula(paste("yw[,",i,"]~1",sep="")))

m.1 <- spMvGLM(frmls, family="poisson", coords=as.matrix(coords),
               starting= list("beta"=beta.starting, "phi"=rep(0.06,q),
                 "A"=A.starting, "w"=0),
               tuning=list("beta"=beta.tuning, "phi"=rep(0.01,q),
                 "A"=rep(0.005,nltr), "w"=0.001),
               priors=list("beta.Flat", "phi.Unif"=list(rep(0.03,q),rep(0.3,q)),
                 "K.IW"=list(q+1, diag(0.1,q))),
               cov.model="exponential",
               n.samples=n.samples,
               verbose=TRUE, n.report=500)

######################################
### Summarize posteriors

## check out the chains
xyplot(m.1$p.beta.theta.samples)

## create new object with just the posterior samples
post=window(m.1$p.beta.theta.samples,start=n.samples*.75)


#m.1$p.samples[,paste("phi_",1:q,sep="")] <- 3/m.1$p.samples[,paste("phi_",1:q,sep="")]
#colnames(m.1$p.beta.theta.samples)
print(summary(mcmc(post)))
beta.hat <- apply(post[,1:q],2,mean)
w.hat <- rowMeans(m.1$p.w.samples)
y.hat <- exp((X%*%beta.hat)+w.hat)
## reshape it for easier plotting
#y.hat=matrix(y.hat,byrow=T,nrow=n)
#w.hat=matrix(w.hat,byrow=T,nrow=n)

## get posterior species correlation matrix
Khat <- matrix(0,q,q)
Khat[lower.tri(Khat,TRUE)] <- colMeans(post[,grep("K",colnames(post))])
Khat <- Khat%*%t(Khat)
cor(Khat) # posterior matrix
cor(K)  # original matrix

### CI on interaction matrix
#for(qn in c(.025,.5,.975)) {
#  Khat <- matrix(0,qn,qn)
#  Khat[lower.tri(Khat,TRUE)]
#  apply(post[,grep("K",colnames(post))],2,function(i) quantile(i,qn))
#}

### Build a single dataframe to hold the observed and predicted
dp=cbind.data.frame(coords,sp=spind,y=y,yhat=y.hat,w=w.hat)

dpl=melt(dp,id.vars=c("x","y","sp"))

xyplot(yhat~y,data=dp,groups=sp)+layer(panel.abline(0,1))

##Take a look at spatial patterns
par(mfrow=c(1,2))
surf <- mba.surf(cbind(coords,yw[,1]),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Observed counts")
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=yw[,1], cex=1, col="blue")


surf <- mba.surf(cbind(coords,yw.hat[,1]),no.X=100, no.Y=100, extend=TRUE)$xyz.est
image(surf, main="Fitted counts")
contour(surf, add=TRUE)
contour(surf, drawlabels=FALSE, add=TRUE)
text(coords, labels=round(yw.hat[,1],0), cex=1, col="blue")
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



