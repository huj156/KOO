library(MSGLasso)
library(R.matlab)
#setwd("/Users/john/Library/CloudStorage/Dropbox/huj/Fujikoshi/model_selection/my_version/Sinica/revision/simulation/simall-20240130/")
setwd("C:/Users/admin/OneDrive/simall-20240130/")



Data <- readMat("XY4lam1G.mat")
X.m<-Data$X
Y.m<-Data$Y

P <- dim(X.m)[2]
Q <- dim(Y.m)[2]
G <- P
R <- 1

gmax <- 1
cmax <- Q*P
GarrStarts <- seq(from=0,to=P-1,by=1)
GarrEnds <- seq(from=0,to=P-1,by=1)
RarrStarts <-c(0)
RarrEnds <- c(Q-1)

tmp <- FindingPQGrps(P, Q, G, R, gmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
PQgrps <- tmp$PQgrps

tmp1 <- Cal_grpWTs(P, Q, G, R, gmax, PQgrps)
grpWTs <- tmp1$grpWTs

tmp2 <- FindingGRGrps(P, Q, G, R, cmax, GarrStarts, GarrEnds, RarrStarts, RarrEnds)
GRgrps <- tmp2$GRgrps

Pen_L <- matrix(rep(1,P*Q),P,Q, byrow=TRUE)
Pen_G <- matrix(rep(1,G*R),G,R, byrow=TRUE)
grp_Norm0 <- matrix(rep(0, G*R), nrow=G, byrow=TRUE)


lam1.v <- seq(1.0, 1.5, length=6)
lamG.v <- seq(0.19, 0.25, length=7)

try.cv<- MSGLasso.cv(X.m, Y.m, grpWTs, Pen_L, Pen_G, PQgrps, GRgrps,
                     lam1.v, lamG.v, fold=5, seed=1)
MSGLassolam1 <- try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam1
MSGLassolamG  <- try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam3
writeMat("lam1G.mat", lam1G=c(MSGLassolam1,MSGLassolamG))
rm(list = ls())
