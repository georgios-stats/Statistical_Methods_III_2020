# Georgios Karagiannis


rm(list=ls())
#install.packages('jpeg')
library(jpeg)

# choose a jpg image
file_name <-'Loutsobos_Greece' ;
file_name <-'monkey-island-burber' ;
file_name <-'monkey-island-ring' ; ## huge reduction
file_name <-'monkey-island-skeleton' ;

myimage <- readJPEG(paste(file_name,'.jpg',sep = ''))
n = dim(myimage)[1]
q = dim(myimage)[2]

# create 3 matrices with red blue green
x.r <- myimage[,,1]
x.b <- myimage[,,2]
x.g <- myimage[,,3]

# compute eigenvectors, and eigen values
x.r.eigen <- eigen(cov(x.r))
x.b.eigen <- eigen(cov(x.b))
x.g.eigen <- eigen(cov(x.g))

# decide how many PC you will use
nview =200
X11()
par(mfrow=c(2,2))
perc = x.r.eigen$values / sum(x.r.eigen$values)
plot(perc,
     type='l',
     xlim=c(1,nview)
)
perc = x.g.eigen$values / sum(x.g.eigen$values)
plot(perc,
     type='l',
     xlim=c(1,nview)
)
perc = x.b.eigen$values / sum(x.b.eigen$values)
plot(perc,
     type='l',
     xlim=c(1,nview)
)

# I cut at 
d = q

# compress
u.r = matrix(NaN,n,d)
for (i in 1:nrow(x.r)) {
  u.r[i,] <- t(x.r.eigen$vectors[,1:d])%*%matrix(x.r[i,]-colMeans(x.r),ncol=1)
}
u.g = matrix(NaN,n,d)
for (i in 1:nrow(x.g)) {
  u.g[i,] <- t(x.g.eigen$vectors[,1:d])%*%matrix(x.g[i,]-colMeans(x.g),ncol=1)
}
u.b = matrix(NaN,n,d)
for (i in 1:nrow(x.b)) {
  u.b[i,] <- t(x.b.eigen$vectors[,1:d])%*%matrix(x.b[i,]-colMeans(x.b),ncol=1)
}


# decompress
x.r.dec = matrix(NaN,n,q)
for (i in 1:nrow(x.r)) {
  x.r.dec[i,] <- colMeans(x.r) + x.r.eigen$vectors[,1:d]%*%matrix(u.r[i,],ncol=1)
}
x.g.dec = matrix(NaN,n,q)
for (i in 1:nrow(x.g)) {
  x.g.dec[i,] <- colMeans(x.g) + x.g.eigen$vectors[,1:d]%*%matrix(u.g[i,],ncol=1)
}
x.b.dec = matrix(NaN,n,q)
for (i in 1:nrow(x.b)) {
  x.b.dec[i,] <- colMeans(x.b) + x.b.eigen$vectors[,1:d]%*%matrix(u.b[i,],ncol=1)
}


# save me as an array
myimage.dec <- array( NaN, dim=c(n,q,3))
myimage.dec[,,1] <-x.r.dec
myimage.dec[,,2] <-x.b.dec
myimage.dec[,,3] <-x.g.dec

# save me as a jpeg
writeJPEG(myimage.dec, paste(file_name, '_dec_',d,'.jpg', sep = ''), quality = 1)


#################################

# size reduction

pcimage <- array( NaN, dim=c(n,q,3))
pcimage[,,1] <- u.r
pcimage[,,2] <- u.b
pcimage[,,3] <- u.g


# the size of the data of the pc's
as.numeric( object.size(pcimage))
# the size of the data of the original picture
as.numeric( object.size(myimage ))
# the size reduction
1- as.numeric( object.size(pcimage)) / as.numeric(object.size(myimage ))


