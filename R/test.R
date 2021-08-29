library(Matrix)
x <- c(1,3,8,11,22)
y <- c(3,9,12,24,63)
z <- c(2,5,7,23,53)
P <- diag(rep(1, 5))
D <- diff(P,diff = 2)
x%*% t(D) %*% D %*% x
t(D) %*% D %*% y
t(D) %*% D %*% z
t(bdiag(D,D,D)) %*% bdiag(D,D,D)
