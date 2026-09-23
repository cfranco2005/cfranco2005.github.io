library(readr)

x <- read_csv("data/image1.csv")
x <- data.matrix(x)

dim(x)
x[1:5, 1:5]

image(x, col = gray.colors(256))

image(t(x[nrow(x):1, ]), col = gray.colors(256), axes = FALSE)

plot_img <- function(M, main = "") {
  image(t(M[nrow(M):1, ]), col = gray.colors(256), axes = FALSE, main = main,
        asp = nrow(M) / ncol(M))
}

plot_img(x, "Image 1")

###################

V <- crossprod(x)                    # V = X'X
e <- eigen(V, symmetric = TRUE)      # V = U Λ U'
U <- e$vectors                       # eigenvectors
lambda <- e$values                   # eigenvalues

dim(V)
dim(U)
round(lambda[1:5], 1)


###################

Z <- x %*% U                         # PC scores (n × p)
dim(Z)

k <- 1
xk <- Z[, 1:k, drop = FALSE] %*% t(U[, 1:k, drop = FALSE])

plot_img(xk, "Rank 1 approximation")


#############

round(xk[1:3, 76:80], 2)



############

plot_img <- function(M, main = "", rng = range(M)) {
  M <- pmin(pmax(M, rng[1]), rng[2])
  image(t(M[nrow(M):1, ]), col = gray.colors(256), axes = FALSE, main = main,
        asp = nrow(M) / ncol(M), zlim = rng)
}

plot_img(xk, "Rank 1 approximation", rng = range(x))


############

pca_compress <- function(X, k) {
  X <- as.matrix(X)
  e <- eigen(crossprod(X), symmetric = TRUE)  # V = X'X = U Λ U'
  U <- e$vectors
  Z <- X %*% U                                # PC scores
  Xk <- Z[, 1:k, drop = FALSE] %*% t(U[, 1:k, drop = FALSE])
  err <- norm(X - Xk, type = "F")             # Frobenius norm
  list(approx = Xk, error = err)
}

###########

out1 <- pca_compress(x, 1)
out1$error

out100 <- pca_compress(x, 100)
out100$error


#######

sqrt(sum(lambda[-1]))


########


par(mfrow = c(3, 4), mar = c(0.5, 0.5, 2, 0.5))

plot_img(x, "Original", rng = range(x))
for (k in 1:10) {
  plot_img(pca_compress(x, k)$approx, paste("k =", k), rng = range(x))
}

par(mfrow = c(1, 1))

####

plot_img(pca_compress(x, 1)$approx, "k = 1", rng = range(x))


#####
p <- ncol(x)
err_all <- sqrt(pmax(sum(lambda) - cumsum(lambda), 0))

plot(1:p, err_all, type = "l", lwd = 2,
     xlab = "k (rank)", ylab = "Frobenius error",
     main = "Approximation error vs. k, image 1")

c(err_all[1], out1$error)


par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

c(err_all[1], out1$error)



plot(1:p, err_all, type = "l", lwd = 2,
     xlab = "k (rank)", ylab = "Frobenius error",
     main = "Approximation error vs. k, image 1")



plot(0:p, c(norm(x, type = "F"), err_all), type = "l", lwd = 2,
     xlab = "k (rank)", ylab = "Frobenius error",
     main = "Approximation error vs. k, image 1")

####
u1 <- U[, 1]
z1 <- Z[, 1]
if (sum(u1) < 0) { u1 <- -u1; z1 <- -z1 }

par(mfrow = c(1, 2))
plot(u1, type = "l", lwd = 2, xlab = "Column index", ylab = "Loading",
     main = "First eigenvector u1")
plot(z1, type = "l", lwd = 2, xlab = "Row index", ylab = "Score",
     main = "First PC score z1")
par(mfrow = c(1, 1))


###############################


x2 <- data.matrix(read_csv("data/image2.csv"))
dim(x2)

plot_img(x2, "Image 2", rng = range(x2))


####

e2 <- eigen(crossprod(x2), symmetric = TRUE)
U2 <- e2$vectors
lambda2 <- e2$values
Z2 <- x2 %*% U2

round(lambda2[1:5], 1)


####

out2 <- pca_compress(x2, 1)
out2$error
norm(x2, type = "F")

par(mfrow = c(1, 2))
plot_img(x2, "Original", rng = range(x2))
plot_img(out2$approx, "k = 1", rng = range(x2))
par(mfrow = c(1, 1))



######

u1b <- U2[, 1]
z1b <- Z2[, 1]
if (sum(u1b) < 0) { u1b <- -u1b; z1b <- -z1b }

par(mfrow = c(1, 2))
plot(u1b, type = "l", lwd = 2, xlab = "Column index", ylab = "Loading",
     main = "First eigenvector u1, image 2")
plot(z1b, type = "l", lwd = 2, xlab = "Row index", ylab = "Score",
     main = "First PC score z1, image 2")
par(mfrow = c(1, 1))



matplot(t(x2[1:5, ]), type = "l", lwd = 2, lty = 1,
        xlab = "Column index", ylab = "Pixel value", main = "First 5 rows of image 2")


#####################


x3 <- data.matrix(read_csv("data/image3.csv", show_col_types = FALSE))
dim(x3)

plot_img(x3, "Image 3", rng = range(x3))


######

pca_parts <- function(X) {
  X <- as.matrix(X)
  e <- eigen(crossprod(X), symmetric = TRUE)
  U <- e$vectors
  s <- sign(colSums(U)); s[s == 0] <- 1   # fix the arbitrary eigenvector signs
  U <- sweep(U, 2, s, "*")
  Z <- X %*% U
  lambda <- pmax(e$values, 0)
  err <- sqrt(pmax(sum(lambda) - cumsum(lambda), 0))
  list(X = X, U = U, Z = Z, lambda = lambda, error = err)
}

approx_k <- function(P, k) {
  P$Z[, 1:k, drop = FALSE] %*% t(P$U[, 1:k, drop = FALSE])
}

P3 <- pca_parts(x3)

c(approx = pca_compress(x3, 5)$error, fast = P3$error[5])

round(P3$lambda[1:10])
round(P3$error[1:10] / norm(x3, type = "F"), 3)


####

par(mfrow = c(3, 4), mar = c(0.5, 0.5, 2, 0.5))
rng3 <- range(x3)
plot_img(x3, "Original", rng = rng3)
for (k in 1:10) plot_img(approx_k(P3, k), paste("k =", k), rng = rng3)
par(mfrow = c(1, 1))

p3 <- ncol(x3)
plot(0:p3, c(norm(x3, type = "F"), P3$error), type = "l", lwd = 2,
     xlab = "k (rank)", ylab = "Frobenius error",
     main = "Approximation error vs. k, image 3")

rel3 <- P3$error / norm(x3, type = "F")
c(k05 = which(rel3 <= 0.05)[1], k02 = which(rel3 <= 0.02)[1])


####

n3 <- nrow(x3)
for (k in c(10, 20, 30, 40, 50, 90)) {
  cat("k =", k,
      " rel error =", round(rel3[k], 3),
      " storage =", round(100 * k * (n3 + p3) / (n3 * p3), 1), "%\n")
}

par(mfrow = c(1, 4), mar = c(0.5, 0.5, 2, 0.5))
plot_img(x3, "Original", rng = rng3)
for (k in c(20, 30, 50)) plot_img(approx_k(P3, k), paste("k =", k), rng = rng3)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


####

par(mfrow = c(1, 2))
plot(P3$U[, 1], type = "l", lwd = 2, xlab = "Column index", ylab = "Loading",
     main = "First eigenvector u1, image 3")
plot(P3$Z[, 1], type = "l", lwd = 2, xlab = "Row index", ylab = "Score",
     main = "First PC score z1, image 3")
par(mfrow = c(1, 1))

########################################

x4 <- data.matrix(read_csv("data/image4.csv", show_col_types = FALSE))
dim(x4)

P4 <- pca_parts(x4)
rng4 <- range(x4)
n4 <- nrow(x4); p4 <- ncol(x4)

par(mfrow = c(3, 4), mar = c(0.5, 0.5, 2, 0.5))
plot_img(x4, "Original", rng = rng4)
for (k in 1:10) plot_img(approx_k(P4, k), paste("k =", k), rng = rng4)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


######
plot(0:p4, c(norm(x4, type = "F"), P4$error), type = "l", lwd = 2,
     xlab = "k (rank)", ylab = "Frobenius error",
     main = "Approximation error vs. k, image 4")

rel4 <- P4$error / norm(x4, type = "F")
for (k in c(10, 20, 30, 40, 50)) {
  cat("k =", k, " rel error =", round(rel4[k], 3),
      " storage =", round(100 * k * (n4 + p4) / (n4 * p4), 1), "%\n")
}

#####

par(mfrow = c(1, 2))
plot(P4$U[, 1], type = "l", lwd = 2, xlab = "Column index", ylab = "Loading",
     main = "First eigenvector u1, image 4")
plot(P4$Z[, 1], type = "l", lwd = 2, xlab = "Row index", ylab = "Score",
     main = "First PC score z1, image 4")
par(mfrow = c(1, 1))

par(mfrow = c(1, 4), mar = c(0.5, 0.5, 2, 0.5))
plot_img(x4, "Original", rng = rng4)
for (k in c(20, 25, 30)) plot_img(approx_k(P4, k), paste("k =", k), rng = rng4)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
