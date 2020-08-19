Ns <<- 4
Nt <<- 8
bet <<- 0.5
N_hit <<- 10
N_meas <<- 100
del <<- 0.1
grid <- array(0,dim=c(Nt*Ns^3*4,4))
i <- sqrt(as.complex(-1))




index <- function(x,mu) {
  x[1] <- x[1]%%Nt
  x[2] <- x[2]%%Ns
  x[3] <- x[3]%%Ns
  x[4] <- x[4]%%Ns
  return(x[1]*Ns^3*4+x[2]*Ns^2*4+x[3]*Ns*4+x[4]*4+mu+1)
}

cold_start <- function() {
  for (i in 0:Nt-1)
    for (j in 0:Ns-1)
      for (k in 0:Ns-1)
        for (l in 0:Ns-1)
          for (m in 0:3) {
            grid[index(c(i,j,k,l),m),1] <<- 1
          }
}

getmat <- function(x) {
  a <- x[1]+x[4]*i
  b <- x[3]+x[2]*i
  M <- matrix(data=c(a,b,-Conj(b),Conj(a)),nrow=2,ncol=2,byrow=TRUE)
  return(M)
}

su2prod <- function(x,y) {
  a <- x[1]+x[4]*i
  b <- x[3]+x[2]*i
  c <- y[1]+y[4]*i
  d <- y[3]+y[2]*i
  u1 <- a*c - b*Conj(d)
  v1 <- a*d + b*Conj(c)
  u <- u1/sqrt(Mod(u1)^2+Mod(v1)^2)
  v <- v1/sqrt(Mod(u1)^2+Mod(v1)^2)
  return(c(Re(u),Im(v),Re(v),Im(u)))
}

random_su2 <- function() {
  alpha <- runif(1,min = 0, max = del)
  u <- runif(1, min = -1, max = 1)
  thet <- runif(1, min = 0, max = 2*pi)
  n1 <- sqrt(1-u^2)*cos(thet)
  n2 <- sqrt(1-u^2)*sin(thet)
  n3 <- u
  return(c(cos(alpha),n1*sin(alpha),n2*sin(alpha),n3*sin(alpha)))
}

get_staple <- function(x,mu) {
  x1 <- x
  x2 <- x
  K <- getmat(rep.int(0,4))
  x1[mu] <- x1[mu] + 1
  for (nu in 0:3) {
    if (nu != mu) {
      x2[nu] <- x2[nu] + 1
      K <- K + getmat(grid[index(x1,nu),])%*%Conj(t(getmat(grid[index(x2,mu),])))%*%Conj(t(getmat(grid[index(x,nu),])))
      x2[nu] <- x2[nu] - 1
    }
  }
  for (nu in 0:3) {
    if (nu != mu) {
      x1[nu] <- x1[nu] - 1
      x2[nu] <- x2[nu] - 1
      K <- K + Conj(t(getmat(grid[index(x1,nu),])))%*%Conj(t(getmat(grid[index(x2,mu),])))%*%getmat(grid[index(x2,nu),])
      x1[nu] <- x1[nu] + 1
      x2[nu] <- x2[nu] + 1
    }
  }
  return(K)
}

sweep <- function() {
  for (t in 0:Nt-1)
    for (x in 0:Ns-1)
      for (y in 0:Ns-1)
        for (z in 0:Ns-1)
          for (mu in 0:3) {
            n <- c(t,x,y,z)
            K <- get_staple(n,mu)
            for (i in 1:N_hit) {
              R <- random_su2()
              U_new <- su2prod(grid[index(n,mu),],R)
              delS <- (bet/2)*Re(sum(diag(getmat(grid[index(n,mu),])%*%getmat(K) - getmat(U_new)%*%K)))
              accept <- delS<0
              if (!accept) {
                accept = runif(1) < exp(-delS)
              }
              if (accept) {
                grid[index(n,mu),] <<- U_new
              }
            }
          }
}

gauge_energy <- function() {
  energy <- 0
  for (t in 0:Nt-1)
    for (x in 0:Ns-1)
      for (y in 0:Ns-1)
        for (z in 0: Ns-1) {
          n <- c(t,x,y,z)
          nplusmu <- n
          nplusnu <- n
          for (mu in 0:2) {
            for (nu in (mu+1):3) {
              nplusmu[mu] <- nplusmu[mu]+1
              nplusnu[nu] <- nplusnu[nu]+1
              energy <- energy + Re(sum(diag(getmat(grid[index(n,mu),])%*%getmat(grid[index(nplusmu,nu),])%*%Conj(t(getmat(grid[index(nplusnu,mu),])))%*%Conj(t(getmat(grid[index(n,nu),]))))))
              nplusmu[mu] <- nplusmu[mu]-1
              nplusnu[nu] <- nplusnu[nu]-1
            }
          }
        }
  return(energy/(Ns^3*Nt*4*2*6))
}

set.seed(100)
start_time <- Sys.time()

cold_start()
energy <- array(NA,dim=N_meas+1)
energy[1] <- gauge_energy()
print(0)
print(energy[1])
for (t in 1:N_meas) {
  st <- Sys.time()
  sweep()
  energy[t+1] <- gauge_energy()
  print(t)
  print(energy[t+1])
  print(Sys.time() - st)
}
print(Sys.time() - start_time)




