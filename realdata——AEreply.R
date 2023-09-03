library(parallel)
num_cores = 15

source('utils/get_pi.R')
source('utils/get_xi.R')
source('utils/get_first_derivatives.R')
source('utils/get_second_derivatives.R')
source('utils/get_C.R')
source('utils/get_MLE.R')
source('utils/get_MSE.R')
source('utils/generate_data.R')

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('utils/Hadamard.cpp')

full_data = as.matrix(read.csv('realdata/flights_clean.csv'))


print('start')
N = nrow(full_data)
model_type = 4
d = 5
J = 4
for (j in 1:J) {
  cat(j, length(which(full_data[,d+1]==j)), '\n')
}

beta0_initial = runif(5, -0.1, 0.1)
beta_initial = list(
  runif(1, -0.1, 0.1),
  runif(1, -0.1, 0.1),
  runif(1, -0.1, 0.1)
)


get_h0 = function(x) {
  return (x)
}

get_h = function(x) {
  return (list(
    c(1),
    c(1),
    c(1)
  ))
}

beta_full = get_MLE(get_h0, get_h, beta0_initial, beta_initial, full_data, model_type)
write.csv(beta_full,'realdata/flights_beta_full.csv')
p = 8
J = 4
H = list()
H[[1]] = full_data[,1:d]
for (j in 2:4) {
  H[[j]] = matrix(1, N, 1)
}

for (i in 1:N) {
  H[[1]][i,] = full_data[i, 1:d]
  H[[2]][i,] = 1
  H[[3]][i,] = 1
  H[[4]][i,] = 1
}
Y <- rep(0,nrow(full_data))
for(i in 1:nrow(full_data)){
  pi = get_pi(get_h0(full_data[i,1:d]), beta_full[[1]], get_h(full_data[i,1:d]), beta_full[[2]], model_type)
  Y[i] <- which(rmultinom(1, 1, pi) == 1)
}

full_data[,d+1] <- Y
get_two_step_estimator_poisson = function (get_h0, get_h, beta0_initial, beta_initial, full_data, n0, ns, model_type, sampling_type, rho=0, H=NULL, r1=0, r2=0, r3=0, flag=T) {
  tryCatch(
    expr = {
      n = ns[1]
      N = nrow(full_data)
      d = ncol(full_data) - 1
      h0 = get_h0(rep(0, d))
      h = get_h(rep(0, d))
      p = length(h0)
      for (k in 1: length(h)) {
        p = p + length(h[[k]])
      }
      J = length(h) + 1
      pilot_idx = (1:N)[rbinom(N, rep(1, N), rep(n0/N, N)) == 1]
      pilot_data = full_data[pilot_idx,]
      tmp = get_MLE(get_h0, get_h, beta0_initial, beta_initial, pilot_data, model_type)
      if (is.null(tmp)) {
        warning('the first step is not convergence')
        return (NULL)
      }
      beta_0_pilot = tmp[[1]]
      beta_pilot = tmp[[2]]
      
      if (sampling_type == 1) {
        M_N = matrix(0, nrow=p, ncol=p)
        for (i in 1:N) {
          x = full_data[i, 1:d]
          y = as.integer(full_data[i, d+1])
          h0 = get_h0(x)
          h = get_h(x)
          M_N = M_N + get_second_derivatives(h0, beta_0_pilot, h, beta_pilot, model_type, y)
        }
        M_N = M_N / N
        M_N_Inverse = solve(M_N)
      } else if (sampling_type == 3) {
        if (flag) {
          Rcpp::sourceCpp('utils/Hadamard.cpp')
        }
        M_N = matrix(0, nrow=p, ncol=p)
        Cs = matrix(0, nrow=N, ncol=J*J)
        for (i in 1:N) {
          x = full_data[i, 1:d]
          y = as.integer(full_data[i, d+1])
          h0 = get_h0(x)
          h = get_h(x)
          Cs[i, ] = as.vector(t(get_C(h0, beta_0_pilot, h, beta_pilot, model_type, y)))
        }
        indices = rep(1, J+1)
        indices[2] = length(h0) + 1
        for (k in 1: (J-1)) {
          indices[k+2] = indices[k+1] + length(h[[k]])
        }
        S = rep(0, N)
        idx = sample(1:N, r1)
        S[idx] = 1
        D = (2 * rbinom(N, 1, 0.5) - 1) / sqrt(N)
        for (k in 1:J) {
          for (l in 1:J) {
            if (l < k) {
              M_N[indices[k]: (indices[k+1]-1), indices[l]: (indices[l+1]-1)] = t(M_N[indices[l]: (indices[l+1]-1), indices[k]: (indices[k+1]-1)])
            }
            else {
              M_N[indices[k]: (indices[k+1]-1), indices[l]: (indices[l+1]-1)] = t(hadamard_cpp(S, D * H[[k]])) %*% (hadamard_cpp(S, D * Cs[,((k-1)*J+l)] * H[[l]]))
            }
          }
        }
        H_candidate = runif(p * r2)
        First_H = which(H_candidate <= 1/6)
        Second_H = which(5/6 < H_candidate)
        H_candidate[First_H] = sqrt(3/r2)
        H_candidate[Second_H] = -sqrt(3/r2)
        H_candidate[-c(First_H, Second_H)] = 0
        T_2 = matrix(H_candidate,p,r2)
        M_N = M_N / N
        M_N_Inverse = t(solve(t(M_N), T_2))
      } else if (sampling_type == 4) {
        M_N = matrix(0, nrow=p, ncol=p)
        Cs = matrix(0, nrow=N, ncol=J*J)
        for (i in 1:N) {
          x = full_data[i, 1:d]
          y = as.integer(full_data[i, d+1])
          h0 = get_h0(x)
          h = get_h(x)
          Cs[i, ] = as.vector(t(get_C(h0, beta_0_pilot, h, beta_pilot, model_type, y)))
        }
        indices = rep(1, J+1)
        indices[2] = length(h0) + 1
        for (k in 1: (J-1)) {
          indices[k+2] = indices[k+1] + length(h[[k]])
        }
        idx = sample(1:N, r3)
        for (k in 1:J) {
          for (l in 1:J) {
            if (l < k) {
              M_N[indices[k]: (indices[k+1]-1), indices[l]: (indices[l+1]-1)] = t(M_N[indices[l]: (indices[l+1]-1), indices[k]: (indices[k+1]-1)])
            }
            else {
              M_N[indices[k]: (indices[k+1]-1), indices[l]: (indices[l+1]-1)] = t(H[[k]][idx, ]) %*% (Cs[,((k-1)*J+l)] * H[[l]])[idx, ]
            }
          }
        }
        H_candidate = runif(p * r2)
        First_H = which(H_candidate <= 1/6)
        Second_H = which(5/6 < H_candidate)
        H_candidate[First_H] = sqrt(3/r2)
        H_candidate[Second_H] = -sqrt(3/r2)
        H_candidate[-c(First_H, Second_H)] = 0
        T_2 = matrix(H_candidate,p,r2)
        M_N = M_N / N
        M_N_Inverse = t(solve(t(M_N), T_2))
      }
      
      numerators = rep(0, N)
      denominator = 0
      for (i in 1:N) {
        x = full_data[i, 1:d]
        y = as.integer(full_data[i, d+1])
        h0 = get_h0(x)
        h = get_h(x)
        if (sampling_type == 1 || sampling_type == 3 || sampling_type == 4) {
          numerators[i] = sqrt(sum((M_N_Inverse %*% get_first_derivatives(h0, beta_0_pilot, h, beta_pilot, model_type, y))^2))
          denominator = denominator + numerators[i]
        }
        else if (sampling_type == 2) {
          numerators[i] = sqrt(sum(get_first_derivatives(h0, beta_0_pilot, h, beta_pilot, model_type, y)^2))
          denominator = denominator + numerators[i]
        }
      }
      
      w = rep(0, N)
      
      for (i in 1:N) {
        x = full_data[i, 1:d]
        y = as.integer(full_data[i, d+1])
        h0 = get_h0(x)
        h = get_h(x)
        w[i] = min((1-rho) * n * numerators[i] / denominator + rho * n / N, 1)
      }
      
      betas = list()
      for (i in 1: length(ns)) {
        idx = (1:N)[rbinom(N, rep(1, N), w * (ns[i] / n)) == 1]
        betas[[i]] = get_MLE(get_h0, get_h, beta0_initial, beta_initial, rbind(pilot_data, full_data[idx, ]), model_type, c(rep((n0+n)/N, nrow(pilot_data)), (n0+n) / n * w[idx]))
      }
      return (betas)
    },
    error = function (e) {
      print(e)
    }
  )
  return (NULL)
}



calc_mse = function(res, beta_full) {
  mse = c()
  len = length(res[[1]])
  for (i in 1:len) {
    tmp = c()
    for (j in 1:length(res)) {
      tmp = c(tmp, get_MSE(beta_full, res[[j]][[i]]))
    }
    mse[i] = mean(tmp)
  }
  return (mse)
}



n0 = 1000
K = 1000
ns = c(4000, 3600, 3200, 2800, 2400, 2000)
rho = 0.2
MSE = matrix(0.01, nrow=5, ncol=length(ns))

start = Sys.time()
print(start)        
res = mclapply(1:K, function(x) {
  if (x %% 10 == 1) {
    cat(x, ' ')
  }
  if (x %% 100 == 1) {
    cur = Sys.time()
    print(cur - start)
  }
  beta0_initial = runif(5, -0.1, 0.1)
  beta_initial = list(
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1)
  )
  return(get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, ns, model_type, 1, rho))
}, mc.cores=num_cores)
save(res, file=paste('realdata/beta/', 'mv1', '.RData',sep=''))
mse = calc_mse(res, beta_full)
for (i in 1:length(mse)) {
  MSE[2, i] = mse[i]
}
print(mse)
write.csv(MSE, 'realdata/po_mv1.csv')



start = Sys.time()
print(start)
res = mclapply(1:K, function(x) {
  if (x %% 10 == 1) {
    cat(x, ' ')
  }
  if (x %% 100 == 1) {
    cur = Sys.time()
    print(cur - start)
  }
  beta0_initial = runif(5, -0.1, 0.1)
  beta_initial = list(
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1)
  )
  return(get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, ns, model_type, 2, rho))
}, mc.cores=num_cores)
save(res, file=paste('realdata/beta/', 'mvc1', '.RData',sep=''))
mse = calc_mse(res, beta_full)
for (i in 1:length(mse)) {
  MSE[3, i] = mse[i]
}
print(mse)
write.csv(MSE, 'realdata/po_mvc1.csv')



start = Sys.time()
print(start)
res = mclapply(1:K, function(x) {
  if (x %% 10 == 1) {
    cat(x, ' ')
  }    
  if (x %% 100 == 1) {
    cur = Sys.time()
    print(cur - start)
  }
  beta0_initial = runif(5, -0.1, 0.1)
  beta_initial = list(
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1)
  )
  return(get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, ns, model_type, 3, rho, H, 20000, 4, flag=F))
}, mc.cores=num_cores)
save(res, file=paste('realdata/beta/', 'fmv_rp1', '.RData',sep=''))
mse = calc_mse(res, beta_full)
for (i in 1:length(mse)) {
  MSE[4, i] = mse[i]
}
print(mse)
write.csv(MSE, 'realdata/fmv_rp1.csv')



start = Sys.time()
print(start) 
res = mclapply(1:K, function(x) {
  if (x %% 10 == 1) {
    cat(x, ' ')
  }
  if (x %% 100 == 1) {
    cur = Sys.time()
    print(cur - start)
  }
  beta0_initial = runif(5, -0.1, 0.1)
  beta_initial = list(
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1)
  )
  return(get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, ns, model_type, 4, rho, H, r2=4, r3=20000))
}, mc.cores=num_cores)
save(res, file=paste('realdata/beta/', 'fmv_rs1', '.RData',sep=''))
mse = calc_mse(res, beta_full)
for (i in 1:length(mse)) {
  MSE[5, i] = mse[i]
}
print(mse)
write.csv(MSE, 'realdata/fmv_rs1.csv')

end = Sys.time()
print(end - start)



##############################result for simple random sampling (need to run this part)############
for(jjjjj in 1:length(ns)){
res = mclapply(1:K, function(x) {
  if (x %% 10 == 1) {
    cat(x, ' ')
  }
  if (x %% 100 == 1) {
    cur = Sys.time()
    print(cur - start)
  }
  beta0_initial = runif(5, -0.1, 0.1)
  beta_initial = list(
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1),
    runif(1, -0.1, 0.1)
  )
  return(get_MLE(get_h0, get_h, beta0_initial, beta_initial, full_data[(1:N)[rbinom(N, 1, rep((n0+ns[jjjjj])/N, N))==1],], model_type))
}, mc.cores=num_cores)

MSE[1, jjjjj] = mse = mean(unlist(lapply(res, get_MSE, beta2=beta_full)))
# for (i in 1:length(mse)) {
#   MSE[1, i] = mse[i]
# }
print(mse)
}
write.csv(MSE, 'realdata/uniform.csv')
#write.csv(MSE, 'realdata/po1.csv')