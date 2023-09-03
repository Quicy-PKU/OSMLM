library(parallel)
num_cores = 15
source('utils/get_pi.R')
source('utils/get_pi_2.R')
source('utils/get_xi.R')
source('utils/get_first_derivatives.R')
source('utils/get_second_derivatives.R')
source('utils/get_C.R')
source('utils/get_MLE.R')
source('utils/get_two_step_estimator.R')
source('utils/get_MSE.R')
source('utils/generate_data.R')

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('utils/Hadamard.cpp')

do_simulation = function(model_type, data_type, rho, csv_file_name, image_file_name, data_file_name) {
  cat('model_type', model_type, 'data_type', data_type, '\n')
  tryCatch(
    expr = {
      N = 2^16
      beta0_true = rep(-0.5, 10)
      beta_true = list(
        rep(0.5, 10),
        rep(1, 10)
      )
      beta0_initial = runif(10, -0.1, 0.1)
      beta_initial = list(
        runif(10, -0.1, 0.1),
        runif(10, -0.1, 0.1)
      )
      get_h0 = function(x) {
        return (x[1:10])
      }
      
      get_h = function(x) {
        return (list(
          x[11:20],
          x[21:30]
        ))
      }
      d = 30
      #dd = 10
      if(data_type==5){
        x_tmp = read.csv('flights_clean.csv')
        x <- matrix(sample(c(as.matrix(x_tmp[,1:5])),N*d,replace = F),nrow = N)
        pi = apply(x,1,get_pi_2, 10, beta0_true, beta_true, 4)
        if(any(is.na(pi))){
          idx <- which(is.na(pi))
          x <- x[-idx,]
          pi <- pi[-idx]
        }
        y = apply(pi,2, function(pi) {which(rmultinom(1, 1, pi) == 1)})
        full_data = cbind(x, y)
      }
      #full_data = generate_data(N, model_type, data_type, beta0_true, beta_true, get_h0, get_h, d)
      for (j in 1:3) {
        print(nrow(full_data[which(full_data[,d+1]==j),]))
      }
      beta_full = get_MLE(get_h0, get_h, beta0_initial, beta_initial, full_data, model_type)
      
      print(get_MSE(list(beta0_true, beta_true), beta_full))
      J = 3
      H = list()
      for (j in 1:3) {
        H[[j]] = matrix(0, N, 10)
      }
      
      for (i in 1:N) {
        H[[1]][i,] = get_h0(full_data[i, 1:d])
        H[[2]][i,] = get_h(full_data[i, 1:d])[[1]]
        H[[3]][i,] = get_h(full_data[i, 1:d])[[2]]
      }
      n0 = 400
      K = 1000
      ns = c(1600, 1400, 1200, 1000, 800, 600)
      MSE = matrix(0.01, nrow=5, ncol=length(ns))
      for (index in 1:length(ns)) {
        start = Sys.time()
        print(start)
        n = ns[index]
        cat(n, ' ')
        tryCatch(
          expr={
            MSE[1, index] = mean(unlist(mclapply(1:K, function(x) {
              beta_uniform = get_MLE(get_h0, get_h, beta0_initial, beta_initial, full_data[(1:N)[rbinom(N, 1, rep((n0+n)/N, N))==1],], model_type)
              return(get_MSE(beta_full, beta_uniform))
            }, mc.cores=num_cores)))
            cat(MSE[1, index], ' ')
            MSE[2, index] = mean(unlist(mclapply(1:K, function(x) {
              beta_MV = get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, n, model_type, 1, rho)
              return(get_MSE(beta_full, beta_MV))
            }, mc.cores=num_cores)))
            cat(MSE[2, index], ' ')
            MSE[3, index] = mean(unlist(mclapply(1:K, function(x) {
              beta_MVc = get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, n, model_type, 2, rho)
              return(get_MSE(beta_full, beta_MVc))
            }, mc.cores=num_cores)))
            cat(MSE[3, index], ' ')
            MSE[4, index] = mean(unlist(mclapply(1:K, function(x) {
              beta_FASA_RP = get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, n, model_type, 3, rho, H, 5000, 10, flag=F)
              return(get_MSE(beta_full, beta_FASA_RP))
            }, mc.cores=num_cores)))
            cat(MSE[4, index], ' ')
            MSE[5, index] = mean(unlist(mclapply(1:K, function(x) {
              beta_FASA_RS = get_two_step_estimator_poisson(get_h0, get_h, beta0_initial, beta_initial, full_data, n0, n, model_type, 4, rho, H, r2=10, r3=5000)
              return(get_MSE(beta_full, beta_FASA_RS))
            }, mc.cores=num_cores)))
            cat(MSE[5, index], ' ')
          },
          error = function(e) {
            print(e)
          }
        )
        end = Sys.time()
        print(end)
        print(end - start)
      }
      write.csv(MSE, csv_file_name)
      print('write csv done!')
    },
    error = function(e) {
      print(e)
    }
  )
}

print('start!')
rho = 0.2
model_type = 4
for (data_type in 5) {
  csv_file_name = paste('csv/model', model_type, 'data', data_type, '.csv', sep='')
  image_file_name = paste('figures/model', model_type, 'data', data_type,'.pdf', sep='')
  data_file_name = paste('data/model', model_type, 'data', data_type, sep='')
  do_simulation(model_type, data_type, rho, csv_file_name, image_file_name, data_file_name)
}
print('end')
