library(parallel)
num_cores = 10
source('utils/get_pi.R')
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

average_time = function(start, end, K) {
    diff = difftime(end, start, units='secs')
    return (diff[[1]] / K)
}

data_type = 1
model_type = 4
N = 2 ^ 13
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

full_data = generate_data(N, model_type, data_type, beta0_true, beta_true, get_h0, get_h, d)
for (j in 1:3) {
    print(nrow(full_data[which(full_data[,d+1]==j),]))
}

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

K = 100

#start = Sys.time()
for (i in 1:K) {
    beta_full = get_MLE(get_h0, get_h, beta0_initial, beta_initial, full_data, model_type) 
    #MSE[i] = get_MSE(beta_full, list(beta0_true,beta_true))
}
end = Sys.time()
print(average_time(start, end, K))

n0 = 400
rho = 0.2
ns = c(1600, 1400, 1200, 1000, 800, 600)
ave_time = matrix(0.01, nrow=5, ncol=length(ns))
full_label = as.integer(full_data[, d+1])

for (index in 1:length(ns)) {
    n = ns[index]
    cat(n, ' ')
    start = Sys.time()
    for (k in 1:K) {
        idx <- (1:N)[rbinom(N, 1, rep((n0+n)/N, N))==1]
        beta_uniform = get_MLE(get_h0, get_h, beta0_initial, beta_initial, full_data[idx,], model_type)
    }
    end = Sys.time()
    ave_time[1, index] = average_time(start, end, K)
    cat(ave_time[1, index], ' ')
    
    start = Sys.time()
    for (k in 1:K) {
        N = nrow(full_data)
        d = ncol(full_data) - 1
        h0 = get_h0(rep(0, d))
        h = get_h(rep(0, d))
        p = length(h0)
        for (k in 1: length(h)) {
            p = p + length(h[[k]])
        }
        J = length(h) + 1
        pilot_data = full_data[(1:N)[rbinom(N, rep(1, N), rep(n0/N, N))==1],]
        tmp = get_MLE(get_h0, get_h, beta0_initial, beta_initial, pilot_data, model_type)
        beta_0_pilot = tmp[[1]]
        beta_pilot = tmp[[2]]
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
        
        xbeta1 <- full_data[,1:20]%*%c(beta_0_pilot,beta_pilot[[1]])
        xbeta2 <- full_data[,c(1:10,21:30)]%*%c(beta_0_pilot,beta_pilot[[2]])
        tmp1 <- cbind(-exp(xbeta1)/(1+exp(xbeta1)),-exp(xbeta2)/(1+exp(xbeta2)))
        G1 <- cbind((1+tmp1[,1])*full_data[,1:10],(1+tmp1[,1])*full_data[,11:20],matrix(0,nrow = N,ncol = 10))
        G2 <- cbind((1+tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],(1+tmp1[,2])*full_data[,21:30])
        G3 <- cbind((tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],tmp1[,2]*full_data[,21:30])
        G1[full_label!=1,] <- 0
        G2[full_label!=2,] <- 0
        G3[full_label!=3,] <- 0
        gradiant1 <- G1+G2+G3
        numerators = sqrt(colSums((M_N_Inverse %*% t(gradiant1))^2))
        w = pmin((1-rho) * n * numerators / sum(numerators) + rho * n / N, 1)
        idx = (1:N)[rbinom(N, rep(1, N), w) == 1]
        beta_MV = get_MLE(get_h0, get_h, beta_0_pilot, beta_pilot, rbind(pilot_data, full_data[idx, ]), model_type, c(rep((n0+n)/N, nrow(pilot_data)), (n0+n) / n * w[idx]))
    }
    end = Sys.time()
    ave_time[2, index] = average_time(start, end, K)
    cat(ave_time[2, index], ' ')
    

    start = Sys.time()
    for (k in 1:K) {
        N = nrow(full_data)
        d = ncol(full_data) - 1
        h0 = get_h0(rep(0, d))
        h = get_h(rep(0, d))
        p = length(h0)
        for (k in 1: length(h)) {
            p = p + length(h[[k]])
        }
        J = length(h) + 1
        pilot_data = full_data[(1:N)[rbinom(N, rep(1, N), rep(n0/N, N))==1],]
        tmp = get_MLE(get_h0, get_h, beta0_initial, beta_initial, pilot_data, model_type)
        beta_0_pilot = tmp[[1]]
        beta_pilot = tmp[[2]]
        xbeta1 <- full_data[,1:20]%*%c(beta_0_pilot,beta_pilot[[1]])
        xbeta2 <- full_data[,c(1:10,21:30)]%*%c(beta_0_pilot,beta_pilot[[2]])
        tmp1 <- cbind(-exp(xbeta1)/(1+exp(xbeta1)),-exp(xbeta2)/(1+exp(xbeta2)))
        G1 <- cbind((1+tmp1[,1])*full_data[,1:10],(1+tmp1[,1])*full_data[,11:20],matrix(0,nrow = N,ncol = 10))
        G2 <- cbind((1+tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],(1+tmp1[,2])*full_data[,21:30])
        G3 <- cbind((tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],tmp1[,2]*full_data[,21:30])
        G1[full_label!=1,] <- 0
        G2[full_label!=2,] <- 0
        G3[full_label!=3,] <- 0
        gradiant1 <- G1+G2+G3
        
        #gradiant <- apply(full_data[,1:31],1,get_first_derivatives1,beta_0_pilot,beta_pilot, model_type)
        numerators = sqrt(colSums((t(gradiant1))^2))
        w = pmin((1-rho) * n * numerators / sum(numerators) + rho * n / N, 1)
        idx = (1:N)[rbinom(N, rep(1, N), w) == 1]
        beta_MVc = get_MLE(get_h0, get_h, beta_0_pilot, beta_pilot, rbind(pilot_data, full_data[idx, ]), model_type, c(rep((n0+n)/N, nrow(pilot_data)), (n0+n) / n * w[idx]))
    }
    end = Sys.time()
    ave_time[3, index] = average_time(start, end, K)
    cat(ave_time[3, index], ' ')
    
    start = Sys.time()
    for (k in 1:K) {
        N = nrow(full_data)
        d = ncol(full_data) - 1
        h0 = get_h0(rep(0, d))
        h = get_h(rep(0, d))
        p = length(h0)
        for (k in 1: length(h)) {
            p = p + length(h[[k]])
        }
        J = length(h) + 1
        pilot_data = full_data[(1:N)[rbinom(N, rep(1, N), rep(n0/N, N))==1],]
        tmp = get_MLE(get_h0, get_h, beta0_initial, beta_initial, pilot_data, model_type)
        beta_0_pilot = tmp[[1]]
        beta_pilot = tmp[[2]]
        
        r1=5000
        r2=10
        M_N = matrix(0, nrow=p, ncol=p)
        #Cs = matrix(0, nrow=N, ncol=J*J)
        #start = Sys.time()
        xbeta1 <- full_data[,1:20]%*%c(beta_0_pilot,beta_pilot[[1]])
        xbeta2 <- full_data[,c(1:10,21:30)]%*%c(beta_0_pilot,beta_pilot[[2]])
        tmp <- cbind(-exp(xbeta1)/(1+exp(xbeta1))^2,-exp(xbeta2)/(1+exp(xbeta2))^2)
        Cs1 <- cbind((tmp[,1]+tmp[,2]),tmp,tmp[,1],tmp[,1],0,tmp[,2],0,tmp[,2])
        Cs2 <- cbind((tmp[,1]),tmp[,1],0,tmp[,1],tmp[,1],0,0,0,0)
        
        Cs1[full_label==1,] <- 0
        Cs2[full_label!=1,] <- 0
        Cs = Cs1+Cs2
        # end = Sys.time()
        # start-end
        # start = Sys.time()
        # for (i in 1:N) {
        #     x = full_data[i, 1:d]
        #     y = as.integer(full_data[i, d+1])
        #     h0 = get_h0(x)
        #     h = get_h(x)
        #     Cs[i, ] = as.vector(t(get_C(h0, beta_0_pilot, h, beta_pilot, model_type, y)))
        # }
        # end = Sys.time()
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
        
        
        tmp1 <- cbind(-exp(xbeta1)/(1+exp(xbeta1)),-exp(xbeta2)/(1+exp(xbeta2)))
        G1 <- cbind((1+tmp1[,1])*full_data[,1:10],(1+tmp1[,1])*full_data[,11:20],matrix(0,nrow = N,ncol = 10))
        G2 <- cbind((1+tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],(1+tmp1[,2])*full_data[,21:30])
        G3 <- cbind((tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],tmp1[,2]*full_data[,21:30])
        G1[full_label!=1,] <- 0
        G2[full_label!=2,] <- 0
        G3[full_label!=3,] <- 0
        gradiant1 <- G1+G2+G3
        #gradiant <- apply(full_data[,1:31],1,get_first_derivatives1,beta_0_pilot,beta_pilot, model_type)
        numerators = sqrt(rowSums((gradiant1%*%t(M_N_Inverse))^2))
        w = pmin((1-rho) * n * numerators / sum(numerators) + rho * n / N, 1)
        idx = (1:N)[rbinom(N, rep(1, N), w) == 1]
        beta_FASA_RP = get_MLE(get_h0, get_h, beta_0_pilot, beta_pilot, rbind(pilot_data, full_data[idx, ]), model_type, c(rep((n0+n)/N, nrow(pilot_data)), (n0+n) / n * w[idx]))
    }
    end = Sys.time()
    ave_time[4, index] = average_time(start, end, K)
    cat(ave_time[4, index], ' ')
    
    start = Sys.time()
    for (k in 1:K) {
        N = nrow(full_data)
        d = ncol(full_data) - 1
        h0 = get_h0(rep(0, d))
        h = get_h(rep(0, d))
        p = length(h0)
        for (k in 1: length(h)) {
            p = p + length(h[[k]])
        }
        J = length(h) + 1
        pilot_data = full_data[(1:N)[rbinom(N, rep(1, N), rep(n0/N, N))==1],]
        tmp = get_MLE(get_h0, get_h, beta0_initial, beta_initial, pilot_data, model_type)
        beta_0_pilot = tmp[[1]]
        beta_pilot = tmp[[2]]
        
        r2=10 
        r3=5000
        M_N = matrix(0, nrow=p, ncol=p)
        # Cs = matrix(0, nrow=N, ncol=J*J)
        # for (i in 1:N) {
        #     x = full_data[i, 1:d]
        #     y = as.integer(full_data[i, d+1])
        #     h0 = get_h0(x)
        #     h = get_h(x)
        #     Cs[i, ] = as.vector(t(get_C(h0, beta_0_pilot, h, beta_pilot, model_type, y)))
        # }
        #M_N = matrix(0, nrow=p, ncol=p)
        #Cs = matrix(0, nrow=N, ncol=J*J)
        #start = Sys.time()
        xbeta1 <- full_data[,1:20]%*%c(beta_0_pilot,beta_pilot[[1]])
        xbeta2 <- full_data[,c(1:10,21:30)]%*%c(beta_0_pilot,beta_pilot[[2]])
        tmp <- cbind(-exp(xbeta1)/(1+exp(xbeta1))^2,-exp(xbeta2)/(1+exp(xbeta2))^2)
        Cs1 <- cbind((tmp[,1]+tmp[,2]),tmp,tmp[,1],tmp[,1],0,tmp[,2],0,tmp[,2])
        Cs2 <- cbind((tmp[,1]),tmp[,1],0,tmp[,1],tmp[,1],0,0,0,0)
        #full_label = as.integer(full_data[, d+1])
        Cs1[full_label==1,] <- 0
        Cs2[full_label!=1,] <- 0
        Cs = Cs1+Cs2
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
        
        tmp1 <- cbind(-exp(xbeta1)/(1+exp(xbeta1)),-exp(xbeta2)/(1+exp(xbeta2)))
        G1 <- cbind((1+tmp1[,1])*full_data[,1:10],(1+tmp1[,1])*full_data[,11:20],matrix(0,nrow = N,ncol = 10))
        G2 <- cbind((1+tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],(1+tmp1[,2])*full_data[,21:30])
        G3 <- cbind((tmp1[,1]+tmp1[,2])*full_data[,1:10],tmp1[,1]*full_data[,11:20],tmp1[,2]*full_data[,21:30])
        G1[full_label!=1,] <- 0
        G2[full_label!=2,] <- 0
        G3[full_label!=3,] <- 0
        gradiant1 <- G1+G2+G3
        numerators = sqrt(rowSums(((gradiant1)%*%t(M_N_Inverse))^2))
        w = pmin((1-rho) * n * numerators / sum(numerators) + rho * n / N, 1)
        idx = (1:N)[rbinom(N, rep(1, N), w) == 1]
        beta_FASA_RS = get_MLE(get_h0, get_h, beta_0_pilot, beta_pilot, rbind(pilot_data, full_data[idx, ]), model_type, c(rep((n0+n)/N, nrow(pilot_data)), (n0+n) / n * w[idx]))
    }
    end = Sys.time()
    ave_time[5, index] = average_time(start, end, K)
    cat(ave_time[5, index], '\n')
}

write.csv(ave_time, 'csv/time1.csv')
print('finish')
