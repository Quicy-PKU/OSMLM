library(MASS)
library(mvtnorm)
generate_data = function(N, model_type, data_type, beta0, beta, get_h0, get_h, d) {
    x = rep(0, d)
    m = length(get_h0(x))
    J = 1
    for (item in get_h(x)) {
        m = m + length(item)
        J = J + 1
    }
    full_data = matrix(0, nrow=N, ncol=m+1)
    sigma = matrix(rep(0.5), m, m) + diag(0.5, m, m)
    if (data_type == 1) {
        for (i in 1:N) {
            x = mvrnorm(1, rep(0, m), sigma)
            h0 = get_h0(c(x))
            h = get_h(c(x))
            pi = get_pi(h0, beta0, h, beta, model_type)
            y = which(rmultinom(1, 1, pi) == 1)
            full_data[i,] = c(x, y)
        }
        return (full_data)
    } else if (data_type == 2) {
        for (i in 1:N) {
            pi = rep(0, J)
            while (any(pi <= 0)) {
                if (rbinom(1, 1, 0.5) == 1) {
                    x = mvrnorm(1, rep(1, m), sigma)
                }
                else {
                    x = mvrnorm(1, rep(-1, m), sigma)
                }
                h0 = get_h0(c(x))
                h = get_h(c(x))
                pi = get_pi(h0, beta0, h, beta, model_type)
            }
            y = which(rmultinom(1, 1, pi) == 1)
            full_data[i,] = c(x, y)
        }
        return (full_data)
    } else if (data_type == 3) {
        for (i in 1:N) {
            x = rmvt(n=1, sigma=diag(1, m, m), df=3)[1,] / 10
            h0 = get_h0(c(x))
            h = get_h(c(x))
            pi = get_pi(h0, beta0, h, beta, model_type)
            y = which(rmultinom(1, 1, pi) == 1)
            full_data[i,] = c(x, y)
        }
        return (full_data)
    } else if (data_type == 4) {
        for (i in 1:N) {
            x = exp(mvrnorm(1, rep(0, m), sigma)) / 10
            h0 = get_h0(c(x))
            h = get_h(c(x))
            pi = get_pi(h0, beta0, h, beta, model_type)
            y = which(rmultinom(1, 1, pi) == 1)
            full_data[i,] = c(x, y)
        }
        return (full_data)
    } else if (data_type == 5){###change 0.995 to 0.99,0.999 and change get_MLE case 5 threshold accordingly.
        for (i in 1:ceiling(N*0.995)) {
            x = mvrnorm(1, rep(0, m), sigma)
            h0 = get_h0(c(x))
            h = get_h(c(x))
            pi = get_pi(h0, beta0, h, beta, model_type)
            y = which(rmultinom(1, 1, pi) == 1)
            full_data[i,] = c(x, y)
        }
        for(i in (ceiling(N*0.995)+1):ceiling(N-floor(0.005*N)/2)){
            x = runif(d,-11,-10)
            h0 = get_h0(c(x))
            h = get_h(c(x))
            #pi = get_pi(h0, beta0, h, beta, model_type)
            y = 2
            full_data[i,] = c(x, y)
        }
        for(i in ceiling(N-floor(0.005*N)/2+1):N){
            x = runif(d,10,11)
            h0 = get_h0(c(x))
            h = get_h(c(x))
            #pi = get_pi(h0, beta0, h, beta, model_type)
            y = 3
            full_data[i,] = c(x, y)
        }
        return (full_data)
    }
}
