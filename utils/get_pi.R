get_pi = function (h0, beta0, h, beta, model_type) {
    J = length(h) + 1
    pi = rep(0, J)
    if (model_type == 1) {
        denominator = 1
        for (k in 1:(J-1)) {
            denominator = denominator + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]])
        }
        for (j in 1:(J-1)) {
            pi[j] = exp(h0 %*% beta0 + h[[j]] %*% beta[[j]]) / denominator
        }
        pi[J] = 1 / denominator
    } else if (model_type == 2) {
        pi[1] = exp(h0 %*% beta0 + h[[1]] %*% beta[[1]]) / (1 + exp(h0 %*% beta0 + h[[1]] %*% beta[[1]]))
        pi[J] = 1 / (1 + exp(h0 %*% beta0 + h[[J-1]] %*% beta[[J-1]]))
        for (k in 2:(J-1)) {
            pi[k] = exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]])) - exp(h0 %*% beta0 + h[[k-1]] %*% beta[[k-1]]) / (1 + exp(h0 %*% beta0 + h[[k-1]] %*% beta[[k-1]]))
        }
    } else if (model_type == 3) {
        tmp = rep(0, J-1)
        tmp[J-1] = exp(h0 %*% beta0 + h[[J-1]] %*% beta[[J-1]])
        denominator = 1 + tmp[J-1]
        for (k in (J-2):1) {
            tmp[k] = tmp[k+1] * exp(h0 %*% beta0 + h[[k]] %*% beta[[k]])
            denominator = denominator + tmp[k]
        }
        for (j in 1:(J-1)) {
            pi[j] = tmp[j] / denominator
        }
        pi[J] = 1 / denominator
    } else if (model_type == 4) {
        tmp = 1
        for (j in 1:(J-1)) {
            tmp = tmp * (1 + exp(h0 %*% beta0 + h[[j]] %*% beta[[j]]))
            pi[j] = exp(h0 %*% beta0 + h[[j]] %*% beta[[j]]) / tmp
        }
        pi[J] = 1 / tmp
    }
    return (pi)
}

