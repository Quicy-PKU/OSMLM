get_second_derivatives = function (h0, beta0, h, beta, model_type, j) {
    J = length(h) + 1
    indices = rep(1, J+1)
    indices[2] = length(h0) + 1
    for (k in 1: (J-1)) {
        indices[k+2] = indices[k+1] + length(h[[k]])
    }
    pi = get_pi(h0, beta0, h, beta, model_type)
    second_derivatives = matrix(0, nrow=indices[length(indices)] - 1, ncol=indices[length(indices)] - 1)
    if (model_type == 1) {
        second_derivatives[indices[1]: (indices[2] - 1), indices[1]: (indices[2] - 1)] = pi[J] * (pi[J] - 1) * h0 %*% t(h0)
        for (k in 1: (J-1)) {
            second_derivatives[indices[1]: (indices[2] - 1), indices[k+1]: (indices[k+2] - 1)] = -pi[J] * pi[k] * h0 %*% t(h[[k]])
            second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[1]: (indices[2] - 1)] = -pi[J] * pi[k] * h[[k]] %*% t(h0)
        }
        for (k in 1:(J-1)) {
            for (l in 1:(J-1)) {
                if (k == l) {
                    second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[k+1]: (indices[k+2] - 1)] = -pi[k] * (1 - pi[k]) * h[[k]] %*% t(h[[k]])
                }
                else {
                    second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[l+1]: (indices[l+2] - 1)] = pi[k] * pi[l] * h[[k]] %*% t(h[[l]])
                }
            }
        }       
    }
    else if (model_type == 2) {
        xi = get_xi(pi)
        if (j == 1) {
            second_derivatives[indices[1]: (indices[2] - 1), indices[1]: (indices[2] - 1)] = (xi[j]^2 - xi[j]) * h0 %*% t(h0)
        } else {
            second_derivatives[indices[1]: (indices[2] - 1), indices[1]: (indices[2] - 1)] = (xi[j]^2 - xi[j] + xi[j-1]^2 - xi[j-1]) * h0 %*% t(h0)
        }
        for (k in 1: (J-1)) {
            if (k == j-1 || k == j) {
                second_derivatives[indices[1]: (indices[2] - 1), indices[k+1]: (indices[k+2] - 1)] = -xi[k] * (1-xi[k]) * h0 %*% t(h[[k]])
                second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[1]: (indices[2] - 1)] = -xi[k] * (1-xi[k]) * h[[k]] %*% t(h0)
            }
        }
        for (k in 1:(J-1)) {
            for (l in 1:(J-1)) {
                if (k == j-1 && l == j-1) {
                    second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[l+1]: (indices[l+2] - 1)] = -(xi[k] * (1 - xi[k]) * (1-2 * xi[k]) / pi[k+1] + xi[k]^2 * (1-xi[k])^2 / pi[k+1]^2) * h[[k]] %*% t(h[[l]])
                }
                else if (k == j-1 && l == j) {
                    second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[l+1]: (indices[l+2] - 1)] = xi[k] * (1 - xi[k]) * xi[l] * (1 - xi[l]) / pi[k+1]^2 * h[[k]] %*% t(h[[l]])
                }
                else if (k == j && l == j-1) {
                    second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[l+1]: (indices[l+2] - 1)] = xi[k] * (1 - xi[k]) * xi[l] * (1 - xi[l]) / pi[k]^2 * h[[k]] %*% t(h[[l]])
                }
                else if (k == j && l == j) {
                    second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[l+1]: (indices[l+2] - 1)] = (xi[k] * (1 - xi[k]) * (1-2 * xi[k]) / pi[k] - xi[k]^2 * (1 - xi[k])^2 / pi[k]^2) * h[[k]] %*% t(h[[l]])
                }
            }
        }
    } else if (model_type == 3) {
        pi = get_pi(h0, beta0, h, beta, model_type)
        first_derivatives = list()
        for (l in 1:J) {
            first_derivatives[[l]] = get_first_derivatives(h0, beta0, h, beta, model_type, l)
        }
        tmp = 0
        for (l in 1:(J-1)){
            tmp = tmp - (J-l) * pi[l] * first_derivatives[[l]]
        }
        second_derivatives[indices[1]: (indices[2] - 1), indices[1]: (indices[2] - 1)] = tmp[indices[1]: (indices[2] - 1)] %*% t(h0)
        for (k in 1:(J-1)) {
            second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[1]: (indices[2] - 1)] = tmp[indices[k+1]: (indices[k+2] - 1)] %*% t(h0)
            second_derivatives[indices[1]: (indices[2] - 1), indices[k+1]: (indices[k+2] - 1)] = t(second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[1]: (indices[2] - 1)])
        }
        for (l in 1:(J-1)) {
            for (k in 1:(J-1)) {
                tmp = 0
                for (ll in 1:k) {
                    tmp = tmp - pi[ll] * first_derivatives[[ll]][indices[l+1]: (indices[l+2] - 1)]
                }
                second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[l+1]: (indices[l+2] - 1)] = h[[k]] %*% t(tmp)
            }
        }
        
    } else if (model_type == 4) {
        tmp = 0
        for (l in 1:min(J-1, j)) {
            tmp = tmp - exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]) / (1 + exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]))^2
        }
        second_derivatives[indices[1]: (indices[2] - 1), indices[1]: (indices[2] - 1)] = tmp[1, 1] * h0 %*% t(h0)
        for (k in 1:(J-1)) {
            if (k <= j) {
                tmp = -exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))^2
                second_derivatives[indices[1]: (indices[2] - 1), indices[k+1]: (indices[k+2] - 1)] = tmp[1, 1] * h0 %*% t(h[[k]])
                second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[1]: (indices[2] - 1)] = tmp[1, 1] * h[[k]] %*% t(h0)
            }
        }
        for (k in 1:(J-1)) {
            if (k <= j) {
                tmp = -exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))^2
                second_derivatives[indices[k+1]: (indices[k+2] - 1), indices[k+1]: (indices[k+2] - 1)] = tmp[1, 1] * h[[k]] %*% t(h[[k]])
            }
        }
    }
    return (second_derivatives)
}

