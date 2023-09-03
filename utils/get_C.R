get_C = function (h0, beta0, h, beta, model_type, j) {
    J = length(h) + 1
    pi = get_pi(h0, beta0, h, beta, model_type)
    C = matrix(0, J, J)
    if (model_type == 1) {
        C[1, 1] = pi[J] * (pi[J] - 1)
        for (k in 1: (J-1)) {
            C[1, k+1] = -pi[J] * pi[k]
            C[k+1, 1] = -pi[J] * pi[k]
        }
        for (k in 1:(J-1)) {
            for (l in 1:(J-1)) {
                if (k == l) {
                    C[k+1, k+1] = -pi[k] * (1 - pi[k])
                }
                else {
                    C[k+1, l+1] = pi[k] * pi[l]
                    C[l+1, k+1] = pi[k] * pi[l]
                }
            }
        }
    }
    else if (model_type == 2) {
        xi = get_xi(pi)
        if (j == 1) {
            C[1, 1] = (xi[j]^2 - xi[j])
        } else {
            C[1, 1] = (xi[j]^2 - xi[j] + xi[j-1]^2 - xi[j-1])
        }
        for (k in 1: (J-1)) {
            if (k == j-1 || k == j) {
                C[1, k+1] = -xi[k] * (1-xi[k])
                C[k+1, 1] = -xi[k] * (1-xi[k])
            }
        }
        for (k in 1:(J-1)) {
            for (l in 1:(J-1)) {
                if (k == j-1 && l == j-1) {
                    C[k+1, l+1] = -(xi[k] * (1 - xi[k]) * (1-2 * xi[k]) / pi[k+1] + xi[k]^2 * (1-xi[k])^2 / pi[k+1]^2)
                }
                else if (k == j-1 && l == j) {
                    C[k+1, l+1] = xi[k] * (1 - xi[k]) * xi[l] * (1 - xi[l]) / pi[k+1]^2
                }
                else if (k == j && l == j-1) {
                    C[k+1, l+1] = xi[k] * (1 - xi[k]) * xi[l] * (1 - xi[l]) / pi[k]^2
                }
                else if (k == j && l == j) {
                    C[k+1, l+1] = xi[k] * (1 - xi[k]) * (1-2 * xi[k]) / pi[k] - xi[k]^2 * (1 - xi[k])^2 / pi[k]^2
                }
            }
        }
    } else if (model_type == 3) {
        indices = rep(1, J+1)
        indices[2] = length(h0) + 1
        for (k in 1: (J-1)) {
            indices[k+2] = indices[k+1] + length(h[[k]])
        }
        first_derivatives = list()
        for (l in 1:J) {
            first_derivatives[[l]] = get_first_derivatives(h0, beta0, h, beta, model_type, l)
        }
        tmp = 0
        for (l in 1:(J-1)){
            tmp = tmp - (J-l) * pi[l] * first_derivatives[[l]]
        }
        C[1, 1] = (tmp[indices[1]: (indices[2] - 1)] / h0)[1]
        for (k in 1:(J-1)) {
            C[k+1, 1] = (tmp[indices[k+1]: (indices[k+2] - 1)] / h[[k]])[1]
            C[1, k+1] = C[k+1, 1]
        }
        for (l in 1:(J-1)) {
            for (k in 1:(J-1)) {
                tmp = 0
                for (ll in 1:k) {
                    tmp = tmp - pi[ll] * first_derivatives[[ll]][indices[l+1]: (indices[l+2] - 1)]
                }
                C[k+1, l+1] = (tmp / h[[l]])[1]
            }
        }
    } else if (model_type == 4) {
        tmp = 0
        for (l in 1:min(J-1, j)) {
            tmp = tmp - exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]) / (1 + exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]))^2
        }
        C[1, 1] = tmp[1, 1]
        for (k in 1:(J-1)) {
            if (k <= j) {
                tmp = -exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))^2
                C[1, k+1] = tmp[1, 1]
                C[k+1, 1] = tmp[1, 1]
            }
        }
        for (k in 1:(J-1)) {
            if (k <= j) {
                tmp = -exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))^2
                C[k+1, k+1] = tmp[1, 1]
            }
        }
    }
    return (C)
}
