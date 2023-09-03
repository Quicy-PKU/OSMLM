get_first_derivatives = function (h0, beta0, h, beta, model_type, j) {
    J = length(h) + 1
    first_derivatives = c()
    pi = get_pi(h0, beta0, h, beta, model_type)
    if (model_type == 1) {
        if (j < J) {
            first_derivatives = c(first_derivatives, pi[J] * h0)
        }
        else {
            first_derivatives = c(first_derivatives, (pi[J] - 1) * h0)
        }
        
        for (k in 1:(J-1)) {
            if (k == j) {
                first_derivatives = c(first_derivatives, (1 - pi[k]) * h[[k]])
            }
            else {
                first_derivatives = c(first_derivatives, -pi[k] * h[[k]])
            }
        }
    } else if (model_type == 2) {
        xi = get_xi(pi)
        if (j == 1) {
            first_derivatives = c(first_derivatives, (1 - xi[j]) * h0)
        }
        else {
            first_derivatives = c(first_derivatives, (1 - xi[j] - xi[j-1]) * h0)
        }
        for (k in 1:(J-1)) {
            if (k == j-1) {
                first_derivatives = c(first_derivatives, -xi[k] * (1-xi[k]) / pi[k+1] * h[[k]])
            }
            else if (k == j){
                first_derivatives = c(first_derivatives, xi[k] * (1-xi[k]) / pi[k] * h[[k]])
            } else {
                first_derivatives = c(first_derivatives, 0 * h[[k]])
            }
        }
    } else if (model_type == 3) {
        xi = get_xi(pi)
        tmp = 0
        for (l in 1:(J-1)) {
            tmp = tmp + xi[l]
        }
        first_derivatives = c(first_derivatives, (J - j - tmp) * h0)
        for (k in 1:(J-1)) {
            if (k < j) {
                first_derivatives = c(first_derivatives, -xi[k] * h[[k]])
            }
            else {
                first_derivatives = c(first_derivatives, (1-xi[k]) * h[[k]])
            }
        }
    } else if (model_type == 4) {
        if (j < J) {
            tmp = 1
            for (l in 1:j) {
                tmp  = tmp - exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]) / (1 + exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]))
            }
        } else {
            tmp = 0
            for (l in 1:(J-1)) {
                tmp  = tmp - exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]) / (1 + exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]))
            }
        }
        
        first_derivatives = c(first_derivatives, tmp[1, 1] * h0)
        for (k in 1:(J-1)) {
            tmp = matrix(0)
            if (k < j) {
                tmp = -exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))
            } else if (k == j) {
                tmp = 1 / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))
                
            }
            first_derivatives = c(first_derivatives, tmp[1, 1] * h[[k]])
        }
    }
    return (first_derivatives)
}

get_first_derivatives1 = function (x, beta0, beta, model_type) {
    h0 = x[1:10]
    h=list(
        x[11:20],
        x[21:30]
    )
    j=x[31]
    J = length(h) + 1
    first_derivatives = c()
    pi = get_pi(h0, beta0, h, beta, model_type)
    if (model_type == 1) {
        if (j < J) {
            first_derivatives = c(first_derivatives, pi[J] * h0)
        }
        else {
            first_derivatives = c(first_derivatives, (pi[J] - 1) * h0)
        }
        
        for (k in 1:(J-1)) {
            if (k == j) {
                first_derivatives = c(first_derivatives, (1 - pi[k]) * h[[k]])
            }
            else {
                first_derivatives = c(first_derivatives, -pi[k] * h[[k]])
            }
        }
    } else if (model_type == 2) {
        xi = get_xi(pi)
        if (j == 1) {
            first_derivatives = c(first_derivatives, (1 - xi[j]) * h0)
        }
        else {
            first_derivatives = c(first_derivatives, (1 - xi[j] - xi[j-1]) * h0)
        }
        for (k in 1:(J-1)) {
            if (k == j-1) {
                first_derivatives = c(first_derivatives, -xi[k] * (1-xi[k]) / pi[k+1] * h[[k]])
            }
            else if (k == j){
                first_derivatives = c(first_derivatives, xi[k] * (1-xi[k]) / pi[k] * h[[k]])
            } else {
                first_derivatives = c(first_derivatives, 0 * h[[k]])
            }
        }
    } else if (model_type == 3) {
        xi = get_xi(pi)
        tmp = 0
        for (l in 1:(J-1)) {
            tmp = tmp + xi[l]
        }
        first_derivatives = c(first_derivatives, (J - j - tmp) * h0)
        for (k in 1:(J-1)) {
            if (k < j) {
                first_derivatives = c(first_derivatives, -xi[k] * h[[k]])
            }
            else {
                first_derivatives = c(first_derivatives, (1-xi[k]) * h[[k]])
            }
        }
    } else if (model_type == 4) {
        if (j < J) {
            tmp = 1
            for (l in 1:j) {
                tmp  = tmp - exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]) / (1 + exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]))
            }
        } else {
            tmp = 0
            for (l in 1:(J-1)) {
                tmp  = tmp - exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]) / (1 + exp(h0 %*% beta0 + h[[l]] %*% beta[[l]]))
            }
        }
        
        first_derivatives = c(first_derivatives, tmp[1, 1] * h0)
        for (k in 1:(J-1)) {
            tmp = matrix(0)
            if (k < j) {
                tmp = -exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]) / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))
            } else if (k == j) {
                tmp = 1 / (1 + exp(h0 %*% beta0 + h[[k]] %*% beta[[k]]))
                
            }
            first_derivatives = c(first_derivatives, tmp[1, 1] * h[[k]])
        }
    }
    return (first_derivatives)
}


