get_MLE = function (get_h0, get_h, beta0_initial, beta_initial, data, model_type, w=NULL, maxiter=1000, tolerance=1e-6) {
    d = ncol(data) - 1
    N = nrow(data)
    h0 = get_h0(rep(0, d))
    h = get_h(rep(0, d))
    J = length(h) + 1
    indices = rep(1, J+1)
    indices[2] = length(h0) + 1
    for (k in 1: (J-1)) {
        indices[k+2] = indices[k+1] + length(h[[k]])
    }
    p = indices[length(indices)] - 1
    beta0 = beta0_initial
    beta = beta_initial
    # if (length(w) != 0) {
    #     w = w / sum(w)
    # }
    error_flag = F
    for (iter in 1:maxiter) {
        first_derivatives = rep(0, p)
        second_derivatives = matrix(0, nrow=p, ncol=p)
        for (i in 1:N) {
            x = data[i, 1:d]
            y = as.integer(data[i, d+1])
            h0 = get_h0(x)
            h = get_h(x)
            if (length(w) == 0) {
                first_derivatives = first_derivatives + get_first_derivatives(h0, beta0, h, beta, model_type, y)
                second_derivatives = second_derivatives + get_second_derivatives(h0, beta0, h, beta, model_type, y)
            }
            else {
                first_derivatives = first_derivatives + get_first_derivatives(h0, beta0, h, beta, model_type, y) / w[i]
                second_derivatives = second_derivatives + get_second_derivatives(h0, beta0, h, beta, model_type, y) / w[i]
            }
        }
        tryCatch(
            expr={
                tmp = solve(second_derivatives / N, first_derivatives / N)
            },
            error= function(e){
                warning(e)
                error_flag = T
            }
        )
        if (error_flag) {
            return (NULL)
        }
        
        if (sum(tmp^2) < tolerance) {
            return (list(beta0, beta))
        }

        beta0 = beta0 - tmp[1: length(beta0)]
        for (k in 1:(J-1)) {
            beta[[k]] = beta[[k]] - tmp[indices[k+1]: (indices[k+2]-1)]
        }
    }
    warning('get MLE not convergence')
    return (NULL)
}
