get_two_step_estimator_poisson = function (get_h0, get_h, beta0_initial, beta_initial, full_data, n0, n, model_type, sampling_type, rho=0, H=NULL, r1=0, r2=0, r3=0, flag=T) {
    tryCatch(
        expr = {
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
            idx = (1:N)[rbinom(N, rep(1, N), w) == 1]
            return (get_MLE(get_h0, get_h, beta0_initial, beta_initial, rbind(pilot_data, full_data[idx, ]), model_type, c(rep((n0+n)/N, nrow(pilot_data)), (n0+n) / n * w[idx])))
        },
        error = function (e) {
            print(e)
        }
    )
    return (NULL)
}
