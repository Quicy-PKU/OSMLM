get_xi = function (pi) {
    J = length(pi)
    xi = rep(0, J)
    xi[1] = pi[1]
    for (k in 2:J) {
        xi[k] = xi[k-1] + pi[k]
    }
    return (xi)
}