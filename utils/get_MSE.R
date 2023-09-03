get_MSE = function (beta1, beta2) {
    if (is.null(beta1) | is.null(beta2)) {
        return (NULL)
    }
    
    res = 0
    res = sum((beta1[[1]] - beta2[[1]])^2)
    for (i in 1:length(beta1[[2]])) {
        res = res + sum((beta1[[2]][[i]] - beta2[[2]][[i]])^2)
    }
    if (is.na(res)) {
        return (NULL)
    }
    return (res)
}

