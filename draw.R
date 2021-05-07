draw = function(filename) {
    x = read.csv(filename, header=F)
    dim = (ncol(x) - 1) / 2
    stopifnot(dim == 2)
    plot(c(0, 1), c(0, 1), type='n')
    densities = rep(0, nrow(x))
    for (i in 1:nrow(x)) {
        densities[i] = x[i, 5] / (x[i, 1] - x[i, 2]) * (x[i, 3] - x[i, 4])
    }
    ranks = rank(densities)

    cols = heat.colors(nrow(x))
    for (i in 1:nrow(x)) {
        rect(x[i, 1], x[i, 3], x[i, 2], x[i, 4], col=cols[ranks[i]])
    }
}
