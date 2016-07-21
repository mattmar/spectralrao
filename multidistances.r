multieuclidean <- function(x) {
    tmp <- lapply(x, function(y) {
        (y[[1]]-y[[2]])^2
    })
    return(sqrt(Reduce(`+`,tmp)))
}

multimanhattan <- function(x) {
    tmp <- lapply(x, function(y) {
        abs(y[[1]]-y[[2]])
    })
    return(sqrt(Reduce(`+`,tmp)))
}

multicanberra <- function(x) {
    tmp <- lapply(x, function(y) {
        abs(y[[1]] - y[[2]]) / (y[[1]] + y[[2]])
    })
    return(Reduce(`+`,tmp))
}
