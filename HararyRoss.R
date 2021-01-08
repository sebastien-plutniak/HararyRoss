# author: Sébastien Plutniak
# date: 2020-01-08

.getC <- function(M, unicliquals, Mrowsums){
            Cgroup.list <- list()
            CgroupIndex.list <- list()
            CgroupPrime.list <- list()
            for(i in 1:nrow(M) ){
                # skip if the point already belongs to a clique:
                if(i %in% unlist(CgroupIndex.list)) next  
                # skip if no unicliqual points:
                if(! i %in% unicliquals ) next
                CgroupIndex <- sort(c(i, which(M[, i] != 0)))
                Cgroup <- sort(colnames(M)[c(i, which(M[, i] != 0))])
                CgroupPrime <- which( Mrowsums == Mrowsums[i] ) 
                CgroupPrime <- CgroupPrime[CgroupPrime %in% CgroupIndex]
                # add to results 
                Cgroup.list <- append(Cgroup.list, list(Cgroup))
                CgroupIndex.list <- append(CgroupIndex.list, list(CgroupIndex))
                CgroupPrime.list <- append(CgroupPrime.list, list(CgroupPrime))
            }
        list(Cgroup.list = Cgroup.list,
             CgroupPrime.list = CgroupPrime.list)
}

.substract.matrix <- function(mat, CgroupPrime.list){
        i <- 1:nrow(mat)
        i <- i[ ! 1:nrow(mat) %in% unlist(CgroupPrime.list) ]
        mat[i, i]
}

.extract.cliques <- function(mat){
    cliques.list <- list() 
    M <- mat * mat %*% t(mat)
    np <- apply(mat, 1, sum) # similar to vertices' degrees
    Mrowsums <- rowSums(M)
    unicliquals <- which( Mrowsums == np * (np - 1) )
    
    if(length(unicliquals) > 0){
        res <- .getC(M, unicliquals, Mrowsums)
        cliques <- res$Cgroup.list
        # get only cliques with at least 3 vertices:
        cliques <- cliques[sapply(cliques, function(x) length(x) > 2 )]
        cliques.list <- append(cliques.list, cliques)
        mat <- .substract.matrix(mat, res$CgroupPrime.list)
    }
    else{
        CgroupIndex <- sort(c(1, which(M[, 1] != 0)))
        Cgroup <- c(1, which(mat[, 1] != 0))
        submat1 <- mat[Cgroup, Cgroup]
        submat2 <- mat[-1, -1]
        mat <- list(submat1, submat2)
    }
    list(mat, cliques.list)
}

haross.cliques <- function(mat){
    # initial tests:
    if( ! is.matrix(mat) ) stop("The argument is not a matrix.") 
    if( ncol(mat) != nrow(mat) ) stop("A square matrix is required.") 
    if( is.null(colnames(mat)) & is.null(rownames(mat)) ){
        colnames(mat) <- 1:ncol(mat)
        rownames(mat) <- 1:nrow(mat)
    }
    # set variables:
    cliques.list.final <- list() 
    mat.list <- list(mat)
    
    repeat{ # repeat while the sum of the matrix values > 0
        # run the main function:
        res <- lapply(mat.list, .extract.cliques)
        # sort results:
        #   1) extract and add the cliques to the list:
        cliques.list.final <- append(cliques.list.final, 
                                     lapply(res, function(x) x[[2]])
                                     )
        #   2) extract the list of matrices:
        mat.list <- lapply(res, function(x) x[[1]] )
        # if the list is too nested, unnest:
        if( is.list(mat.list[[1]]) & length(mat.list[[1]]) > 1 ) {
            mat.list <- unlist(mat.list, recursive = F)
        }
        # keep only the matrices with more than 2 points:
        mat.list <- mat.list[ sapply(mat.list, function(x) sum(x) > 2 ) ]
        # if there is no more matrices, break
        if( length(mat.list) == 0 ) break
    }
    
    unlist(cliques.list.final, recursive = F)
}
 


