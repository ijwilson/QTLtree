
 

"bifurcation_split" <- function(d, quiet = TRUE) {
    n <- nrow(d)
    maxedges <- 4 * (n - 1)
    
    a <- .C("bsplit", 
            haps = as.integer(t(d)), 
            n    = as.integer(n), 
            nSNP = as.integer(ncol(d)), 
            mat  = as.integer(numeric(4*maxedges)), 
            leaves = as.integer(numeric(1)), 
            PACKAGE = "snptree")
    
    len <- 2*(a$leaves-1)
    
    yyy <- matrix(a$mat[4*len], ncol=4, byrow=FALSE)
    print(yyy)


    edge <- yyy[, 1:2]

    bb <- list(edge = edge, Nnode = a$leaves - 1, edge.length = rep(1, 2 * (a$leaves - 1)), tip.label = 1:n)
    class(bb) <- c("phylo")
    

    
    b <- list(tree = bb, counts = yyy[,3], positions=yyy[,4])
    
    class(b) <- c("bsplit")
    b
}

