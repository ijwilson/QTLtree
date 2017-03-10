
 

"bifurcation_split" <- function(d) {
    n <- nrow(d)
    maxedges <- 4 * (n - 1)
    
    a <- .C("bsplit", 
            haps = as.integer(t(d)), 
            n    = as.integer(n), 
            nSNP = as.integer(ncol(d)), 
            mat  = as.integer(numeric(4*maxedges)), 
            leaves = as.integer(numeric(1)), 
            PACKAGE = "snptree"
            )
    
    nedges <- 2*(a$leaves-1)
    print(nedges)
    cat("length a$mat = ",length(a$mat), "we have ", nedges, " edges, and", a$leaves,  "leaves\n")
    yyy <- matrix(a$mat[1:(4*nedges)], ncol=4, byrow=FALSE)
    bb <- list(edge = yyy[, 1:2], Nnode = a$leaves - 1, edge.length = yyy[, 4], tip.label = 1:a$leaves, counts = yyy[,3])
    class(bb) <- c("phylo")

    bb
}




"bifurcation_left_right" <- function(d, first_right) {
  ## split after the 
  n <- nrow(d)
  nSNPs <- ncol(d)
  if (first_right<2 | first_right>nSNPS)
    stop("first_right must be between 2 and the number of SNPs")
  
  dleft <- d[, (first_right-1):1]
  dright <- d[, first_right:nSNPs]
  
  left_tree <- bifurcation_split(dleft)
  right_tree <- bifurcation_split(dright)
  
  
}

