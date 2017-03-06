
"split_simple" <- function(d, SplitPositions, positions, quiet = TRUE) {
    if (missing(SplitPositions)) 
        SplitPositions <- 1:ncol(d)
    if (missing(positions)) 
        positions <- 1:ncol(d)
    if (length(positions) != ncol(d)) 
        stop("length of positions does not match columns")
    if (min(SplitPositions) < 1 || max(SplitPositions) > ncol(d)) 
        stop("values in Splitpositions should lie between 1 and ncol(d)")
    
    n <- nrow(d)
    maxedges <- 4 * (n - 1)
    
    a <- .C("split_simple", 
            as.integer(t(d)), 
            as.integer(n), 
            as.integer(ncol(d)), 
            as.integer(SplitPositions - 1), as.integer(length(SplitPositions)), 
            edge = as.integer(numeric(maxedges)), leaves = as.integer(numeric(1)), 
            leafcount = as.integer(numeric(n)), 
            labels = as.integer(numeric(n)), 
            nodepos = as.integer(numeric(2 * n)), 
            PACKAGE = "snptree")
    
    leaves <- a$leaves
    nedge <- 4 * (leaves - 1)
    
    edge <- matrix(a$edge[1:(4 * (leaves - 1))], ncol = 2, byrow = FALSE)
    labs <- tapply(a$labels, rep(1:leaves, a$leafcount[1:a$leaves]), c)
    
    bb <- list(edge = edge, Nnode = leaves - 1, edge.length = rep(1, 2 * (leaves - 1)), tip.label = labs)
    class(bb) <- c("phylo")
    
    haplotype_string <- apply(d, 1, paste, collapse = "")
    firstlabel <- sapply(labs, function(x) x[1])  ## haplotypes at the 
    tip_haplotypes <- haplotype_string[firstlabel]
    
    
    if (!quiet) 
        cat(leaves, " leaves on tree\n")
    
    b <- list(tree = bb, nodepos = a$nodepos[1:(leaves - 1)] + 1, n = n, labels = labs, tip.haplotypes = tip_haplotypes)
    
    class(b) <- c("split")
    b
}

