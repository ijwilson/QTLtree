

"split_qtrait" <- function(d, qtrait, SplitPositions, positions, quiet = TRUE) {
    if (missing(SplitPositions)) 
        SplitPositions <- 1:ncol(d)
    if (missing(positions)) 
        positions <- 1:ncol(d)
    if (length(positions) != ncol(d)) 
        stop("length of positions does not match columns")
    if (min(SplitPositions) < 1 || max(SplitPositions) > ncol(d)) 
        stop("values in Splitpositions should lie between 1 and ncol(d)")
    
    n <- nrow(d)
    if (length(qtrait) != n) 
        stop("We need length of quantitative trait to be the same as the number of samples")
    
    maxedges <- 4 * (n - 1)
    
    a <- .C("GetQTLSplit", 
            as.integer(t(d)), 
            as.integer(n), 
            as.integer(ncol(d)), 
            as.integer(SplitPositions - 1), 
            as.integer(length(SplitPositions)), 
            as.double(qtrait), 
            edge = as.integer(numeric(maxedges)), 
            leaves = as.integer(numeric(1)),  # number of leaves
            leafcount = as.integer(numeric(n)),
            labels = as.integer(numeric(n)),  # labels at the leaves
            nodepos = as.integer(numeric(2 * n)),  # stuff for ape
            tippos = as.integer(numeric(2 * n)), 
            qtlnode = as.double(numeric(n)), 
            qtlleaf = as.double(numeric(n)), 
             PACKAGE = "snptree")
    
    
    leaves <- a$leaves
    nedge <- 4 * (leaves - 1)
    
    edge <- matrix(a$edge[1:(4 * (leaves - 1))], ncol = 2, byrow = FALSE)
    labs <- tapply(a$labels, rep(1:leaves, a$leafcount[1:leaves]), c)
    
    nodepos <- matrix(a$nodepos[1:(2 * (leaves - 1))], ncol = 2, byrow = T) + 1
    NodePos <- list(left = nodepos[, 1], right = nodepos[, 2])
    
    tippos <- matrix(a$tippos[1:(2 * leaves)], ncol = 2, byrow = T) + 1
    TipPos <- list(left = tippos[, 1], right = tippos[, 2])
    
    bb <- list(edge = edge, Nnode = leaves - 1, edge.length = rep(1, 2 * (leaves - 1)), tip.label = labs)
    class(bb) <- c("phylo")
    
    if (!quiet) 
        cat(leaves, " leaves on tree\n")
    b <- list(tree = bb, 
              nodepos = NodePos, 
              tippos = TipPos, 
              leafcount = a$leafcount[1:leaves], 
              nodecount = a$leafcount[(leaves +  1):(2 * leaves)], 
              node_val = a$qtlnode[1:(leaves - 1)], 
              leaf_val = a$qtlleaf[1:leaves], 
              n = n, 
              labels = labs)
    class(b) <- "splitqtrait"
    return(b)
}


"splitQTLTest" <- function(d, qtl, positions, maxk = 4, reps = 1000, pickStat = "A") {
    if (missing(positions)) {
        positions <- 1:ncol(d)
    }
    if (length(positions) != ncol(d)) {
        stop("length of positions does not match columns")
    }
    n <- nrow(d)
    if (length(qtl) != n) {
        stop("number of QTLs does not match number of haplotypes")
    }
    
    if (!(pickStat %in% c("A", "Z", "P", "N"))) {
        stop("Error, statPick should be one of A,P,Z,N")
    }
    
    a <- .C("splittestQTL", 
            as.integer(t(d)), 
            as.integer(n), 
            as.integer(ncol(d)), 
            as.integer(positions - 1), 
            as.integer(length(positions)), 
            as.double(qtl), 
            as.integer(reps), 
            as.integer(maxk), 
            teststat = as.double(numeric(maxk)), 
            randteststats = as.double(numeric(maxk * reps)), 
            leaves = as.integer(numeric(1)), 
            as.character(pickStat), 
            PACKAGE = "snptree")
    
    rs <- matrix(a$randteststats, ncol = maxk, byrow = TRUE)
    
    p_ranks <- apply(rbind(a$teststat, rs), 2, function(x) (rank(x)[1]))
    p <- 1 - (p_ranks - 0.5)/(reps + 1)
    
    list(testStat = a$teststat, randTestStats = rs, leaves = a$leaves, p.value = p)
    
}





80
60
40
