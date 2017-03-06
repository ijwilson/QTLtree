split_category <- function(d, categories, SplitPositions, positions, quiet = TRUE) {
    n <- nrow(d)
    if (missing(SplitPositions)) 
        SplitPositions <- 1:ncol(d)
    if (missing(positions)) 
        positions <- 1:ncol(d)
    if (length(positions) != ncol(d)) 
        stop("length of positions does not match columns")
    if (min(SplitPositions) < 1 || max(SplitPositions) > ncol(d)) 
        stop("values in Splitpositions should lie between 1 and ncol(d)")
    if (length(categories) != n) 
        stop("We need length of categories to be the same as the number of samples")
    
    maxedges <- 4 * (n - 1)
    ncat <- max(categories)
    
    a <- .C("GetCategorySplit", 
            as.integer(t(d)), 
            as.integer(n), 
            as.integer(ncol(d)), 
            as.integer(SplitPositions - 1), 
            as.integer(length(SplitPositions)), 
            as.integer(categories), 
            edge = as.integer(numeric(maxedges)), 
            leaves = as.integer(numeric(1)),                               ## number of leaves
            leafcount = as.integer(numeric(n)), 
            labels = as.integer(numeric(n)),                               ## labels at the leaves
            nodepos = as.integer(numeric(2 * n)),                          ## stuff for ape
            tippos = as.integer(numeric(2 * n)), 
            nodecat = as.integer(numeric(ncat * n)), 
            leafcat = as.integer(numeric(ncat * n)), 
            as.integer(ncat), 
            PACKAGE = "snptree")
    
    leaves <- a$leaves
    nedge <- 4 * (leaves - 1)
    
    edge <- matrix(a$edge[1:(4 * (leaves - 1))], ncol = 2, byrow = FALSE)
    
    # print(length(a$labels)) print(leaves) print(a$leafcount[1:leaves])
    labs <- tapply(a$labels, rep(1:leaves, a$leafcount[1:leaves]), c)
    
    nodepos <- matrix(a$nodepos[1:(2 * (leaves - 1))], ncol = 2, byrow = T) + 1
    NodePos <- list(left = nodepos[, 1], right = nodepos[, 2])
    
    tippos <- matrix(a$tippos[1:(2 * leaves)], ncol = 2, byrow = T) + 1
    TipPos <- list(left = tippos[, 1], right = tippos[, 2])
    
    bb <- list(edge = edge, Nnode = leaves - 1, edge.length = rep(1, 2 * (leaves - 1)), tip.label = labs)
    class(bb) <- c("phylo")
    
    if (!quiet) 
        cat(leaves, " leaves on tree\n")
    b <- list(tree = bb, nodepos = NodePos, tippos = TipPos, 
              nodecat = matrix(a$nodecat, ncol = ncat, byrow = T)[1:(leaves - 1), ], 
              leafcat = matrix(a$leafcat, ncol = ncat, byrow = T)[1:leaves, ], n = n, labels = labs)
    class(b) <- "splitcategory"
    b
}

