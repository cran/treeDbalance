#' Useful extensions to the phylo format
#'
#' \code{extendPhylo} - Extends a tree in phylo or phylo3D format, i.e., adds 
#' or updates several useful attributes of the tree that facilitate various
#' computations and allow it to be independent of a certain node enumeration.
#' These are: information on descendants, ancestors, and depths 
#' as well as on the centroids of all pending subtrees. The latter is 
#' only computed if the tree is in phylo3D format and as such contains the
#' attributes 'node.coord' and 'edge.weight'.
#'
#' @author Sophie Kersting
#' @param tree A rooted tree in phylo3D format (no special node enumeration 
#' required, except that nodes are numbered from 1 to |V| = the total number of
#' nodes). There must be at least 2 nodes, i.e., one edge. The attributes 
#' 'node.coord' and 'edge.weight' are strictly required.
#' @return \code{extendPhylo} Tree in extended phylo(3D) format, 
#' i.e., phylo(3D) format with further attributes.
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' tree <- treeDbalance::example3Dtrees$bean09
#' ext_tree <- extendPhylo(tree)
extendPhylo <- function(tree){
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  tree$node.descs <- getDescs(tree)
  tree$node.ancs <- getAncs(tree)
  tree$node.depth <- getNodeDepths(tree)
  if(sum(c("node.coord","edge.weight") %in% attributes(tree)$names)==2){
    tree$node.subtrCentr <- getSubtrCentr(tree)
  }
  return(tree)
}
#' Useful extensions to the phylo format
#'
#' \code{getDescs} - Creates a matrix with two rows, the second contains 
#' in the \eqn{i}-th entry the index of the first row in which the 
#' descendants of node \eqn{i} start. Use the function \code{getChildren} 
#' to quickly retrieve the direct descendants of each node.
#'
#' @return \code{getDescs} Numeric matrix with 2 rows.
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' getDescs(tree)
getDescs <- function(tree){
  n <- length(tree$tip.label)
  m <- tree$Nnode
  descs <- rep(NA,m+n) #for the m+n-1 descendants (the root is no descendant)+NA
  edge_to_desc <- rep(NA,m+n) #for the m+n-1 edges to descendants (as above) +NA
  descs_index <- rep(1,m+n) #where the nodes' descendants start
  for(i in 1:nrow(tree$edge)){ # count number of descendants
    source <- tree$edge[i,1]
    descs_index <- descs_index + c(rep(0,source),rep(1,m+n-source))
  }
  cur_index <- rep(0,m+n)
  for(i in 1:nrow(tree$edge)){ # transfer every row (edge) to descs
    source <- tree$edge[i,1]
    descs[descs_index[source]+cur_index[source]] <- tree$edge[i,2]
    edge_to_desc[descs_index[source]+cur_index[source]] <- i
    cur_index[source] <- cur_index[source]+1
  }
  # Add one column to facilitate getChildren.
  return(cbind(rbind(descs,edge_to_desc,descs_index), c(NA,NA,m+n))) 
}
#' Useful extensions to the phylo format
#'
#' \code{getChildren} - Creates a vector containing the direct children of a 
#' node. If the method indicates that also edges should be returned, this 
#' function will also return the number (identifier) of the incoming edge of 
#' each child.
#'
#' @param node Numeric/integer value representing a node of the tree.
#' @param method A string specifying if only descending nodes or also
#' descending edges should be returned. Can be one of 'onlyNodes' or 
#' 'alsoEdges'.
#' @return \code{getChildren} Depending on the method the function either
#' returns an integer vector containing the direct descendants of a node
#' or an integer matrix with two rows, the first containing the direct
#' descendants and the second the corresponding incoming edges.
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' getChildren(ext_tree, 3, method="alsoEdges")
getChildren <- function(tree, node, method="onlyNodes"){
  if(!"node.descs" %in% attributes(tree)$names){
    comment(paste("This may take longer as the attribute 'node.descs'",
            "does not exist and has to be calculated first."))
    tree$node.descs <- getDescs(tree)
  }
  desc <- tree$node.descs
  
  if(method=="onlyNodes"){
    if(desc["descs_index",node]<desc["descs_index",node+1]){
      return(unname(desc["descs",
                      desc["descs_index",node]:(desc["descs_index",node+1]-1)]))
    } else {
      return(NULL)
    }
    
  }else if(method=="alsoEdges"){
    if(desc["descs_index",node]<desc["descs_index",node+1]){
      return(as.matrix(unname(desc[c("descs","edge_to_desc"),
                     desc["descs_index",node]:(desc["descs_index",node+1]-1)])))
    } else {
      return(NULL)
    }
    
  }else{
    stop("Unknown method to getChildren of a node.")
  }
}
#' Useful extensions to the phylo format
#'
#' \code{getDescendants} - Creates a vector containing all descendants of a 
#' node.
#' 
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' getDescendants(ext_tree,3)
getDescendants <- function(tree, node){
  if(!"node.descs" %in% attributes(tree)$names){
    comment(paste("This may take longer as the attribute 'node.descs'",
                  "does not exist and has to be calculated first."))
    tree$node.descs <- getDescs(tree)
  }
  children <- getChildren(tree, node, method = "onlyNodes")
  if (is.null(children)) {
    return(NULL)
  } else {
    return(c(children, 
             unlist(sapply(children, 
                           function(x) getDescendants(tree = tree, node = x)))))
  }
}
#' Useful extensions to the phylo format
#'
#' \code{getAncs} - Creates a matrix that contains the parent (direct ancestor) 
#' of node \eqn{i} as well as the corresponding edge number in column \eqn{i}.
#'
#' @return \code{getAncs} Integer matrix with 2 rows. The first 
#' row contains the direct ancestor of each node, the second row the incoming
#' edge of this node, i.e., the edge that leads to its ancestor.
#' @rdname usefulPhyloExt
#' @export
#' @examples
#' getAncs(tree)
getAncs <- function(tree){
  n <- length(tree$tip.label)
  ancestor <- rep(NA,tree$Nnode+n)
  inc_edge <- rep(NA,tree$Nnode+n)
  for(i in 1:nrow(tree$edge)){ # transfer every row (edge) to ancestor
    ancestor[tree$edge[i,2]] <- tree$edge[i,1]
    inc_edge[tree$edge[i,2]] <- i
  }
  return(rbind(ancestor,inc_edge))
}
#' Useful extensions to the phylo format
#'
#' \code{getNodeDepths} - Creates a matrix with three rows:
#' The first contains the nodes ordered by increasing depth. The second
#' contains the indices at which the next depth starts in the first row, i.e.,
#' these first two rows are similar to the output matrix of \code{getDescs}.
#' The last row contains the depth of each node.
#'
#' @return \code{getNodeDepths} Numeric matrix with 3 rows.
#' @rdname usefulPhyloExt
#' @export
#' @examples
#' getNodeDepths(tree)
getNodeDepths <- function(tree){
  n <- length(tree$tip.label) #number of leaves
  m <- tree$Nnode #number of inner nodes
  if(sum(c("node.descs","node.ancs") %in% attributes(tree)$names)<2){
    comment(paste("This may take longer as at least one of the attributes,",
                  "'node.descs' and 'node.ancs' does not exist ",
                  "and has to be calculated first."))
  }
  if(!"node.descs" %in% attributes(tree)$names){
    tree$node.descs <- getDescs(tree)
  }
  if(!"node.ancs" %in% attributes(tree)$names){
    tree$node.ancs <- getAncs(tree)
  }
  ancs <- tree$node.ancs["ancestor",]
  orderByIncrDepth <- NULL # the m+n nodes ordered by depth
  depth_number <- rep(0,m+n) # number of nodes at this depth (begins with 0)
  depthOfNodes <- rep(NA, m+n) # depth of nodes
  
  #----------- go through all nodes top-down
  lastNodes <- which(is.na(ancs)) #start with the root
  current_depth <- 0
  while(length(lastNodes)>0){ 
    depthOfNodes[lastNodes] <- current_depth
    depth_number[current_depth+1] <- length(lastNodes)
    orderByIncrDepth <- c(orderByIncrDepth,lastNodes) # append nodes
    # get descendants, i.e., nodes of next depth
    lastNodes <- do.call(c, lapply(lastNodes, 
                                  function(x) getChildren(tree,x) ))
    current_depth <- current_depth +1
  }
  
  depth_index <- utils::head(cumsum(append(1,depth_number)),-1)
  return(rbind(orderByIncrDepth,depth_index,depthOfNodes))
}
#' Useful extensions to the phylo format
#'
#' \code{getNodesAtDepth} - Creates a vector containing the nodes at a certain 
#' depth.
#'
#' @param depth An integer value representing the depth of interest in the tree.
#' A depth of 0 indicates the root layer, 1 the layer of its children, and so 
#' forth.
#' @return \code{getNodesAtDepth} Integer/numeric vector 
#' containing all nodes at the desired depth.
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' getNodesAtDepth(tree,4)
getNodesAtDepth <- function(tree,depth){
  if(!"node.depth" %in% attributes(tree)$names){
    comment(paste("This may take longer as the attribute node.depth",
            "does not exist and has to be calculated first."))
    tree$node.depth <- getNodeDepths(tree)
  }
  depthmatr <- tree$node.depth
  
  if(depthmatr["depth_index",depth+1]<depthmatr["depth_index",depth+2]){
    return(unname(depthmatr["orderByIncrDepth",
                            depthmatr["depth_index",depth+1]:
                              (depthmatr["depth_index",depth+2]-1)]))
  } else {
    return(NULL)
  }
}
#' Useful extensions to the phylo format
#'
#' \code{getLeaves} - Creates a logical vector that indicates if the \eqn{i}-th
#' node is a leaf.
#'
#' @return \code{getLeaves} Creates a logical vector that indicates if the 
#' \eqn{i}-th node is a leaf, TRUE for leaf and FALSE for interior node.
#' @rdname usefulPhyloExt
#' @export
#' @examples
#' getLeaves(tree)
getLeaves <- function(tree){
  ancs <- NULL
  if("node.ancs" %in% attributes(tree)$names){
    ancs <- tree$node.ancs
  } else {
    ancs <- getAncs(tree)
  }
  is_leaf <- rep(TRUE, tree$Nnode+length(tree$tip.label))
  is_leaf[ancs[1,]] <- FALSE
  return(is_leaf)
}
#' Useful extensions to the phylo format
#'
#' \code{getSubtrCentr} - Calculates the centroid of each pending subtree. 
#' Returns a matrix containing the 3D coordinates (3 columns) where row \eqn{i}
#' gives the position of the centroid of \eqn{T_i}, the pending subtree rooted
#' in node \eqn{i}.
#' 
#' @return \code{getSubtrCentr} Numeric matrix with 2 columns.
#' @rdname usefulPhyloExt
#' @export
#' @examples
#' getSubtrCentr(ext_tree)
getSubtrCentr <- function(tree){
  if(sum(c("node.coord","edge.weight") %in% attributes(tree)$names)<2){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "at least one of these attributes:'node.coord','edge.weight'."))
  }
  if(sum(c("node.descs","node.ancs","node.depth") %in% 
         attributes(tree)$names)<3){
    comment(paste("This may take longer as at least one of the attributes,",
                  "'node.descs', 'node.ancs' and 'node.depth' does not exist ",
                  "and has to be calculated first."))
  }
  if(!"node.descs" %in% attributes(tree)$names){
    tree$node.descs <- getDescs(tree)
  }
  if(!"node.ancs" %in% attributes(tree)$names){
    tree$node.ancs <- getAncs(tree)
  }
  if(!"node.depth" %in% attributes(tree)$names){
    tree$node.depth <- getNodeDepths(tree)
  }
  desc <- tree$node.descs
  node_order <- rev(tree$node.depth["orderByIncrDepth",]) # order "bottom-up"
  subtrcentr <- matrix(NA, nrow = length(node_order), ncol = 3)
  subtrweight <- rep(NA, length(node_order))
  for(i in node_order){
    descs_i <- getChildren(tree,i,method = "alsoEdges")
    if(is.null(descs_i)){
      # If the node is a leaf, it is also its subtree's centroid with weight 0.
      subtrweight[i] <- 0
      subtrcentr[i,] <- tree$node.coord[i,] #c(T_l) is the leaf l itself
    }else{
      # If it is an internal node, obtain the centroid from a weighted mean of 
      # the childrens' centroids and the outgoing edges.
      subtrweight[i] <- sum(subtrweight[descs_i[1,]])+
                        sum(tree$edge.weight[descs_i[2,]])
      subtrcentr[i,] <- (colSums(subtrweight[descs_i[1,]]*
                                   matrix(subtrcentr[descs_i[1,],],ncol=3))+
        colSums(tree$edge.weight[descs_i[2,]]*
                  matrix(tree$node.coord[rep(i,length(descs_i[1,])),]+
                   tree$node.coord[descs_i[1,],],ncol=3)/2))/
        (sum(subtrweight[descs_i[1,]])+sum(tree$edge.weight[descs_i[2,]]))
    }
  }
  colnames(subtrcentr) <- c("centr_x","centr_y","centr_z")
  return(cbind(subtrcentr,subtrweight))
}
#' Useful extensions to the phylo format
#'
#' \code{getDistFromRoot} - Creates a vector containing the length of the path 
#' from the node to the root, i.e., the sum of the corresponding edge lengths.
#'
#' @return \code{getDistFromRoot} Integer/numeric vector containing the length 
#' of the path from each node to the root.
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' getDistFromRoot(ext_tree)
getDistFromRoot <- function(tree){
  if(!"node.coord" %in% attributes(tree)$names){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "this attribute:'node.coord'."))
  }
  if(sum(c("node.descs","node.ancs","node.depth") %in% 
         attributes(tree)$names)<3){
    comment(paste("This may take longer as at least one of the attributes,",
                  "'node.descs', 'node.ancs' and 'node.depth' does not exist ",
                  "and has to be calculated first."))
  }
  if(!"node.descs" %in% attributes(tree)$names){
    tree$node.descs <- getDescs(tree)
  }
  if(!"node.ancs" %in% attributes(tree)$names){
    tree$node.ancs <- getAncs(tree)
  }
  if(!"node.depth" %in% attributes(tree)$names){
    tree$node.depth <- getNodeDepths(tree)
  }
  node_order <- tree$node.depth["orderByIncrDepth",] # order "top-down"
  distFrRoot <- rep(NA, length(node_order))
  distFrRoot[node_order[1]] <- 0 # set 0 for root
  for(i in node_order[-1]){
    distFrRoot[i] <- distFrRoot[tree$node.ancs["ancestor",][i]]+
      sqrt(sum((tree$node.coord[i,]-
                  tree$node.coord[tree$node.ancs["ancestor",][i],])^2))
  }
  return(distFrRoot)
}
#' Useful extensions to the phylo format
#'
#' \code{getDistFromLeaf} - Creates a vector containing the length of the path 
#' from the node to the nearest descendant leaf, i.e., the sum of the 
#' corresponding edge lengths.
#'
#' @return \code{getDistFromLeaf} Integer/numeric vector containing the length 
#' of the path from each node to its nearest descendant leaf.
#' @export
#' @rdname usefulPhyloExt
#' @examples
#' getDistFromLeaf(ext_tree)
getDistFromLeaf <- function(tree){
  if(!"node.coord" %in% attributes(tree)$names){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "this attribute:'node.coord'."))
  }
  if(sum(c("node.descs","node.ancs","node.depth") %in% 
         attributes(tree)$names)<3){
    comment(paste("This may take longer as at least one of the attributes,",
                  "'node.descs', 'node.ancs' and 'node.depth' does not exist ",
                  "and has to be calculated first."))
  }
  if(!"node.descs" %in% attributes(tree)$names){
    tree$node.descs <- getDescs(tree)
  }
  if(!"node.ancs" %in% attributes(tree)$names){
    tree$node.ancs <- getAncs(tree)
  }
  if(!"node.depth" %in% attributes(tree)$names){
    tree$node.depth <- getNodeDepths(tree)
  }
  node_order <- rev(tree$node.depth["orderByIncrDepth",]) # order "bottom-up"
  distFrLeaf <- rep(NA, length(node_order))
  for(i in node_order){
    descs_i <- getChildren(tree,i, method = "alsoEdges")
    if(is.null(descs_i)){
      # If the node is a leaf, it gets 0.
      distFrLeaf[i] <- 0
    }else{
      # If it is an internal node take the minimum of: descendants' value plus
      # edge length.
      distFrLeaf[i] <- min(distFrLeaf[descs_i[1,]]+
                    sqrt(rowSums(
                      matrix((tree$node.coord[rep(i,length(descs_i[1,])),]-
                              tree$node.coord[descs_i[1,],])^2,ncol=3))))
    }
  }
  return(distFrLeaf)
}
#' Useful extensions to the phylo format
#'
#' \code{getIncEdgeLens} - Returns the length of the incoming edge of every
#' node as a numeric vector.
#' 
#' @return \code{getIncEdgeLens} Numeric vector containing the 
#' length of the incoming edge of each node, i.e., the length of the
#' edge from its direct ancestor to the node itself.
#'
#' @export
#' @rdname usefulPhyloExt
#'
#' @examples
#' getIncEdgeLens(tree)
getIncEdgeLens <- function(tree){
  if(!"node.coord" %in% attributes(tree)$names){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "this attribute:'node.coord'."))
  }
  n_edges <- nrow(tree$edge)
  inc_EL <- rep(NA,n_edges)
  for(i in 1:n_edges){
    inc_EL[tree$edge[i,2]] <- stats::dist(tree$node.coord[tree$edge[i,],])
  }
  return(inc_EL)
}
#' Useful extensions to the phylo format
#'
#' \code{getIncEdgeWeights} - Returns the weight of the incoming edge of every
#' node as a numeric vector.
#' 
#' @return \code{getIncEdgeWeights} Numeric vector containing the 
#' weight of the incoming edge of each node, i.e., the weight of the
#' edge from its direct ancestor to the node itself.
#'
#' @export
#' @rdname usefulPhyloExt
#'
#' @examples
#' getIncEdgeWeights(tree)
getIncEdgeWeights <- function(tree){
  if(!"edge.weight" %in% attributes(tree)$names){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "this attribute:'edge.weight'."))
  }
  n_edges <- nrow(tree$edge)
  inc_EW <- rep(NA,n_edges)
  for(i in 1:n_edges){
    inc_EW[tree$edge[i,2]] <- tree$edge.weight[i]
  }
  return(inc_EW)
}
