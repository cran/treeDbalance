#' Calculation of 3D imbalance profiles
#'
#' \code{imbalProfile} - Calculates the database for a 3D imbalance profile of  
#' a 3D tree in phylo3D format for any of the four node imbalance approaches: 
#' relative centroid distance, expanded relative centroid distance, centroid 
#' angle, or minimal centroid angle. It is also used as a basis to visualize the 
#' imbalance in a 3D plot, i.e., functions \code{plotImbalPhylo3D} and 
#' \code{addImbalPhylo3D}.\cr
#' The profile is computed with regards to the z-coordinate (height), path 
#' length to the root, and path length to the nearest descendant leaf of the 
#' nodes.\cr
#' The accuracy can be determined by defining the maximal section length
#' \code{max.seclen}, which means that an edge with length \eqn{l} will be
#' subdivided ceiling(\eqn{l}/\code{max.seclen})-1 times into parts of equal 
#' length and all subdividing nodes will be analyzed for their node imbalance.
#' For example, an edge \eqn{(p,v)} of length 3 with a maximal section length 
#' of 1 would be subdivided twice and would then be evaluated for three points
#' in total: for \eqn{v} itself and for the two subdivisions at \eqn{1/3} and 
#' \eqn{2/3} of the edge length.
#' 
#' @author Sophie Kersting
#' 
#' @param tree A rooted tree in phylo3D format (no special node enumeration 
#' required, except that nodes are numbered from 1 to |V| = the total number of
#' nodes). There must be at least 2 nodes, i.e., one edge. The attributes 
#' 'node.coord' and 'edge.weight' are strictly required.
#' @param imbal_type Specifies which node imbalance measurement should be 
#' used. Available are:\cr
#' "A"     - centroid angle\cr
#' "alpha" - minimal centroid angle\cr
#' "M"     - expanded relative centroid distance\cr
#' "mu"    - relative centroid distance
#' @param max.seclen Numeric value >0 that specifies the maximal section length.
#' 
#' @return \code{imbalProfile} Numeric matrix with five columns. The rows each
#' represent the values of a single tree node or edge subdivision. The first
#' column contains the z-coordinate (height), the second the root path length, 
#' the third the (nearest) descendant leaf path length, and the fourth the 
#' imbalance value. The fifth column stores the number of the corresponding 
#' edge.
#'
#' @export
#' @rdname profiles
#' 
#' @examples
#' tree <- treeDbalance::extendPhylo(treeDbalance::example3Dtrees$bean09)
#' imbalProfile(tree, imbal_type="mu", max.seclen=1)
imbalProfile <- function(tree, imbal_type, max.seclen){
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  if(!imbal_type %in% c("A","alpha", "M", "mu")){
    stop(paste("Unknown node imbalance measurement. Available are:","'A'", 
               "'alpha'","'M'", "'mu'."))
  }
  if(max.seclen<=0){
    stop(paste("Invalid 'max.seclen'. Has to be a positive numeric value."))
  }
  if(sum(c("node.coord","edge.weight") %in% attributes(tree)$names)!=2){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "at least one of these attributes: 'node.coord','edge.weight'."))
  }
  if(sum(c("node.descs", "node.ancs","node.depth", "node.subtrCentr") %in% 
         attributes(tree)$names)<4){
    comment(paste("This may take longer as at least one of the attributes,",
                  "'node.descs', 'node.ancs',",
                  "'node.depth' and 'node.subtrCentr' does not",
                  "exist and has to be calculated first."))
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
  if(!"node.subtrCentr" %in% attributes(tree)$names){
    tree$node.subtrCentr <- getSubtrCentr(tree)
  }
  node_order <- tree$node.depth["orderByIncrDepth",] #order "top-down"
  root_dists <- getDistFromRoot(tree)
  leaf_dists <- getDistFromLeaf(tree)
  #-----------
  # Create the imbalance profile data.
  is_leaf <- getLeaves(tree)
  imbal_profile <- NULL
  for(i in node_order[-1]){ # leave out the root (has no ancestor)
    # Calculate the desired edge subdivisions.
    edge_len <- stats::dist(tree$node.coord[c(i,
                                              tree$node.ancs["ancestor",][i]),])
    edge_ht <- tree$node.coord[tree$node.ancs["ancestor",][i],3]-
                               tree$node.coord[i,3]
    subdiv_numb <- ceiling(edge_len/max.seclen)-1
    xs_subdiv <- c(0,seq(from=1/(subdiv_numb+1),by=1/(subdiv_numb+1),
                     length.out=subdiv_numb))
    edge_profile <- matrix(NA, nrow = subdiv_numb+1, ncol = 5)
    edge_profile[,1] <- tree$node.coord[i,3]+xs_subdiv*edge_ht
    edge_profile[,2] <- root_dists[i]-xs_subdiv*edge_len
    edge_profile[,3] <- leaf_dists[i]+xs_subdiv*edge_len
    edge_profile[,5] <- tree$node.ancs["inc_edge",][i]
    if(is_leaf[i]){ #if leaf edge
      edge_profile[,4] <- 0
    }else{ #if internal edge
      edge_profile[,4] <- imbalProfile_e(xs=xs_subdiv,
        p=tree$node.coord[tree$node.ancs["ancestor",][i],],
        v=tree$node.coord[i,],
        centr_v=tree$node.subtrCentr[i,1:3],
        centr_v_weight=tree$node.subtrCentr[i,4],
        edge_weight=tree$edge.weight[tree$node.ancs["inc_edge",][i]],
        imbal_type=imbal_type)
    }
    imbal_profile <- rbind(imbal_profile,edge_profile)
  }
  colnames(imbal_profile) <- c("z-coordinate","root-path-length",
                               "desc-leaf-path-length", 
                               paste("imbalance value", imbal_type), 
                               "edge")
  return(imbal_profile)
}