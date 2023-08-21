#' Calculation of combined 3D imbalance indices
#'
#' \code{combined3DIndex} - Calculates either the pure root imbalance value
#' with regard to a specified vertical axis or the combined 3D imbalance 
#' index value of a 3D tree in phylo3D format. 
#' The latter is a weighted mean of the integral-based 3D imbalance index 
#' value (i.e., \code{A_Index}, \code{alpha_Index}, \code{M_Index}, or 
#' \code{mu_Index} with edge length or edge weight based weighting) as well as 
#' the root imbalance value.
#' 
#' @author Sophie Kersting
#' 
#' @param tree A rooted tree in phylo3D format (no special node enumeration 
#' required, except that nodes are numbered from 1 to |V| = the total number of
#' nodes). There must be at least 2 nodes, i.e., one edge. The attributes 
#' 'node.coord' and 'edge.weight' are strictly required.
#' @param imbal_type Specifies which node imbalance measurement should be used.
#' Available are:\cr
#' "A"     - centroid angle\cr
#' "alpha" - minimal centroid angle\cr
#' "M"     - expanded relative centroid distance\cr
#' "mu"    - relative centroid distance
#' @param weight Specifies how the node imbalance values should be weighted.
#' Available weighting methods are: \cr
#' "edge_weight" (default) -> Imbalance with regards to the total edge 
#' weight.\cr
#' "edge_length" -> Imbalance with regards to the total edge length.
#' @param root_weight_factor Numeric value >0 (default 1), which specifies the
#' weight of the root imbalance value in the weighted mean. For example, a 
#' weight of 1 means that both the imbalance index value of the whole tree as 
#' well as the root imbalance value contribute equally, i.e., the unweighted 
#' mean of the two values is returned. For a larger value, the influence of 
#' the root imbalance value increases.\cr
#' If the weight is set to infinity (\code{Inf}), then the pure root
#' imbalance value is returned.
#' @param vertical_axis Numeric vector of length 3 (default (0,0,1)), which
#' specifies the given vertical axis for the given tree model. For example,
#' use the default (0,0,1) for models that grow straight upwards (e.g., trees) 
#' and (0,0,-1) for models that grow downwards (e.g., roots). The vector (0,0,0)
#' is not allowed.
#' 
#' @return \code{combined3DIndex} Numeric value in the interval between 0 
#' (included) and 1 (excluded). A value near 1 indicates a higher degree and 
#' near 0 a lower degree of asymmetry.
#'
#' @export
#' @rdname combinedIndices
#' 
#' @examples
#' tree <- treeDbalance::extendPhylo(treeDbalance::example3Dtrees$bean09)
#' combined3DIndex(tree, imbal_type = "A", weight = "edge_weight", 
#'                 root_weight_factor = 2, vertical_axis = c(0,0,1))
#' combined3DIndex(tree, imbal_type = "A", root_weight_factor = Inf, 
#'                 vertical_axis = c(0,0,1))
combined3DIndex <- function(tree, imbal_type, weight = "edge_weight", 
                            root_weight_factor = 1, vertical_axis = c(0,0,1)){
  if (!weight %in% c("edge_weight", "edge_length")) {
    warning(paste(
      "Unknown weighting method. The default 'edge_weight' is used",
      "instead."))
    weight <- "edge_weight"
  }
  if (!imbal_type %in% c("mu", "M", "A", "alpha")) {
    stop(paste("Unknown node imbalance measurement '", imbal_type,
               "'. Available are: 'A', 'alpha', 'M' and 'mu'.", sep = ""))
  }
  if (!is.numeric(root_weight_factor) || root_weight_factor<=0){
    stop("The root_weight_factor must be >0.")
  }
  if (!is.numeric(vertical_axis) || sum(vertical_axis==c(0,0,0))==3){
    stop(paste("The vertical_axis must be a numerical vector of length 3 ",
               "and not (0,0,0).", sep = ""))
  }
  tree <- treeDbalance::extendPhylo(tree)
  # Compute the root imbalance value.
  old_root <- which(is.na(tree$node.ancs["ancestor",]))
  root_imbal <- get(paste("imbalSubdiv_",imbal_type,sep = ""))(
    x=0, 
    p=tree$node.coord[old_root,]-vertical_axis, 
    v=tree$node.coord[old_root,],
    centr_v=tree$node.subtrCentr[old_root,1:3],
    centr_v_weight=tree$node.subtrCentr[old_root,4], 
    edge_weight=1) #edge weight irrelevant for x=0
  if(is.infinite(root_weight_factor)){
    return(root_imbal)
  } else {
    # Compute the imbalance value of the whole tree.
    tree_imbal <- get(paste(imbal_type,"_Index",sep = ""))(tree, weight)
    # Compute the weighted mean.
    weight_T <- sum(tree$edge.weight)
    combIndex <- (root_weight_factor*weight_T*root_imbal+
                    weight_T*tree_imbal)/((1+root_weight_factor)*weight_T)
    return(combIndex)
  }
}
