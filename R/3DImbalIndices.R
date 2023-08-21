#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{all3DImbalIndices} - This serves as a wrapper function to calculate 
#' a 3D imbalance index value of a 3D tree in phylo3D format according to the 
#' specified imbalance measurement and weighting scheme.\cr
#' If problems occur with the estimation of the integrals, try to increase 
#' the optional parameters  \cr\code{rel.tolerance} and \code{max.subdiv}.
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
#' @param rel.tolerance (Optional) Numeric value which specifies the relative 
#' tolerance which should be used for estimating the integral using 
#' stats::integrate. Set to 1e-10 by default (the stats::integrate default 
#' value is approx 3e-16).
#' @param max.subdiv (Optional) Integer value which specifies the maximal number 
#' of interval subdivisions for estimating the integral using stats::integrate. 
#' Set to 200 by default (stats::integrate default value 100).
#' 
#' @return \code{all3DImbalIndices} Numeric value indicating the internal 3D
#' imbalance according to the chosen method.
#'
#' @export
#' @rdname integralImbalIndices
#' 
#' @examples
#' tree <- treeDbalance::extendPhylo(treeDbalance::example3Dtrees$bean09)
#' all3DImbalIndices(tree, imbal_type = "A", weight="edge_length")
all3DImbalIndices <- function(tree, imbal_type, weight = "edge_weight", 
                  rel.tolerance = 1e-10, max.subdiv = 200L){
  if (imbal_type == "mu") {
    mu_Index(tree = tree, weight = weight, 
             rel.tolerance = rel.tolerance, max.subdiv = max.subdiv)
  } else if (imbal_type == "M") {
    M_Index(tree = tree, weight = weight, 
             rel.tolerance = rel.tolerance, max.subdiv = max.subdiv)
  } else if (imbal_type == "A") {
    A_Index(tree = tree, weight = weight, 
             rel.tolerance = rel.tolerance, max.subdiv = max.subdiv)
  } else if (imbal_type == "alpha") {
    alpha_Index(tree = tree, weight = weight, 
             rel.tolerance = rel.tolerance, max.subdiv = max.subdiv)
  } else {
    stop(paste("Unknown node imbalance measurement '", imbal_type,
               "'. Available are: 'A', 'alpha', 'M' and 'mu'.",
               sep = ""))
  }
}
#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{A_Index} - Calculates the 3D imbalance index "weighted integral-based 
#' centroid angle" of a 3D tree in phylo3D format using either the 
#' edge weights or the edge lengths as weights.\cr
#' If problems occur with the estimation of the integrals, try to increase 
#' the optional parameters  \cr\code{rel.tolerance} and \code{max.subdiv}.
#'
#' @return \code{A_Index} Numeric value in the interval 
#' between 0 (included) and \eqn{\pi} (excluded). A value near \eqn{\pi} 
#' indicates a higher degree and near 0 a lower degree of asymmetry.
#'
#' @export
#' @rdname integralImbalIndices
#' 
#' @examples
#' A_Index(tree, weight="edge_weight")
A_Index <- function(tree, weight = "edge_weight", rel.tolerance = 1e-10,
                    max.subdiv = 200L) {
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  if (!weight %in% c("edge_weight", "edge_length")) {
    warning(paste(
      "Unknown weighting method. The default 'edge_weight' is used",
      "instead."))
    weight <- "edge_weight"
  }
  if (sum(c("node.coord", "edge.weight") %in% attributes(tree)$names) != 2) {
    stop(paste(
      "Cannot calculate subtree centroids. The tree is missing at least one",
      "of these attributes:'node.coord', 'edge.weight'."))
  }
  if (sum(c("node.descs", "node.ancs", "node.depth", "node.subtrCentr") %in%
          attributes(tree)$names) < 4) {
    comment(paste(
      "This may take longer as at least one of the attributes,",
      "'node.descs', 'node.ancs',", "'node.depth' and 'node.subtrCentr'",
      "does not exist and has to be calculated first."))
  }
  if (!"node.descs" %in% attributes(tree)$names) {
    tree$node.descs <- getDescs(tree)
  }
  if (!"node.ancs" %in% attributes(tree)$names) {
    tree$node.ancs <- getAncs(tree)
  }
  if (!"node.depth" %in% attributes(tree)$names) {
    tree$node.depth <- getNodeDepths(tree)
  }
  if (!"node.subtrCentr" %in% attributes(tree)$names) {
    tree$node.subtrCentr <- getSubtrCentr(tree)
  }
  node_order <- tree$node.depth["orderByIncrDepth", ] # order "top-down"
  #-----------
  wei_v <- NULL # Calculate weights for the edge imbalance values.
  inc_edgeweights <- getIncEdgeWeights(tree)
  if (weight == "edge_weight") {
    wei_v <- inc_edgeweights
  } else {
    wei_v <- getIncEdgeLens(tree)
    if(sum(inc_edgeweights[node_order[-1]]==0)>0){
      warning(paste("It might not be advisable to use edge_length as",
                    "weighting method if there are edges with zero weight.\n"))
    }
  }
  #-----------
  is_leaf <- getLeaves(tree)
  centangles_e <- rep(NA, length(node_order))
  for (i in node_order[-1]) { # leave out the root (has no ancestor)
    if (is_leaf[i]) { # if leaf edge
      centangles_e[i] <- 0
    } else { # if internal edge
      if (imbalSubdiv_A(
        x = 0, p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
        v = tree$node.coord[i, ],
        centr_v = tree$node.subtrCentr[i, 1:3],
        centr_v_weight = tree$node.subtrCentr[i, 4],
        edge_weight = inc_edgeweights[i]
      ) == 0) {
        centangles_e[i] <- 0 # integral is zero iff v's node imbalance is zero
      } else { # v's node imbalance is not zero
        centangles_e[i] <- imbalInt_e(
          p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
          v = tree$node.coord[i, ],
          centr_v = tree$node.subtrCentr[i, 1:3],
          centr_v_weight = tree$node.subtrCentr[i, 4],
          edge_weight = inc_edgeweights[i],
          imbal_type = "A",
          rel.tolerance = rel.tolerance, max.subdiv = max.subdiv
        )
      }
    }
  }
  A <- sum(centangles_e * wei_v, na.rm = TRUE) / sum(wei_v, na.rm = TRUE)
  return(A)
}
#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{alpha_Index} - Calculates the 3D imbalance index "weighted 
#' integral-based minimal centroid angle" of a 3D tree in phylo3D format using 
#' either the edge weights or the edge lengths as weights.\cr
#' If problems occur with the estimation of the integrals, try to increase 
#' the optional parameters  \cr\code{rel.tolerance} and \code{max.subdiv}.
#'
#' @return \code{alpha_Index} Numeric value in the interval 
#' between 0 (included) and \eqn{\pi/2} (excluded). A value near \eqn{\pi/2} 
#' indicates a higher degree and near 0 a lower degree of asymmetry.
#'
#' @export
#' @rdname integralImbalIndices
#' 
#' @examples
#' alpha_Index(tree)
alpha_Index <- function(tree, weight = "edge_weight", rel.tolerance = 1e-10,
                        max.subdiv = 200L) {
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  if (!weight %in% c("edge_weight", "edge_length")) {
    warning(paste(
      "Unknown weighting method. The default 'edge_weight' is used",
      "instead."))
    weight <- "edge_weight"
  }
  if (sum(c("node.coord", "edge.weight") %in% attributes(tree)$names) != 2) {
    stop(paste(
      "Cannot calculate subtree centroids. The tree is missing at least one",
      "of these attributes:'node.coord', 'edge.weight'."))
  }
  if (sum(c("node.descs", "node.ancs", "node.depth", "node.subtrCentr") %in%
          attributes(tree)$names) < 4) {
    comment(paste(
      "This may take longer as at least one of the attributes,",
      "'node.descs', 'node.ancs',",
      "'node.depth' and 'node.subtrCentr' does not",
      "exist and has to be calculated first."))
  }
  if (!"node.descs" %in% attributes(tree)$names) {
    tree$node.descs <- getDescs(tree)
  }
  if (!"node.ancs" %in% attributes(tree)$names) {
    tree$node.ancs <- getAncs(tree)
  }
  if (!"node.depth" %in% attributes(tree)$names) {
    tree$node.depth <- getNodeDepths(tree)
  }
  if (!"node.subtrCentr" %in% attributes(tree)$names) {
    tree$node.subtrCentr <- getSubtrCentr(tree)
  }
  node_order <- tree$node.depth["orderByIncrDepth", ] # order "top-down"
  #-----------
  wei_v <- NULL # Calculate weights for the edge imbalance values.
  inc_edgeweights <- getIncEdgeWeights(tree)
  if (weight == "edge_weight") {
    wei_v <- inc_edgeweights
  } else {
    wei_v <- getIncEdgeLens(tree)
    if(sum(inc_edgeweights[node_order[-1]]==0)>0){
      warning(paste("It might not be advisable to use edge_length as",
                    "weighting method if there are edges with zero weight.\n"))
    }
  }
  #-----------
  is_leaf <- getLeaves(tree)
  mcentangles_e <- rep(NA, length(node_order))
  for (i in node_order[-1]) { # leave out the root (has no ancestor)
    if (is_leaf[i]) { # if leaf edge
      mcentangles_e[i] <- 0
    } else { # if internal edge
      if (imbalSubdiv_alpha(
        x = 0, p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
        v = tree$node.coord[i, ],
        centr_v = tree$node.subtrCentr[i, 1:3],
        centr_v_weight = tree$node.subtrCentr[i, 4],
        edge_weight = inc_edgeweights[i]
      ) == 0) {
        mcentangles_e[i] <- 0 # integral is zero iff v's node imbalance is zero
      } else { # v's node imbalance is not zero
        mcentangles_e[i] <- imbalInt_e(
          p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
          v = tree$node.coord[i, ],
          centr_v = tree$node.subtrCentr[i, 1:3],
          centr_v_weight = tree$node.subtrCentr[i, 4],
          edge_weight = inc_edgeweights[i],
          imbal_type = "alpha",
          rel.tolerance = rel.tolerance, max.subdiv = max.subdiv
        )
      }
    }
  }
  alpha <- sum(mcentangles_e * wei_v, na.rm = TRUE) / sum(wei_v, na.rm = TRUE)
  return(alpha)
}
#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{M_Index} - Calculates the 3D imbalance index "weighted integral-based 
#' expanded relative centroid distance" of a 3D tree in phylo3D format using 
#' either the edge weights or the edge lengths as weights.\cr
#' If problems occur with the estimation of the integrals, try to increase 
#' the optional parameters  \cr\code{rel.tolerance} and \code{max.subdiv}.
#' 
#' @return \code{M_Index} Numeric value in the interval 
#' between 0 (included) and 1 (excluded). A value near 1 indicates a higher 
#' degree and near 0 a lower degree of asymmetry.
#'
#' @export
#' @rdname integralImbalIndices
#' 
#' @examples
#' M_Index(tree, weight="edge_length")
M_Index <- function(tree, weight = "edge_weight", rel.tolerance = 1e-10,
                  max.subdiv = 200L) {
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  if (!weight %in% c("edge_weight", "edge_length")) {
    warning(paste(
      "Unknown weighting method. The default 'edge_weight' is used",
      "instead."))
    weight <- "edge_weight"
  }
  if (sum(c("node.coord", "edge.weight") %in% attributes(tree)$names) != 2) {
    stop(paste(
      "Cannot calculate subtree centroids. The tree is missing at least one",
      "of these attributes:'node.coord', 'edge.weight'."))
  }
  if (sum(c("node.descs", "node.ancs", "node.depth", "node.subtrCentr") %in%
          attributes(tree)$names) < 4) {
    comment(paste(
      "This may take longer as at least one of the attributes,",
      "'node.descs', 'node.ancs',",
      "'node.depth' and 'node.subtrCentr' does not",
      "exist and has to be calculated first."))
  }
  if (!"node.descs" %in% attributes(tree)$names) {
    tree$node.descs <- getDescs(tree)
  }
  if (!"node.ancs" %in% attributes(tree)$names) {
    tree$node.ancs <- getAncs(tree)
  }
  if (!"node.depth" %in% attributes(tree)$names) {
    tree$node.depth <- getNodeDepths(tree)
  }
  if (!"node.subtrCentr" %in% attributes(tree)$names) {
    tree$node.subtrCentr <- getSubtrCentr(tree)
  }
  node_order <- tree$node.depth["orderByIncrDepth", ] # order "top-down"
  #-----------
  wei_v <- NULL # Calculate weights for the edge imbalance values.
  inc_edgeweights <- getIncEdgeWeights(tree)
  if (weight == "edge_weight") {
    wei_v <- inc_edgeweights
  } else {
    wei_v <- getIncEdgeLens(tree)
    if(sum(inc_edgeweights[node_order[-1]]==0)>0){
      warning(paste("It might not be advisable to use edge_length as",
                    "weighting method if there are edges with zero weight.\n"))
    }
  }
  #-----------
  is_leaf <- getLeaves(tree)
  centdist_e <- rep(NA, length(node_order))
  for (i in node_order[-1]) { # leave out the root (has no ancestor)
    if (is_leaf[i]) { # if leaf edge
      centdist_e[i] <- 0
    } else { # if internal edge
      if (imbalSubdiv_M(
        x = 0,
        p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
        v = tree$node.coord[i, ],
        centr_v = tree$node.subtrCentr[i, 1:3],
        centr_v_weight = tree$node.subtrCentr[i, 4],
        edge_weight = inc_edgeweights[i]
      ) == 0) {
        centdist_e[i] <- 0 # integral is zero iff v's node imbalance is zero
      } else { # v's node imbalance is not zero
        centdist_e[i] <- imbalInt_e(
          p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
          v = tree$node.coord[i, ],
          centr_v = tree$node.subtrCentr[i, 1:3],
          centr_v_weight = tree$node.subtrCentr[i, 4],
          edge_weight = inc_edgeweights[i],
          imbal_type = "M",
          rel.tolerance = rel.tolerance, max.subdiv = max.subdiv
        )
      }
    }
  }
  M <- sum(centdist_e * wei_v, na.rm = TRUE) / sum(wei_v, na.rm = TRUE)
  return(M)
}
#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{mu_Int} - Calculates the 3D imbalance index "weighted integral-based 
#' relative centroid distance" of a 3D tree in phylo3D format using either the 
#' edge weights or the edge lengths as weights.\cr
#' If problems occur with the estimation of the integrals, try to increase 
#' the optional parameters  \cr\code{rel.tolerance} and \code{max.subdiv}. 
#' 
#' @return \code{mu_Index} Numeric value in the interval between 0 (included) 
#' and 1 (excluded). A value near 1 indicates a higher degree and near 0 a 
#' lower degree of asymmetry.
#'
#' @export
#' @rdname integralImbalIndices
#' 
#' @examples
#' mu_Index(tree, weight="edge_length")
mu_Index <- function(tree, weight = "edge_weight", rel.tolerance = 1e-10,
                     max.subdiv = 200L) {
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  if (!weight %in% c("edge_weight", "edge_length")) {
    warning(paste(
      "Unknown weighting method. The default 'edge_weight' is used",
      "instead."))
    weight <- "edge_weight"
  }
  if (sum(c("node.coord", "edge.weight") %in% attributes(tree)$names) != 2) {
    stop(paste(
      "Cannot calculate subtree centroids. The tree is missing at least one",
      "of these attributes:'node.coord', 'edge.weight'."))
  }
  if (sum(c("node.descs", "node.ancs", "node.depth", "node.subtrCentr") %in%
          attributes(tree)$names) < 4) {
    comment(paste(
      "This may take longer as at least one of the attributes,",
      "'node.descs', 'node.ancs',",
      "'node.depth' and 'node.subtrCentr' does not",
      "exist and has to be calculated first."))
  }
  if (!"node.descs" %in% attributes(tree)$names) {
    tree$node.descs <- getDescs(tree)
  }
  if (!"node.ancs" %in% attributes(tree)$names) {
    tree$node.ancs <- getAncs(tree)
  }
  if (!"node.depth" %in% attributes(tree)$names) {
    tree$node.depth <- getNodeDepths(tree)
  }
  if (!"node.subtrCentr" %in% attributes(tree)$names) {
    tree$node.subtrCentr <- getSubtrCentr(tree)
  }
  node_order <- tree$node.depth["orderByIncrDepth", ] # order "top-down"
  #-----------
  wei_v <- NULL # Calculate weights for the edge imbalance values.
  inc_edgeweights <- getIncEdgeWeights(tree)
  if (weight == "edge_weight") {
    wei_v <- inc_edgeweights
  } else {
    wei_v <- getIncEdgeLens(tree)
    if(sum(inc_edgeweights[node_order[-1]]==0)>0){
      warning(paste("It might not be advisable to use edge_length as",
                    "weighting method if there are edges with zero weight.\n"))
    }
  }
  #-----------
  is_leaf <- getLeaves(tree)
  centdist_e <- rep(NA, length(node_order))
  for (i in node_order[-1]) { # leave out the root (has no ancestor)
    if (is_leaf[i]) { # if leaf edge
      centdist_e[i] <- 0
    } else { # if internal edge
      if (imbalSubdiv_mu(
        x = 0, p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
        v = tree$node.coord[i, ],
        centr_v = tree$node.subtrCentr[i, 1:3],
        centr_v_weight = tree$node.subtrCentr[i, 4],
        edge_weight = inc_edgeweights[i]
      ) == 0) {
        centdist_e[i] <- 0 # integral is zero iff v's node imbalance is zero
      } else { # v's node imbalance is not zero
        centdist_e[i] <- imbalInt_e(
          p = tree$node.coord[tree$node.ancs["ancestor", ][i], ],
          v = tree$node.coord[i, ],
          centr_v = tree$node.subtrCentr[i, 1:3],
          centr_v_weight = tree$node.subtrCentr[i, 4],
          edge_weight = inc_edgeweights[i],
          imbal_type = "mu",
          rel.tolerance = rel.tolerance, max.subdiv = max.subdiv
        )
      }
    }
  }
  mu <- sum(centdist_e * wei_v, na.rm = TRUE) / sum(wei_v, na.rm = TRUE)
  return(mu)
}
#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{imbalInt_e} - Calculates the integral of the node imbalance values
#' over all possible subdividing nodes on an edge.
#'
#' @param v Numeric vector of size 3 (3D coordinates of node \eqn{v}).
#' @param p Numeric vector of size 3 (3D coordinates of parent node \eqn{p}).
#' @param centr_v Numeric vector of size 3 (3D coordinates of the centroid of
#' the pending subtree of node \eqn{v}).
#' @param centr_v_weight Numeric value >=0 (weight of the pending subtree of 
#' node \eqn{v}).
#' @param edge_weight Numeric value >=0 (weight of the edge \eqn{(p,v)}).
#'
#' @return \code{imbalInt_e} Numeric value (0 minimal value, higher values
#'indicate a higher degree of asymmetry).
#'
#' @export
#' @rdname integralImbalIndices
#'
#' @examples
#' imbalInt_e(p=c(0,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
#'     centr_v_weight=1,edge_weight=1,imbal_type="mu")
imbalInt_e <- function(p, v, centr_v, centr_v_weight, edge_weight,
                       imbal_type, rel.tolerance = 1e-8, max.subdiv = 200L) {
  if (sum(is.na(c(p, v, centr_v, centr_v_weight, edge_weight))) > 0 ||
      length(p) != 3 || length(v) != 3 || length(centr_v) != 3 ||
      length(centr_v_weight) != 1 || length(edge_weight) != 1 ||
      sum(p == v) == 3 ||
      edge_weight < 0 || centr_v_weight < 0 || (edge_weight+centr_v_weight)==0){
    stop(paste("Invalid input",
               ": p=", paste(p, collapse = " "),
               "; v=", paste(v, collapse = " "),
               "; centr_v=", paste(centr_v, collapse = " "),
               "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
               "; edge_weight=", paste(edge_weight, collapse = " "),
               sep = ""))
  }
  if (!imbal_type %in% c("mu", "M", "A", "alpha")) {
    stop(paste("Unknown node imbalance measurement '", imbal_type,
               "'. Available are: 'A', 'alpha', 'M' and 'mu'.",
               sep = ""))
  }
  node_imbal_fun <- function(x) {
    imbalProfile_e(
      xs = x, p, v, centr_v, centr_v_weight,
      edge_weight, imbal_type
    )
  }
  result <- tryCatch(
    stats::integrate(node_imbal_fun,
                     lower = 0, upper = 1,
                     subdivisions = max.subdiv, rel.tol = rel.tolerance
    ),
    error = function(cond) {
      message(paste("stats::integrate has severe troubles integrating. Input",
                    ": p=", paste(p, collapse = " "),
                    "; v=", paste(v, collapse = " "),
                    "; centr_v=", paste(centr_v, collapse = " "),
                    "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
                    "; edge_weight=", paste(edge_weight, collapse = " "),
                    "; imbal_type=", paste(imbal_type, collapse = " "),
                    "; rel.tolerance=", paste(rel.tolerance, collapse = " "),
                    "; max.subdiv=", paste(max.subdiv, collapse = " "),
                    ". This mostly occurs for integrals that are approx. zero.",
                    "Thus, zero is returned.",
                    sep = ""))
      message("The original error message:")
      message(cond)
      return(0)
    },
    warning = function(cond) {
      message(paste("stats::integrate has some troubles integrating. Input",
                    ": p=", paste(p, collapse = " "),
                    "; v=", paste(v, collapse = " "),
                    "; centr_v=", paste(centr_v, collapse = " "),
                    "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
                    "; edge_weight=", paste(edge_weight, collapse = " "),
                    "; imbal_type=", paste(imbal_type, collapse = " "),
                    "; rel.tolerance=", paste(rel.tolerance, collapse = " "),
                    "; max.subdiv=", paste(max.subdiv, collapse = " "),
                    ". This mostly occurs for integrals that are approx. zero.",
                    "Thus, zero is returned.",
                    sep = ""))
      message("The original warning message:")
      message(cond)
      return(0)
    }
  )
  if (is.na(result$value)) {
    stop(paste("A problem occured when trying to integrate using this input",
               ": p=", paste(p, collapse = " "),
               "; v=", paste(v, collapse = " "),
               "; centr_v=", paste(centr_v, collapse = " "),
               "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
               "; edge_weight=", paste(edge_weight, collapse = " "),
               "; imbal_type=", paste(imbal_type, collapse = " "),
               sep = ""))
  } else {
    return(result$value)
  }
}
#' Calculation of integral-based distance- and angle-based 3D imbalance indices
#'
#' \code{imbalProfile_e} - Calculates the node imbalance values
#' for a given set of subdivisions of an edge.
#'
#' @param xs Numeric vector with values between 0 (included) and 1 (excluded).
#' Set of edge subdivisions.
#'
#' @return \code{imbalProfile_e} Numeric vector of imbalance values (0 minimal 
#' value, higher values indicate a higher degree of asymmetry) for the edge
#' subdivisions indicated by input \code{xs}.
#'
#' @export
#' @rdname integralImbalIndices
#'
#' @examples
#' imbalProfile_e(xs=c(0,0.2,0.4),p=c(1,1,0),v=c(0,0,0),centr_v=c(0.5,0,0),
#'     centr_v_weight=1,edge_weight=1,imbal_type="A")
imbalProfile_e <- function(xs, p, v, centr_v, centr_v_weight, edge_weight,
                           imbal_type) {
  if (sum(is.na(c(p, v, centr_v, centr_v_weight, edge_weight))) > 0 ||
      length(p) != 3 || length(v) != 3 || length(centr_v) != 3 ||
      length(centr_v_weight) != 1 || length(edge_weight) != 1 ||
      sum(p == v) == 3 || min(xs) < 0 || max(xs) >= 1 ||
      edge_weight < 0 || centr_v_weight < 0 || (edge_weight+centr_v_weight)==0){
    stop(paste("Invalid input",
               ": xs=", paste(xs, collapse = " "),
               "; p=", paste(p, collapse = " "),
               "; v=", paste(v, collapse = " "),
               "; centr_v=", paste(centr_v, collapse = " "),
               "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
               "; edge_weight=", paste(edge_weight, collapse = " "),
               sep = ""))
  }
  node_imbal_vals <- NULL
  if (imbal_type == "mu") {
    node_imbal_vals <- sapply(xs, function(my_x) {
      imbalSubdiv_mu(my_x, p, v, centr_v, centr_v_weight, edge_weight)
    })
  } else if (imbal_type == "M") {
    node_imbal_vals <- sapply(xs, function(my_x) {
      imbalSubdiv_M(my_x, p, v, centr_v, centr_v_weight, edge_weight)
    })
  } else if (imbal_type == "A") {
    node_imbal_vals <- sapply(xs, function(my_x) {
      imbalSubdiv_A(my_x, p, v, centr_v, centr_v_weight, edge_weight)
    })
  } else if (imbal_type == "alpha") {
    node_imbal_vals <- sapply(xs, function(my_x) {
      imbalSubdiv_alpha(my_x, p, v, centr_v, centr_v_weight, edge_weight)
    })
  } else {
    stop(paste("Unknown node imbalance measurement '", imbal_type,
               "'. Available are: 'A', 'alpha', 'M' and 'mu'.",
               sep = ""))
  }
  if (sum(is.na(node_imbal_vals)) > 0) {
    stop(paste("A problem occured when trying to integrate using this input",
               ": p=", paste(p, collapse = " "),
               "; v=", paste(v, collapse = " "),
               "; centr_v=", paste(centr_v, collapse = " "),
               "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
               "; edge_weight=", paste(edge_weight, collapse = " "),
               "; imbal_type=", paste(imbal_type, collapse = " "),
               sep = ""))
  } else {
    return(node_imbal_vals)
  }
}
