#' Turn a rooted 3D tree into one of its balanced versions
#'
#' \code{makePhylo3DBalanced} - Creates a balanced version of a rooted 3D tree
#' in phylo3D format. From the leaves of lowest depth to the root, each node
#' is turned into a balanced node by rotating its pending subtree such that
#' it is in line with its incoming edge. The edge and subtree lengths and 
#' weights as well as the angles between the outgoing edges of a single node 
#' (the sister edges) are left intact. \cr
#' Note that this function yields only one of many possible balanced 
#' version of the given rooted 3D tree (most often the minimal tree under the 
#' aforementioned requirements is not unique).
#'
#' @author Sophie Kersting, Luise Kühn
#' @param tree A rooted tree in phylo3D format (no special node enumeration 
#' required, except that nodes are numbered from 1 to |V| = the total number of
#' nodes). There must be at least 2 nodes, i.e., one edge. The attributes 
#' 'node.coord' and 'edge.weight' are strictly required.
#' @return \code{makePhylo3DBalanced} Tree in phylo3D format which is balanced
#' with regards to all four node imbalance measurements and their 
#' corresponding imbalance indices.
#' @export
#' @rdname makeBalanced
#' @examples
#' tree <- treeDbalance::extendPhylo(treeDbalance::example3Dtrees$bean22)
#' tree_bal <- makePhylo3DBalanced(tree)
makePhylo3DBalanced <- function(tree) {
  if (!inherits(tree, "phylo") && !inherits(tree, "phylo3D")) {
    stop("The input tree must have class phylo or phylo3D.")
  }
  if (sum(c("node.coord", "edge.weight") %in% attributes(tree)$names) != 2) {
    stop(paste(
      "Cannot calculate subtree centroids. The tree is missing",
      "these attributes:'node.coord', 'edge.weight'."
    ))
  }
  if (sum(c("node.descs", "node.ancs", "node.depth", "node.subtrCentr") %in%
          attributes(tree)$names) < 4) {
    comment(paste(
      "This may take longer as at least one of the attributes,",
      "'node.descs', 'node.ancs',",
      "'node.depth' and 'node.subtrCentr' does not",
      "exist and has to be calculated first."
    ))
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
  if (!"edge.length" %in% attributes(tree)$names) {
    tree$edge.length <- sapply(1:nrow(tree$edge), function(x) {
      stats::dist(tree$node.coord[tree$edge[x, ], ])
    })
  }
  node_order <- rev(tree$node.depth["orderByIncrDepth", ]) # order "bottom-up"
  for (i in node_order) {
    descs_i <- getDescendants(tree = tree, node = i)
    anc_i <- tree$node.ancs["ancestor", i]
    # Ignore leaves and the root:
    if (!is.null(descs_i) && !is.na(anc_i)) {
      centr_i <- tree$node.subtrCentr[i, 1:3]
      # Check if the vertex is not already balanced:
      if (imbalSubdiv_A(
        x = 0, p = tree$node.coord[anc_i, ],
        v = tree$node.coord[i, ],
        centr_v = centr_i,
        centr_v_weight = tree$node.subtrCentr[i, 4],
        edge_weight = tree$node.ancs["inc_edge", i]
      ) != 0) {
        # Calculate where a balanced centroid should be:
        centr_dist <- unname(stats::dist(rbind(
          tree$node.coord[i, ],
          centr_i)))
        bal_centr_i <- tree$node.coord[i, ] -
          (tree$node.coord[anc_i, ] - tree$node.coord[i, ]) * centr_dist /
          tree$edge.length[tree$node.ancs["inc_edge", i]]
        # Calculate the rotation axis (orthogonal to vectors from node i to
        # its old and new centroid).
        axis_i <- tryCatch(
          solve(
            rbind(
              c(bal_centr_i[2], bal_centr_i[3]),
              c(centr_i[2], centr_i[3])),
            c(-bal_centr_i[1], -centr_i[1])),
          # An error is thrown if there is no computational difference between
          # the old and new centroid because the linear equation system is not
          # uniquely solvable.
          error = function(cond) {
            message(paste(
              "Node", i, "could not be made completely balanced",
              "because it was - computationally - already ",
              "perfectly balanced." ))
            message("The original error message:")
            message(cond)
            return(NA)
          },
          warning = function(cond) {
            message(paste(
              "Node", i, "could not be made completely balanced",
              "because it was - computationally - already ",
              "perfectly balanced."))
            message("The original warning message:")
            message(cond)
            return(NA)
          }
        )
        if (is.na(axis_i[1])) { # skip to the next if there was an error
          next()
        }
        axis_i <- c(1, axis_i)
        axis_i <- axis_i / sqrt(sum(axis_i^2))
        # Calculate the rotation angle.
        angle_i <- angle3dVec(
          bal_centr_i - tree$node.coord[i, ],
          centr_i - tree$node.coord[i, ]
        )
        # Test signum of angle_i:
        if (sum((rotate3dVec(angle = angle_i, axis = axis_i,
                             vec = centr_i - tree$node.coord[i, ]) +
                 tree$node.coord[i, ] - bal_centr_i)^2) >
            sum((rotate3dVec(angle = -angle_i, axis = axis_i,
                             vec = centr_i - tree$node.coord[i, ]) +
                 tree$node.coord[i, ] - bal_centr_i)^2)) {
          angle_i <- -angle_i
        }
        # Calculate the new coordinates for all descendants.
        for (dec in descs_i) {
          tree$node.coord[dec, ] <- rotate3dVec(angle = angle_i, axis = axis_i,
                         vec = tree$node.coord[dec, ] - tree$node.coord[i, ]) +
            tree$node.coord[i, ]
        }
        # Compute the new subtree centroids.
        tree$node.subtrCentr <- getSubtrCentr(tree)
      }
    }
  }
  return(tree)
}
#' Turn a rooted 3D tree into one of its balanced versions
#'
#' \code{rotate3dVec} - Rotates a vector in 3D space for a given angle and 
#' rotation axis.
#'
#' @param angle Angle for the rotation.
#' @param axis Rotation axis.
#' @param vec Numeric vector of size 3 (3D coordinates of the
#' vector that shall be rotated).
#' @return \code{rotate3dVec} Numeric vector of size 3 (3D coordinates of the
#' rotated vector).
#' @export
#' @rdname makeBalanced
#' @examples
#' rotate3dVec(angle = pi/2, axis = c(0,-1,0), vec = c(5,0,0)) # approx. (0,0,5)
#' round(rotate3dVec(angle = pi/2, axis = c(0,-1,0), vec = c(5,0,0)),15)
rotate3dVec <- function(angle, axis, vec) {
  length_axis <- sqrt(sum(axis^2))
  if (sum(is.na(c(angle, axis, vec))) > 0 || 
      length(axis) != 3 || length(vec) != 3 || 
      length(angle) != 1 || length_axis == 0) {
    stop(paste("Invalid input: angle=", paste(angle, collapse = " "),
               "; axis=", paste(axis, collapse = " "),
               "; vec=", paste(vec, collapse = " "),
               sep = ""))
  }
  if (length_axis != 1) { # make the axis a unit vector if necessary
    axis <- axis / length_axis
  }
  vec_rot <- axis * sum(axis * vec) +
    cos(angle) * cross3d_prod(cross3d_prod(axis, vec), axis) +
    sin(angle) * cross3d_prod(axis, vec)
  return(vec_rot)
}
