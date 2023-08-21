#' Calculation of the centroid distances
#'
#' \code{imbalSubdiv_mu} - Calculates the node imbalance value "relative
#' centroid distance" of a vertex which subdivides the edge \eqn{(p,v)} at 
#' \eqn{v+x \cdot (p-v)}{v+x*(p-v)} with \eqn{x \in [0,1]}. For example, we 
#' can obtain the node imbalance value of \eqn{v} if \eqn{x=0}, and 
#' \eqn{x=0.5} would indicate a subdividing node exactly in the middle of 
#' \eqn{v} and \eqn{p}. \cr
#' Attention: If \eqn{x=1}, this function will not calculate the node imbalance
#' value of \eqn{p} with respect to its own incoming edge but with respect to 
#' the edge \eqn{(p,v)} itself. This enables us to estimate the 
#' node imbalance integrals over the entire edge length.
#'
#' @param x Numeric value \eqn{\in [0,1]} which indicates where on the 
#' edge \eqn{(p,v)} the subdivision takes place: 
#' \eqn{v+x \cdot (p-v)}{v+x*(p-v)}.
#' @param v Numeric vector of size 3 (3D coordinates of node \eqn{v}).
#' @param p Numeric vector of size 3 (3D coordinates of parent node \eqn{p}).
#' @param centr_v Numeric vector of size 3 (3D coordinates of the centroid of
#' the pending subtree of node \eqn{v}).
#' @param centr_v_weight Numeric value >=0 (weight of the pending subtree of 
#' node \eqn{v}).
#' @param edge_weight Numeric value >=0 (weight of the edge \eqn{(p,v)}).
#' 
#' @return \code{imbalSubdiv_mu} Numeric value \eqn{\in [0,1]} (higher values
#' indicate a higher degree of asymmetry).
#'
#' @author Sophie Kersting, Luise Kühn
#' 
#' @export
#' @rdname centroidDistances
#'
#' @examples
#' imbalSubdiv_mu(x=0.5,p=c(1,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
#' centr_v_weight=1,edge_weight=1)
imbalSubdiv_mu <- function(x, p, v, centr_v, centr_v_weight, edge_weight) {
  if (sum(is.na(c(x, p, v, centr_v, centr_v_weight, edge_weight))) > 0 ||
      length(p) != 3 || length(v) != 3 || length(centr_v) != 3 ||
      length(centr_v_weight) != 1 || length(edge_weight) != 1 || 
      length(x) != 1 || x < 0 || x > 1 || sum(p == v) == 3|| 
      edge_weight < 0 || centr_v_weight < 0 || (edge_weight+centr_v_weight)==0){
    stop(paste("Invalid input",
               ": x=", paste(x, collapse = " "),
               "; p=", paste(p, collapse = " "),
               "; v=", paste(v, collapse = " "),
               "; centr_v=", paste(centr_v, collapse = " "),
               "; centr_v_weight=", paste(centr_v_weight, collapse = " "),
               "; edge_weight=", paste(edge_weight, collapse = " "),
               sep = ""))
  }
  # Determine the coordinates of the subdividing node x_subdiv and the
  # centroid of its pending subtree.
  x_subdiv <- v + x * (p - v)
  centr_x <- NULL # Check if leaf edge, i.e. if zero weight of T_v.
  if(centr_v_weight == 0){ # is leaf edge
    return(0)
  } else { # is inner edge
    centr_x <- (centr_v * centr_v_weight + (v + x/2 *(p-v))* x *edge_weight) /
      (centr_v_weight + x * edge_weight)
  }
  # Calculate the respective node imbalance value.
  maxcentdist_x <- sqrt(sum((x_subdiv - centr_x)^2))
  if (maxcentdist_x == 0) {
    return(0)
  } else {
    centdist_x <- dist3dToLine(point = centr_x, a = p, b = v)
    return(centdist_x / maxcentdist_x)
  }
}
#' Calculation of the centroid distances
#'
#' \code{imbalSubdiv_M} - Calculates the node imbalance value "expanded
#' relative centroid distance" of a vertex which subdivides the edge \eqn{(p,v)} 
#' at \eqn{v+x \cdot (p-v)}{v+x*(p-v)} with \eqn{x \in [0,1]}. For example, we 
#' can obtain the node imbalance value of \eqn{v} if \eqn{x=0}, and \eqn{x=0.5} 
#' would indicate a subdividing node exactly in the middle of \eqn{v} and 
#' \eqn{p}. \cr
#' Attention: If \eqn{x=1}, this function will not calculate the node imbalance
#' value of \eqn{p} with respect to its own incoming edge but with respect to 
#' the edge \eqn{(p,v)} itself. This enables us to estimate the 
#' node imbalance integrals over the entire edge length.
#' 
#' @return \code{imbalSubdiv_M} Numeric value \eqn{\in [0,2]} (higher values
#' indicate a higher degree of asymmetry).
#'
#' @export
#' @rdname centroidDistances
#'
#' @examples
#' imbalSubdiv_M(x=0.5,p=c(1,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
#' centr_v_weight=1,edge_weight=1)
imbalSubdiv_M <- function(x, p, v, centr_v, centr_v_weight, edge_weight) {
  centangle_x <- imbalSubdiv_A(x, p, v, centr_v, centr_v_weight, edge_weight)
  centreldist_x <- imbalSubdiv_mu(x, p, v, centr_v, centr_v_weight, edge_weight)
  if (centangle_x > (pi / 2)) {
    return(2 - centreldist_x)
  } else {
    return(centreldist_x)
  }
}
#' Calculation of the centroid distances
#'
#' \code{dist3dToLine} - Calculates the distance of a \eqn{point} 
#' to the infinite line between two points \eqn{a} and \eqn{b} in 3D space.
#'
#' @param point Numeric vector of size 3 (e.g. 3D coordinates).
#' @param a Numeric vector of size 3 (e.g. 3D coordinates).
#' @param b Numeric vector of size 3 (e.g. 3D coordinates).
#'
#' @return \code{dist3dToLine} Numeric value. 
#'
#' @export
#' @rdname centroidDistances
#'
#' @examples
#' dist3dToLine(point=c(1,1,1),a=c(0,0,0),b=c(1,2,2)) # 0.47140...
dist3dToLine <- function(point, a, b) {
  if (sum(is.na(c(point, a, b))) > 0 || length(point) != 3 || length(a) != 3 ||
      length(b) != 3) {
    stop(paste("Invalid 3D coordinates: point=",
               paste(point, collapse = " "), "; a=",
               paste(a, collapse = " "), "; b=",
               paste(b, collapse = " "),
               sep = "" ))
  }
  # Calculate distance to infinite line through a and b.
  dist_to_line <- sqrt(sum(cross3d_prod(point - a, b - a)^2)) /
    sqrt(sum((b - a)^2))
  return(dist_to_line)
}
#' Calculation of the centroid distances
#'
#' \code{cross3d_prod} - Calculates the cross-product of two 3D vectors.
#'
#' @return \code{cross3d_prod} Numeric vector of size 3.
#'
#' @export
#' @rdname centroidDistances
#'
#' @examples
#' cross3d_prod(a=c(1,-1,1),b=c(1,2,2)) # c(-4, -1, 3)
cross3d_prod <- function(a, b) {
  return(c( a[2] * b[3] - a[3] * b[2],
            a[3] * b[1] - a[1] * b[3],
            a[1] * b[2] - a[2] * b[1]))
}