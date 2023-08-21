#' Calculation of the centroid angles
#'
#' \code{imbalSubdiv_A} - Calculates the node imbalance value "centroid
#' angle" of a vertex which subdivides the edge \eqn{(p,v)} at 
#' \eqn{v+x \cdot (p-v)}{v+x*(p-v)} with \eqn{x \in [0,1]}. For example, 
#' we can obtain the node imbalance value of \eqn{v} if \eqn{x=0}, and 
#' \eqn{x=0.5} would indicate a subdividing node exactly in the middle of 
#' \eqn{v} and \eqn{p}. \cr
#' Attention: If \eqn{x=1}, this function will not calculate the node imbalance
#' value of \eqn{p} with respect to its incoming edge but with respect to the
#' edge \eqn{(p,v)} itself. This enables us to estimate the 
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
#' @return \code{imbalSubdiv_A} Numeric value \eqn{\in [0,\pi]} (higher values
#' indicate a higher degree of asymmetry).
#'
#' @author Sophie Kersting, Luise Kühn
#'
#' @export
#' @rdname centroidAngles
#'
#' @examples
#' imbalSubdiv_A(x=0.5,p=c(1,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
#' centr_v_weight=1,edge_weight=1)
imbalSubdiv_A <- function(x, p, v, centr_v, centr_v_weight, edge_weight) {
  if (sum(is.na(c(x, p, v, centr_v, centr_v_weight, edge_weight))) > 0 ||
      length(p) != 3 || length(v) != 3 || length(centr_v) != 3 ||
      length(centr_v_weight) != 1 || length(edge_weight) != 1 || 
      length(x) != 1 || x < 0 || x > 1 || sum(p == v) == 3 ||
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
  centr_x <- (centr_v * centr_v_weight + (v + x/2 * (p-v)) * x * edge_weight) /
    (centr_v_weight + x * edge_weight)
  # Calculate the respective node imbalance value.
  centangle_x <- angle3dVec(a = v - p, b = centr_x - x_subdiv)
  return(centangle_x)
}
#' Calculation of the centroid angles
#'
#' \code{imbalSubdiv_alpha} - Calculates the node imbalance value "minimal
#' centroid angle" of a vertex which subdivides the edge \eqn{(p,v)} at 
#' \eqn{v+x \cdot (p-v)}{v+x*(p-v)} with \eqn{x \in [0,1]}. For example, 
#' we can obtain the node imbalance value of \eqn{v} if \eqn{x=0}, and 
#' \eqn{x=0.5} would indicate a subdividing node exactly in the middle of 
#' \eqn{v} and \eqn{p}. \cr
#' Attention: If \eqn{x=1}, this function will not calculate the node imbalance
#' value of \eqn{p} with respect to its incoming edge but with respect to the
#' edge \eqn{(p,v)} itself. This enables us to estimate the 
#' node imbalance integrals over the entire edge length.
#' 
#' @return \code{imbalSubdiv_alpha} Numeric value \eqn{\in [0,\pi/2]} 
#' (higher values indicate a higher degree of asymmetry).
#'
#' @export
#' @rdname centroidAngles
#'
#' @examples
#' imbalSubdiv_alpha(x=0.5,p=c(1,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
#' centr_v_weight=1,edge_weight=1)
imbalSubdiv_alpha <- function(x, p, v, centr_v, centr_v_weight, edge_weight) {
  centangle_x <- imbalSubdiv_A(x, p, v, centr_v, centr_v_weight, edge_weight)
  if (centangle_x > (pi / 2)) {
    return(pi - centangle_x)
  } else {
    return(centangle_x)
  }
}
#' Calculation of the centroid angles
#'
#' \code{angle3dVec} - Calculates the angle in the interval \eqn{[0,\pi]} 
#' between two 3D vectors \eqn{a} and \eqn{b}.
#' Note that the function returns 0 if one entry vector is \eqn{(0,0,0)}.
#'
#' @param a Numeric vector of size 3 (e.g., 3D coordinates).
#' @param b Numeric vector of size 3 (e.g., 3D coordinates).
#'
#' @return \code{angle3dVec} Numeric value in \eqn{[0,\pi]}.
#'
#' @export
#' @rdname centroidAngles
#'
#' @examples
#' angle3dVec(a=c(1,0,0),b=c(0,1,0)) # right angle = pi/2 = 1.5707...
angle3dVec <- function(a, b) {
  if (sum(is.na(c(a, b))) > 0 || length(a) != 3 ||
      length(b) != 3) {
    stop(paste("Invalid 3D coordinates: a=",
               paste(a, collapse = " "), "; b=",
               paste(b, collapse = " "),
               sep = ""))
  }
  if (sum(a == c(0, 0, 0)) == 3 || sum(b == c(0, 0, 0)) == 3) {
    return(0)
  }
  # Calculate angle between the two vectors using arccos.
  angle <- acos(sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))))
  if (is.na(angle)) {
    stop(paste("Cannot calculate angle for a=",
               paste(a, collapse = " "), "; b=",
               paste(b, collapse = " "),
               sep = ""))
  }
  return(angle)
}