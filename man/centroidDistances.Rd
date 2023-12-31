% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centroidDistances.R
\name{imbalSubdiv_mu}
\alias{imbalSubdiv_mu}
\alias{imbalSubdiv_M}
\alias{dist3dToLine}
\alias{cross3d_prod}
\title{Calculation of the centroid distances}
\usage{
imbalSubdiv_mu(x, p, v, centr_v, centr_v_weight, edge_weight)

imbalSubdiv_M(x, p, v, centr_v, centr_v_weight, edge_weight)

dist3dToLine(point, a, b)

cross3d_prod(a, b)
}
\arguments{
\item{x}{Numeric value \eqn{\in [0,1]} which indicates where on the 
edge \eqn{(p,v)} the subdivision takes place: 
\eqn{v+x \cdot (p-v)}{v+x*(p-v)}.}

\item{p}{Numeric vector of size 3 (3D coordinates of parent node \eqn{p}).}

\item{v}{Numeric vector of size 3 (3D coordinates of node \eqn{v}).}

\item{centr_v}{Numeric vector of size 3 (3D coordinates of the centroid of
the pending subtree of node \eqn{v}).}

\item{centr_v_weight}{Numeric value >=0 (weight of the pending subtree of 
node \eqn{v}).}

\item{edge_weight}{Numeric value >=0 (weight of the edge \eqn{(p,v)}).}

\item{point}{Numeric vector of size 3 (e.g. 3D coordinates).}

\item{a}{Numeric vector of size 3 (e.g. 3D coordinates).}

\item{b}{Numeric vector of size 3 (e.g. 3D coordinates).}
}
\value{
\code{imbalSubdiv_mu} Numeric value \eqn{\in [0,1]} (higher values
indicate a higher degree of asymmetry).

\code{imbalSubdiv_M} Numeric value \eqn{\in [0,2]} (higher values
indicate a higher degree of asymmetry).

\code{dist3dToLine} Numeric value.

\code{cross3d_prod} Numeric vector of size 3.
}
\description{
\code{imbalSubdiv_mu} - Calculates the node imbalance value "relative
centroid distance" of a vertex which subdivides the edge \eqn{(p,v)} at 
\eqn{v+x \cdot (p-v)}{v+x*(p-v)} with \eqn{x \in [0,1]}. For example, we 
can obtain the node imbalance value of \eqn{v} if \eqn{x=0}, and 
\eqn{x=0.5} would indicate a subdividing node exactly in the middle of 
\eqn{v} and \eqn{p}. \cr
Attention: If \eqn{x=1}, this function will not calculate the node imbalance
value of \eqn{p} with respect to its own incoming edge but with respect to 
the edge \eqn{(p,v)} itself. This enables us to estimate the 
node imbalance integrals over the entire edge length.

\code{imbalSubdiv_M} - Calculates the node imbalance value "expanded
relative centroid distance" of a vertex which subdivides the edge \eqn{(p,v)} 
at \eqn{v+x \cdot (p-v)}{v+x*(p-v)} with \eqn{x \in [0,1]}. For example, we 
can obtain the node imbalance value of \eqn{v} if \eqn{x=0}, and \eqn{x=0.5} 
would indicate a subdividing node exactly in the middle of \eqn{v} and 
\eqn{p}. \cr
Attention: If \eqn{x=1}, this function will not calculate the node imbalance
value of \eqn{p} with respect to its own incoming edge but with respect to 
the edge \eqn{(p,v)} itself. This enables us to estimate the 
node imbalance integrals over the entire edge length.

\code{dist3dToLine} - Calculates the distance of a \eqn{point} 
to the infinite line between two points \eqn{a} and \eqn{b} in 3D space.

\code{cross3d_prod} - Calculates the cross-product of two 3D vectors.
}
\examples{
imbalSubdiv_mu(x=0.5,p=c(1,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
centr_v_weight=1,edge_weight=1)
imbalSubdiv_M(x=0.5,p=c(1,0,1),v=c(0,0,0),centr_v=c(0.5,0,0),
centr_v_weight=1,edge_weight=1)
dist3dToLine(point=c(1,1,1),a=c(0,0,0),b=c(1,2,2)) # 0.47140...
cross3d_prod(a=c(1,-1,1),b=c(1,2,2)) # c(-4, -1, 3)
}
\author{
Sophie Kersting, Luise Kühn
}
