% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPhylo3D.R
\name{plotPhylo3D}
\alias{plotPhylo3D}
\alias{addPhylo3D}
\alias{plotImbalPhylo3D}
\alias{addImbalPhylo3D}
\title{Plot a phylo3D object}
\usage{
plotPhylo3D(tree, show_node_enum = FALSE)

addPhylo3D(tree, offset = c(0, 0, 0), show_node_enum = FALSE)

plotImbalPhylo3D(
  tree,
  imbal_type,
  max.seclen,
  color.imbal = TRUE,
  show.gradient = FALSE
)

addImbalPhylo3D(
  tree,
  offset = c(0, 0, 0),
  imbal_type,
  max.seclen,
  color.imbal = TRUE,
  show.gradient = FALSE
)
}
\arguments{
\item{tree}{A rooted tree in phylo3D format (no special node enumeration 
required, except that nodes are numbered from 1 to |V| = the total number of
nodes). There must be at least 2 nodes, i.e., one edge. The attributes 
'node.coord' and 'edge.weight' are strictly required.\cr
Optional: Add attribute 'edge.type' and/or 'edge.color' (character vectors 
e.g., c("d","b","l","l") and c("red","blue", "red", "gray30"))
to the pylo3D object to plot the edges differently.\cr
The \eqn{i}-th entry of 'edge.type' has to contain the type of the 
\eqn{i}-th edge. The following types are available: "d"= default cylinder, 
"b" = bud, and "l"=leaf. If 'edge.type' is not given, every edge is by
default visualized as a cylinder.\cr
Similar for the attribute 'edge.color'. If it is not given, every edge is 
depicted as "gray30" by default.\cr
Furthermore, the attribute 'edge.diam' can also be added optionally to
set the edge diameters. Otherwise, the edge diameters will be calculated 
based on the edge lengths and the edge weights (treated as volume).}

\item{show_node_enum}{A boolean value (default FALSE). If true, each node
of the visualized phylo3D object is marked with its number. This helps to
identify specific nodes and edges.}

\item{offset}{Numeric vector of length 3, contains 3D coordinates by which 
the phylo object should be shifted (default = c(0,0,0), i.e., no shift).}

\item{imbal_type}{Specifies which node imbalance measurement should be used.
Available are:\cr
"A"     - centroid angle\cr
"alpha" - minimal centroid angle\cr
"M"     - expanded relative centroid distance\cr
"mu"    - relative centroid distance}

\item{max.seclen}{Numeric value >0 that specifies the maximal section length.}

\item{color.imbal}{Boolean value (default TRUE). If true, colors are 
used to depict the imbalance. Otherwise, a grayscale image is produced.}

\item{show.gradient}{Boolean value (default FALSE). If true 
the color or grayscale gradient is depicted.}
}
\value{
\code{plotPhylo3D} No return value, called for side effects 
(plotting).

\code{addPhylo3D} No return value, called for side effects 
(plotting).

\code{plotImbalPhylo3D} No return value, called for side effects 
(plotting).

\code{addImbalPhylo3D} No return value, called for side effects 
(plotting).
}
\description{
\code{plotPhylo3D} - Plots a phylo3D object using functions of the 
package 'rgl'.

\code{addPhylo3D} - This function plots a phylo3D object without any 
coordinate axis or adds the tree to an existing 
plot (e.g., for plotPhylo3D).

\code{plotImbalPhylo3D} - Plots a phylo3D object using functions of the 
package 'rgl'. Moreover, it uses either brightness or a color scale to
indicate the imbalance.\cr
Edge sections are shown darker or red with higher degree of imbalance and
brighter or cyan if they are balanced. 
This function does not use the parameter \code{edge.color}\cr
Attention: Edges of type 'bud' or 'leaf' will always be depicted as 
balanced, because they should represent leaf edges that are by definition
always balanced.

\code{addImbalPhylo3D} - This function plots a phylo3D object without any 
coordinate axis or adds the tree to an existing plot (e.g., for 
\code{plotImbalPhylo3D}). Moreover, it uses either brightness or a color 
scale to indicate the imbalance.\cr
Edge sections are shown darker or red with higher degree of imbalance and
brighter or cyan if they are balanced. 
This function does not use the parameter \code{edge.color}.\cr
Attention: Edges of type 'bud' or 'leaf' will always be depicted as 
balanced because they should represent leaf edges that are always 
balanced.
}
\examples{
tree <- list(edge = matrix(c(1,2, 2,3), byrow = TRUE, ncol = 2),
             tip.label = "", Nnode = 2,
             node.coord = matrix(c(0,0,0, 1,2,2, 3,1,3), byrow = TRUE, 
                                 ncol = 3),
             edge.weight = c(0.5, 0.25), edge.type = c("d","l"))
class(tree) <- "phylo3D"
# Alternatively try: tree <- treeDbalance::example3Dtrees$bean09
plotPhylo3D(tree, show_node_enum = FALSE)
addPhylo3D(tree, offset = c(1,1,0)) 
plotImbalPhylo3D(tree, imbal_type="mu", max.seclen=0.5, color.imbal=TRUE,
    show.gradient=FALSE)
addImbalPhylo3D(tree, imbal_type="mu", offset = c(1,0,0), max.seclen=0.5, 
                color.imbal=FALSE, show.gradient = FALSE) 
}
\author{
Sophie Kersting
}
