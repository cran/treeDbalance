#' Plot a phylo3D object
#'
#' \code{plotPhylo3D} - Plots a phylo3D object using functions of the 
#' package 'rgl'.
#'
#' @author Sophie Kersting
#' @param tree A rooted tree in phylo3D format (no special node enumeration 
#' required, except that nodes are numbered from 1 to |V| = the total number of
#' nodes). There must be at least 2 nodes, i.e., one edge. The attributes 
#' 'node.coord' and 'edge.weight' are strictly required.\cr
#' Optional: Add attribute 'edge.type' and/or 'edge.color' (character vectors 
#' e.g., c("d","b","l","l") and c("red","blue", "red", "gray30"))
#' to the pylo3D object to plot the edges differently.\cr
#' The \eqn{i}-th entry of 'edge.type' has to contain the type of the 
#' \eqn{i}-th edge. The following types are available: "d"= default cylinder, 
#' "b" = bud, and "l"=leaf. If 'edge.type' is not given, every edge is by
#' default visualized as a cylinder.\cr
#' Similar for the attribute 'edge.color'. If it is not given, every edge is 
#' depicted as "gray30" by default.\cr
#' Furthermore, the attribute 'edge.diam' can also be added optionally to
#' set the edge diameters. Otherwise, the edge diameters will be calculated 
#' based on the edge lengths and the edge weights (treated as volume).
#' @param show_node_enum A boolean value (default FALSE). If true, each node
#' of the visualized phylo3D object is marked with its number. This helps to
#' identify specific nodes and edges.
#' 
#' @return \code{plotPhylo3D} No return value, called for side effects 
#' (plotting).
#' 
#' @export
#' @rdname plotPhylo3D
#' @examples
#' tree <- list(edge = matrix(c(1,2, 2,3), byrow = TRUE, ncol = 2),
#'              tip.label = "", Nnode = 2,
#'              node.coord = matrix(c(0,0,0, 1,2,2, 3,1,3), byrow = TRUE, 
#'                                  ncol = 3),
#'              edge.weight = c(0.5, 0.25), edge.type = c("d","l"))
#' class(tree) <- "phylo3D"
#' # Alternatively try: tree <- treeDbalance::example3Dtrees$bean09
#' plotPhylo3D(tree, show_node_enum = FALSE)
plotPhylo3D <- function(tree, show_node_enum = FALSE){
  # Measure the tree to create the box for the plot.
  square_width <- max(max(tree$node.coord[,1]) - min(tree$node.coord[,1]),
                      max(tree$node.coord[,2]) - min(tree$node.coord[,2]),
                      max(tree$node.coord[,3]) - min(tree$node.coord[,3]))
  my_epsilon <- square_width/20
  square_mids <- c((max(tree$node.coord[,1]) + min(tree$node.coord[,1]))/2,
                   (max(tree$node.coord[,2]) + min(tree$node.coord[,2]))/2,
                   (max(tree$node.coord[,3]) + min(tree$node.coord[,3]))/2)
  # Construct basic plot.
  rgl::plot3d(1,1,1,type = "n", 
              xlab = "x", ylab = "y", zlab = "z",
              xlim = c(square_mids[1]-square_width/2-my_epsilon,
                       square_mids[1]+square_width/2+my_epsilon), 
              ylim = c(square_mids[2]-square_width/2-my_epsilon,
                       square_mids[2]+square_width/2+my_epsilon), 
              zlim = c(square_mids[3]-square_width/2-my_epsilon,
                       square_mids[3]+square_width/2+my_epsilon),
              surface=FALSE, axis.scales = FALSE, aspect = TRUE)
  addPhylo3D(tree, offset = c(0,0,0), show_node_enum)
}
#' Plot a phylo3D object
#'
#' \code{addPhylo3D} - This function plots a phylo3D object without any 
#' coordinate axis or adds the tree to an existing 
#' plot (e.g., for plotPhylo3D).
#'
#' @param offset Numeric vector of length 3, contains 3D coordinates by which 
#' the phylo object should be shifted (default = c(0,0,0), i.e., no shift).
#' 
#' @return \code{addPhylo3D} No return value, called for side effects 
#' (plotting).
#' 
#' @export
#' @rdname plotPhylo3D
#' @examples
#' addPhylo3D(tree, offset = c(1,1,0)) 
addPhylo3D <- function(tree, offset = c(0,0,0), show_node_enum = FALSE){
  if(sum(c("node.coord") %in% attributes(tree)$names)!=1){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "this attribute: 'node.coord'."))
  }
  tree$node.coord <- tree$node.coord +  # shift all coordinates by offset
    outer(rep.int(1L, nrow(tree$node.coord)), offset)
  
  if (!"edge.diam" %in% attributes(tree)$names){
    edge.len <- sapply(1:nrow(tree$edge), function(x) {
      sqrt(sum((tree$node.coord[tree$edge[x,1],]-
                  tree$node.coord[tree$edge[x,2],])^2))})
    tree$edge.diam <- 2*sqrt(tree$edge.weight/pi/edge.len)
    message("Edge diameters were computed based on edge length and weight.\n")
  }
  
  if ("edge.type" %in% attributes(tree)$names){
    etype <- tree$edge.type
  } else {
    etype <- rep("d",nrow(tree$edge))
  }
  if ("edge.color" %in% attributes(tree)$names){
    ecol <- tree$edge.color
  } else {
    ecol <- rep("gray30",nrow(tree$edge))
  }
  edges_to_draw <- which(tree$edge.diam>0)
  for(i in edges_to_draw){
    if(etype[i] == "d"){ # default edge
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center = tree$node.coord[tree$edge[i,1:2],], 
                        radius = tree$edge.diam[i]/2)), 
        col = ecol[i], alpha=0.8)
    } else if (etype[i] == "b") { # bud approximated with 4 cylinders
      startnode <- tree$node.coord[tree$edge[i,1],1:3]
      endnode <- tree$node.coord[tree$edge[i,2],1:3]
      bud_centers <- rbind(startnode,
                           startnode+0.25*(endnode-startnode),
                           startnode+0.6*(endnode-startnode),
                           startnode+0.85*(endnode-startnode),
                           endnode)
      bud_diams <- c(tree$edge.diam[i]*0.7,tree$edge.diam[i],
                     tree$edge.diam[i]*0.7,tree$edge.diam[i]*0.4,
                     tree$edge.diam[i]*0.1)/2
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center=bud_centers, 
                        radius = bud_diams)),
        col = ecol[i],alpha=0.8)
    } else if (etype[i] == "l") { # leaf depiction as hexagon
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center=tree$node.coord[tree$edge[i,1:2],], 
                        radius = c(tree$edge.diam[i],tree$edge.diam[i]*0.4)/2)), 
        col = ecol[i],alpha=0.8)
      startnode <- tree$node.coord[tree$edge[i,1],1:3]
      endnode <- tree$node.coord[tree$edge[i,2],1:3]
      # The spanning vector is orthogonal to the z-axis and the leaf edge.
      span_vec <- NULL
      # Check if leaf edge is parallel to the z-axis.
      if(startnode[1]!=endnode[1] && startnode[2]!=endnode[2]){
        span_vec <- c((startnode[2]-endnode[2])/(endnode[1]-startnode[1]),1,0)
      } else if(startnode[1]!=endnode[1]){ # then s[2]=e[2]
        span_vec <- c(0,1,0)
      } else if(startnode[2]!=endnode[2]){ # then s[1]=e[1]
        span_vec <- c(1,0,0)
      } else { # here s[1]=e[1] and s[2]=e[2]
        span_vec <- c(1,1,0)
      }
      span_vec <- span_vec/stats::dist(rbind(span_vec, rep(0,3)), #normalizing
                                       method ="euclidean")
      # Make leaf width proportional to leaf length.
      leaf_length <- sqrt(sum((endnode-startnode)^2))
      leaf_centers <- rbind(startnode,
                            startnode+0.25*(endnode-startnode)+ 
                              0.2*span_vec*leaf_length,
                            startnode+0.6*(endnode-startnode)+ 
                              0.17*span_vec*leaf_length,
                            startnode+0.85*(endnode-startnode)+ 
                              0.07*span_vec*leaf_length,
                            endnode)
      leaf_diams <- c(tree$edge.diam[i]*0.9,tree$edge.diam[i]*0.8,
                      tree$edge.diam[i]*0.7,tree$edge.diam[i]*0.6,
                      tree$edge.diam[i]*0.3)/2
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center = leaf_centers, 
                        radius = leaf_diams)), 
        col = ecol[i],alpha=0.8)
      leaf_centers <- rbind(startnode,
                            startnode+0.25*(endnode-startnode)- 
                              0.2*span_vec*leaf_length,
                            startnode+0.6*(endnode-startnode)- 
                              0.17*span_vec*leaf_length,
                            startnode+0.85*(endnode-startnode)- 
                              0.07*span_vec*leaf_length,
                            endnode)
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center = leaf_centers, 
                        radius = leaf_diams)), 
        col = ecol[i],alpha=0.8)
    } else {
      warning("Unknown edge type for edge ", i,".\n")
    }
  }
  if(show_node_enum){
    for(v in 1:nrow(tree$node.coord)){
      rgl::text3d(x = tree$node.coord[v,1], 
             y = tree$node.coord[v,2], 
             z = tree$node.coord[v,3], texts=paste(v), pos = 5, offset = 1)
    }
  }
}
#' Plot a phylo3D object
#'
#' \code{plotImbalPhylo3D} - Plots a phylo3D object using functions of the 
#' package 'rgl'. Moreover, it uses either brightness or a color scale to
#' indicate the imbalance.\cr
#' Edge sections are shown darker or red with higher degree of imbalance and
#' brighter or cyan if they are balanced. 
#' This function does not use the parameter \code{edge.color}\cr
#' Attention: Edges of type 'bud' or 'leaf' will always be depicted as 
#' balanced, because they should represent leaf edges that are by definition
#' always balanced.
#' 
#' @param imbal_type Specifies which node imbalance measurement should be used.
#' Available are:\cr
#' "A"     - centroid angle\cr
#' "alpha" - minimal centroid angle\cr
#' "M"     - expanded relative centroid distance\cr
#' "mu"    - relative centroid distance
#' @param max.seclen Numeric value >0 that specifies the maximal section length.
#' @param color.imbal Boolean value (default TRUE). If true, colors are 
#' used to depict the imbalance. Otherwise, a grayscale image is produced. 
#' @param show.gradient Boolean value (default FALSE). If true 
#' the color or grayscale gradient is depicted. 
#' 
#' @return \code{plotImbalPhylo3D} No return value, called for side effects 
#' (plotting).
#' 
#' @export
#' @rdname plotPhylo3D
#' @examples
#' plotImbalPhylo3D(tree, imbal_type="mu", max.seclen=0.5, color.imbal=TRUE,
#'     show.gradient=FALSE)
plotImbalPhylo3D <- function(tree, imbal_type, max.seclen,
                             color.imbal=TRUE, show.gradient=FALSE){
  # Measure the tree to create the box for the plot.
  square_width <- max(max(tree$node.coord[,1]) - min(tree$node.coord[,1]),
                      max(tree$node.coord[,2]) - min(tree$node.coord[,2]),
                      max(tree$node.coord[,3]) - min(tree$node.coord[,3]))
  my_epsilon <- square_width/20
  square_mids <- c((max(tree$node.coord[,1]) + min(tree$node.coord[,1]))/2,
                   (max(tree$node.coord[,2]) + min(tree$node.coord[,2]))/2,
                   (max(tree$node.coord[,3]) + min(tree$node.coord[,3]))/2)
  # Construct basic plot
  rgl::plot3d(1,1,1,type = "n", 
              xlab = "x", ylab = "y", zlab = "z",
              xlim = c(square_mids[1]-square_width/2-my_epsilon,
                       square_mids[1]+square_width/2+my_epsilon), 
              ylim = c(square_mids[2]-square_width/2-my_epsilon,
                       square_mids[2]+square_width/2+my_epsilon), 
              zlim = c(square_mids[3]-square_width/2-my_epsilon,
                       square_mids[3]+square_width/2+my_epsilon),
              surface=FALSE, axis.scales = FALSE, aspect = TRUE)
  addImbalPhylo3D(tree, offset = c(0,0,0), imbal_type, max.seclen,
                  color.imbal, show.gradient)
}
#' Plot a phylo3D object
#'
#' \code{addImbalPhylo3D} - This function plots a phylo3D object without any 
#' coordinate axis or adds the tree to an existing plot (e.g., for 
#' \code{plotImbalPhylo3D}). Moreover, it uses either brightness or a color 
#' scale to indicate the imbalance.\cr
#' Edge sections are shown darker or red with higher degree of imbalance and
#' brighter or cyan if they are balanced. 
#' This function does not use the parameter \code{edge.color}.\cr
#' Attention: Edges of type 'bud' or 'leaf' will always be depicted as 
#' balanced because they should represent leaf edges that are always 
#' balanced.
#'
#' @return \code{addImbalPhylo3D} No return value, called for side effects 
#' (plotting).
#' 
#' @export
#' @rdname plotPhylo3D
#' @examples
#' addImbalPhylo3D(tree, imbal_type="mu", offset = c(1,0,0), max.seclen=0.5, 
#'                 color.imbal=FALSE, show.gradient = FALSE) 
addImbalPhylo3D <- function(tree, offset = c(0,0,0), imbal_type, max.seclen,
                            color.imbal=TRUE, show.gradient=FALSE){
  # Ensure not to change the user's options.
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  # Check input.
  if(sum(c("node.coord") %in% attributes(tree)$names)!=1){
    stop(paste("Cannot calculate subtree centroids. The tree is missing",
               "this attribute:'node.coord'."))
  }
  if(!"node.descs" %in% attributes(tree)$names){
    tree$node.descs <- getDescs(tree)
  }
  # First, gather the imbalance data and create a coloring scheme.
  e_profile <- imbalProfile(tree, imbal_type=imbal_type, 
                            max.seclen=max.seclen)[,c(4,5)]
  e_brightness <- NULL
  if(imbal_type=="mu"){
    e_brightness <- cbind(round(e_profile[,1]*80)+1,
                          e_profile[,2])
  }else if(imbal_type=="M"){
    e_brightness <- cbind(round(e_profile[,1]/2*80)+1,
                          e_profile[,2])
  }else if(imbal_type=="A"){
    e_brightness <- cbind(round(e_profile[,1]/pi*80)+1,
                          e_profile[,2])
  }else if(imbal_type=="alpha"){
    e_brightness <- cbind(round(e_profile[,1]/(pi/2)*80)+1,
                          e_profile[,2])
  }else{
    stop(paste("Unknown imbalance type: ",imbal_type))
  }
  imbal_palette <- NULL
  if(color.imbal){
    colfunc <- grDevices::colorRampPalette(c("cyan","yellow","red",
                                             "deeppink4","black"))
    imbal_palette <- colfunc(81)
  }else{# gray: BRIGHTEST <- "gray81"; DARKEST <- "gray1";
    imbal_palette <- paste("gray",81:1,sep = "")
  }
  # Then proceed with the other preparations for plotting the tree
  tree$node.coord <- tree$node.coord +  # shift all coordinates by offset
    outer(rep.int(1L, nrow(tree$node.coord)), offset)
  
  if (!"edge.diam" %in% attributes(tree)$names){
    edge.len <- sapply(1:nrow(tree$edge), function(x) {
      sqrt(sum((tree$node.coord[tree$edge[x,1],]-
                  tree$node.coord[tree$edge[x,2],])^2))})
    tree$edge.diam <- 2*sqrt(tree$edge.weight/pi/edge.len)
    message("Edge diameters were computed based on edge length and weight.\n")
  }
  
  if ("edge.type" %in% attributes(tree)$names){
    etype <- tree$edge.type
  } else {
    etype <- rep("d",nrow(tree$edge))
  }
  
  # Now, start adding the edges to the plot.
  edges_to_draw <- which(tree$edge.diam>0)
  for(i in edges_to_draw){ 
    if(etype[i] == "d"){ # default edge
      e_bright_rows <- which(e_brightness[,2]==i)
      subdiv_numb <- length(e_bright_rows)
      xs_subdiv <- c(seq(from=0,by=1/subdiv_numb,
                           length.out=subdiv_numb),1)
      p <- tree$node.coord[tree$edge[i,1],1:3]
      v <- tree$node.coord[tree$edge[i,2],1:3]
      for(section in 1:subdiv_numb){ 
        rgl::shade3d(rgl::addNormals(
          rgl::cylinder3d(center = rbind(v+xs_subdiv[section]*(p-v),
                                         v+xs_subdiv[section+1]*(p-v)), 
                          radius = tree$edge.diam[i]/2)), 
          col = imbal_palette[e_brightness[e_bright_rows[section],1]])
      }
    } else if (etype[i] == "b") { # bud approximated with 4 cylinders
      if(!is.null(getChildren(tree,tree$edge[i,2])[1])){
        warning(paste("Edge",i,"has type 'bud' but is not a leaf edge."))
      }
      startnode <- tree$node.coord[tree$edge[i,1],1:3]
      endnode <- tree$node.coord[tree$edge[i,2],1:3]
      bud_centers <- rbind(startnode,
                           startnode+0.25*(endnode-startnode),
                           startnode+0.6*(endnode-startnode),
                           startnode+0.85*(endnode-startnode),
                           endnode)
      bud_diams <- c(tree$edge.diam[i]*0.7,tree$edge.diam[i],
                     tree$edge.diam[i]*0.7,tree$edge.diam[i]*0.4,
                     tree$edge.diam[i]*0.1)/2
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center=bud_centers, 
                        radius = bud_diams)),
        col = imbal_palette[1])
    } else if (etype[i] == "l") { # leaf depiction as hexagon
      if(!is.null(getChildren(tree,tree$edge[i,2])[1])){
        warning(paste("Edge",i,"has type 'leaf' but is not a leaf edge."))
      }
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center=tree$node.coord[tree$edge[i,1:2],], 
                        radius = c(tree$edge.diam[i],tree$edge.diam[i]*0.4)/2)), 
        col = imbal_palette[1],alpha=0.8)
      startnode <- tree$node.coord[tree$edge[i,1],1:3]
      endnode <- tree$node.coord[tree$edge[i,2],1:3]
      # The spanning vector is orthogonal to the z-axis and the leaf edge.
      span_vec <- NULL
      # Check if leaf edge is parallel to the z-axis.
      if(startnode[1]!=endnode[1] && startnode[2]!=endnode[2]){
        span_vec <- c((startnode[2]-endnode[2])/(endnode[1]-startnode[1]),1,0)
      } else if(startnode[1]!=endnode[1]){ # then s[2]=e[2]
        span_vec <- c(0,1,0)
      } else if(startnode[2]!=endnode[2]){ # then s[1]=e[1]
        span_vec <- c(1,0,0)
      } else { # here s[1]=e[1] and s[2]=e[2]
        span_vec <- c(1,1,0)
      }
      span_vec <- span_vec/stats::dist(rbind(span_vec, rep(0,3)), #normalizing
                                       method ="euclidean")
      # Make leaf width proportional to leaf length.
      leaf_length <- sqrt(sum((endnode-startnode)^2))
      leaf_centers <- rbind(startnode,
                            startnode+0.25*(endnode-startnode)+ 
                              0.2*span_vec*leaf_length,
                            startnode+0.6*(endnode-startnode)+ 
                              0.17*span_vec*leaf_length,
                            startnode+0.85*(endnode-startnode)+ 
                              0.07*span_vec*leaf_length,
                            endnode)
      leaf_diams <- c(tree$edge.diam[i]*0.9,tree$edge.diam[i]*0.8,
                      tree$edge.diam[i]*0.7,tree$edge.diam[i]*0.6,
                      tree$edge.diam[i]*0.3)/2
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center = leaf_centers, 
                        radius = leaf_diams)), 
        col = imbal_palette[1])
      leaf_centers <- rbind(startnode,
                            startnode+0.25*(endnode-startnode)- 
                              0.2*span_vec*leaf_length,
                            startnode+0.6*(endnode-startnode)- 
                              0.17*span_vec*leaf_length,
                            startnode+0.85*(endnode-startnode)- 
                              0.07*span_vec*leaf_length,
                            endnode)
      rgl::shade3d(rgl::addNormals(
        rgl::cylinder3d(center = leaf_centers, 
                        radius = leaf_diams)), 
        col = imbal_palette[1])
    } else {
      warning("Unknown edge type for edge ", i,".\n")
    }
  }
  if(show.gradient){ 
    if(color.imbal){
      rgl::show2d({
        graphics::par(mar=c(0,0,0,0))
        colfunc <- grDevices::colorRampPalette(rev(c("cyan","yellow","red",
                                                     "deeppink4","black")))
        legend_image <- grDevices::as.raster(matrix(colfunc(20), ncol=1))
        plot(c(-5,5),c(0,1.2),type = 'n', axes = F,xlab = '', ylab = '')
        graphics::text(x=1.5, y = 1.1, family="serif", labels = expression(mu),cex=1.5)
        graphics::text(x=2.5, y = 1.094, family="serif", labels = "\U1D4DC", cex=2.2)
        graphics::text(x=3.5, y = 1.104, family="serif", labels = expression(alpha),
             cex=1.5)
        graphics::text(x=4.5, y = 1.094, family="serif", labels = "\U1D49C", cex=2.2)
        graphics::text(x=1.5, y = seq(0,1,l=5), family="serif", labels = seq(0,1,l=5))
        graphics::text(x=2.5, y = seq(0,1,l=5), family="serif", labels = seq(0,2,l=5))
        graphics::text(x=3.5, y = seq(0,1,l=5), family="serif", 
             labels = c(0,expression(paste(0.25,pi)), expression(paste(0.5,pi)),
                        expression(paste(0.75,pi)), expression(pi)))
        graphics::text(x=4.5, y = seq(0,1,l=5), family="serif", 
             labels = c(0,expression(paste(0.125,pi)), 
                        expression(paste(0.25,pi)),
                        expression(paste(0.375,pi)), expression(paste(0.5,pi))))
        graphics::rasterImage(legend_image, 0, 0, 1,1)}, face = "ys")
    }else{
      rgl::show2d({
        graphics::par(mar=c(0,0,0,0))
        colfunc <- grDevices::colorRampPalette(rev(c("gray81","gray1")))
        legend_image <- grDevices::as.raster(matrix(colfunc(20), ncol=1))
        plot(c(-5,5),c(0,1.2),type = 'n', axes = F,xlab = '', ylab = '')
        graphics::text(x=1.5, y = 1.1, family="serif", labels = expression(mu),cex=1.5)
        graphics::text(x=2.5, y = 1.094, family="serif", labels = "\U1D4DC", cex=2.2)
        graphics::text(x=3.5, y = 1.104, family="serif", labels = expression(alpha),cex=1.5)
        graphics::text(x=4.5, y = 1.094, family="serif", labels = "\U1D49C", cex=2.2)
        graphics::text(x=1.5, y = seq(0,1,l=5), family="serif", labels = seq(0,1,l=5))
        graphics::text(x=2.5, y = seq(0,1,l=5), family="serif", labels = seq(0,2,l=5))
        graphics::text(x=3.5, y = seq(0,1,l=5), family="serif", 
             labels = c(0,expression(paste(0.25,pi)), expression(paste(0.5,pi)),
                        expression(paste(0.75,pi)), expression(pi)))
        graphics::text(x=4.5, y = seq(0,1,l=5), family="serif", 
             labels = c(0,expression(paste(0.125,pi)), expression(paste(0.25,pi)),
                        expression(paste(0.375,pi)), expression(paste(0.5,pi))))
        graphics::rasterImage(legend_image, 0, 0, 1,1)}, face = "ys")
    }
  }
}
