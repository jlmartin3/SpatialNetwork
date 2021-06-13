
# code for running a cost based network analysis on GIS data
# requires points and lines, where points are origin/desinations
#    and lines represent the network (roads, paths, ect.)
# points are snaped to the nearest point on the network
# the shortest path between two given nodes is calcuated
#    and the distance in meters is returned


options(scipen=999) # take sci notation off

setwd("C:/Users/James/Desktop/huaq_network")
library(sp)
library(spatstat)
library(igraph)
library(shp2graph)
library(rgdal)
library(rgeos)
library(maptools)

# Object ID of origin and destination features 
origin_point      <- 76
destination_point <- 33

#_______________ Imports ______________________

# projection used for this example 
utm_17s <- '+proj=utm +zone=17 +south 
+datum=WGS84 +units=m +no_defs 
+ellps=WGS84 +towgs84=0,0,0 '

# import road network shapefile
roads.shp <- readOGR("small_shapes/roads.shp")
roads.shp <-spTransform(roads.shp, utm_17s)

# import buildings shapefile
builds.shp <- readOGR("small_shapes/buildings.shp")
builds.shp <-spTransform(builds.shp, utm_17s)

# take a look 
par(mar = c(0,0,0,0))
plot(builds.shp)
lines(roads.shp, col = "red", lwd = 3)


#__________ convert building shapes to points ________

# In this example we're using buildings as 
#   destinations, they are polygons but we need points
# We get the center points of each polygon

# get centroids and make a new spatial points df
builds.pnt <- SpatialPointsDataFrame(
               gCentroid(builds.shp, byid=TRUE), 
                builds.shp@data, match.ID=FALSE)
# extract coordinates as new df
builds.coord <- data.frame(builds.pnt@coords)
coordinates(builds.coord) <- c("x", "y")
proj4string(builds.coord) <- CRS(utm_17s)
# plot point numbers at polygon centroids 
with(builds.pnt, 
     text(coordinates(builds.pnt), 
          labels = as.character(1:length(builds.pnt))))
#points(builds.pnt)

#__________ make nodes at intersections _______________

# the first 'nodes' are made where lines intersect
# really what we're doing in this step is cutting 
# the roads shapefile for later

# convert to psp for processing
#   this also drops the original attribute data 
roads.psp <- as.psp.SpatialLinesDataFrame(roads.shp)

# cut lines where intersections occur so
#    later we can count these as nodes
roads.cut   <- selfcut.psp(roads.psp)

# convert back to sldf, add some dummy meta data so it works
roads.shp2 <- as.SpatialLines.psp(roads.cut)
roads.shp2 <- SpatialLinesDataFrame(roads.shp2,
                   data.frame(c(1:length(roads.shp2@lines))))
proj4string(roads.shp2) <- utm_17s


#_____ make nodes on network closest to points _______

# In this step we cut lines again. 
# This time calculate which line is closest to 
#   a given point that isnt on the network. 
# Then we put a point on that line where its closest 
#   to the point.
# Then we split the line where this point is.

# for a given point, which line is nearest?
line_nearest <- c()
for(i in 1:length(builds.coord)){ # for each point
  # get distances to all lines in shapefile
  line_dist <- c(gDistance(builds.coord[i,], 
                           roads.shp2, byid=TRUE, 
                           hausdorff=FALSE, densifyFrac = NULL))
  # get the number of the closest line
  line_nearest[i] <- which(line_dist == min(line_dist))[1]
  
}
 # number of lines that should be split
length(unique(line_nearest)) 

# how many times is should each line be split?/ 
# how many points are on each line(that has points to begin with)
splits_reps <- c()
for(i in 1:length(unique(line_nearest))){ # for each line with a point
  splits_reps[i] <- length(which( # how many points are on it
    line_nearest == unique(line_nearest)[i]))
}

# new list to store new lines
all_new_lines <- vector(mode = "list", length(unique(line_nearest)))
for(i in 1:length(unique(line_nearest))){ # for each line to be split
  # get the points associated with the line
  point.num <- which(line_nearest == unique(line_nearest)[i])
  # make a data frame to store coordinates and ids
  point.coords <- data.frame("x" = NA, "y" = NA, "id" = NA)
  # add the end points of the line at the begining of the df, their coords and ids
  point.coords[1,] <- c(coordinates(roads.shp2@lines[[unique(line_nearest)[i]]])[[1]][1,],
                        paste0("line",unique(line_nearest)[i],"start"))
  point.coords[2,] <- c(coordinates(roads.shp2@lines[[unique(line_nearest)[i]]])[[1]][2,],
                        paste0("line",unique(line_nearest)[i],"end"))
  # add the coords and ids of each point as a row after the line endpoints
  for(j in 1:length(point.num)){
    point.coords[j+2,] <- c(coordinates(builds.coord[point.num[j]]),
                            paste0("point", point.num[j]))
  }
  # make sure formating is correct to run geometry functions 
  point.coords[[1]] <- as.numeric(point.coords[[1]])
  point.coords[[2]] <- as.numeric(point.coords[[2]])
  point.coords2 <- point.coords
  coordinates(point.coords) <- c("x", "y")
  proj4string(point.coords) <- CRS(utm_17s)
  
  # we swap out the coordinates for each point, with the coordinates 
  #   of the point on the line closest to it
  #   This is where we "snap" points to the network
  for(j in 3:length(point.coords)){
    point.coords2[j,] <- c(nearestPointOnSegment(
      rbind(coordinates(point.coords[1,]),coordinates(point.coords[2,])),
      coordinates(point.coords[j,]))[1:2],  point.coords@data[j,])
  }
  # format again, this time to sort
  point.coords2[[1]] <- as.numeric(point.coords2[[1]])
  point.coords2[[2]] <- as.numeric(point.coords2[[2]])
  point.coords3 <-  point.coords2
  coordinates(point.coords2) <- c("x", "y")
  proj4string(point.coords2) <- CRS(utm_17s)
  
  # we're sorting the rows of the data frame based on how close each 
  #   point is to the line start point. This way we can make lines in 
  #   sequence
  dist4sort <- c()
  for(j in 1:length(point.coords2)){
    dist4sort[j] <- gDistance(point.coords2[1,], point.coords2[j,])
  }
  # reorder and reset nrow names
  point.coords3            <-  point.coords3[order(dist4sort),]
  row.names(point.coords3) <-  c(1:nrow(point.coords3))
  
  # list to store new lines, the result of all the spliting
  new_lines <- vector(mode = "list", length = (nrow(point.coords3)-1))
  
  # using the coordinates in the df, we make lines iteratively
  #   connecting all the 'nodes' on the old line
  for(j in 1:(nrow(point.coords3)-1)){
    
    new_line <-   rbind( point.coords3[j,   1:2] ,
                         point.coords3[j+1, 1:2]) 
    
    new_lines[[j]]     <- Lines(sp::Line(new_line), 
                      paste0(point.coords3[j, 3] ,"-", point.coords3[j+1, 3]) )
  }
  # add to the lines list of lists
  all_new_lines[[i]] <- new_lines
  
}

# turn a list of lists into just a list 
all_new_lines <- unlist(all_new_lines, recursive = FALSE)

# make a new spatial lines df containing all the lines
#   that were not split (the lines that aren't having nodes placed on them)
roads.shp2.1 <- roads.shp2[-unique(line_nearest),]
for(i in 1:length(roads.shp2.1)){ # assign IDs
  roads.shp2.1@lines[[i]]@ID <- paste0("line",  roads.shp2.1@lines[[i]]@ID,"start", "-",
                                   "line",  roads.shp2.1@lines[[i]]@ID,"end")
}

# combine the unsplit and split lines
roads.shp3 <- c(roads.shp2.1@lines, all_new_lines)
roads.shp3 <- SpatialLines(roads.shp3, proj4string = CRS(utm_17s))
# extract IDs and add them again when making spatial lines df
IDs <- c()
for(i in 1:length(roads.shp3)){
  IDs[i] <- roads.shp3@lines[[i]]@ID
}
roads.shp3 <- SpatialLinesDataFrame(roads.shp3,
                                    data.frame(IDs),
                                    match.ID = FALSE)
# when plotted, roads.shp3 should look exactly like roads.shp and roads.shp2
plot(roads.shp3)


#____________ converting to a network object ____________________

# for distance weighting, we calcuate the lengths in meters of 
#   each line segment in the sptail lines df
edge_lengths <- c()
for(i in 1:length(roads.shp3@lines)){
  edge_lengths[i] <- LinesLength(roads.shp3@lines[[i]])
}

# extract a formated nodes and edgelist to make the network
roads.nw <- readshpnw(roads.shp3,
                    ELComputed=TRUE, longlat=FALSE,
                    Detailed=FALSE, 
                    ea.prop= c(1))

roads_nodes <- roads.nw[[2]]
roads_edges <- roads.nw[[3]]

# create network using igraph, weight by lengths
roads.network <- nel2igraph(roads_nodes,
                            roads_edges,
                  weight = edge_lengths)
# remove self connections ("loops")
roads.network <-simplify(roads.network)

# plot
plot(roads.network, 
     vertex.size = 1, 
     # vertex.label = NA,  
     vertex.color = "red")

# check to see if everything is connected 
components(roads.network)


#_______ relate original point features to their nodes ___

# Through many conversions we can't easly tell which node
#    in the network represents a given feature
# In this part we connect the node IDs back to the point IDs

# Split the edge id string into 2 strings, 
#   one for each endpoint
og_edge_ids <- strsplit( as.character(roads.nw[[5]][,2]), "-")

# create a dataframe, 
# 'ID' column is the original ID of the feature
# 'node' is the network node that each feature is associated with 
node.ID <- data.frame(
  ID = c(do.call(rbind, og_edge_ids)[,1], 
         do.call(rbind, og_edge_ids)[,2]),
  node = c(roads_edges[,2], roads_edges[,3])
)

# create a list, that for each node,
#   contains the ID number of each point associated with it
node.relates <- vector(mode = "list", length(unique(node.ID$node)))
for(i in 1:length(unique(node.ID$node))){
  # retrieve all IDs associated with a node
  id.locs <- which(node.ID$node == i)
  ids     <- unique(as.character(node.ID$ID[id.locs]))
  ids     <-  ids[grepl("point", ids)] # drop lines
  if(length(ids) == 0){
    ids <- NA
  } # convert IDs from a character to a number 
  node.relates[[i]] <- as.numeric(gsub("point", "", ids ))
}

# get the node numbers for the specified origin and destination points
origin_node <- which(
  unlist(lapply(node.relates, function(x) any(x == origin_point))))

destination_node <- which(
       unlist(lapply(node.relates, function(x) any(x == destination_point))))


#_________________ results and plotting ______________________________

d_path <- distances(roads.network, v = V(roads.network)[origin_node],
          to = V(roads.network)[destination_node])[,1]

s_path <- shortest_paths(roads.network, from = V(roads.network)[origin_node],
                    to = V(roads.network)[destination_node])
E(roads.network)$color <- "grey"
E(roads.network, path=unlist(s_path[[1]]))$color <- "red"

par(mar = c(0,3,3,0))
plot(roads.network, 
     vertex.size = 1, 
     vertex.color = "red", 
     main = paste( round(d_path),
                  "meters"))
s_path[[1]]
d_path

















