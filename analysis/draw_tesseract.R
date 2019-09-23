library(ggraph)
library(igraph)
library(binaryLogic)
library(tidygraph)
library(tidyr)
library(dplyr)
library(graphlayouts)
library(stringr)
library(gganimate)
library(readr)

# Adjust this to choose file
filename <- "pop_sizes.csv"

# Adjust this to choose timepoint for snapshot
time_point <- 10000

# Adjust these to determine time range for animation
lower_time_bound <- 10000
upper_time_bound <- 20000

# Helper for converting node ids to labels
make_label <- function(x) {paste(as.character(as.binary(x-1, n=4)), collapse="")}
# Helper for finding y coordinate based on node id
get_y <- function(x){sum(as.vector(as.binary(x-1)))}

# Deal with ugly ggraph interal stuff
make_layout <- function() {
  adjmat <- read.csv("adj_mat.csv", header = FALSE)
  adjmat <- as.matrix(adjmat)
  landscape <- graph_from_adjacency_matrix(adjmat)
  layout <- create_layout(landscape, "stress")
  
  layout$y <- sapply(seq(16), get_y)
  layout$x <- c(4,6,5,7,3,6,3,6,2,2,5,5,1,3,2,4)
  layout$label <- sapply(layout$.ggraph.index, make_label) 
  return(layout)
}

# Sneakily add supplemental data to ggraph object
add_data_to_graph <- function(layout, df) {
  spread_pop_df <- df %>% spread(generation, pop)
  new_layout <- layout %>% full_join(spread_pop_df) %>% gather(-name, -circular, -label, -.ggraph.orig_index, -.ggraph.index, -x, -y, key="time", value="pop")    
  attr(new_layout, "graph") <- attr(layout, "graph")
  attr(new_layout, "circular") <- attr(layout, "circular")
  attr(new_layout, "class") <- attr(layout, "class")
  new_layout$time <- as.numeric(new_layout$time)
  return(new_layout)
}

# Read in file output by C++ code and convert it to correct format
prepare_file <- function(filename){
  pop_sizes <- read_csv(filename)
  df <- pop_sizes %>% gather(-generation, key="name", value="pop")
  df$name <- str_replace(df$name, "pop", "V")
  df$name <- str_replace(df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
  return(df)
}

# Make an animation over time
make_temporal_animation <- function(filename, lower_time_bound, upper_time_bound, frame_freq = 10) {
  df <- prepare_file(filename)
  layout <- make_layout()
  data_layout <- add_data_to_graph(layout, df %>% filter(generation > lower_time_bound) %>% filter(generation < upper_time_bound) %>% filter(generation %% frame_freq == 0))
  ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=sqrt(pop)/100), fill="white") + theme_graph(background = "white") + geom_node_text(aes(label=label)) + transition_time(time)
}

# Make a figure of a single snapshot in time
make_snapshot <- function(filename, time_point) {
  df <- prepare_file(filename)
  layout <- make_layout()
  data_layout <- add_data_to_graph(layout, df %>% filter(generation == time_point))
  ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=sqrt(pop)/100), fill="white") + theme_graph(background = "white") + geom_node_text(aes(label=label))  
}

# Uncomment this line to make a temporal animation
# make_temporal_animation(filename, lower_time_bound, upper_time_bound)

# This line makes a snapshot
make_snapshot(filename, time_point)
