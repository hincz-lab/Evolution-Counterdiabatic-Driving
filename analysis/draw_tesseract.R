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
filename1 <- "../pop_sizes.csv"
filename2 <- "../pop_sizes.csv"

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
  ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75), fill="lightblue", color="lightblue") + theme_graph(background = "white") + geom_node_text(aes(label=label)) + transition_time(time)
}

# Make a figure of a single snapshot in time
make_snapshot <- function(filename, time_point) {
  df <- prepare_file(filename)
  layout <- make_layout()
  data_layout <- add_data_to_graph(layout, df %>% filter(generation == time_point))
  ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75), fill="lightblue", color="lightblue") + theme_graph(background = "white") + geom_node_text(aes(label=label))  
}
  

make_circle_overlay <- function(filename1, filename2, time_point) {
  layout <- make_layout()
  df1 <- prepare_file(filename1)
  df2 <- prepare_file(filename2)
  data_layout <- add_data_to_graph(layout, df1 %>% filter(generation == time_point))
  other_data_layout <- add_data_to_graph(layout, df2 %>% filter(generation == time_point))
  ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75), fill="blue", color="blue", alpha=.5, linetype="blank") + geom_node_circle(data= other_data_layout, aes(r=(log10(pop+1))^2/75), fill="red", color="red",alpha=.5, linetype="blank") + theme_graph(background = "white") + geom_node_text(aes(label=label))
}

draw_box_and_whiskers_animation <- function(filenameglob, lower_time_bound, upper_time_bound, frame_freq = 10) {
  fileNames <- Sys.glob(filenameglob)
  df <- data.frame()
  for (file in fileNames) {
    df <- rbind(df, read_csv(file))
  }
  
  median_df <- df %>% gather(-generation, key="name", value="popsize") %>% group_by(generation, name) %>% summarise(pop=median(popsize))
  innerring_df <- df %>% gather(-generation, key="name", value="popsize") %>% group_by(generation, name) %>% summarise(pop=quantile(popsize)[2])
  outerring_df <- df %>% gather(-generation, key="name", value="popsize") %>% group_by(generation, name) %>% summarise(pop=quantile(popsize)[4])
  
  
  median_df$name <- str_replace(median_df$name, "pop", "V")
  median_df$name <- str_replace(median_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
  
  innerring_df$name <- str_replace(innerring_df$name, "pop", "V")
  innerring_df$name <- str_replace(innerring_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
  
  outerring_df$name <- str_replace(outerring_df$name, "pop", "V")
  outerring_df$name <- str_replace(outerring_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
  
  median_layout <- add_data_to_graph(layout, median_df %>% filter(generation > lower_time_bound) %>% filter(generation < upper_time_bound) %>% filter(generation  %% frame_freq == 0))
  innerring_layout <- add_data_to_graph(layout, innerring_df%>% filter(generation > lower_time_bound) %>% filter(generation < upper_time_bound) %>% filter(generation  %% frame_freq == 0))
  outerring_layout <- add_data_to_graph(layout, outerring_df%>% filter(generation > lower_time_bound) %>% filter(generation < upper_time_bound) %>% filter(generation  %% frame_freq == 0))
  
  ggraph(median_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(data=outerring_layout, aes(r=(log10(pop+1))^2/75), fill="white", color="blue") + geom_node_circle(data= median_layout, aes(r=(log10(pop+1))^2/75),color="white")+ geom_node_circle(data= innerring_layout, aes(r=(log10(pop+1))^2/75), color="red") + theme_graph(background = "white") + geom_node_text(aes(label=label)) + transition_time(time)
  }

# Uncomment this line to make a temporal animation
# make_temporal_animation(filename1, lower_time_bound, upper_time_bound)

# This line makes a snapshot
#make_snapshot(filename1, time_point)

#draw_box_and_whiskers_animation("../data/90*/pop_sizes.csv", lower_time_bound, upper_time_bound)

make_circle_overlay(filename1, filename2, time_point)
