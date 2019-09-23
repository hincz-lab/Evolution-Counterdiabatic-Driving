library(ggraph)
library(igraph)
library(binaryLogic)
library(tidygraph)
library(tidyr)
library(dplyr)
library(graphlayouts)
library(stringr)
library(gganimate)

# Helper for converting node ids to labels
make_label <- function(x) {paste(as.character(as.binary(x-1, n=4)), collapse="")}
# Helper for finding y coordinate based on node id
get_y <- function(x){sum(as.vector(as.binary(x-1)))}

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

add_data_to_graph <- function(layout, df) {
  spread_pop_df <- df %>% spread(generation, pop)
  new_layout <- layout %>% full_join(spread_pop_df) %>% gather(-name, -circular, -label, -.ggraph.orig_index, -.ggraph.index, -x, -y, key="time", value="pop")    
  attr(new_layout, "graph") <- attr(layout, "graph")
  attr(new_layout, "circular") <- attr(layout, "circular")
  attr(new_layout, "class") <- attr(layout, "class")
  new_layout$time <- as.numeric(new_layout$time)
  return(new_layout)
}

make_graph <- function() {
  layout <- make_layout()  
  ggraph(layout) + geom_edge_link() + geom_node_label(aes(label=label))
}

# Edit this file name to change the input file
pop_sizes <- read_csv("../data/1/pop_sizes.csv")
df <- pop_sizes %>% gather(-generation, key="name", value="pop")

# Adjust the numbers in the calls to filter to adjust lower and upper bounds
add_data_to_graph(layout, df %>% filter(generation > 10000) %>% filter(generation < 11000))
df$name <- str_replace(df$name, "pop", "V")
df$name <- str_replace(df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})

layout <- make_layout()
data_layout <- add_data_to_graph(layout, df %>% filter(generation > 10000) %>% filter(generation < 11000) %>% filter(generation %% 10 == 0))
ggraph(new_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=sqrt(pop)/100), fill="white") + theme_graph(background = "white") + geom_node_text(aes(label=label)) + transition_time(time)


