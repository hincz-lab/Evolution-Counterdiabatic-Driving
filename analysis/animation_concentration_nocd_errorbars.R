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

source("draw_tesseract.R")

concentrations <- read_csv("/home/emily/repos/evolution_speed/concentrations.csv")

cd_data <- import_data_set("/home/emily/Downloads/final_data/accelerated_equilib/nocd*/pop_sizes.csv")

theoretical <- prepare_file("/home/emily/Downloads/xeqi.csv")

layout <- make_layout()
other_data_layout <- add_data_to_graph(layout, theoretical %>% filter(generation %% 100 == 0))




grouped_df <- cd_data %>% gather(-generation, key="name", value="popsize") %>% group_by(generation, name)
median_df <- grouped_df %>% summarise(pop=median(popsize))
innerring_df <- grouped_df %>% summarise(pop=quantile(popsize)[2])
outerring_df <- grouped_df %>% summarise(pop=quantile(popsize)[4])
innerbar_df <- grouped_df %>% summarise(pop=quantile(popsize)[1])
outerbar_df <- grouped_df %>% summarise(pop=quantile(popsize)[5])

median_df$name <- str_replace(median_df$name, "pop", "V")
median_df$name <- str_replace(median_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
innerring_df$name <- str_replace(innerring_df$name, "pop", "V")
innerring_df$name <- str_replace(innerring_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
outerring_df$name <- str_replace(outerring_df$name, "pop", "V")
outerring_df$name <- str_replace(outerring_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
innerbar_df$name <- str_replace(innerbar_df$name, "pop", "V")
innerbar_df$name <- str_replace(innerbar_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})
outerbar_df$name <- str_replace(outerbar_df$name, "pop", "V")
outerbar_df$name <- str_replace(outerbar_df$name, "[:digit:]+", function(x){return(as.numeric(x)+1)})

median_layout <- add_data_to_graph(layout, median_df %>% filter(generation %% 100 == 0))
innerring_layout <- add_data_to_graph(layout, innerring_df %>% filter(generation %% 100 == 0))
outerring_layout <- add_data_to_graph(layout, outerring_df %>% filter(generation %% 100 == 0))
innerbar_layout <- add_data_to_graph(layout, innerbar_df%>% filter(generation %% 100 == 0))
outerbar_layout <- add_data_to_graph(layout, outerbar_df %>% filter(generation %% 100 == 0))
median_layout$min <- innerbar_layout$pop
median_layout$max <- outerbar_layout$pop
median_layout$quart2 <- innerring_layout$pop
median_layout$quart4 <- outerring_layout$pop

#a_conc <- ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75, fill="Equilibrium"), alpha=.5, linetype="blank") + geom_node_circle(data= other_data_layout, aes(r=(log10(pop+1))^2/75, fill="Observed"), color="red",alpha=.5, linetype="blank") + theme_graph(background = "white") + geom_node_text(aes(label=label)) + geom_text(aes(x=1.5, y=0, label=paste0("Time: ", time)), size = 6) + scale_fill_manual("", values=c("red", "blue")) + geom_text(data=concentrations, aes(x=6.2, y=.1, label=paste0("Concentration: \n", concentration)), size=6) + transition_time(time)
#a_conc <- ggraph(data_layout) + geom_edge_link(start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75, fill="Equilibrium"), alpha=.5, linetype="blank") + geom_node_circle(data= other_data_layout, aes(r=(log10(pop+1))^2/75, fill="Observed"), color="red",alpha=.5, linetype="blank") + theme_graph(background = "white") + geom_node_text(aes(label=label)) + scale_fill_manual("", values=c("red", "blue")) + transition_time(time) +labs(title = 'Time: {frame_time}', subtitle='Drug Concentration: {concentrations[frame_time+1,]["concentration"]}') + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#a_conc <- ggraph(data_layout) + geom_edge_link(start_cap = circle(15, 'mm'), end_cap = circle(15, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75, fill="Equilibrium"), alpha=.5, linetype="blank") + geom_node_circle(data= other_data_layout, aes(r=(log10(pop+1))^2/75, fill="Observed"), color="red",alpha=.5, linetype="blank") + theme_graph(background = "white") + geom_node_text(aes(label=label), size=12) + scale_fill_manual("", values=c("red", "blue")) + labs(title = 'Time: {as.integer(frame_time)}', subtitle='Concentration: {concentrations[frame_time+1,]["concentration"]}') + theme(plot.title = element_text(hjust = 0.5, size=30), plot.subtitle = element_text(hjust = 0.5, size = 30), legend.text = element_text(hjust = 0.5, size = 30)) +transition_time(time)
a_conc <- ggraph(median_layout) + geom_edge_link(start_cap = circle(15, 'mm'), end_cap = circle(15, 'mm')) + geom_node_circle(aes(r=(log10(pop+1))^2/75, fill="Observed"), alpha=.5, linetype="blank") + geom_errorbarh(aes(y=y, xmin=x-(log10(min+1))^2/75, xmax=x-(log10(max+1))^2/75), color="darkgreen", height=.1) + geom_node_circle(data= other_data_layout, aes(r=(log10(pop+1))^2/75, fill="Equilibrium"), color="red",alpha=.5, linetype="blank") + theme_graph(background = "white") + geom_node_text(aes(label=label), size=12) + scale_fill_manual("", values=c("red", "blue")) + labs(title = 'Time: {as.integer(frame_time)}', subtitle='Concentration: {concentrations[frame_time+1,]["concentration"]}') + theme(plot.title = element_text(hjust = 0.5, size=30), plot.subtitle = element_text(hjust = 0.5, size = 30), legend.text = element_text(hjust = 0.5, size = 30)) +transition_time(time)

anim_conc <- animate(a_conc, nframes = 1000, fps = 100, width = 1000, height = 1000)


anim_save("conc_anim_nocd_errorbars.gif", anim_conc)