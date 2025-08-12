library(plotly)

library(readr)
esbl <- read_csv("esbl.csv")
View(esbl)

install.packages("networkD3")
library(networkD3) 
library(viridis) # 
data <- esbl
nodes <- unique(c(data$Type, data$Source, data$Element.symbol))
nodes <- data.frame(name = nodes)


links <- rbind(
  data.frame(source = match(data$Type, nodes$name) - 1, target = match(data$Source, nodes$name) - 1, value = 1),
  data.frame(source = match(data$Source, nodes$name) - 1, target = match(data$Element.symbol, nodes$name) - 1, value = 1)
)

# Create Sankey plot without legend and with color mapping
sankey <- sankeyNetwork(
  Links = links, 
  Nodes = nodes, 
  Source = "source", 
  Target = "target", 
  Value = "value", 
  NodeID = "name", 
  units = "Count", 
  fontSize = 16, 
  nodeWidth = 8, 
  fontFamily = 'Times New Roman' 
 # colourScale = JS(paste0("d3.scaleOrdinal().domain(", jsonlite::toJSON(names(spectrum_colors)), ").range(", jsonlite::toJSON(spectrum_colors), ")")),
  #showLegend = FALSE # Removes the legend
)

# Plot the Sankey diagram

sankey
