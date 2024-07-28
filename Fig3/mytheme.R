library(ggplot2)

theme_blue <- theme(
  plot.title = element_text(
    size = 13,
    face = "bold",
    color = "darkred",
    hjust = 0,
    lineheight = 1.2
  ), 
  plot.subtitle = element_text(
    size = 13,
    face = "bold",
    color = "grey30",
    lineheight = 1.2,
    hjust = 0
  ), 
  panel.background = element_rect(fill = "#F6F4D2"), 
  panel.grid.major = element_line(
    colour = "gray80",
    size = 1,
    linetype = "dashed"
  ),
  panel.grid.minor = element_blank(), 
  axis.line.x = element_blank(), 
  axis.line.y = element_blank(),
  axis.title.x = element_text(
    vjust = 1,
    face = "bold",
    size = 14,
    color = "darkred"
  ), 
  axis.title.y = element_text(
    size = 14,
    face = "bold",
    color = "darkred"
  ), 
  axis.text.x = element_text(
    size = 12, colour = "black"
  ), 
  legend.title = element_text(size = 12, colour = "black"), 
  axis.text.y = element_text(size = 12, colour = "black"), 
  legend.text = element_text(size = 12, colour = "black"),
  panel.border = element_rect(color = "#80B199", fill = NA, size = 2), 
  legend.key = element_blank()
)
