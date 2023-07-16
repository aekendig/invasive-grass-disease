# theme
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_blank(),
        legend.margin = margin(-0.3, 0, -0.1, 0, unit = "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7,
                                    margin = margin(0, 0, 0.1, 0, unit = "cm")),
        strip.text.y = element_text(size = 7,
                                    margin = margin(0, 0.1, 0, 0, unit = "cm")),
        strip.placement = "outside",
        plot.title = element_text(size = 8, hjust = 0.5, face = "bold"))

# colors
col_pal = c("black", "#238A8D")

# dodge size
dodge_width <- 0.5

# text size
textSize = 2.5
