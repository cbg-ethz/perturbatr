#!/usr/bin/env Rscript

library(hexSticker)
library(ggplot2)
library(ggsci)
library(here)

data(iris)

outfile <- paste0(here::here(), "/inst/figure/sticker.png")

p <- ggplot(aes(x = Sepal.Length, y = -Sepal.Width, shape = Species, color = Species), data = iris) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE, size=.75) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#2D2D2D")) +
  ggsci::scale_color_tron()  +
  guides(shape=FALSE, color=FALSE)

sticker(p, package="", p_size=8, s_x=1, s_y=1,1, s_width=1.75, s_height=1.0,
        h_fill="#2D2D2D", h_color="#2D2D2D",
        filename=outfile,
        url="perturbR", u_color="white", u_size=4)
