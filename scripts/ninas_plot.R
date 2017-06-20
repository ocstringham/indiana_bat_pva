# ninas plot

rm(list=ls())

library(ggplot2)
library(magrittr)

setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs")
df = read.csv("ninas_data.csv") %>% t() %>% as.data.frame(stringsAsFactors = F)
colnames(df)[1] = df$V1[1]
colnames(df)[2] = df$V2[1]
df = df[2:26 ,]
df$Year = 2011:2035
df$Linear = df$Linear %>% as.numeric()
df$Log = df$Log %>% as.numeric()

#plot
bat_plot = function(y_max)
{
p = ggplot() +
  
  geom_line(data = df, aes(x=Year, y=Linear), linetype=2) +
  
  geom_line(data = df, aes(x=Year, y=Log)) +
  
  
  scale_colour_manual(labels = c("Linear", "Log"),
                      values=c("black", "black"))  +
  
  theme(
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    panel.grid.major.y = element_line(colour="#d9d9d9", size=0.25),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(colour="#d9d9d9", size=0.25),
    legend.position = c(0.8,0.8)#,
    #axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  # scale_colour_manual("", 
  #                     breaks = c("Linear", "Log")) +
  scale_y_continuous("Population size",
                     expand = c(0, 0), 
                     limits = c(0,y_max), 
                     breaks = seq(0, y_max, 100))

X11()
print(p)
return(p)
}


p = bat_plot(600)


setwd("C:/Users/oliver/Google Drive/PhD/Research/Indiana bat/figs/")
tempname = "ninas_plot.pdf"
ggsave(tempname, plot = p, device = NULL, path = NULL,
       scale = 1, width = NA, height = NA,
       units = c("cm"), dpi = 300)
tempname = "ninas_plot.png"
ggsave(tempname, plot = p, device = NULL, path = NULL,
       scale = 1, width = NA, height = NA,
       units = c("cm"), dpi = 300)

