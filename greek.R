
library(derivmkts)
library(ggplot2)
library(lemon)
s <- 36; r <- 0.05; v <- 0.4; d <- 0
k <- seq(.1, 60, by=.1)
period = 12
tt = (1:12)/period
y = list()

for(i in 1:period){
  y[[i]] = greeks(bscall(s, k, v, r, tt[i], d), long=TRUE)
}

colors = c('1' = 'lightcoral', '2' = 'indianred4', '3' = 'lightskyblue4',
           '4' = 'midnightblue', '5' = 'olivedrab', '6' = 'mediumaquamarine',
           '7' = 'orangered3', '8' = 'palegreen4', '9' = 'red',
           '10' = 'paleturquoise4', '11' = 'peru', '12' = 'violetred4')


pdf(file = "C:\\Users\\sikde\\OneDrive\\Studium\\M.Sc. Statistics\\5_semester\\greek.pdf",width = 6.5, height = 4)
greek_plot = ggplot()+
  geom_line(data = y[[1]], aes(x = k, y = value, color = '1')) + 
  geom_line(data = y[[2]], aes(x = k, y = value, color = '2')) + 
  geom_line(data = y[[3]], aes(x = k, y = value, color = '3')) + 
  geom_line(data = y[[4]], aes(x = k, y = value, color = '4')) + 
  geom_line(data = y[[5]], aes(x = k, y = value, color = '5')) + 
  geom_line(data = y[[6]], aes(x = k, y = value, color = '6')) + 
  geom_line(data = y[[7]], aes(x = k, y = value, color = '7')) + 
  geom_line(data = y[[8]], aes(x = k, y = value, color = '8')) + 
  geom_line(data = y[[9]], aes(x = k, y = value, color = '9')) + 
  geom_line(data = y[[10]], aes(x = k, y = value, color = '10')) + 
  geom_line(data = y[[11]], aes(x = k, y = value, color = '11')) + 
  geom_line(data = y[[12]], aes(x = k, y = value, color = '12')) + 
  labs(x = 'Strike price', color = 't (month)')+
  facet_wrap(~ greek, scales = 'free_y') + 
  scale_color_manual(values = colors) +
  theme(legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.position = 'bottom',
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.45, 'cm'))+ 
  guides(color=guide_legend(ncol=3)) + scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(panel.background = element_blank())
  
reposition_legend(greek_plot, 'center', panel = 'panel-3-3')
 dev.off()
