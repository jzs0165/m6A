library(ggplot2)
theme_set(theme_classic())

# Plot
g <- ggplot(m6A, aes(Location))
g + geom_density(aes(fill=factor(Frequency)), alpha=0.8) + 
  labs(title="Density plot", 
       subtitle="m6A modifications on PAN RNA",
       caption="Source: ",
       x="Location",
       fill="Frequency")

