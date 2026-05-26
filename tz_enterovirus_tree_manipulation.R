# all dependency installs now done using r2u on ubuntu

#remotes::install_github("ms609/TreeDist")

library(ggplot2)
library(ape)
library(dplyr)
library(phangorn)
library(treeio)
library(tidytree)
library(ggtree)
library(rio)
library(curl)
library(remotes)
library(TreeDist)
library(lubridate)
library(stringr)

correct_date <- function(date_string, month="06", day="15") {
  if (str_detect(date_string, "^\\d+$")[1]) {
    date_string <- paste(date_string, month, sep="-")
  }
  if (str_detect(date_string, regex("^\\d+-\\d+$"))[1]) {
    date_string <- paste(date_string, day, sep="-")
  }
  date_string
}

# load full genome and VP1 region trees

# from TZ Red Eyes 5 history
vp1_tree_filename <- 'Galaxy39-[IQ-TREE on data 30 and data 32_ Tree labeled with dates].nex'
vp1_tree = read.nexus(paste("data", vp1_tree_filename, sep="/"))
metadata = import("data/all_microreact_metadata.csv") %>% 
  mutate(Collection_Date = ymd(mapply(correct_date, Collection_Date)))
vp1_tree_annotated = full_join(vp1_tree, metadata, join_by(label == Accession))

# vp1_tree_annotated$edge.length <- vp1_tree_annotated$edge.length * 3270


# plot VP1 tree with clades of interest highlighted
tree_plot <- ggtree(vp1_tree_annotated, mrsd="2024-02-01") + 
  theme_tree2() +
  geom_tippoint(aes(color=Country)) +
  geom_highlight(node=382, fill="blue") + # highlight and label clade A
  geom_cladelab(382, "A", offset=-35, offset.text=-18) + 
  # geom_strip("MF419263.1", "MK989722.1", offset=-35, label="A1") +
  geom_highlight(node=323, fill="pink") + # highlight and label clade B
  geom_cladelab(323, "B", offset=-34.65, offset.text=-18, extend=-2) +
  geom_treescale(x=2010, y=20, width=10, offset=5) + # add scale bar
  vexpand(0.01, direction=-1) # give extra space at bottom
tree_plot

full_genome_tree_filename <- 'Galaxy100-[IQ-TREE on data 8 and data 93_ Tree labeled with dates].nex'
full_genome_tree = read.nexus(paste("data", full_genome_tree_filename, sep="/"))
full_genome_tree_annotated <- full_join(full_genome_tree, metadata, join_by(label == Accession))

# compute subsect of VP1 tree that matches full genome tree

full_tree <- ggtree(vp1_tree_annotated, mrsd="2024-02-01") + 
  theme_tree2() +
  geom_tippoint(aes(color=Country), size=2) +
  theme_tree2()
full_genome_tips <- full_genome_tree$tip.label

not_matching_tip_labels <- Filter(function(x) !(x %in% full_genome_tips), vp1_tree$tip.label)

# compute tanglegram of VP1 subset and full genome tree

vp1_tree_matching_full_tree <- drop.tip(vp1_tree, not_matching_tip_labels)

full_genome_ggtree = ggtree(midpoint(full_genome_tree))
vp1_ggtree = ggtree(midpoint(vp1_tree_matching_full_tree)) 

d1 <- full_genome_ggtree$data
d2 <- vp1_ggtree$data

d1$tree <- 'full_genome_ggtree'
d2$tree <- 'vp1_ggtree'

d2$x <- max(d2$x) - d2$x + max(d1$x) + max(d1$x)*0.3
pp <- full_genome_ggtree + geom_tree(data=d2)
pp

dd <- bind_rows(d1, d2) %>% 
  filter(isTip == TRUE)
dd1 <- as.data.frame(dd) 

pp + geom_line(aes(x, y, group=label), data=dd1, color='#009E73')
pp

vp1_phylo <- as.phylo(vp1_tree)
full_phylo <- as.phylo(full_genome_tree)

# compute RobinsonFoulds distance - not currently working
# VisualizeMatching(InfoRobinsonFoulds, full_phylo, vp1_tree_matching_full_tree, Plot=TreeDistPlot)

blue_clade <- tree_subset(vp1_tree_annotated, 386)
ggtree(blue_clade, mrsd="2024-02-01") + 
  geom_tiplab(aes(subset=label %in% c("28AO", "26HB"), label=label), fontface="bold", size=2, offset=.1) +
  geom_tiplab(size=2, offset=.1) +
  theme_tree2() +
  geom_tippoint(aes(color=Country))
export(blue_clade %>% select(label, Country, year), "blue_clade.csv")

pink_clade <- tree_subset(vp1_tree_annotated, 328)
ggtree(pink_clade, mrsd="2024-02-01") + 
  theme_tree2() +
  geom_tiplab(size=1, offset=.1) +
  geom_tippoint(aes(color=Country))
export(pink_clade %>% select(label, Country, year), "pink_clade.csv")
