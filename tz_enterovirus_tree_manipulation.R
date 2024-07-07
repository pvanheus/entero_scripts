if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install(c("ggtree", "treeio", "ggtreeExtra", "phangorn", "TreeDist"))

library('librarian')
librarian::shelf(ggplot2, ape, dplyr, phangorn, treeio, tidytree, ggtree, ggtreeExtra, rio, TreeDist)

# load full genome and VP1 region trees

vp1_tree = read.nexus("vp1_tree_with_dates.nex")
full_genome_tree = read.nexus("full_length_tree_with_dates.nex")
metadata = import("all_microreact_format_metadata.csv")
vp1_tree_annotated = full_join(vp1_tree, metadata, join_by(label == Accession))

full_tree <- ggtree(vp1_tree_annotated, mrsd="2024-02-01") + 
  theme_tree2() +
  geom_tippoint(aes(color=Country))

# plot VP1 tree with clades of interest highlighted
tree_plot <- ggtree(vp1_tree_annotated, mrsd="2024-02-01") + 
  theme_tree2() +
  geom_tippoint(aes(color=Country)) +
  geom_highlight(node=382, fill="blue") +
  geom_cladelab(382, "A", offset=-35, offset.text=-18) +
  # geom_strip("MF419263.1", "MK989722.1", offset=-35, label="A1") +
  geom_highlight(node=323, fill="pink") +
  geom_cladelab(323, "B", offset=-34.65, offset.text=-18, extend=-2) +
  geom_treescale(x=2010, y=20, width=10, offset=5)
tree_plot

# compute subsect of VP1 tree that matches full genome tree

full_genome_tips <- d1 %>% filter(isTip == TRUE) %>% select(label)

not_matching_tips_nodes <- d2 %>% filter(isTip == TRUE) %>% filter(!label %in% full_genome_tips$label) %>% select(node)

# compute tanglegram of VP1 subset and full genome tree

vp1_tree_matching_full_tree <- drop.tip(vp1_tree, not_matching_tips_nodes$node)

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
