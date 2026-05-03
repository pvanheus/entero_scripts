if (!require('pacman', character.only=TRUE)) {
  install.packages('pacman', dependencies = TRUE)
  library('pacman', character.only = TRUE)
}

pacman::p_load(ape)
vp1_tree = ape::read.tree(file = "data/vp1_tree.nwk")
ape::write.nexus(vp1_tree, file = "data/vp1_tree.nex")

full_genome_tree = ape::read.tree(file = "data/full_genome_tree.nwk")
ape::write.nexus(full_genome_tree, file = "data/full_genome_tree.nex")
