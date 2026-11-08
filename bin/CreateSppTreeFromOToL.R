# Load required package
library(ape)
library(dplyr) # For data manipulation

root_dir <- "~/projects/TaxonSampling/development/ParseOToL/"
setwd(root_dir)
set.seed(42)

# Functions

# sanitize terminal OTUs (tip labels)
sanitize_ottid_labels <- function(label) {
  label <- gsub("[']", "", label)  # Remove single quotes (e.g. 'Solms-laubachia linearis ott816127')
  label <- gsub("_ott[0-9]+$", "", label)  # Remove "_ott" followed by digits at the end of the label
  label <- gsub(" ott[0-9]+$", "", label)  # Remove " ott" followed by digits at the end of the label (e.g. 'Schizomeria_gorumensis_ott6106743')
  return(label)
}

# identify specific internal labels and sanitize/remove them
sanitize_mrca_labels <- function(label) {
  label <- gsub("[']", "", label)  # Remove single quotes (e.g. 'Starmera/Candida clade ott5674058')
  
  # Remove labels matching "mrcaott<number>ott1<number>"
  label <- ifelse(grepl("^mrcaott[0-9]+ott[0-9]+$", label), NA, label)
  return(label)
}

# mapping function of names to TaxID
map_taxid <- function(taxa_name, mapping_df) {
  taxid <- mapping_df %>% filter(Name == taxa_name) %>% pull(TaxID)
  if (length(taxid) > 0) {
    return(paste0(taxa_name, "[TaxID:", taxid, "]"))
  } else {
    return(taxa_name) # Return the name unchanged if no TaxID is found
  }
}


# Read the original Newick tree
# Downloaded from https://files.opentreeoflife.org/synthesis/opentree15.1/opentree15.1tree.tgz
# untar, expand and point to this directory

tree_file <- "data/opentree15.1_tree/labelled_supertree/labelled_supertree_ottnames.tre"  # Replace with your file path
tree <- read.tree(tree_file)

tip_labels <- tree$tip.label

i <- 1
for (label in tip_labels) {
  print(i)
  label <- sanitize_ottid_labels(label)
#  print(label)
  i <- i + 1
}

# Rename terminal nodes (tips) by cleaning up OTTIDs
tree$tip.label <- unname(sapply(tree$tip.label, sanitize_ottid_labels))

# Identify and remove duplicate terminal labels
table(!duplicated(tree$tip.label))
unique_tips <- !duplicated(tree$tip.label)
tree <- drop.tip(tree, tree$tip.label[!unique_tips])

# Remove branch lengths, as these are not needed (set to NULL)
tree$edge.length <- NULL

# sanitize node labels as well
tree$node.label <- unname(sapply(tree$node.label, sanitize_ottid_labels))

# sanitize a little bit more to remove internal OToL IDs matching "mrcaott<number>ott1<number>
if (!is.null(tree$node.label)) {
  tree$node.label <- unname(sapply(tree$node.label, sanitize_mrca_labels))
}


# Save the intermediate cleaned tree to a new file
output_file <- "results/proc_trees/OToL_taxon_names.tre" 
write.tree(tree, file = output_file)
# load, if needed
#tree <- read.tree(tree, file = output_file)


##Add NCBI TaxIDs

# Step 2: Load NCBI Taxonomy data (names.dmp)
# Ensure names.dmp is in the working directory
names_dmp <- read.delim("data/NCBI_Taxonomy/names.dmp", header = FALSE, sep = "|", stringsAsFactors = FALSE, strip.white = TRUE)

# Clean up column names for readability
colnames(names_dmp) <- c("TaxID", "Name", "Unique_Name", "Name_Class", "Extra")

# Filter for scientific names
scientific_names <- names_dmp %>% 
  filter(Name_Class == "scientific name") %>% 
  select(TaxID, Name)

scientific_names$Name <- gsub(" ", "_", scientific_names$Name)

tree_labels <- unique(c(tree$node.label, tree$tip.label))

#subset to reduce search space
tmp <- scientific_names[scientific_names$Name %in% tree_labels,]

# Step 3: Map taxa names in the tree to TaxIDs

# Apply the mapping to all taxa in the tree
tree$tip.label <- unname(sapply(tree$tip.label, map_taxid, mapping_df = tmp))
tree$node.label <- unname(sapply(tree$node.label, map_taxid, mapping_df = tmp))

# Step 4: Write the updated tree to a new file
output_file <- "results/proc_trees/OToL_taxon_names_with_TaxIDs.tre"
write.tree(tree, file = output_file)

# load, if needed
tree <- read.tree(file = output_file)

tip_labels <- tree$tip.label

# Step 2: Extract TaxIDs from the tip labels
extract_taxid <- function(label) {
  if (grepl("\\[TaxID-", label)) {
    label <- sub(".*\\[TaxID-([0-9]+)\\].*", "\\1", label)
    label
  } else {
    NA # Tips with no TaxID
  }
}

# Apply the extraction function to all tip labels
taxids <- sapply(tree$tip.label, extract_taxid)

# Step 3: Create a recursive function to count species for each node
get_species_count <- function(node, tree) {
  if (node <= length(tree$tip.label)) {
    # It's a tip (leaf)
    return(1)
  } else {
    # It's an internal node, sum counts of all descendants
    children <- tree$edge[tree$edge[, 1] == node, 2]
    return(sum(sapply(children, get_species_count, tree = tree)))
  }
}

# Step 4: Calculate species counts for all nodes
node_counts <- sapply(1:(length(tree$tip.label) + tree$Nnode), get_species_count, tree = tree)

# Step 5: Create a data frame of TaxIDs and their species counts
taxid_species_counts <- data.frame(
  Node = 1:(length(tree$tip.label) + tree$Nnode),
  TaxID = c(tree$taxid, rep(NA, tree$Nnode)), # NA for internal nodes without TaxIDs
  SpeciesCount = node_counts
)

# Step 6: Save the results to a CSV file
output_file <- "Tree_TaxID_species_counts.csv"
write.csv(taxid_species_counts, output_file, row.names = FALSE)

#cat("Species counts based on the tree structure saved to:", output_file, "\n")
