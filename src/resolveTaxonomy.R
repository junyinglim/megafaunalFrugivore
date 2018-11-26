## TAXONOMY TOOLS
library(plyr)
library(stringr)

frug.dir <- "~/Dropbox/Projects/2019/palms/data/frugivores/"
main.dir <- "~/Dropbox/Projects/2019/palms/projects/megafaunalFrugivore/data/"
iucnTaxonRef <- read.csv(file.path(frug.dir, "iucnMammalTaxonList.csv"))
iucnTaxonRef$SpecName <- paste(iucnTaxonRef$Genus, iucnTaxonRef$Species)

createSynTable <- function(x){
  if(x$Synonyms == ""){
    name <- x$SpecName
    accName <- x$SpecName
  } else {
    name <- c(x$SpecName, str_split(x$Synonyms, pattern = "\\|")[[1]])
    #name <- str_trim(name)
    accName <- x$SpecName
  }
  data.frame("name" = name, "accName" = accName)
}
iucnSynTab <- ddply(.data = iucnTaxonRef, .variables = .(SpecName), .fun = createSynTable)
iucnSynTab$name <- gsub(iucnSynTab$name, pattern = "  ", replacement = "")
#saveRDS(iucnSynTab, file = file.path(data.dir, "iucnSynTab.rds"))




resolveTaxonomy <- function(x, ref, syn.col, acc.col){
  # Resolve taxonomy based on a reference data frame
  #
  # Args:
  #   x = string or vector of species name for evaluation
  #   ref = data frame containing taxonomic references
  #   syn.col = column name in ref dataframe with synonyms
  #   acc.col = column name in ref dataframe with accepted names
  #
  # Returns:
  #   vector with accepted names. Names with no match are returned as NA
  #
  # Example:
  # y <- data.frame(syn = c("test1", "test2"), acc = c("Test", "Test"))
  # x <- c("test1", "test2", "test3")
  # x = x; ref = y; syn.col = "syn"; acc.col = "acc"
  return(as.vector(ref[[acc.col]][match(x, ref[[syn.col]])]))
}