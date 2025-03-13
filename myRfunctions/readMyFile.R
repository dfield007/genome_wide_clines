readMyFile <- function(folderStructure,folderLocation,fileName) {
  # folderLocation <- "folder_base_antSpec"
  # dataSet <- 
  # fileName <- "Antspec2009NEW.csv"
  # folderFind <- "~/Dropbox/work/PostdocIST/projects_Rmd/SNPs/HZanalyses/output/pedigree/results/"
  folderFind <- as.vector(folderStructure[folderStructure[,"FolderName"]==folderLocation,"Location"])
  fileToRead <- read.csv(paste(c(folderFind,fileName),collapse=""),header=TRUE,fill=FALSE,na.strings = "?",stringsAsFactors = F)
  return(fileToRead)
}