flowerPhenoSingle <- function(red,yellow) {
  # if two readings, randomly sample
  # red <- Red
  # yellow <- Yellow
  red <- red[which(!is.na(red))]
  yellow <- yellow[which(!is.na(yellow))]
  pheno <- 0
  if (length(red)>1) {
    red <- as.numeric(sample(as.character(red),1))
  }
  if (length(yellow)>1) {
    yellow <- as.numeric(sample(as.character(yellow),1))
  }
  if (red=="NA" | yellow=="NA" | red=="noData" | yellow=="noData" | is.na(red) | is.na(yellow) | red==-9 | yellow==-9) {
    red <- -9
    yellow <- -9
    return(-9)
  }
  if (red!=-9 & yellow!=-9) {
    # Yellow
    if ((0 <= red & red < 1.5) & (2.0 <= yellow & yellow <= 4.0)) {
      pheno <- "Y"
    }
    # White phenotype
    if ((0 <= red & red < 1.5) & (0 <= yellow & yellow <= 1.5)) {
      pheno <- "W"
    }
    # Weak Orange
    if ((1.5 <= red & red <= 2.75) & (2.0 <= yellow & yellow <= 4.0)) {
      pheno <- "WO"
    }
    # Full Orange - note updated 2019 changed yellow > 1.5 
    if ((3.0 <= red & red <= 5.0) & (1.5 < yellow & yellow <= 5.0)) {
      pheno <- "FO"
    }
    # Weak Red 
    if ((1.5 <= red & red <= 2.75) & (0 <= yellow & yellow <= 1.5)) {
      pheno <- "WR"
    }
    # Full Red 
    if ((3.0 <= red & red <= 5.0) & (0 <= yellow & yellow <= 1.5)) {
      pheno <-"FR"
    }
    # Missing
    if ((red == "") | (yellow == "")) {
      pheno <- -9
    }
  }
  return(pheno)
}