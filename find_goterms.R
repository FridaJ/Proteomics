if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install() #To install latest version if multiple versions are installed
BiocManager::install(c("dplyr", "tidyverse", "readxl", "stringr"))

library(tidyverse)
library(readxl)
library(dplyr)
library(stringr)

setwd("<dir>")



#####----- USER INPUT -----#####

inp_file <- "my_data.xlsx"
rows_to_skip <- 0
#patterns <- c("antibacterial", "defense response", "antimicrobial", "immune response", "bacterium", "neutrophil degranulation")
patterns <- "neutrophil degranulation"
col_name <- my_data$`Gene ontology (biological process)`

out_file <- "my_data_nd_raw.xlsx"



#####----- READ FILE -----#####

my_data_go <- read_xlsx(inp_file, skip = rows_to_skip)


#####----- Method for finding substrings in GO-terms -----#####

# Function that looks through a list of string lists (col_name), searching for ontology terms (patterns) and
# returns a vector of unique row numbers corresponding to the proteins with the ontology terms.

find_goterms <- function(col_name, patterns) {
  
  proteins <- list()
  for (pattern in patterns) {
    occurences <- grep(pattern, col_name)
    proteins <- c(proteins, occurences)
    proteins <- unique(proteins)
    proteins <- unlist(proteins)
  }
  return(proteins)
}

#####----- FIND ontology terms and save proteins to file -----#####

GO_rows <- find_goterms(col_name = col_name, patterns = patterns)
my_proteins <- my_data_go[GO_rows,]
xlsx::write.xlsx(my_proteins, file = out_file)

# Go do analysis using the proteins in the dataframe!

