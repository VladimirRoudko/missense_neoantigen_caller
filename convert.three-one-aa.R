#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
convert <- function(l) {
  
  map <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  names(map) <- c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln",
                  "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
                  "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
  
  sapply(strsplit(l, "(?<=[A-Za-z]{3})", perl = TRUE),
         function(x) paste(map[x], collapse = ""))
}
convert(args[1])

