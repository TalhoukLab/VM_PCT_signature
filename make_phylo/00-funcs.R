## this file will have all functions needed to do aim 1 analysis 
library(stringr)
library(phyloseq)
library(dplyr)
library(phangorn)

generate_otu_table <- function(otus){
  otus <- as.data.frame(otus)
  otus[] <- lapply(otus, as.numeric)
  rownames(otus) <- sapply(rownames(otus), function(x) {
    if (x %in% c("EC_S72", "NTC_S73", "EC_S73", "PC_S72")) {
      return(x)
    }
    base <- sub("_.*", "", x)
    num <- sub("^X", "", base)  
    parts <- strsplit(num, "\\.")[[1]]
    num_fmt <- sprintf("%03d", as.integer(parts[1]))
    dec <- if (length(parts) > 1) parts[2] else "0"
    paste0("OVRST", num_fmt, "-", dec)
  })
  colnames(otus) <- sub("^X", "", colnames(otus))  
  return(otus)
}

generate_tax_table <- function(tax_s, filtered){
  tax_s <- data.frame(tax_s[-1, ])
  tax_s <- tax_s %>%
            mutate(across(everything(), ~ str_replace_all(.x, ";\\s+", ";")))
  tax_s <- data.frame(tax_s %>%
              tidyr::separate_wider_delim(
              cols = Taxon,
              delim = ";",
              names = c("domain", "phylum", "class", "order", "family", "genus", "species"),
              too_few = "align_start"))
  rownames(tax_s) <- tax_s[,1]  
  tax_s <- tax_s[,-1] 
  tax_s <- tax_s %>%
                  mutate(across(everything(), as.character)) %>%  
                  mutate(across(everything(), ~ str_remove_all(.x, "[a-z]__"))) %>%     
                  mutate(across(everything(), ~ 
                                  ifelse(.x == "" | .x == "-" |
                                          str_detect(tolower(.x), "uncultured|unassigned|unidentified|metagenome"), NA, .x)))

  tax_s <- tax_s[ , !(colnames(tax_s) %in% "Confidence") ]
  if(filtered == TRUE){
    tax_s <- tax_s[!is.na(tax_s$phylum), ]
  }
  return(tax_s)
}

generate_phy_tree <- function(phylo_tree){
  rooted_tree <- phangorn::midpoint(phylo_tree)
  tree_phy <- phy_tree(rooted_tree)
  return(tree_phy)
}

plot_library_sizes <- function(phylo_obj){
  df <- as.data.frame(sample_data(phylo_obj)) 
  df$LibrarySize <- sample_sums(phylo_obj)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=control)) + geom_point()
}
