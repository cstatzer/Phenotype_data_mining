set.seed(42)
library(tidyverse)
library(ggalluvial)
library(cowplot)
library(UpSetR)
library(homologene)
library(eulerr)
library(broom)
library(ggforce)
library(igraph)
library(tidygraph)
library(ggraph)
library(STRINGdb)
library(xlsx)
library(viridis)
library(openxlsx)
library(ggplotify)
library(conflicted)
library(ggpubr)
library(ggrepel)
conflict_prefer("select",winner = "dplyr")
conflict_prefer("filter",winner = "dplyr")
conflict_prefer("group_by",winner = "dplyr")
conflict_prefer("saveWorkbook",winner = "openxlsx")

# Updated homology database:
if(FALSE) {
  homologene_db_20200327 = updateHomologene()
  write_rds(path = "./data/homologene_db_20200327.rds",x = homologene_db_20200327)
}
homologene_db_20200327 <- read_rds("./data/homologene_db_20200327.rds")

# General settings
apply_grouping <- TRUE
dir_output <- "./output"
remap_diopt <- FALSE
reimport_string <- FALSE
resave_string <- FALSE
tax_human <- 9606
tax_mouse <- 10090
tax_worm <- 6239
tax_fly <- 7227
tax_fish <- 7955


# Plot settings
pal_manual <- list(
  "pal_Cyril1" = c("#000080","#b03060","#006400","#00ff00","#ff00ff","#ff0000","#ffd700","#00bfff","#acb54c","#b57b4c","#4cb5af","#724cb5"),
  "pal_Cyril2" = c("#4363d8","#3cb44b","#e6194b","#911eb4","#9a6324","#f58231","#000075","#808000","#ffe119","#bcf60c","#46f0f0","#f032e6"),
  "pal_subway" = c("#EE473E","#A9BD38","#0172BE","#C71F65","#A26337","#03B498","#8D7EB9","#D0A95F","#04B0EF","#F7A527","#C7BFB2","#1AB268"),
  "colorblind" = c("#8F0104","#0368D3","#0A4C4B","#380686","#874B10","#8ACEC8","#BF66F5","#E06A11","#BF89A0","#73B2FC","#2DF628","#EBBAD7","#AEDCFA","#F9FF93"),
  "deuteranopia" = c("#C8433A","#4B7AB3","#2F5D53","#9B251C","#1F4C92","#213C39","#7BA73D","#CD7EBD","#64A996","#E1C59A","#84E9E0","#DDB7C2","#D4F7E1","#D7FEFA","#FAF3F5"),
  "protanopia" = c("#CF402B","#697BBC","#2F635F","#273517","#374992","#6B212F","#ED6C52","#64A6C0","#58A296","#EFAB98","#D9D8E7","#78D2C3","#F9E8DB","#F2F0FB","#B2FEF3"),
  "tritanopia" = c("#B02C54","#1D5399","#DA5E41","#3D8C45","#5E5782","#EB6545","#79ABD3","#977EDD","#EEAE99","#BFEE74","#C9C5E2","#F4FA42","#F1F5FE","#EEF4CE"),
  "viridis" = viridis::viridis(12)
)
col_pal <- pal_manual[["viridis"]]


species_color = c("hub" = "black", 
                  "Danio rerio" = "#3E2369", 
                  "Homo sapiens" = "#225F7F", 
                  "Mus musculus" = "#02987E", 
                  "Mus musculus domesticus" = "#02987E", 
                  "Rattus norvegicus" = "#AED942", 
                  "Drosophila melanogaster" = "#3495eb",
                  "Caenorhabditis elegans" = "#872771", 
                  "Canis lupus familiaris" = "#60b347",
                  "Gallus gallus" = "#63327d",
                  "Bos taurus" = "#7d9970",
                  "Saccharomyces cerevisiae" = "#cf9013",
                  "shared" = "#000000"  #"#DD3A32"
                  )

species_abbreviations <- tibble(
  subject_taxon_label = c("Homo sapiens","Mus musculus","Rattus norvegicus","Drosophila melanogaster","Canis lupus familiaris", "Bos taurus","Saccharomyces cerevisiae","Danio rerio","Caenorhabditis elegans","Gallus gallus"),
  species_abbreviation_latin = c("hs","mm","rn","dm","clf","bt","sc","dr","ce","gg"),
  species_abbreviation = c("human","mouse","rat","fly","dog","cow","yeast","fish","worm","chicken")
  
)



species_order <- c("Homo sapiens", "Mus musculus", "Bos taurus","Canis lupus familiaris", "Rattus norvegicus",  "Mus musculus domesticus","Danio rerio", "Drosophila melanogaster", "Caenorhabditis elegans","Saccharomyces cerevisiae S288C", 
                   "Bos grunniens", "Bubalus bubalis", "Capra hircus","Coturnix japonica", "Equus caballus","Felis catus","Gallus gallus",              
                   "Macaca mulatta","Mammuthus primigenius","Meleagris gallopavo","Mesocricetus auratus", "Mustela putorius furo","Neovison vison", "Numida meleagris",
                   "Oryctolagus cuniculus","Oryzias latipes", "Ovis aries","Panthera tigris tigris","Peromyscus polionotus","Phoenicopterus ruber","Sciurus carolinensis",
                   "Sus scrofa","Ursus americanus","Vulpes lagopus","Vulpes vulpes",
                   "NCBITaxon:27706","NCBITaxon:32536","NCBITaxon:38666","NCBITaxon:43597","NCBITaxon:61402","NCBITaxon:68352","NCBITaxon:70340","NCBITaxon:87177","NCBITaxon:9690"
                   )


main_color <- "#f45042"
# matrisome colors
matrisome_color <- list()
matrisome_color$division <- list("core" = "#1336F4",
                                 "associated" = "#DB8530")
matrisome_color$division_lighted <- list("core" = "#8594EF",
                                         "associated" = "#DEBB9B")
matrisome_color$division_all <- list("Core matrisome" = "#1336F4",
                                     "Matrisome-associated" = "#DB8530",
                                     "Nematode-specific core matrisome" = "#8b9aef",
                                     "Nematode-specific matrisome-associated" = "#ce9e6d")
matrisome_color$category <- list("Collagens" = "#4BACFB", 
                                 "ECM Glycoproteins" = "#432E69", 
                                 "Proteoglycans" = "#8AF8ED",
                                 "ECM Regulators" = "#BFAD9A",
                                 "ECM-affiliated Proteins" = "#F4B076",
                                 "Secreted Factors" = "#F9CFEC",
                                 "Cuticlin" = "#D9F5D9")

matrisome_color$category_lighted <- list("Collagens" = "#4BACFB", 
                                         "ECM Glycoproteins" = "#845dcc", 
                                         "Proteoglycans" = "#8AF8ED",
                                         "ECM Regulators" = "#BFAD9A",
                                         "ECM-affiliated" = "#F4B076",
                                         "Secreted Factors" = "#F9CFEC",
                                         "Cuticlin" = "#D9F5D9")

palette_category_viridis <- c(
  "Collagens" = "#440154FF",
  "ECM Glycoproteins" = "#424186FF",
  "ECM Regulators" = "#2B758EFF",
  "ECM-affiliated Proteins" = "#20A386FF",
  "Proteoglycans" = "#6CCD5AFF",
  "Secreted Factors" = "#EBE51AFF"
)

palette_division_inferno <- c("Core matrisome" = "#7f068d",
                              "Matrisome-associated" = "#d24e71")



# matrisome_color$division <-  lapply(matrisome_color$division, function(x) adjustcolor(col = x,alpha.f = 0.5))
matrisome_color$category <-  lapply(matrisome_color$category, function(x) adjustcolor(col = x,alpha.f = 0.5))




expanded_letters <- c("",LETTERS) %>% map(~paste0(.x,LETTERS)) %>% flatten()





# Taxonomy replacements:
taxonomy_replacements <- c(
  "NCBITaxon:32536" = "Acinonyx jubatus",
  "NCBITaxon:38666" = "Chaetodipus intermedius",
  "NCBITaxon:43597" = "Holbrookia maculata",
  "NCBITaxon:61402" = "Puma yagouaroundi",
  "NCBITaxon:68352" = "Aspidoscelis inornatus",
  "NCBITaxon:70340" = "Anser caerulescens caerulescens",
  "NCBITaxon:87177" = "Coereba flaveola",
  "NCBITaxon:9690"  = "Panthera onca"
)
ecm_categories_ordered <- c("Collagens","ECM Glycoproteins","ECM Regulators","ECM-affiliated Proteins","Proteoglycans","Secreted Factors")



