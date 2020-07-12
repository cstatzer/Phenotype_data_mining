# # Brainstorming:
# DONE:  Tables about how many phenotypes per species and so on.
#
# Phenotype perspective:
# DONE: - What are the most frequent phenotypes for each species --> which has the most genes associated with it. Illustrate this as a histogram --> "This is the phenome space"
# DONE: - In the same histogram color each bar in which at least X% of the genes are matrisome associated. --> "These are the phenotypic matrisome footprint"
# DONE: - Make a histogram with "% Matrisome involvement" on the y axis to illustrate which phenotypes are most controlled by the matrisome.
# 
# Gene perspective:
# - See which genes are involved in the most phenotypes overall --> histogram
# - See which matrisome genes are involved in the most phenotypes overall
# - Group the matrisome by category and division to see which has the highest impact.
# - Combine the species by homology mapping and make a stacked histogram for the matrisome genes to see which genes are important across species.
# 
# Homology to human
#  - Make a network where a human gene is in the middle, show the connections to species it has a homolog in and then show the phenotypes connected to the species this gene is involved in.




# Lessons learned: for stacked barplots use a regular axis with zoom and not a log axis since it distorts the ratio.


source("./config.R")
source("./helper.R")
# detach("package:biomaRt", unload = TRUE)

##################################################################
####### Import ##########
##################################################################
##### import monarch data
gene_pheno <- read_tsv(file = "./data/gene_phenotype.all.tsv") %>% separate(col = subject,into = c("prefix","geneid"),sep = ":") %>% separate(col = subject_taxon,into = c("prefix_taxon","taxon"),sep = ":")
gene_pheno <- gene_pheno %>% 
  mutate(phenotype_original = TRUE, phenotype_as_imported = object_label) %>%
  mutate(subject_taxon_label = ifelse(subject_taxon_label == "Saccharomyces cerevisiae S288C","Saccharomyces cerevisiae",subject_taxon_label)) %>%
  mutate(object_label = ifelse(object_label == "phenotype","Uncharacterized phenotype",object_label)) %>%
  mutate(unique_id = paste(subject_taxon_label,subject_label,phenotype_as_imported)) %>%
  filter(subject_taxon_label != "organism") %>%
  filter(!duplicated(unique_id))
dim(gene_pheno)
# thresholding required by reviewers:

gene_pheno <- gene_pheno %>% 
  filter(relation_label == "has phenotype") %>% # "causes condition", "contributes to condition" and "pathogenic_for_condition" are excluded. They make up only 20 rows of the table.
  filter(!grepl(pattern = "without traceable support",x = evidence_label)) %>% # remove evidence labels which do not provide us with enough confidence
  filter( !is.na(source) | !is.na(evidence_label) )
dim(gene_pheno)


if(apply_grouping){
  imported_grouping_import <- read_csv(file = "./data/200423_phenotype_grouping.csv")
  out <- wrangle_grouping(imported_grouping_file = imported_grouping_import, gene_pheno_df = gene_pheno)
  gene_pheno <- out$output_grouped
  out$fail_to_group %>% write_csv(path = paste0(dir_output,"/phenotype_grouping/fail_to_group.csv",collapse = ""))
  out$grouping_statistics %>% write_csv(path = paste0(dir_output,"/phenotype_grouping/grouping_statistics.csv",collapse = ""))
  out$grouping_statistics %>% rename(`Assigned simplified phenotype group description` = object_label, `Original phenotype` = phenotype_as_imported, `Number of occurences of original phenotype (across taxa)` = n) %>% write_csv(path = paste0(dir_output,"/Paper_figures/Supplement_tables/Sx2.csv"))
  grouping_mentions_in_text <- list(n_species_gene_phenotype_associations = nrow(out$output_grouped), n_grouped_phenotypes = out$grouping_statistics %>% pull(object_label) %>% n_distinct(na.rm = TRUE))
  # Manual fusion of categories:
  gene_pheno <- gene_pheno %>% mutate(object_label = case_when(object_label %in% c("Short stature","Altered body size","abnormal(ly) decreased length whole organism") ~ "Altered body size", 
                                                               object_label %in% c("clear") ~ "Increased tissue transparency", 
                                                               
                                                               TRUE ~ object_label))
}

write_csv(x = gene_pheno %>% select(-prefix_taxon,-relation,-relation_label,-evidence,-is_defined_by,-qualifier,-phenotype_original,-unique_id,-prefix) %>% 
            rename (`Gene identifier` = geneid, `Gene name` = subject_label, `Taxon identifier` = taxon, `Taxon name` = subject_taxon_label, `Phenotype identifier` = object, `Simplified phenotype group` = object_label, `Association evidence` = evidence_label, `Association source (if directly provided)` = source,`Phenotype name` = phenotype_as_imported) %>% 
            select(`Taxon name`, `Taxon identifier`, `Gene name`,`Gene identifier`,`Phenotype name`, `Phenotype identifier`, `Association source (if directly provided)`, `Association evidence`, `Simplified phenotype group`) %>% arrange(desc(`Taxon name`,`Gene name`,`Phenotype name`)),paste0(dir_output,"/Paper_figures/Supplement_tables/Sx1.csv"))
gene_pheno <- gene_pheno %>% mutate(taxon = ifelse(taxon == "559292","4932",taxon)) # Saccharomyces cerevisiae S288C is not found in homologene thus use "bakers yeast" in general: 4932


# split by taxon
all_pheno <- gene_pheno %>% split(x = .,f = .$subject_taxon_label)
all_species <- gene_pheno %>% pull(subject_taxon_label) %>% unique()
# assemble matrisome by species
ce_mat <- read_csv(file = "./data/ce_matrisome.csv") %>% select(method:systematic_wikigene,combined_orthology) %>% rename(geneid = WormBaseID) %>% mutate(gene_symbol = wikigene) %>% mutate(division_nematode = division, division = ifelse(division %in% c("Core matrisome","Nematode-specific core matrisome"),"Core matrisome","Matrisome-associated"))
mm_mat_import <- read_csv(file = "./data/mm_matrisome.csv") %>% select(division:Refseq_IDs) %>% separate(col = MGI_IDs,into = c("prefix","geneid"),sep = ":") %>% mutate(gene_symbol = `Gene Symbol`) 
mm_mat <- wrangle_hom_mouse_to_human(df = mm_mat_import)
hs_mat <- read_csv2(file = "./data/hs_matrisome.csv") %>% select(division:Refseq_IDs) %>% rename(geneid = HGNC_IDs) %>% mutate(geneid = as.character(geneid)) %>% mutate(gene_symbol = `Gene Symbol`) %>% mutate(combined_orthology = `Gene Symbol`)
fish_mat <- read_csv(file = "./data/fish_matrisome.csv") %>% select(division:`UniProt IDs`,combined_orthology)  %>% mutate(gene_symbol = `Gene Symbol`) 
fly_mat <- read_csv(file = "./data/fly_matrisome.csv") %>% select(division:`UniProt ID`,combined_orthology)  %>% mutate(gene_symbol = `Gene Name`) 
ecm_pheno_worm <- full_join(x = ce_mat,y = gene_pheno, by = "geneid") %>% mutate(observed_in = case_when(is.na(division) ~ "phenome_only", is.na(subject_taxon_label) ~ "matrisome_only", !is.na(division) & !is.na(subject_taxon_label) ~ "both", TRUE ~ "other")) %>% filter(subject_taxon_label == "Caenorhabditis elegans" | !is.na(division)) %>% mutate(subject_taxon_label = "Caenorhabditis elegans")
ecm_pheno_mouse <- full_join(x = mm_mat,y = gene_pheno, by = "geneid") %>% mutate(observed_in = case_when(is.na(division) ~ "phenome_only", is.na(subject_taxon_label) ~ "matrisome_only", !is.na(division) & !is.na(subject_taxon_label) ~ "both", TRUE ~ "other")) %>% filter(subject_taxon_label == "Mus musculus" | !is.na(division))  %>% mutate(subject_taxon_label = "Mus musculus")
ecm_pheno_human <- full_join(x = hs_mat,y = gene_pheno, by = "geneid") %>% mutate(observed_in = case_when(is.na(division) ~ "phenome_only", is.na(subject_taxon_label) ~ "matrisome_only", !is.na(division) & !is.na(subject_taxon_label) ~ "both", TRUE ~ "other")) %>% filter(subject_taxon_label == "Homo sapiens" | !is.na(division))  %>% mutate(subject_taxon_label = "Homo sapiens")
ecm_pheno_fish <- full_join(x = fish_mat,y = gene_pheno, by = "geneid") %>% mutate(observed_in = case_when(is.na(division) ~ "phenome_only", is.na(subject_taxon_label) ~ "matrisome_only", !is.na(division) & !is.na(subject_taxon_label) ~ "both", TRUE ~ "other")) %>% filter(subject_taxon_label == "Danio rerio" | !is.na(division))  %>% mutate(subject_taxon_label = "Danio rerio")
ecm_pheno_fly <- full_join(x = fly_mat,y = gene_pheno, by = "geneid") %>% mutate(observed_in = case_when(is.na(division) ~ "phenome_only", is.na(subject_taxon_label) ~ "matrisome_only", !is.na(division) & !is.na(subject_taxon_label) ~ "both", TRUE ~ "other")) %>% filter(subject_taxon_label == "Drosophila melanogaster" | !is.na(division))  %>% mutate(subject_taxon_label = "Drosophila melanogaster")
ecm_pheno <- list(worm = ecm_pheno_worm, mouse = ecm_pheno_mouse, human = ecm_pheno_human, fish = ecm_pheno_fish, fly = ecm_pheno_fly)


# Direct Homolog Step 1#
# Matrisome ortholog conversion: df with species, gene, homolog, +(origin_species_gene_symbol, gene_symbol)
mat_species <- bind_rows(ce_mat,mm_mat,fish_mat,fly_mat,hs_mat) %>% 
  filter(!is.na(geneid)) %>% 
  filter(!geneid %in% c("-")) %>%
  filter(!duplicated(geneid)) %>%
  select(geneid, gene_symbol,combined_orthology)
ecm_homolog_conversion_table <- mat_species %>% 
  rowwise() %>% 
  mutate(homolog = str_split(combined_orthology, ";", simplify = TRUE)%>% str_squish() %>% list()) %>%
  unnest(cols = homolog) %>%
  select(-combined_orthology) %>%
  filter(!homolog %in% c("","-","0",NA)) %>%
  mutate(origin_species_gene_symbol = gene_symbol)
# Direct Homolog Step 1#

# ## Test grouping
# pheno_grouping_df <- read_csv(file = "./data/200118_phenotype_grouping.csv") %>% 
#   mutate(matched_substrings = matched_substrings %>% tolower()) %>%
#   rowwise() %>%
#   mutate(matched_substrings_split = map(matched_substrings, ~str_split(string = matched_substrings,pattern = ",")),
#          targets_trim = map(matched_substrings_split, ~.x %>% str_trim() %>% wrangle_remove_empty_strings()))
# 
# output_grouping <- wrangle_phenotype_grouping(pheno_df = gene_pheno, group_df = pheno_grouping_df)
# 
# output_grouping$grouping_statistics %>% select(object_label, `replacement table`) %>% unnest(cols = c(`replacement table`)) %>% 
#   write_csv(path = paste0(dir_output,"/phenotype_grouping/grouping_statistics.csv",collapse = ""))
# 
# fail_to_group <- output_grouping$pheno_grouped %>% filter(phenotype_original) %>% select(object_label,phenotype_as_imported,phenotype_original) 
# fail_to_group %>% write_csv(path = paste0(dir_output,"/phenotype_grouping/fail_to_group.csv",collapse = ""))
# ###### ONLY Phenotype grouping


##################################################################
####### Phenotype_perspective for each species ##########
##################################################################
top_n_all_phenotypes <- 50
all_pheno_occ <- all_pheno %>%
  map(~.x %>% group_by(object_label) %>% 
        summarize(count = n(),
                  species = subject_taxon_label %>% unique()) %>%
        arrange(desc(count)) %>%
        filter(!is.na(object_label), count >= 1) %>%
        head(top_n_all_phenotypes))

nrow_elements <- all_pheno_occ %>% map_dbl(~.x$object_label %>% n_distinct())
all_pheno_occ <- all_pheno_occ[order(-nrow_elements)]
all_pheno_occ_plot <- all_pheno_occ %>% 
  map(~.x %>% 
        ggplot(aes(x = reorder(object_label, count), y = count, fill = count)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(begin = 0.2, end=0.8,direction = 1) +
        # scale_y_log10() +
        labs(title = ifelse(unique(.x[["species"]]) %in% names(taxonomy_replacements), taxonomy_replacements[[unique(.x[["species"]])]] ,unique(.x[["species"]])), y = "# occurances", x = "Phenotype") +
        theme_classic() +
        coord_flip() +
        theme(axis.text.y = element_text(face = "italic",size = 20, angle = 0, hjust = 1),plot.title = element_text(face = "bold.italic", size = 28)))
all_pheno_occ_plot_names <- all_pheno_occ %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Phenotype_perspective/",.,"_all","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(all_pheno_occ_plot,all_pheno_occ_plot_names, ~ggsave(.x,filename = .y,height = 10, width = 7))
ggarrange(plotlist = all_pheno_occ_plot,labels = "AUTO",font.label = list(size = 36),
         ncol = 3,nrow = 4,align = "hv") %>%
  ggexport(
         filename = paste0(dir_output,"/Paper_figures/Supp_fig_n3/","Grid_of_phenotype_perspective","_grouping_",apply_grouping,".pdf",collapse = ""),width = 50, 
         height = 70, 
         units = "mm"
         )

top_n_ecm_phenotypes <- 30
ecm_pheno_occ <- ecm_pheno %>%
  map(~.x %>% 
        filter(observed_in == "both") %>%
        filter(!object_label %in% c("Uncharacterized phenotype","no abnormal phenotype detected")) %>%
        mutate(object_label = str_remove_all(string = object_label, pattern = "abnormal\\(ly\\) "),
               object_label = str_replace_all(string = object_label,pattern = ", complete penetrance",replacement = ""),
               object_label = case_when(object_label == "increased or absent threshold for auditory brainstem response" ~ "auditory brainstem response threshold",
                                        grepl(pattern = "first instar larval cuticle phenotype",object_label) ~ "larval cuticle phenotype",
                                        TRUE ~ object_label),
               object_label = object_label %>% str_to_sentence()) %>%
        group_by(object_label) %>% 
        summarize(count = n(),
                  species = subject_taxon_label %>% unique()) %>%
        arrange(desc(count)) %>%
        filter(!is.na(object_label), count > 1) %>%
        head(top_n_ecm_phenotypes))
pheno_ecm_names <- ecm_pheno_occ %>% map(~.x$object_label) %>% unlist() %>% unname() %>% unique()
pheno_ecm_palette <- setNames(rep("grey",length(pheno_ecm_names)),pheno_ecm_names)

if(apply_grouping){
  pheno_ecm_palette[c("Aging-related phenotype","Resistance phenotype")] <- col_pal[[1]]
  pheno_ecm_palette[c("Short stature","Altered body size")] <- col_pal[[2]]
  pheno_ecm_palette[c("Poor viability","Systemic health decrease")] <- col_pal[[3]]
  pheno_ecm_palette[c("Developmental phenotype","molt defect")] <- col_pal[[4]]
  pheno_ecm_palette[c("Sterility / Reduced reproduction potential","Reproductive system phenotype")] <- col_pal[[5]]
  pheno_ecm_palette[c("Bone phenotype","Cranial deformations","Facial deformations")] <- col_pal[[6]]
  pheno_ecm_palette[c("Angiogenetic and blood vessels formation","Blood phenotype","Cardiovascular impairments","Edema occurence")] <- col_pal[[7]]
  pheno_ecm_palette[c("Immune system phenotype","Inflamation phenotype")] <- col_pal[[8]]
  pheno_ecm_palette[c("Nervous and neural tissues affected","Brain phenotype","Complex brain phenotype")] <- col_pal[[9]]
  pheno_ecm_palette[c("Muscle tissue phenotype")] <- col_pal[[10]]
  pheno_ecm_palette[c("Kinestetic derived phenotypes","Limbs related phenotype","Major articulations","Podalic, digitalis phenotype")] <- col_pal[[11]]
  pheno_ecm_palette[c("Skin phenotype", "Connective tissue phenotype")] <- col_pal[[12]]
}
ecm_pheno_occ_plot <- ecm_pheno_occ %>% 
  map(~.x %>% 
         ggplot(aes(x = reorder(object_label, count), y = count,fill = object_label)) +
         geom_bar(stat = "identity") +
         # scale_y_log10() +
        scale_fill_manual(values = pheno_ecm_palette) +
         labs(title = paste(unique(.x[["species"]])), y = "# occurances (log 10)", x = "Phenotype") +
         theme_classic() +
        coord_flip() +
         theme(axis.text.y = element_text(face = "italic",size = 22, angle = 0, hjust = 1),
               axis.title = element_blank(), 
               legend.position = "none", 
               plot.title = element_text(face = "bold.italic", size = 27))
  )
ecm_pheno_occ_plot_names <- ecm_pheno_occ %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Phenotype_perspective/",.,"_ecm","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(ecm_pheno_occ_plot,ecm_pheno_occ_plot_names, ~ggsave(.x,filename = .y,height = 10, width = 7))

ggarrange(plotlist = ecm_pheno_occ_plot[c(3,2,4,6,5,1)],labels = c("A","B","C","","E","F"),font.label = list(size = 30),
          ncol = 2,nrow = 3,align = "hv") %>%
ggexport(
  filename = paste0(dir_output,"/Paper_figures/ecm_pheno_occ_plot/","Grid_of_ecm_phenotypes","_grouping_",apply_grouping,".pdf",collapse = ""),
  width = 27, 
  height = 29, 
  units = "mm"
)

### Export this figure als facetted plot
# plot_df <- ecm_pheno_occ %>% 
#   reduce(bind_rows) %>% 
#   mutate(species = factor(species,c("Homo sapiens","Mus musculus","Danio rerio","Drosophila melanogaster","Caenorhabditis elegans"))) %>%
#   group_by(species) %>%
#   arrange(count) %>%
#   mutate(order = row_number()) %>%
#   ungroup()
# ecm_pheno_occ_plot_export <- plot_df %>%
#         ggplot(aes(x = order, y = count, fill = object_label)) +
#         geom_bar(stat = "identity") +
#   facet_wrap(~species,scales = "free",drop = TRUE, nrow = 1) +
#   coord_flip() +
#   theme(axis.text.y = element_text(face = "italic",size = 5, angle = 0, hjust = 1)) +
#         # scale_y_log10() +
#         scale_fill_manual(values = pheno_ecm_palette) +
#         theme_classic() +
#   theme(axis.line.y = element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         legend.position = "bottom") 
#  ggsave(filename = paste(dir_output,"/Paper_figures/Fig1/ecm_pheno_occ_plot_export","_grouping_",apply_grouping,".pdf",sep = ""),plot = ecm_pheno_occ_plot_export,width = 14,height = 16)


##################################################################
####### Matrisome involvement in all phenotypes ##########
##################################################################
top_n_phenotypes <- 25
matrisome_fraction <- ecm_pheno %>%
  map(~.x %>% summarize(species = unique(subject_taxon_label),
                          genes_total = n_distinct(geneid),
                          genes_matrisome_and_phenome = n_distinct(geneid[observed_in == "both"]),
                          genes_only_matrisome = n_distinct(geneid[observed_in == "matrisome_only"]),
                          genes_only_phenome = n_distinct(geneid[observed_in == "phenome_only"]),
                          matrisome_backround = genes_matrisome_and_phenome / (genes_matrisome_and_phenome + genes_total)
                           ))

ecm_pheno_occ <- ecm_pheno %>%
  map(~.x %>% 
        group_by(object_label) %>% 
        summarize(count_overall = n(),
                  matrisome = sum(observed_in == "both"),
                  species = subject_taxon_label %>% unique(),
                  other_genes = count_overall - matrisome,
                  matrisome_fraction = matrisome / count_overall
                  ) %>%
        arrange(desc(count_overall)) %>%
        filter(!is.na(object_label), count_overall >= 1),
      ungroup())

pheno_occ_ecm_contribution <- map2(ecm_pheno_occ,matrisome_fraction, ~.x %>% 
                                     mutate(count_expected_by_background = count_overall * .y[["matrisome_backround"]],
                                            matrisome_enrichment = matrisome / count_expected_by_background) %>%
                                     arrange(desc(matrisome_enrichment)))
# plot
pheno_occ_ecm_contribution_plot <- pheno_occ_ecm_contribution %>% 
  map(~.x  %>%
        head(top_n_phenotypes) %>% 
        ggplot(aes(x = reorder(object_label, -count_overall), y = count_overall, fill = matrisome_enrichment)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_c(begin = 0, end=1,direction = 1) +
        scale_y_log10() +
        labs(title = paste("Most prevalent phenotypes in ", unique(.x[["species"]]), " and the contribution of the matrisome"),subtitle = "Phenotype mining using Monarch", y = "# occurances (log 10)", x = "Phenotype") +
        theme_classic() +
        theme(axis.text.x = element_text(face = "italic",size = 4, angle = 45, hjust = 1)))
pheno_occ_ecm_contribution_plot_names <- pheno_occ_ecm_contribution %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Matrisome_involvement/",.,"_ecm_contribution","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(pheno_occ_ecm_contribution_plot,pheno_occ_ecm_contribution_plot_names, ~ggsave(.x,filename = .y,height = 8, width = 13))

# plot
pheno_occ_ecm_absolute_top_plot <- pheno_occ_ecm_contribution %>% 
  map(~.x %>% 
        filter(count_overall > 10) %>%
        mutate(n_matrisome = matrisome) %>%
        arrange(desc(matrisome),desc(count_overall)) %>%
        head(top_n_phenotypes) %>% 
        mutate(object_label = factor(object_label,levels = unique(object_label))) %>%
        gather(key = "Matrisome membership",value = "value",-object_label,-count_overall,-species,-matrisome_enrichment,-count_expected_by_background,-n_matrisome) %>%
        mutate(`Matrisome membership` = factor(`Matrisome membership`,levels = c("other_genes","matrisome"))) %>%
        ggplot(aes(x = object_label, y = value, fill = `Matrisome membership`)) +
        geom_bar(stat = "identity", color = "black", size = 0.1) +
        scale_fill_manual(values = c("matrisome" = species_color[[unique(.x[["species"]])]],"other_genes" = "#e4f0e7")) +
        labs(title = paste("Most prevalent phenotypes in ", unique(.x[["species"]]), " and the contribution of the matrisome"),subtitle = "Phenotype mining using Monarch", y = "# Genes involved", x = "Phenotype") +
        theme_classic() +
        theme(axis.text.x = element_text(face = "italic",size = 8, angle = 60, hjust = 1), legend.position = "none") +
        facet_zoom(ylim = c(0, 100),zoom.size = 4))
pheno_occ_ecm_absolute_top_plot_names <- pheno_occ_ecm_contribution %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Matrisome_involvement/",.,"_ecm_most_contribution","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(pheno_occ_ecm_absolute_top_plot,pheno_occ_ecm_contribution_plot_names, ~ggsave(.x,filename = .y,height = 8, width = 13))

# plot
pheno_occ_ecm_mostenriched_top_plot <- pheno_occ_ecm_contribution %>% 
  map(~.x %>% 
        mutate(phenotype_magnitude = case_when(count_overall > 100 ~ "systemic phenotype",
                                               count_overall > 20 ~ "medium-sized phenotype",
                                               count_overall > 5 ~ "small phenotype",
                                               TRUE ~ "very small")) %>%
        mutate(n_matrisome = matrisome,
               phenotype_magnitude = factor(phenotype_magnitude,levels = c("very small","small phenotype","medium-sized phenotype","systemic phenotype"))) %>%
        group_by(phenotype_magnitude) %>%
        arrange(desc(matrisome_fraction),desc(count_overall)) %>%
        slice(1:top_n_phenotypes) %>% 
        ungroup() %>%
        mutate(object_label = factor(object_label,levels = unique(object_label))) %>%       
        ggplot(aes(x = object_label, y = matrisome_fraction, fill = count_overall)) +
        geom_bar(stat = "identity", color = species_color[[unique(.x[["species"]])]], size = 0.1) +
        scale_fill_viridis_c(begin = 0.5, end=1,direction = 1) +
        facet_wrap(. ~ phenotype_magnitude,scales = "free_x",ncol = 4) +
        labs(title = paste("Phenotypes with the highest (relative) ECM contribution ", unique(.x[["species"]]), ""),subtitle = "Phenotype mining using Monarch", y = "Degree of matrisome genes involvement in the phenotype [%]", x = "Phenotype") +
        theme_classic() +
        theme(axis.text.x = element_text(face = "italic",size = 4, angle = 45, hjust = 1)))
pheno_occ_ecm_mostenriched_top_plot_names <- pheno_occ_ecm_contribution %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Matrisome_involvement/",.,"_ecm_most_contribution","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(pheno_occ_ecm_mostenriched_top_plot,pheno_occ_ecm_mostenriched_top_plot_names, ~ggsave(.x,filename = .y,height = 8, width = 13))

# plot
pheno_occ_ecm_fraction_all_plot <- pheno_occ_ecm_contribution %>% 
  map(~.x %>% 
        filter(count_overall > 10,!object_label %in% c("Uncharacterized phenotype","no abnormal phenotype detected")) %>%
        mutate(n_overall= count_overall) %>%
        mutate(object_label = str_remove_all(string = object_label, pattern = "abnormal\\(ly\\) "),
               object_label = str_replace_all(string = object_label,pattern = ", complete penetrance",replacement = ""),
               object_label = case_when(object_label == "increased or absent threshold for auditory brainstem response" ~ "auditory brainstem response threshold",
                                        grepl(pattern = "first instar larval cuticle phenotype",object_label) ~ "larval cuticle phenotype",
                                        TRUE ~ object_label),
               object_label = object_label %>% str_to_sentence()) %>%
        arrange(desc(count_overall)
                ,desc(matrisome)
                ) %>%
        head(top_n_phenotypes) %>%
        mutate(object_label = factor(object_label,levels = unique(object_label))) %>%
        gather(key = "Matrisome membership",value = "value",-object_label,-count_overall,-species,-matrisome_enrichment,-count_expected_by_background,-n_overall) %>%
        mutate(`Matrisome membership` = factor(`Matrisome membership`,levels = c("other_genes","matrisome"))) %>%
        arrange(object_label) %>%
        ggplot(aes(x = object_label, y = value, fill = `Matrisome membership`)) +
        geom_bar(stat = "identity", color = "black", size = 0.05) +
        scale_fill_manual(values = c("matrisome" = species_color[[unique(.x[["species"]])]],"other_genes" = "#e8e8e8")) +
        labs(title = paste(unique(.x[["species"]])), y = "# Unique genes involved", x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(face = "italic",size = 12, angle = 50, hjust = 1), legend.position = "none",plot.title = element_text(face = "bold.italic", size = 16)) +
  facet_zoom(ylim = c(0, 150),zoom.size = 4))
pheno_occ_ecm_fraction_all_plot_names <- pheno_occ_ecm_contribution %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Matrisome_involvement/",.,"_all_most_contribution","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(pheno_occ_ecm_fraction_all_plot,pheno_occ_ecm_fraction_all_plot_names, ~ggsave(.x,filename = .y,height = 6, width = 10))
# figure: 
ggarrange(plotlist = pheno_occ_ecm_fraction_all_plot[c(3,3,2,4,5,1)],labels = "AUTO",font.label = list(size = 20),
          ncol = 2,nrow = 3,align = "hv") %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/pheno_occ_ecm_fraction_all_plot/","Grid_of_phenotype_perspective","_grouping_",apply_grouping,".pdf",collapse = ""),width = 17, 
    height = 17, 
    units = "mm"
  )

degree_matrisome_involvement <- pheno_occ_ecm_contribution %>%
  map(~.x %>% summarise(sum_count_overall = sum(count_overall), 
                        sum_matrisome = sum(matrisome), 
                        sum_other_genes = sum(other_genes),
                        matrisome_percentage = sum_matrisome/sum_count_overall * 100,
                        sum_mat_othergenes = sum_matrisome + sum_other_genes
                        ))



##################################################################
####### Overview table ##########
##################################################################
tab_monarch_overview <- gene_pheno %>% 
  group_by(subject_taxon_label) %>%
  summarise(Taxon = paste(unique(taxon),collapse = ", "),
            `# observations (total)` = n(),
            `# unique genes covered` = n_distinct(geneid),
            `# unique phenotypes reported` = n_distinct(object_label),
            `# combined sources (publications and others)` = n_distinct(source),
            `Modes of information assembly` = paste(unique(qualifier),collapse = ", ")
            ) %>%
  arrange(desc(`# observations (total)`))
write_csv(x = tab_monarch_overview,path = paste0(dir_output,"/Overview/tab_monarch_overview","_grouping_",apply_grouping,".csv",collapse = ""))




##################################################################
########## Gene & phenotype overview  #############
#########################################################

# - See which genes are involved in the most phenotypes overall --> histogram

pheno_by_taxa <- gene_pheno %>% 
  group_by(subject_taxon_label) %>% 
  summarise(distinct_phenotypes = n_distinct(object_label),
            distinct_genes = n_distinct(geneid)) %>% 
  arrange(desc(distinct_genes)) 

pheno_by_taxa_genes <- pheno_by_taxa %>% 
  arrange(desc(distinct_genes)) %>%
  mutate(subject_taxon_label = factor(subject_taxon_label,levels = subject_taxon_label)) %>%
  ggplot(aes(x = subject_taxon_label, y = distinct_genes)) +
  scale_y_log10() +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = "Gene perspective",subtitle = "Number of unique genes for which phenotypes were described for each species",y = "# of genes", x = "Species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,face = "italic"))
ggsave(pheno_by_taxa_genes, filename = paste0(dir_output,"/Overview/pheno_by_taxa_genes","_grouping_",apply_grouping,".pdf"),width = 7,height = 7)
ggsave(pheno_by_taxa_genes, filename = paste0(dir_output,"/Paper_figures/Fig_supp_n1/pheno_by_taxa_genes","_grouping_",apply_grouping,".pdf"),width = 7,height = 7)

pheno_by_taxa_pheno <- pheno_by_taxa %>% 
  arrange(desc(distinct_phenotypes)) %>%  
  mutate(subject_taxon_label = factor(subject_taxon_label,levels = subject_taxon_label)) %>%
  ggplot(aes(x = subject_taxon_label, y = distinct_phenotypes)) +
  scale_y_log10() +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = "Phenotype perspective",subtitle = "Number of unique phenotypes which were observed for each species",y = "# of phenotypes", x = "Species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,face = "italic"))
ggsave(pheno_by_taxa_pheno, filename = paste0(dir_output,"/Overview/pheno_by_taxa_pheno","_grouping_",apply_grouping,".pdf"),width = 7,height = 7)
ggsave(pheno_by_taxa_pheno, filename = paste0(dir_output,"/Paper_figures/Fig_supp_n1/pheno_by_taxa_pheno","_grouping_",apply_grouping,".pdf"),width = 7,height = 7)


cor_pheno_and_genes <- pheno_by_taxa %>% 
  ggplot(aes(x = distinct_genes, y = distinct_phenotypes)) +
  ggrepel::geom_label_repel(data = pheno_by_taxa[1:10,],aes(label = subject_taxon_label),fontface = "italic" , alpha = 0.5) +
  geom_point() +
  # coord_equal() +
  theme_classic() +
  labs(title = "Combined gene and phenotype perspectives",subtitle = "Relationship between the number of studied genes and phenotypes",y = "# of phenotypes", x = "# of genes")
ggsave(cor_pheno_and_genes, filename = paste0(dir_output,"/Overview/cor_pheno_and_genes","_grouping_",apply_grouping,".pdf"),width = 7,height = 7)
ggsave(cor_pheno_and_genes, filename = paste0(dir_output,"/Paper_figures/pheno_occ_ecm_fraction_all_plot/cor_pheno_and_genes","_grouping_",apply_grouping,".pdf"),width = 7,height = 5)



##################################################################
########## Gene perspective #############
##################################################################
# See which genes are involved in the most phenotypes overall --> histogram
top_n_genes <- 40
matrisome_fraction_genes <- ecm_pheno %>%
  map(~.x %>% summarize(species = unique(subject_taxon_label),
                        pheno_total = n_distinct(object_label),
                        pheno_matrisome_and_phenome = n_distinct(object_label[observed_in == "both"]),
                        pheno_only_matrisome = n_distinct(object_label[observed_in == "matrisome_only"]),
                        pheno_only_phenome = n_distinct(object_label[observed_in == "phenome_only"]),
                        matrisome_backround = pheno_matrisome_and_phenome / (pheno_matrisome_and_phenome + pheno_total)
  ))

ecm_genes_occ <- ecm_pheno %>%
  map(~.x %>% 
        group_by(gene_symbol) %>% 
        summarize(count_overall = n_distinct(object_label),
                  division = unique(division),
                  category = unique(category),
                  matrisome = sum(observed_in == "both"),
                  species = subject_taxon_label %>% unique(),
                  other_genes = count_overall - matrisome,
                  matrisome_fraction = matrisome / count_overall
        ) %>%
        arrange(desc(count_overall)) %>%
        filter(!is.na(gene_symbol), count_overall >= 1),
      ungroup())

genes_occ_ecm_contribution <- map2(ecm_genes_occ,matrisome_fraction_genes, ~.x %>% 
                                     mutate(count_expected_by_background = count_overall * .y[["matrisome_backround"]],
                                            matrisome_enrichment = matrisome / count_expected_by_background) %>%
                                     arrange(desc(count_overall),gene_symbol))


# plot
genes_occ_ecm_divis_plot <- genes_occ_ecm_contribution %>% 
  map(~.x  %>%
        arrange(desc(count_overall)) %>%
        head(top_n_genes) %>% 
        mutate(gene_symbol = factor(gene_symbol, levels = unique(gene_symbol))) %>% 
        ggplot(aes(x = gene_symbol, y = count_overall, fill = division)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d(begin = 0.3, end=0.8,direction = 1) +
        labs(title = paste("Matrisome genes which are most involved phenotypes in", unique(.x[["species"]]), ""),subtitle = "Phenotype mining using Monarch", y = "# distinct phenotypes this gene is involved in", x = "Gene") +
        theme_classic() +
        theme(axis.text.x = element_text(face = "italic",size = 13, angle = 60, hjust = 1)))
genes_occ_ecm_divis_plot_names <- genes_occ_ecm_contribution %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Gene_perspective/",.,"_ecm_divis_contribution","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(genes_occ_ecm_divis_plot,genes_occ_ecm_divis_plot_names, ~ggsave(.x,filename = .y,height = 8, width = 13))

genes_occ_ecm_cat_plot <- genes_occ_ecm_contribution %>% 
  map(~.x  %>%
        arrange(desc(count_overall)) %>%
        head(top_n_genes) %>% 
        mutate(gene_symbol = factor(gene_symbol, levels = unique(gene_symbol))) %>% 
        ggplot(aes(x = gene_symbol, y = count_overall, fill = category)) +
        geom_bar(stat = "identity") +
        scale_fill_viridis_d(begin = 0, end=1,direction = 1) +
        labs(title = paste(unique(.x[["species"]])), y = "", x = "Gene") +
        theme_classic() +
        theme(axis.title.x = element_blank(),axis.text.x = element_text(face = "italic",size = 18, angle = 90, hjust = 1, vjust = 0.5),axis.text.y = element_text(size = 16),plot.title = element_text(face = "bold.italic", size = 20)))
genes_occ_ecm_cat_plot_names <- genes_occ_ecm_contribution %>% map(~.x %>% pull(species) %>% unique() %>% str_replace_all(pattern = "\\.",replacement = "") %>% paste(dir_output,"/Gene_perspective/",.,"_ecm_cat_contribution","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(genes_occ_ecm_cat_plot,genes_occ_ecm_cat_plot_names, ~ggsave(.x,filename = .y,height = 8, width = 13))
# export figure
ggarrange(plotlist = genes_occ_ecm_cat_plot[c(3,2,4,5,1)],labels = "AUTO",font.label = list(size = 20),
          ncol = 2,nrow = 3,align = "hv",common.legend = TRUE,legend = "right") %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/gene_perspective_ecm_cat/","Grid_of_gene_perspective_ecm_cat","_grouping_",apply_grouping,".pdf",collapse = ""),width = 22, 
    height = 22, 
    units = "mm"
  )


##################################################################
########## Homology to taxa of interest ############# --> make this modular for every taxa, the whole section
##################################################################
# Direct Homolog Step 2#
# Here we use the provided homology information in the matrisomes to map all species matrisome genes to the corresponding human genes
gene_pheno_matrisomespecies <- gene_pheno %>% filter(subject_taxon_label %in% c("Drosophila melanogaster","Mus musculus","Caenorhabditis elegans","Homo sapiens","Danio rerio"))
homology_for_human <- left_join(gene_pheno_matrisomespecies,ecm_homolog_conversion_table,"geneid") %>% mutate(origin_species_gene_symbol = gene_symbol)
# Homolog is equal to gene_symbol (human). origin_species_gene_symbol refers to the gene_symbol the gene had in its original species
# Direct Homolog Step 2#

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

# Indirect inferred Homolog, Step 3#
# querry for diopt orthologs the remaining species that have no matrisome attached and keep the highest confidence one (e.g. rank "high" and "moderate")
diopt_taxa <- list("10116", #rat,
                   "4932" # Saccharomyces cerevisiae
                   )
if(remap_diopt){
  diopt_genes <- diopt_taxa %>% map(~gene_pheno %>% filter(taxon %in% .x) %>% filter(!is.na(subject_label)) %>% pull(subject_label) %>% unique())
  diopt_translated <- map2(.x = diopt_taxa,.y = diopt_genes,.f = ~diopt(.y,inTax = as.numeric(.x), outTax =9606))
  diopt_translation_table_export <- map2_df(.x = diopt_translated,.y = diopt_taxa,.f = ~.x %>% as_tibble() %>% mutate(taxon = .y,`Human GeneID` = `Human GeneID` %>% as.character()) %>% 
                                              select(taxon,`Search Term`,`Human Symbol`,`Human GeneID`,`DIOPT Score`,`Weighted Score`, Rank,`Best Score`,`Prediction Derived From`)
  )
  write_csv(path = paste0("./data/diopt_orthology/","diopt_translation_table_export_grouping_",apply_grouping,".csv"),x = diopt_translation_table_export)
}
diopt_translation_table <- read_csv(file = paste0("./data/diopt_orthology/","diopt_translation_table_export_grouping_",apply_grouping,".csv"))
diopt_translation_filt <- diopt_translation_table %>%
  filter(Rank %in% c("high","moderate")) %>%
  mutate(homolog = `Human Symbol`,origin_species_gene_symbol = `Search Term`, gene_symbol = `Human Symbol`, diopt_merge_id = paste(taxon,origin_species_gene_symbol)) %>%
  select(gene_symbol,homolog,origin_species_gene_symbol,diopt_merge_id)

gene_pheno_diopt_taxa <- gene_pheno %>% filter(taxon %in% diopt_taxa) %>% mutate(diopt_merge_id = paste(taxon,subject_label))
homology_diopt <- left_join(gene_pheno_diopt_taxa,diopt_translation_filt,"diopt_merge_id") 
# Indirect inferred Homolog, Step 3#

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

# Indirect inferred Homolog, Step 4#
already_covered_taxons <- c("10090","7227","9606","7955","6239",   "10116","4932")
blacklisted_taxons <- c("10092") #Mus musculus domesticus
homology_for_remaining <- gene_pheno %>%
  filter(!taxon %in% c(already_covered_taxons,blacklisted_taxons)) %>%
  group_by(taxon) %>%
  do(wrangle_homolog(.,target_taxa = "9606")) %>%
  ungroup() %>%
  mutate(homolog = as.character(homolog),
         homolog = case_when(taxon == 9606 ~ subject_label, TRUE ~ homolog),
         gene_symbol = homolog)


# Indirect inferred Homolog, Step 4#

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

# Assemble homology, Step 5#
combined_homology_human <- bind_rows(homology_for_human,homology_for_remaining,homology_diopt)
df_hom_for_human <- left_join(x = hs_mat %>% mutate(homolog = gene_symbol) %>% 
                                select(c(-geneid,-gene_symbol,-combined_orthology)),y = combined_homology_human %>% select(-gene_symbol), by = "homolog") %>%
  filter(division != "Retired") %>%
  mutate(gene_symbol = homolog)

homology_ecm_huborder <- df_hom_for_human %>% wrangle_homology_ecm_huborder()
hub_genes <- homology_ecm_huborder %>% filter(in_n_species >= 3) %>% arrange(desc(in_n_species),gene_symbol) %>% pull(gene_symbol)
# Assemble homology, Step 5#

df_hom_for_human %>% 
  filter(!is.na(subject_taxon_label)) %>%
  # filter(subject_taxon_label %in% c("Bos taurus", "Canis lupus familiaris","Rattus norvegicus","Saccharomyces cerevisiae S288C")) %>%
  rename(`Gene identifier in species of origin` = geneid, `Predicted human ortholog` = homolog,`Gene name in species of origin` = origin_species_gene_symbol,
         `Extended definition URL` = is_defined_by,Qualifier = qualifier,`Gene name` = subject_label, `Taxon identifier` = taxon, `Taxon name` = subject_taxon_label, 
         `Phenotype identifier` = object, `Simplified phenotype group` = object_label, `Association evidence label` = evidence_label,`Association evidence identifier` = evidence, `Association source (if directly provided)` = source,`Phenotype name` = phenotype_as_imported) %>%
  select(`Taxon name`, `Taxon identifier`,division,category, `Predicted human ortholog`,`Gene name in species of origin`,`Gene identifier in species of origin`,`Phenotype name`, `Phenotype identifier`, `Association source (if directly provided)`, `Association evidence label`,`Association evidence identifier`, `Simplified phenotype group`,`Extended definition URL`,Qualifier) %>%
  mutate(`Taxon identifier` = ifelse(`Taxon identifier` == "4932","559292 (4932 used for mapping)",`Taxon identifier`)) %>% 
  arrange(`Taxon name`,`Predicted human ortholog`,`Phenotype name`) %>%
  write_csv(x = .,path = paste0(dir_output,"/Paper_figures/Supplement_tables/Sx3.csv"),na = "")
  
##################################################### 
#### pipeline for all phenotypes mapped to human #### 
##################################################### 
homology_ecm <- df_hom_for_human %>% filter(!is.na(division),!is.na(gene_symbol),!is.na(object_label)) 
tax <- 9606
hub <- hub_genes

# single-gene plots
homology_ecm_linebreak <- homology_ecm %>%
  rowwise() %>%
  mutate(object_label = wrangle_insert_linebreaks(text = object_label,max_char = 48)) %>%
  ungroup()

net <- hub %>% 
  map(~.x %>% wrangle_generate_hub_connections(show_max_n_additional_per_species = 10,homology_ecm = homology_ecm_linebreak,species_abbreviations = species_abbreviations) %>% .[["g_unfolded"]])
dend_lin <- net %>%
  map(~.x %>% ggraph(layout = 'dendrogram') +
        geom_edge_diagonal() +
        geom_node_text(aes( label=name, filter=leaf, color = origin) , angle=90 , hjust=1, nudge_y = -0.04, size = 2.5) +
        geom_node_point(aes(filter=leaf) , alpha=0.6) +
        ylim(-.5, NA) +
        theme_void())
dend_lin_names <- net %>% map(~.x %>% activate(nodes) %>% filter(origin == "hub") %>% pull(name) %>% str_replace_all(pattern = "[[:punct:]]",replacement = "") %>% paste(dir_output,"/Homolog_taxon",tax,"/linear_dendrogram_",.,"","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
# walk2(dend_lin,dend_lin_names, ~ggsave(.x,filename = .y,height = 10, width = 20))
dend_circ <- net %>%
  map(~ draw_circular_dendrogram(.x,species_color = species_color,species_size = 3,leaf_size = 3,hub_size = 8,point_size_range = c(1.2,3.5)))
dend_circ_names <- net %>% map(~.x %>% activate(nodes) %>% filter(origin == "hub") %>% pull(name) %>% str_replace_all(pattern = "[[:punct:]]",replacement = "") %>% paste(dir_output,"/Homolog_taxon",tax,"/circular_dendrogram_",.,"","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
# walk2(dend_circ,dend_circ_names, ~ggsave(.x,filename = .y,height = 10, width = 10))
dend_circ[[1]]
ggarrange(plotlist = dend_circ,labels = "AUTO",font.label = list(size = 50),
          ncol = 2,nrow = 3,align = "hv") %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/dend_circ_supplement_grid/","Grid_of_dend_circ","_grouping_",apply_grouping,".pdf",collapse = ""),width = 25, # width = 25, height = 35
    height = 35, 
    units = "mm"
  )

# Figure 4
homology_ecm_noyeast <- homology_ecm %>% filter(subject_taxon_label != "Saccharomyces cerevisiae")
end_circ_fig4 <- c("CTSD") %>% #c("MSTN","CTSD","LAMB2","HSPG2","COL11A2")
  map(~.x %>% wrangle_generate_hub_connections(show_max_n_additional_per_species = 10,homology_ecm = homology_ecm_noyeast,species_abbreviations = species_abbreviations) %>% .[["g_unfolded"]]) %>%
  map(~ draw_circular_dendrogram(.x,species_color = species_color,species_size = 5,leaf_size = 4,hub_size = 8,point_size_range = c(1.2,3.5)) + theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+ theme(plot.margin = unit(c(0,0,0,0), "lines")))
dend_circ_fig4_names <- c("CTSD") %>% map(~.x  %>% str_replace_all(pattern = "[[:punct:]]",replacement = "") %>% paste(dir_output,"/Paper_figures/Fig4/hub_",.,"_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
walk2(end_circ_fig4,dend_circ_fig4_names, ~ggsave(.x,filename = .y,height = 10, width = 10))

networkname <- "Matrisome phenotype across species"
plot <- homology_ecm %>% 
  # filter(!subject_taxon_label %in% c("Homo sapiens", "Mus musculus")) %>%
  mutate(object_label = case_when(object_label %in% c("Brain phenotype","Complex brain phenotype") ~ "Brain \n phenotype", 
                                  object_label %in% c("Cardiovascular impairments","Cardiovascular phenotype","Heart function impairment") ~ "Cardio- \n vasculature", 
                                  object_label %in% c("Muscle weakness phenotype","Muscle tissue phenotype") ~ "Muscle \n phenotype", 
                                  # object_label %in% c("Diabetes related phenotype","Glucose metabolism related") ~ "Glucose \n metabolism", 
                                  object_label %in% c("Inflammation phenotype","Autoimmune diseases","Immune system phenotype","Immunodeficiency","Immunologic hypersensitivity") ~ "Immune \n system", 
                                  object_label %in% c("Cancer","Aberrant proliferation") ~ "Proliferation \n phenotype", 
                                  object_label %in% "Skin phenotype" ~ "Skin \n phenotype",
                                  TRUE ~ object_label)) %>%
  filter(grepl(pattern = "\n",x = object_label)) %>%
     gene_pheno_net_subset(df = .,n_phenotypes = 10,
                        genes = c("MSTN","HSPG2","LAMA5","LAMA2","BMP2","WNT1","GPC4","SHH","MMP9","AGRN","COL7A1")
                        # phenotypes = c("Brain \n phenotype","Cardio- \n vasculature","Muscle \n phenotype","Immune \n system","Proliferation \n phenotype","Skin \n phenotype") 
                        ) %>%
  wrangle_generate_any_network(species_order = species_order) %>%
  gene_pheno_net_plot(title = networkname,node_color_fill = division,species_color = species_color,nrow = 3,edge_widths = c(0.15,0.4,0.7,1.1,1.5,2),edge_alphas = c(1,1,0.7,0.7,0.7,0.7),size_label = 4.5,size_text = 4.0)
export   <- plot + theme(legend.position = "None", plot.subtitle = element_blank(), strip.text.x = element_text(size = 16)) + theme(panel.spacing.x=unit(1, "lines"))
export %>% ggsave(filename = "./output/Paper_figures/Fig6/disease_relevant_network.pdf",plot = ., width = 12, height = 16)

networkname <- "Matrisome phenotype across species: body size"
plot <- homology_ecm %>% 
  # filter(category %in% c("ECM Regulators")) %>%
  gene_pheno_net_subset(df = .,
                        n_genes = 20,
                        phenotypes = c("Altered body size","Bone phenotype","Connective tissue phenotype", "Cranial deformations","Facial deformations","Short stature","Spinal chord phenotype"),
                        n_phenotypes = 5
  ) %>%
  wrangle_generate_any_network(species_order = species_order) %>%
  gene_pheno_net_plot(title = networkname,node_color_fill = division,species_color = species_color)
network_name <- paste(dir_output,"/Homolog_taxon",tax,"/network_",networkname,"","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " ")
ggsave(filename = network_name,plot = plot, width = 12, height = 7)


# Figure 5 is assembled here: each category
categories <- homology_ecm %>% filter(category!="n/a") %>% pull(category) %>% unique()
networkname <-  categories %>% map(~.x %>% paste0("",.))

homology_ecm_toplot <- homology_ecm %>% 
  mutate(object_label = case_when(object_label %in% c("Blood phenotype") ~ "Blood \n phenotype", 
                                  object_label %in% c("Cardiovascular impairments") ~ "Cardiovascular \n impairments", 
                                  object_label %in% c("Altered body size") ~ "Altered \n body size", 
                                  object_label %in% c("Sterility / reproductive potential") ~ "Sterility/ \n reproductive potential", 
                                  object_label %in% c("Morphology altered") ~ "Morphology \n altered", 
                                  object_label %in% c("Muscle tissue phenotype") ~ "Muscle tissue \n phenotype", 
                                  object_label %in% c("Eye phenotype") ~ "Eye \n phenotype", 
                                  object_label %in% c("Developmental phenotype") ~ "Developmental \n phenotype", 
                                  object_label %in% c("Connective tissue phenotype") ~ "Connective tissue \n phenotype", 
                                  object_label %in% c("Bone phenotype") ~ "Bone \n phenotype", 
                                  TRUE ~ object_label))

plot <- map2(.x = categories,.y = networkname,{~
    homology_ecm_toplot %>% 
    mutate(subject_taxon_label = factor(subject_taxon_label,levels = species_order)) %>%
    filter(category %in% .x) %>%
    filter(subject_taxon_label %in% c("Homo sapiens","Mus musculus","Danio rerio")) %>% # Filter for human and mouse
    gene_pheno_net_subset(n_genes = 5,
                          n_phenotypes = 5
    ) %>%
    wrangle_generate_any_network(species_order = species_order) %>%
    gene_pheno_net_plot(title = .y,node_color_fill = division,species_color = species_color,margin = 15,border = FALSE,size_label = 4.5,size_text = 4.0)
    })
plot <- plot %>% map(~.x + theme(panel.spacing.x=unit(1.5, "lines"))) # reduce space between facets
keep <- !map_lgl(plot,is.null)
network_name <-  networkname %>% map(~.x %>%  paste(dir_output,"/Homolog_taxon",tax,"/network_",.,"","_grouping_",apply_grouping,".pdf",collapse = "") %>% str_remove_all(pattern = " "))
map2(.x = network_name[keep],.y = plot[keep], ~ ggsave(filename = .x,plot = .y, width = 10, height = 6))
ggarrange(plotlist = plot,
          ncol = 1,nrow = 6,align = "hv",common.legend = TRUE) %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/Fig5/","all_combined_categories","_grouping_",apply_grouping,".pdf",collapse = ""),width = 12, 
    height = 16, 
    units = "mm"
  )
ggarrange(plotlist = plot[c(1,2,3)],
          ncol = 1,nrow = 3,align = "hv",common.legend = TRUE) %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/Fig5/","structural_combined_categories","_grouping_",apply_grouping,".pdf",collapse = ""),width = 12, 
    height = 16, 
    units = "mm"
  )



# All categories for all species --> supp. fig 8
categories <- homology_ecm %>% filter(category!="n/a") %>% pull(category) %>% unique()
networkname <-  categories %>% map(~.x %>% paste0("",.))
plot <- map2(.x = categories,.y = networkname,{~
    homology_ecm %>% 
    filter(object_label != "Uncharacterized phenotype", subject_taxon_label != "Saccharomyces cerevisiae") %>%
    mutate(subject_taxon_label = factor(subject_taxon_label,levels = species_order)) %>%
    filter(category %in% .x) %>%
    gene_pheno_net_subset(n_genes = 5,
                          n_phenotypes = 5
    ) %>%
    wrangle_generate_any_network(species_order = species_order) %>%
    gene_pheno_net_plot(title = .y,node_color_fill = division,species_color = species_color,margin = 15,border = TRUE,nrow = 1) 
})
plot <- plot %>% map(~.x + theme(panel.spacing.x=unit(1, "lines")))
# plot <- plot %>% map(~.x + theme(legend.position = "none"))
keep <- !map_lgl(plot,is.null)
network_name <-  networkname %>% map(~.x %>%  paste(dir_output,"/Homolog_taxon",tax,"/network_",.,"","_grouping_",apply_grouping,"Supp_n8.pdf",collapse = "") %>% str_remove_all(pattern = " "))
map2(.x = network_name[keep],.y = plot[keep], ~ ggsave(filename = .x,plot = .y, width = 16, height = 16))
ggarrange(plotlist = plot[c(1,5,4,6,2,3)], # 
          ncol = 1,common.legend = TRUE) %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/Supp_n8/","all_combined_categories","_grouping_",apply_grouping,".pdf",collapse = ""),width = 18, 
    height = 26, 
    units = "mm"
  )

# All categories for all species --> supp. fig 9
divisions <- c("Core matrisome","Matrisome-associated")
networkname <-  divisions %>% map(~.x %>% paste0("",.))
plot <- map2(.x = divisions,.y = networkname,{~
    homology_ecm %>%
    filter(object_label != "Uncharacterized phenotype", subject_taxon_label != "Saccharomyces cerevisiae") %>%
    mutate(subject_taxon_label = factor(subject_taxon_label,levels = species_order)) %>%
    filter(division %in% .x) %>%
    gene_pheno_net_subset(n_genes = 7,
                          n_phenotypes = 7
    ) %>%
    wrangle_generate_any_network(species_order = species_order) %>%
    gene_pheno_net_plot(title = .y,node_color_fill = division,species_color = species_color,margin = 15,border = TRUE,nrow = 2) 
})
# plot <- plot %>% map(~.x + theme(legend.position = "none"))
keep <- !map_lgl(plot,is.null)
network_name <-  networkname %>% map(~.x %>%  paste(dir_output,"/Homolog_taxon",tax,"/network_",.,"","_grouping_",apply_grouping,"Supp_n8.pdf",collapse = "") %>% str_remove_all(pattern = " "))
map2(.x = network_name[keep],.y = plot[keep], ~ ggsave(filename = .x,plot = .y, width = 16, height = 16))
ggarrange(plotlist = plot, # 
          ncol = 1,common.legend = TRUE) %>%
  ggexport(
    filename = paste0(dir_output,"/Paper_figures/Supp_n9/","all_combined_divisions","_grouping_",apply_grouping,".pdf",collapse = ""),width = 17, 
    height = 24, 
    units = "mm"
  )


# highlight a specific catgory:
plot <- homology_ecm %>%  
  filter(category %in% "ECM Glycoproteins") %>%
    mutate(subject_taxon_label = factor(subject_taxon_label,levels = species_order)) %>%
    gene_pheno_net_subset(genes = c("LAMB1","LAMBA1","LAMC1","FBN1","FN1", "LAMA5","CRIM1","LAMB2","BMPER","EFEMP2"),phenotypes = c("Bone phenotype","Eye phenotype","Sterility / Reduced reproduction potential","Developmental phenotype","Cardiovascular impairments","Altered body size")
    ) %>%
  mutate(object_label = case_when(object_label == "Bone phenotype" ~ "Bone \n phenotype", 
                                  object_label == "Sterility / Reduced reproduction potential" ~ "Sterility / Reduced \n reproduction potential", 
                                  object_label == "Eye phenotype" ~ "Eye \n phenotype", 
                                  object_label == "Developmental phenotype" ~ "Developmental \n phenotype",
                                  object_label ==  "Cardiovascular impairments" ~ "Cardiovascular \n impairments",
                                  object_label ==  "Altered body size" ~ "Altered \n body size",
                                  TRUE ~ object_label)) %>%
    wrangle_generate_any_network(species_order = species_order) %>%
    gene_pheno_net_plot(title = "ECM Glycoproteins highlight",node_color_fill = division,species_color = species_color)
export   <- plot + theme(legend.position = "None", plot.subtitle = element_blank()) 
export %>% ggsave(filename = "./output/Paper_figures/Fig5/Highlight_ECM_Glycoproteins.pdf",plot = ., width = 15, height = 12)

############################## Export tables ############################## 
exp <- homology_ecm %>% 
  filter(!is.na(subject_taxon_label) & !is.na(object_label)) %>%
  select(division,category,gene_symbol,subject_taxon_label,object_label) %>% 
  mutate(subject_taxon_label = factor(subject_taxon_label,levels = species_order)) %>%
  group_by(gene_symbol, subject_taxon_label,division, category) %>% 
  arrange(object_label) %>%
  summarise(phenotypes = paste0(object_label,collapse = "\r\n ")) %>%
  ungroup() %>%
  spread(key = subject_taxon_label,value = phenotypes) %>%
  arrange(division,category)

col_width <- 60
sheet_name <- paste0("Taxon ", tax)
wb <- openxlsx::createWorkbook() 
openxlsx::addWorksheet(wb, sheet_name)
openxlsx::addStyle(wb, sheet = sheet_name, 
                   createStyle(fontSize = 14, fontColour = "black", halign = "left",wrapText = TRUE), 
                   rows = 1:nrow(exp), cols = 1:ncol(exp), gridExpand = TRUE)
writeData(wb, sheet_name, exp)
setColWidths(wb, sheet_name , cols = 1:ncol(exp), widths = c(rep(col_width,ncol(exp)-1),100))
saveWorkbook(wb, paste0(dir_output,"/Homolog_taxon",tax,"/gene_phenotype_table_taxon",tax,"grouping_",apply_grouping,".xlsx",collapse = ""),overwrite = TRUE)

# All phenotypes:
phenotype_table <- gene_pheno %>% 
  group_by(object_label) %>%
  mutate(`Involved in species` = paste(subject_taxon_label %>% unique(),collapse = ", "), n_species = n_distinct(subject_taxon_label)) %>%
  select(n_species,object_label,`Involved in species`) %>%
  filter(!duplicated(object_label)) %>%
  arrange(desc(n_species))
write_csv(phenotype_table,paste0(dir_output,"/phenotype_grouping/phenotypes_for_grouping.csv",collapse = ""))
##########################################################################################################
########################################## End of core analysis ########################################## 
##########################################################################################################






if(apply_grouping){
  # ##################################################################
  # ####### String analysis ##########
  # ##################################################################
  # #### Run once to establish network
  if(reimport_string){
    string_db <- STRINGdb$new(version="10", species=9606, score_threshold=400, input_directory="" )
    graph_full <- string_db$get_graph()
    graph_full_backup <- graph_full
    if(resave_string) write_rds(x = graph_full_backup,path = "./data/string_graph_full_backup.rds")
    homology_ecm_string <- string_db$map( homology_ecm %>% as.data.frame(), "gene_symbol", removeUnmappedRows = TRUE ) %>% as_tibble() %>% filter(STRING_id %in% names(igraph::degree(graph_full)))
    graph_fullecm <- induced_subgraph(graph_full, homology_ecm_string %>% pull(STRING_id))
    graph_core <- induced_subgraph(graph_full, homology_ecm_string %>% filter(division == "Core matrisome") %>% pull(STRING_id))
    graph_assoc <- induced_subgraph(graph_full, homology_ecm_string %>% filter(division == "Matrisome-associated") %>% pull(STRING_id))
    phenotypes_to_validate <- homology_ecm %>% group_by(object_label) %>% summarise(n_distinct_genes = n_distinct(gene_symbol), n_distinct_original_phenotypes = n_distinct(phenotype_as_imported)) %>% arrange(desc(n_distinct_genes),desc(n_distinct_original_phenotypes))
    # set vertex names
    v_name <- V(graph_fullecm)[[]] %>% enframe() %>% select(name) %>% rename(STRING_id = name)
    homology_ecm_string_uniq <- homology_ecm_string %>% filter(!duplicated(STRING_id))
    v_attr <- left_join(x = v_name, y = homology_ecm_string_uniq,by = "STRING_id")
    graph_fullecm <- set.vertex.attribute(graph_fullecm, "name", value=v_attr$gene_symbol)
    fullecm <- graph_fullecm %>% as_tbl_graph()
    if(resave_string) write_rds(x = fullecm,path = "./data/fullecm.rds")
  }
  # Loop over phenotypes
  n_times_resample <- 100
  universe <- homology_ecm %>% pull(`Gene Symbol`) %>% unique()
  phenotypes_to_study <- homology_ecm %>% filter(!is.na(object_label),!duplicated(object_label)) %>% pull(object_label)
  network_stats <- phenotypes_to_study %>% map_df(~.x %>% wrapper_network_statistics(analyze_phenotype = .,df_gene_pheno = homology_ecm,n_times_resample = n_times_resample,universe = universe))
  network_stats[is.na(network_stats)] <- NA
  # mark grouped phenotypes
  grouped_pheno_names <- homology_ecm %>% select(object_label, phenotype_as_imported) %>% mutate(object_label_lower = object_label %>% tolower(),
                                                                                                 phenotype_as_imported_lower = phenotype_as_imported %>% tolower(),
                                                                                                 manually_changed_phenotype = object_label_lower != phenotype_as_imported_lower) %>% 
    filter(manually_changed_phenotype) %>% pull(object_label) %>% unique(na.rm = TRUE)
  network_stats <- network_stats %>% filter(n_genes >= 10) %>% 
    mutate(`Grouped phenotype` = ifelse(phenotype %in% grouped_pheno_names,"grouped","original")) 
  # export table
  network_stats %>% 
    arrange(desc(n_genes),desc(`node connectivity relative to randomly sampled [%]`,phenotype,type)) %>%
    write_csv(x = .,path = paste0(dir_output,"/Paper_figures/Supplement_tables/Sx4.csv"),na = "")
  # Answer questions in manuscript in later section
  network_stats_string_analysis <- network_stats
  
  
  #####################################################################################
  ###### Investigate if matrisome categories are enriched for certain phenotypes ######
  #####################################################################################
  ##################################################################
  ### TABLE 5: Supplementary table: ENTIRE phenome
  n_genes_in_pheno_cutoff <- 5
  stats_species <- c("Homo sapiens","Mus musculus", "Danio rerio", "Drosophila melanogaster","Caenorhabditis elegans")
  stats_taxa_dataset <- list("Homo sapiens" = ecm_pheno_human,
                             "Mus musculus" = ecm_pheno_mouse, 
                             "Danio rerio" = ecm_pheno_fish, 
                             "Drosophila melanogaster" = ecm_pheno_fly,
                             "Caenorhabditis elegans" = ecm_pheno_worm) %>% 
    map(~ .x %>% filter(observed_in %in% c("both", "phenome_only")) %>%
          mutate(
            category = ifelse(is.na(category), "Non-matrisome", category),
            category_overall = ifelse(observed_in == "both", "Matrisome", "Non-matrisome")
          ))
  stats_pheno <- stats_taxa_dataset %>%
    map( ~{.x %>%
        filter(!is.na(object_label)) %>%
        group_by(object_label) %>%
        summarise(n_distinct_genes = n_distinct(subject_label, na.rm = TRUE)) %>%
        filter(n_distinct_genes >= n_genes_in_pheno_cutoff) %>%
        pull(object_label)
    })
  pmap_input <- tibble(pheno_of_interest = stats_pheno, taxon_of_interest = stats_species, df_gene_pheno = stats_taxa_dataset,ecm_and_nonecm_categories = list(c("Non-matrisome",ecm_categories_ordered)))
  stats_list <- pmap(.l = pmap_input,wrangle_overall_phenotype_enrichment) # Very slow
  stats_df <- stats_list  %>% bind_rows() %>% select(species, everything()) %>% group_by(species,phenotype) %>% 
    mutate(`Matrisome contribution` = wrangle_decorate_table_ecmcontribution(category = category,pval = pval)) %>% 
    ungroup() %>%
    mutate(`Matrisome contribution` = factor(`Matrisome contribution`,levels = c("No matrisome signature","Matrisome enriched (1 category)","Matrisome enriched (2 categories)","Matrisome enriched (3 categories)","Matrisome enriched (4 categories)","Matrisome enriched (5 categories)","Matrisome enriched (6 categories)")),
           species = factor(species,levels = species_order)) %>%
    arrange(species,desc(`Matrisome contribution`),phenotype) %>% rename(`n genes in category` = `n genes in matrisome category`)
  write_csv(x = stats_df,path = paste0(dir_output,"/Paper_figures/Supplement_tables/Sx5.csv"))
  ### TABLE 6: Supplementary table: HUMAN ORTHOLOG phenome
  n_genes_in_pheno_cutoff <- 5
  stats_species <- c("Homo sapiens","Mus musculus", "Danio rerio", "Drosophila melanogaster","Caenorhabditis elegans", "Rattus norvegicus")
  
  stats_pheno <- stats_species %>% map(~{
    homology_ecm %>%
      filter(subject_taxon_label == .x) %>%
      filter(!is.na(object_label)) %>%
      filter(!object_label %in% c("Uncharacterized phenotype")) %>%
      group_by(object_label) %>%
      summarise(n_distinct_genes = n_distinct(gene_symbol, na.rm = TRUE)) %>%
      filter(n_distinct_genes >= n_genes_in_pheno_cutoff) %>%
      pull(object_label)
  })
  stats_subdf<- map2(.x = stats_species,.y = stats_pheno,~{
    homology_ecm %>%
      filter(subject_taxon_label == .x) %>%
      filter(!is.na(object_label)) %>%
      filter(object_label %in% .y)
  })
  pmap_input <- tibble(pheno_of_interest = stats_pheno, taxon_of_interest = stats_species, df_gene_pheno = stats_subdf,ecm_and_nonecm_categories = list(c(ecm_categories_ordered)))
  stats_list <- pmap(.l = pmap_input,wrangle_overall_phenotype_enrichment) # Very slow
  stats_df <- stats_list  %>% bind_rows() %>% select(species, everything()) %>% group_by(species,phenotype) %>% 
    mutate(`Matrisome contribution` = wrangle_decorate_table_ecmcontribution(category = category,pval = pval)) %>% 
    ungroup() %>%
    mutate(`Matrisome contribution` = factor(`Matrisome contribution`,levels = c("No matrisome signature","Matrisome enriched (1 category)","Matrisome enriched (2 categories)","Matrisome enriched (3 categories)","Matrisome enriched (4 categories)","Matrisome enriched (5 categories)","Matrisome enriched (6 categories)")),
           species = factor(species,levels = species_order)) %>%
    arrange(species,desc(`Matrisome contribution`),phenotype) %>% rename(`n genes in category` = `n genes in matrisome category`)
  write_csv(x = stats_df,path = paste0(dir_output,"/Paper_figures/Supplement_tables/Sx6.csv"))
  ##################################################################
  ##### FIGURE Fig 5:  Human, Mouse, Fish - phenotypes shown in main figure
  path_fig5_supp <- paste0("./output/Paper_figures/ecmcat_enriched/Fig5_supplement/separated/")
  fig5_species <- c("Homo sapiens","Mus musculus","Danio rerio", "Drosophila melanogaster","Caenorhabditis elegans"
                    )
  # fig5_pheno <- list("Group1" = c("Morphology altered", "Altered body size", "Sterility / reproductive potential"),"Group2" = c("Blood phenotype", "Cardiovascular impairments" , "Sterility / reproductive potential"))
  fig5_pheno <- list(
    # "ECM Glycoproteins" = c("Altered body size", "Blood phenotype", "Cardiovascular impairments" ,"Morphology altered", "Sterility / reproductive potential"),
    #                  "Collagens" = c("Altered body size","Developmental phenotype", "Eye phenotype","Morphology altered","Muscle tissue phenotype"),
    #                  "Proteoglycans" = c("Blood phenotype","Bone phenotype","Connective tissue phenotype","Eye phenotype","Morphology altered"),
    "Figure 5 follow up:" = c("Altered body size","Bone phenotype","Cardiovascular impairments","Connective tissue phenotype","Muscle tissue phenotype", "Sterility / reproductive potential")
    # "Figure 6 follow up" = c("Skin phenotype", "Immune system phenotype", "Aging-related phenotype", "Aplasia/Hypoplasia of bone", "Aplasia/Hypoplasia of skin","atrophic scars")
  )
  map(fig5_species,~ map(fig5_pheno, ~{
    print(.x)
    print(.y)
    filename <- wrangle_file_path(path_fig5_supp,species = .y,phenotype = paste0(.x,collapse = "") %>%  str_remove_all(string = .,pattern = " ") )
    print(filename)
    wrangle_ecmcategory_enriched_phenotypes(df_gene_pheno = homology_ecm,
                                            taxon_of_interest = .y,
                                            venn_fontsize = 16,
                                            pheno_of_interest = .x,
                                            ecm_categories = ecm_categories_ordered) %>% 
      
      # plot_ecmcategory_enriched_phenotypes(output_list = ., path_output = filename,plot_hight = 22,plot_width = 19) # all on one page
      plot_ecmcategory_enriched_phenotypes_separate(output_list = ., path_output = filename,upset_height = 8,upset_width = 19) # all plots separate
    
  },.y=.x)) 
  ##################################################################
  ### FIGURE Fig 6:  Human, Mouse, Fish - phenotypes shown in main figure
  # Assemble larger meta group as in figure 6
  path_fig6_supp <- paste0("./output/Paper_figures/ecmcat_enriched/Fig6_supplement/")
  fig6_species <- c("Homo sapiens","Mus musculus", "Danio rerio", "Drosophila melanogaster","Caenorhabditis elegans")
  fig6_pheno <- list("Group2" = c("Brain phenotype","Muscle tissue phenotype","Immune system phenotype","Cardiovascular phenotype"),
                     "Group1" = c("Skin phenotype","Immune system phenotype","Cardiovascular phenotype"),
                     "Testgroup3" = c("Skin phenotype","Immune system phenotype","Cardiovascular phenotype","Brain phenotype","Muscle tissue phenotype","Immune system phenotype","Cardiovascular phenotype")
  )
  map(fig6_species,~ map(fig6_pheno, ~{
    print(.x)
    print(.y)
    filename <- wrangle_file_path(path_fig6_supp,species = .y,phenotype = paste0(.x,collapse = "") %>%  str_remove_all(string = .,pattern = " ") )
    print(filename)
    try({
      wrangle_ecmcategory_enriched_phenotypes(df_gene_pheno = homology_ecm,
                                              taxon_of_interest = .y,
                                              venn_fontsize = 12,
                                              pheno_of_interest = .x,
                                              ecm_categories = ecm_categories_ordered) %>% 
        plot_ecmcategory_enriched_phenotypes(output_list = ., path_output = filename)
    },silent = TRUE)},.y=.x)) 
  
  
  ######################################
  ########### Section on MSTN ##########
  ######################################
  df_hom_for_human %>% 
    filter(!is.na(subject_taxon_label)) %>%
    # filter(subject_taxon_label %in% c("Bos taurus", "Canis lupus familiaris","Rattus norvegicus","Saccharomyces cerevisiae S288C")) %>%
    rename(`Gene identifier in species of origin` = geneid, `Predicted human ortholog` = homolog,`Gene name in species of origin` = origin_species_gene_symbol,
           `Extended definition URL` = is_defined_by,Qualifier = qualifier,`Gene name` = subject_label, `Taxon identifier` = taxon, `Taxon name` = subject_taxon_label, 
           `Phenotype identifier` = object, `Simplified phenotype group` = object_label, `Association evidence label` = evidence_label,`Association evidence identifier` = evidence, `Association source (if directly provided)` = source,`Phenotype name` = phenotype_as_imported) %>%
    select(`Taxon name`, `Taxon identifier`,division,category, `Predicted human ortholog`,`Gene name in species of origin`,`Gene identifier in species of origin`,`Phenotype name`, `Phenotype identifier`, `Association source (if directly provided)`, `Association evidence label`,`Association evidence identifier`, `Simplified phenotype group`,`Extended definition URL`,Qualifier) %>%
    mutate(`Taxon identifier` = ifelse(`Taxon identifier` == "4932","559292 (4932 used for mapping)",`Taxon identifier`)) %>% 
    arrange(`Taxon name`,`Predicted human ortholog`,`Phenotype name`)  %>%
    filter(`Predicted human ortholog` == "MSTN") %>% 
    mutate(muscle_associated = ifelse(grepl("Muscl",`Simplified phenotype group`),1,0),
           adipose_associated = ifelse(grepl("Adip",`Simplified phenotype group`),1,0),
           heart_associated = ifelse(grepl("Cardi|Heart",`Simplified phenotype group`),1,0),
           liver_associated = ifelse(grepl("Liver",`Simplified phenotype group`),1,0),
           combined_association = muscle_associated  + adipose_associated + heart_associated + liver_associated) %>% 
    arrange(desc(combined_association),`Taxon name`,desc(muscle_associated),desc(adipose_associated)) %>%
    select(`Predicted human ortholog`,`Gene identifier in species of origin`,`Taxon name`,muscle_associated,adipose_associated,heart_associated,liver_associated,`Gene name in species of origin`, `Simplified phenotype group`, everything()) %>%
    select(c(-division,-category,-`Taxon identifier`)) %>% 
    write_csv(paste0(dir_output,"/Paper_figures/MSTN_info.csv"))
  
  homology_ecm %>%
    filter(!is.na(subject_taxon_label)) %>%
    filter(gene_symbol == "MSTN") %>%
    mutate(muscle_associated = ifelse(grepl("Muscl",object_label),1,0),
           adipose_associated = ifelse(grepl("Adip",object_label),1,0),
           heart_associated = ifelse(grepl("Cardi|Heart",object_label),1,0),
           liver_associated = ifelse(grepl("Liver",object_label),1,0),
           combined_association = muscle_associated  + adipose_associated + heart_associated + liver_associated) %>%
    arrange(desc(combined_association),subject_taxon_label,desc(muscle_associated),desc(adipose_associated)) %>%
    select(`Gene Symbol`,subject_taxon_label,muscle_associated,adipose_associated,heart_associated,liver_associated,subject_label, object_label, everything()) %>%
    select(c(-division:-Refseq_IDs,-prefix,-prefix_taxon,-taxon)) 
}


