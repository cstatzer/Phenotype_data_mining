wrangle_homolog <- function(df,target_taxa){
  # designed to be used with group by on taxon
  print(paste("Analyzing homologs of:",target_taxa))
  orignal_taxa <- unique(df %>% pull(taxon)) %>% .[]
  orignal_taxa_label <- unique(df %>% pull(subject_taxon_label)) %>% .[]
  print(paste("original taxa: ",orignal_taxa, " name: ",paste(orignal_taxa_label,collapse = ",")))
  df$homolog <- NA
  if(orignal_taxa == target_taxa) return(df)
  
  covered_genes <- df %>% pull(subject_label) %>% unique()
  trans_tab_raw <- homologene(covered_genes, inTax = orignal_taxa, outTax = target_taxa,db = homologene_db_20200327) # now using and updated database since 27. März
  colnames(trans_tab_raw)[colnames(trans_tab_raw) == orignal_taxa] <- "origin_species"
  colnames(trans_tab_raw)[colnames(trans_tab_raw) == target_taxa] <- "target_species"
  non_dup <- !duplicated(trans_tab_raw["origin_species"])
  trans_tab <- trans_tab_raw[non_dup,]
  # pb <- txtProgressBar(min = 1, max = nrow(df), style = 3)
  for (row in 1:nrow(df)) {
    old_gene <- df[row,"subject_label"] %>% unlist() %>% unname()
    translated <- trans_tab %>% filter(origin_species == old_gene) %>% pull(target_species) %>% .[]
    df[row,"homolog"] <- ifelse(length(translated) == 0,NA,translated)
    # setTxtProgressBar(pb, as.numeric(row))
  }
  df %>% mutate(origin_species_gene_symbol = subject_label)
}

wrangle_insert_linebreaks <- function(text, max_char = 20){
  str_len <- stringr::str_length(text)
  if(str_len <= max_char) return(text)
  chars <- unlist(strsplit(text, ""))
  str_whitespace_idx <- str_which(string = chars,pattern = " ")
  pieces <- ceiling(length(chars) / max_char)
  potential_breakpoints <- c(1:(pieces-1)) * (str_len / pieces)
  if((pieces-1) <= length(potential_breakpoints)) {
    absolute_breakpoints <- map_int(potential_breakpoints, ~str_whitespace_idx[which(abs(str_whitespace_idx-.x)==min(abs(str_whitespace_idx-.x)))[1]])
    chars[absolute_breakpoints] <- "\n"
  } else {
    absolute_breakpoints <- round(potential_breakpoints,0)
    chars[absolute_breakpoints] <- paste0(chars[absolute_breakpoints],"\n")
  }
  str_return <- paste(chars,collapse = "")
  return(str_return)
}

distinct <- function(x){
  x <- x %>% n_distinct()
  unique(x)
}

wrangle_is_dataframe_full <- function(df, check = c("cols","rows"), remove_NA_rows = TRUE){
  if(remove_NA_rows){
    all_na <- apply(df, 1, function(x) is.na(x) %>% all())
    df <- df[!all_na,]
  }
  dims <- list(rows = nrow(df), cols = ncol(df))
  dims <- dims[names(dims) %in% check]
  full <- map_lgl(dims,~.x != 0) %>% all()
  full
}



wrangle_generate_hub_connections <- function(show_hub, homology_ecm,species_abbreviations,show_max_n_additional_per_species = 10){
  # randomize to be able to randomly show up n additional phenotypes per taxon
  homology_ecm <- homology_ecm[sample(nrow(homology_ecm),replace = FALSE),]
  phenotypes_conserved <- homology_ecm %>% 
    filter(gene_symbol %in% show_hub) %>%
    group_by(object_label) %>%
    summarise(n_species = subject_taxon_label %>% n_distinct()) %>%
    filter(n_species > 1) %>% pull(object_label)
  
  df_conserved <- homology_ecm %>% filter(gene_symbol %in% show_hub) %>% filter(object_label %in% phenotypes_conserved)
  
  df_additional <- homology_ecm %>% 
    filter(gene_symbol %in% show_hub) %>%
    filter(!object_label %in% phenotypes_conserved) %>%
    group_by(subject_taxon_label) %>%
    slice(1:show_max_n_additional_per_species)
  
  df_to_plot <- bind_rows(df_conserved,df_additional)
  
  # edges
  relations_homology <- df_to_plot %>% rename(from = homolog, to = subject_taxon_label) %>% mutate(interaction = "homolog") %>% filter(!is.na(from) & !is.na(to)) %>% select(from,to,interaction) %>% dplyr::distinct()
  relations_pheno <- df_to_plot %>% rename(from = subject_taxon_label, to = object_label) %>% mutate(interaction = from) %>% filter(!is.na(from) & !is.na(to)) %>% select(from,to,interaction) %>% dplyr::distinct()
  relations <- bind_rows(relations_homology,relations_pheno) %>% ungroup()
  relations <- relations %>% group_by(to) %>% mutate(n_incoming = n()) %>% 
    rowwise() %>%
    mutate(names_unfolded = ifelse(n_incoming > 1,paste0(to,", ", species_abbreviations %>% filter(subject_taxon_label == interaction) %>% pull(species_abbreviation),""),to)) %>% 
    ungroup()
  # nodes
  vertices <- df_to_plot %>% 
    select(homolog,subject_taxon_label,object_label) %>% 
    mutate(species = subject_taxon_label)
  vertices <- vertices %>% 
    gather(key = "category",value = "name",-species) 
  vertices <- vertices %>%
    group_by(name) %>% 
    dplyr::distinct() %>%
    mutate(n_occurences = n(), 
           shared_regulation = ifelse(n_occurences > 1 & !(name %in% show_hub), "shared","not-shared"), 
           origin = case_when(name %in% show_hub ~ "hub", n_occurences > 1 ~ "shared", TRUE ~ species),
           origin_label = ifelse(name %in% unique(origin),name,"other"),
           root = ifelse(name %in% show_hub,"hub","no-hub")) %>%
    ungroup() %>% 
    dplyr::distinct() %>%
    rowwise() %>%
    mutate(names_unfolded = ifelse(n_occurences > 1 & root != "hub",paste0(name,", ", species_abbreviations %>% filter(subject_taxon_label == species) %>% pull(species_abbreviation),""),name)) %>%
    ungroup() %>%
    select(name, origin, n_occurences, shared_regulation, origin_label, category, root,species,names_unfolded)
  g_folded <- graph_from_data_frame(relations, directed=TRUE, vertices=vertices %>% filter(!duplicated(name))) %>% as_tbl_graph()
  g_unfolded <- graph_from_data_frame(relations %>% mutate(to = names_unfolded), directed=TRUE, vertices=vertices %>% mutate(name = names_unfolded) %>% filter(!duplicated(name))) %>% as_tbl_graph()
  g_unfolded <- g_unfolded %>% activate(nodes) %>% mutate(fontface = ifelse(shared_regulation == "shared","italic","italic")) # change font if shared phenotypes should be highlighted
  return(list(g_folded = g_folded, g_unfolded = g_unfolded))
}




draw_circular_dendrogram <- function(graph,
                                     leaf_size = 2.5,
                                     hub_size = 6,
                                     species_size = 1.5,
                                     species_repel = TRUE,
                                     legend_position = "none",
                                     point_size_range = c(1.2,3.5),
                                     species_color){
  graph %>%
    ggraph(layout = 'dendrogram', circular = TRUE) + 
    geom_edge_diagonal(color = "grey") +
    geom_node_text(aes(x = x*1.05, y=y*1.05, 
                       label = name,
                       filter=leaf,
                       alpha = shared_regulation,
                       fontface = fontface,
                       color = origin,
                       angle = atan(y/x)*360/(2*pi), hjust='outward'), 
                   size = leaf_size) +
    geom_node_label(mapping = aes(filter=origin=="hub", label = name), size = hub_size) +
    geom_node_label(mapping = aes(filter=name%in%all_species, label = name,color=origin), size = species_size, alpha = 0.8, repel = species_repel,fontface = "italic") +
    scale_edge_colour_distiller(palette = "RdPu") +
    geom_node_point(aes(filter = leaf, x = x*1.0, y=y*1.0, color=origin, size=n_occurences, alpha = shared_regulation)) +
    scale_size_continuous( range = point_size_range ) +
    scale_alpha_manual(values = c("shared" = 1, "not-shared" = 0.6)) +
    scale_color_manual(values = species_color) +
    theme_void() +
    guides(color = guide_legend(title = "Species",override.aes = list(size = 5)),
           alpha = "none",
           size = "none") +
    theme(
      legend.position= legend_position,
      plot.margin=unit(c(0,0,0,0),"cm")
    ) +
    expand_limits(x = c(-2, 2), y = c(-2, 2))
}


wrangle_generate_ecm_network <- function(tab_ecm){
  tab_ecm <- tab_ecm %>% ungroup() %>% filter(!is.na(gene_symbol) & !is.na(object_label) & !is.na(subject_taxon_label)) %>% rename(species = subject_taxon_label)
  # edges
  relations <- tab_ecm %>% rename(from = homolog, to = object_label) %>% mutate(interaction = "homolog") %>% filter(!is.na(from) & !is.na(to)) %>% select(from,to,species,interaction,involved_n_species) %>% dplyr::distinct()
  relations <- relations %>% group_by(to) %>% mutate(n_genes_per_pheno = n()) %>% ungroup() %>% 
    group_by(from) %>% mutate(n_affected_pheno_per_gene = n()) %>% 
    ungroup() %>% mutate(n_affected_pheno_per_gene_order = order(-n_affected_pheno_per_gene), n_genes_per_pheno_order = order(-n_genes_per_pheno))
  # vertices
  vertices_gene <- tab_ecm %>% select(gene_symbol, everything()) %>% mutate(type = "gene") %>% rename(name = gene_symbol) %>% select(name:Refseq_IDs,type,gene_in_n_species:involved_n_species)
  vertices_pheno <- tab_ecm %>% select(object_label, everything()) %>% mutate(type = "phenotype") %>% rename(name = object_label) %>% select(name, type,gene_in_n_species:involved_n_species)
  vertices <- bind_rows(vertices_gene,vertices_pheno) %>% select(name,type, everything())
  vertices <- vertices %>% select(name:Refseq_IDs,gene_in_n_species:involved_n_species)  %>% filter(!duplicated(name))
  # make network
  net <- graph_from_data_frame(relations, directed=TRUE, vertices=vertices) %>% as_tbl_graph()
  net <- net %>% 
    activate(edges) %>%
    mutate(centrality_bet = centrality_edge_betweenness()) %>%
    activate(nodes) %>%
    mutate(centrality_auth = centrality_authority(),
           # community_infomap = group_infomap() %>% as.character(),
           centrality_page = centrality_pagerank()
           )
  return(net)
}

wrangle_generate_any_network <- function(tab, species_order){
  if(nrow(tab) == 0) return(NULL)
  tab <- tab %>% ungroup() %>% filter(!is.na(gene_symbol) & !is.na(object_label) & !is.na(subject_taxon_label)) %>% rename(species = subject_taxon_label)
  # edges
  relations <- tab %>% rename(from = homolog, to = object_label) %>% mutate(interaction = "homolog") %>% filter(!is.na(from) & !is.na(to)) %>% select(from,to,species,interaction,involved_n_species) %>% dplyr::distinct()
  relations <- relations %>% group_by(to) %>% mutate(n_genes_per_pheno = n()) %>% ungroup() %>% 
    group_by(from) %>% mutate(n_affected_pheno_per_gene = n()) %>% 
    ungroup() %>% mutate(n_affected_pheno_per_gene_order = order(-n_affected_pheno_per_gene), n_genes_per_pheno_order = order(-n_genes_per_pheno))
  # vertices
  vertices_gene <- tab %>% select(gene_symbol, everything()) %>% mutate(type = "gene") %>% rename(name = gene_symbol) 
  vertices_pheno <- tab %>% select(object_label, everything()) %>% mutate(type = "phenotype") %>% rename(name = object_label) 
  vertices <- bind_rows(vertices_gene,vertices_pheno) %>% select(name,type, everything())
  vertices <- vertices %>% filter(!duplicated(name))
  # make network
  net <- graph_from_data_frame(relations, directed=TRUE, vertices=vertices) %>% as_tbl_graph()
  net <- net %>% 
    activate(edges) %>%
    mutate(species = species %>% factor(.,levels = species_order)) %>%
    mutate(centrality_bet = centrality_edge_betweenness()) %>%
    activate(nodes) %>%
    mutate(species = species %>% factor(.,levels = species_order)) %>%
    mutate(centrality_auth = centrality_authority(),
           # community_infomap = group_infomap() %>% as.character(),
           centrality_page = centrality_pagerank()
    )
  return(net)
}

wrangle_hub_characteristic <- function(df, grouping_variable){
  grouping_variable <- enquo(grouping_variable)
  df %>% 
    group_by(!!grouping_variable) %>% 
    summarize(
      in_n_species = n_distinct(taxon),
      in_n_phenotypes = n_distinct(object_label),
      in_n_genes = n_distinct(gene_symbol)
    ) %>% arrange(desc(in_n_species)) %>% 
    ungroup() %>% tbl_df
}


gene_pheno_net_subset <- function(df,genes = NULL,phenotypes = NULL,n_genes = NULL,n_phenotypes = NULL,filter_multi_species = FALSE){
  # if genes are provided directly they superceed the n_genes input
  sub <- df
  sub_top_genes <- sub %>% wrangle_hub_characteristic(grouping_variable = gene_symbol) %>% rename(gene_in_n_species = in_n_species, gene_in_n_phenotypes = in_n_phenotypes) %>% select(gene_symbol,gene_in_n_species,gene_in_n_phenotypes)
  sub_top_phenotypes <- sub %>% wrangle_hub_characteristic(grouping_variable = object_label) %>% rename(phenotype_in_n_species = in_n_species, phenotype_in_n_genes = in_n_genes) %>% select(object_label,phenotype_in_n_species,phenotype_in_n_genes)
  sub_keyhubs <- left_join(sub,sub_top_genes, by = "gene_symbol") %>% left_join(.,sub_top_phenotypes, by = "object_label")
  if(length(genes) > 0){sub_keyhubs <- sub_keyhubs %>% filter(gene_symbol %in% genes)
  }else{
    message("Computing top ", n_genes, " genes")
    top_genes <- sub_keyhubs %>% arrange(desc(gene_in_n_phenotypes)) %>% pull(gene_symbol) %>% unique() %>% .[1:n_genes]
    sub_keyhubs <- sub_keyhubs %>% filter(gene_symbol %in% top_genes)
  }
  if(length(phenotypes) > 0){sub_keyhubs <- sub_keyhubs %>% filter(object_label %in% phenotypes) 
  }else{
    message("Computing top ", n_phenotypes, " phenotypes")
    top_pheno <- sub_keyhubs %>% arrange(desc(phenotype_in_n_genes)) %>% pull(object_label) %>% unique() %>% .[1:n_phenotypes]
    sub_keyhubs <- sub_keyhubs %>% filter(object_label %in% top_pheno)
  }  
  involved_n_species <- sub_keyhubs %>% select(gene_symbol,object_label,subject_taxon_label) %>% group_by(gene_symbol,object_label) %>% summarise(involved_n_species = n_distinct(subject_taxon_label,na.rm = TRUE)) %>% ungroup()
  out <- left_join(sub_keyhubs,involved_n_species, by = c("gene_symbol","object_label"))
  if(filter_multi_species) out <- out %>% filter(involved_n_species > 1)
  return(out)
}

gene_pheno_net_plot <- function(network,title = "network", node_color_fill,species_color,border = TRUE,margin = 30, nrow = 1,size_text = 2.3,size_label = 2.5,edge_widths = c(0.3,0.8,1.6,2.5,3,3.6), edge_alphas = c(1,1,1,1,1,1,1,1)){
  if(is.null(network)) return(NULL)
  node_color_fill <- enquo(node_color_fill)
  # requires ggraph network input, generate with wrangle_generate_*_network function family
  plot <- network %>% 
    ggraph(layout = 'nicely') + 
    geom_edge_link(aes(color = species, edge_alpha = as.factor(involved_n_species), edge_width = as.factor(involved_n_species)),arrow = arrow(length = unit(1.5, 'mm')), end_cap = circle(1.7, 'lines'), start_cap = circle(1, 'lines')) + 
    geom_node_text(aes(label = name,filter=type=="phenotype"), alpha = 1, size = size_text, fontface = "bold", color = "black") +
    geom_node_label(aes(label = name,filter=type=="gene",  color = !!node_color_fill, fill = !!node_color_fill), alpha = 0.03, size = size_label) +
    labs(title = title) +
    scale_edge_color_manual(values = species_color) +
    scale_edge_width_manual(values = edge_widths) +
    scale_edge_alpha_manual(values = edge_alphas) +
    scale_color_manual(values = c("Core matrisome" = "#5E2C83","Matrisome-associated" = "#005FB8","Retired" = "#32AB43")) +
    scale_fill_manual(values = c("Core matrisome" = "#5E2C83","Matrisome-associated" = "#005FB8","Retired" = "#32AB43")) +
    facet_edges(. ~ species,nrow = nrow) +
    theme_graph(base_family = 'Helvetica',border = border,foreground = "grey40",plot_margin = margin(rep(margin,4))
    ) + theme(strip.background = element_blank(),strip.text.x  = element_text(face = "bold.italic", size = 12),panel.spacing.x=unit(6, "lines"))
  return(plot)
}



wrangle_homology_ecm_huborder <- function(df){
  df %>% 
    group_by(gene_symbol) %>% 
    summarize(
      division = unique(division),
      category = unique(category),
      in_n_species = n_distinct(taxon),
      in_n_phenotypes = n_distinct(object_label),
      in_n_genes = n_distinct(gene_symbol),
      species = paste(unique(subject_taxon_label), collapse = ", "),
      cross_species_phenome = paste(unique(object_label), collapse = ", ")
    ) %>% arrange(desc(in_n_species)) %>% 
    select(division,category,in_n_species,in_n_phenotypes,in_n_genes,gene_symbol,species,cross_species_phenome) %>% 
    ungroup() %>% tbl_df
}


wrangle_remove_empty_strings <- function(vec){
  vec[str_squish(vec) != ""]
}

wrangle_phenotype_grouping <- function(pheno_df, group_df){
  pheno_grouped <- pheno_df %>% mutate(phenotype_original = TRUE, object_label = object_label %>% tolower())
  for(idx in seq_len(nrow(group_df))){
    name <- group_df[idx,] %>% pull(group_term)
    message(paste("Processing phenotype group",name))
    target <-  group_df[idx,] %>% pull(targets_trim) %>% unlist() %>% str_c(.,collapse = "|")
    target <- target %>% gsub(pattern = "rrrrr",replacement = "(?!\\S)",fixed = TRUE)
    occurance <- pheno_grouped$object_label %>% grepl(pattern = target,x = .,perl = TRUE)
    pheno_grouped[occurance & pheno_grouped$phenotype_original,"object_label"] <- name
    pheno_grouped[occurance & pheno_grouped$phenotype_original,"phenotype_original"] <- FALSE
  }

  grouping_table <- pheno_grouped %>% 
    filter(phenotype_original == FALSE) %>% group_by(object_label) %>% nest() %>% ungroup() %>%
    mutate(`# occurences` = map_dbl(data, ~.x %>% nrow()),
           `replacement table` = map(data, ~.x %>% dplyr::select(phenotype_as_imported) %>% group_by(phenotype_as_imported) %>% summarise(n = n()) %>% arrange(desc(n))),
           `Phenotypes replaced` = map_chr(`replacement table`, ~.x$phenotype_as_imported %>% paste(collapse = ", "))
    ) %>%
    arrange(desc(`# occurences`))
  
  grouping_information <- list(frac_captured_in_groups = round((1-sum(pheno_grouped$phenotype_original)/nrow(pheno_grouped)) * 100,3),
                               grouping_statistics = grouping_table)
  message(paste("--> ",grouping_information$frac_captured_in_groups, "% phenotypes captured in groups"))
  return(list(pheno_grouped = pheno_grouped, grouping_statistics = grouping_information$grouping_statistics, frac_captured_in_groups = grouping_information$frac_captured_in_groups))
}

wrangle_grouping <- function(imported_grouping_file, gene_pheno_df){
  pheno_grouping_df <- imported_grouping_file %>% 
    mutate(matched_substrings = matched_substrings %>% tolower()) %>%
    rowwise() %>%
    mutate(matched_substrings_split = map(matched_substrings, ~str_split(string = matched_substrings,pattern = ",")),
           targets_trim = map(matched_substrings_split, ~.x %>% str_trim() %>% wrangle_remove_empty_strings()))
  # return objects
  output_grouping <- wrangle_phenotype_grouping(pheno_df = gene_pheno_df, group_df = pheno_grouping_df)
  grouping_statistics <-  output_grouping$grouping_statistics %>% select(object_label, `replacement table`) %>% unnest(cols = c(`replacement table`)) 
  fail_to_group <- output_grouping$pheno_grouped %>% filter(phenotype_original) %>% select(object_label,phenotype_as_imported,phenotype_original) 
  return(list(grouping_statistics = grouping_statistics, fail_to_group = fail_to_group, output_grouped = output_grouping$pheno_grouped))
}


# R DIOPT function mentioned in homologene package (code from here: https://rdrr.io/github/oganm/homologene/src/R/diopt.R)
diopt = function(genes, inTax, outTax, delay = 10){
  # rtxt = robotstxt::robotstxt(domain = "flyrnai.org")
  # delay = rtxt$crawl_delay %>% filter(useragent =='*') %$% value %>% as.integer()
  session = rvest::html_session('https://www.flyrnai.org/cgi-bin/DRSC_orthologs.pl')
  form = rvest::html_form(session)[[1]]
  
  acceptableInTax= form$fields$input_species$options
  acceptableOutTax = form$fields$output_species$options
  
  assertthat::assert_that(inTax %in% acceptableInTax)
  assertthat::assert_that(outTax %in% acceptableOutTax)
  
  form = rvest::set_values(form,
                           input_species = inTax,
                           output_species = outTax,
                           gene_list = paste(genes,collapse = '\n\r'))
  
  additional_filters = which(names(form$fields) == 'additional_filter')
  
  additional_filter_names = form$fields[additional_filters] %>% purrr::map_chr('value')
  
  form$fields[additional_filters][additional_filter_names %in% 'None'][[1]]$checked = 'checked'
  form$fields[additional_filters][additional_filter_names %in% 'NoLow'][[1]]$checked = NULL
  
  Sys.sleep(delay)
  
  response = rvest::submit_form(session,form)
  
  # writeLines(as.char(response$response),'hede.html')
  # utils::browseURL('hede.html')
  
  output = response %>% 
    rvest::html_node('#results') %>% 
    rvest::html_table() %>% 
    dplyr::select(-`Gene2FunctionDetails`,-`Feedback`,-`Alignment & Scores`)
  return(output)
}
# use: diopt(c("emb-9","col-99","cdc-42","col-10"),inTax = 6239, outTax = 9606)


wrangle_hom_mouse_to_human <- function(df,orignal_taxa = 10090,target_taxa = 9606){
  trans_df_raw <- homologene(df$`Gene Symbol`, inTax = orignal_taxa, outTax = target_taxa) %>% as_tibble() #,db = homologene_db_20200327) # now using and updated database since 27. März
  trans_df <- trans_df_raw %>% group_by(`10090`) %>% summarise(combined_orthology = paste(`9606`,collapse = "; ")) %>% rename(`Gene Symbol` = `10090`)
  df_out <- left_join(df,trans_df,by = "Gene Symbol")
  return(df_out)
}

wrangle_random_sampling <- function(universe,n_genes){
  sample(x = universe,size = n_genes,replace = FALSE)
}
wrangle_string_network_stats <- function(genes_of_interest,universe,network,fold){
  subset_network <- network %>% activate(nodes) %>% filter(name %in% genes_of_interest)
  subset_stats <- subset_network %>% activate(nodes) %>% mutate(n_edges = graph_size(), n_nodes = graph_order(), n_edges_per_node = n_edges/n_nodes,
                                                                graph_adhesion = graph_adhesion(), mean_dist_btw_all_nodepairs = graph_mean_dist(),
                                                                graph_girth = graph_girth()) %>% as_tibble()
  
  out <- subset_stats %>% select(-name) %>% mutate(fold = fold) %>% select(fold,everything()) %>% .[1,]
  out
}

wrapper_network_statistics <- function(analyze_phenotype,df_gene_pheno,n_times_resample,universe){
  # Phenotype
  genes_in_pheno_group <- df_gene_pheno %>% filter(object_label %in% c(analyze_phenotype)) %>% pull(`Gene Symbol`) %>% unique()
  pheno_stats_summary <- genes_in_pheno_group %>% wrangle_string_network_stats(genes_of_interest = .,universe = universe,network = fullecm, fold = "none") %>%
    summarise(n_edge_per_node_mean = mean(n_edges_per_node),n_edge_per_node_sd = sd(n_edges_per_node), dist_btw_all_nodepairs_mean = mean(mean_dist_btw_all_nodepairs), dist_btw_all_nodepairs_sd = sd(mean_dist_btw_all_nodepairs)) %>%
    mutate(type = "phenotype genes",phenotype = paste(analyze_phenotype,collapse = ", "), n_genes = length(genes_in_pheno_group))
  # Random Sample
  random_stats <-c(1:n_times_resample) %>% map_df(~wrangle_random_sampling(universe = universe,n_genes = length(genes_in_pheno_group)) %>% wrangle_string_network_stats(genes_of_interest = .,universe = universe,network = fullecm, fold = .x))
  random_stats_summary <- random_stats %>% summarise(n_edge_per_node_mean = mean(n_edges_per_node,na.rm = TRUE),n_edge_per_node_sd = sd(n_edges_per_node,na.rm = TRUE), dist_btw_all_nodepairs_mean = mean(mean_dist_btw_all_nodepairs,na.rm = TRUE), dist_btw_all_nodepairs_sd = sd(mean_dist_btw_all_nodepairs,na.rm = TRUE))%>%
    mutate(type = "randomly sampled genes",phenotype = paste(analyze_phenotype,collapse = ", "), n_genes = length(genes_in_pheno_group))
  # Prepare for export
  pheno_tab_enr <- bind_rows(pheno_stats_summary,random_stats_summary) 
  increase_pernode <- round((pheno_tab_enr[1,1]/pheno_tab_enr[2,1])*100,0)
  edge_per_node_increase <- c(increase_pernode,increase_pernode) %>% unname() %>% unlist() %>% enframe(value = "node connectivity relative to randomly sampled [%]") %>% select(`node connectivity relative to randomly sampled [%]`)
  increase_meandist <- round((pheno_tab_enr[1,"dist_btw_all_nodepairs_mean"]/pheno_tab_enr[2,"dist_btw_all_nodepairs_mean"])*100,0)
  dist_btw_all_nodepairs_increase <- c(increase_meandist,increase_meandist) %>% unname() %>% unlist() %>% enframe(value = "mean distance relative to randomly sampled [%]") %>% select(`mean distance relative to randomly sampled [%]`)
  out <- bind_cols(pheno_tab_enr,edge_per_node_increase) %>% bind_cols(.,dist_btw_all_nodepairs_increase) %>% mutate(`# resampling` = n_times_resample, `Genes in phenotype group` = paste(genes_in_pheno_group,collapse = ", ")) %>% select(phenotype, type, n_genes, everything())
  return(out)
}

wrangle_bool_empty_set <- function(x){
  ifelse(n_distinct(x,na.rm = TRUE) > 0,TRUE,FALSE)
}
wrangle_ndistinct_in_cell <- function(x){
  n_distinct(x,na.rm = TRUE)
}

# Stats for enrichment
test_overrepresentation <- function(n_background, n_A, n_B, n_overlap){
  mat <- matrix(c(n_overlap, n_A-n_overlap, n_B-n_overlap, n_background-n_A-n_B+n_overlap), 2, 2)
  fisher.test(mat, alternative='greater')
}

wrangle_stats_overlap <- function(df,ecm_categories,pheno_of_interest){
  # within each category
  A_B_each <- df %>% select(-category,-gene) %>% colSums() %>% enframe() 
  n_A <- A_B_each %>% filter(name %in% pheno_of_interest) %>% pull(value)# n_A: number of unique genes associated with phenotype (e.g. Brain phenotype) (A and B can be exchanged)
  n_B <- A_B_each  %>% filter(name %in% ecm_categories) %>% pull(value)  # n_B: number of unique genes belonging to a matrisome category (e.g. Collagens) (A and B can be exchanged)
  overlap_bool <- df %>% select(-category,-gene) %>% rowSums() == 2 # n_overlap: number of genes BOTH in matrisome category and associated with the phenotype
  n_overlap <- sum(overlap_bool)
  overlapping_genes <-  df[overlap_bool,] %>% arrange(gene) %>% pull(gene) %>% paste(.,collapse = ", ")
  n_background_genes <- df  %>% filter(!duplicated(gene)) %>% pull(gene) %>% n_distinct(na.rm = TRUE) # n_background: total number of unique gene-phenotype pairs --> this then answers the question what the odds are that the drawn genes (e.g. a matrisome category) all also belong to a certain phenotype (e.g. Brain phenotype)
  overrep_test <- test_overrepresentation(n_background = n_background_genes,n_A = n_A,n_B = n_B,n_overlap = n_overlap)
  ecm_cat <- A_B_each %>% filter(name %in% ecm_categories) %>% pull(name)
  phenotype <- A_B_each %>% filter(name %in% pheno_of_interest) %>% pull(name)
  tibble(phenotype = phenotype, category = ecm_cat, pval = overrep_test$p.value, oddsratio = overrep_test$estimate, alternative = overrep_test$alternative, `n genes in phenotype` = n_A, `n genes in matrisome category` = n_B, `n genes shared` = n_overlap, `genes shared` = overlapping_genes, `n genes in background` = n_background_genes)
}


wrangle_ecmcategory_enriched_phenotypes <- function(df_gene_pheno,taxon_of_interest,pheno_of_interest,ecm_categories, show_venn_legend = TRUE, venn_fontsize = 6, do_not_plot = FALSE){
  # Generate overall parameters
  message(paste("processing:",pheno_of_interest))
  covered_pheno <- df_gene_pheno %>% filter(subject_taxon_label == taxon_of_interest) %>% pull(object_label) %>% unique()
  pheno_of_interest <- pheno_of_interest[pheno_of_interest %in% covered_pheno]
  palette_phenotypes_grey <- grey.colors(n = length(pheno_of_interest),start = 0.45,end = 0.9,alpha = 1) %>% setNames(pheno_of_interest) #grey.colors
  palette_phenotypes_grey_transparent <- grey.colors(n = length(pheno_of_interest),start = 0.45,end = 0.9,alpha = 0.3) %>% setNames(pheno_of_interest) #grey.colors
  palette_phenotypes_with_multimap_grey <- grey.colors(n = length(pheno_of_interest)+1,start = 0.45,end = 0.9,alpha = 1)%>% setNames(c(pheno_of_interest,"Involved in \n multiple phenotypes"))    
  palette_categories_venn <- viridis(length(ecm_categories),alpha = 0.4)
  ########################################################################
  ###############################   Venn   ###############################
  ########################################################################
  p_venns <- list()
  stats_venn <- list()
  for (pheno_single in pheno_of_interest) {
    palette_phenotypes_single <- palette_phenotypes_grey_transparent[names(palette_phenotypes_grey_transparent) %in% pheno_single]
    message(pheno_single)
    occurance <- df_gene_pheno %>% 
      filter(subject_taxon_label == taxon_of_interest) %>%
      select(category, object_label,gene_symbol) %>%
      mutate(gene = gene_symbol) %>%
      pivot_wider(names_from = object_label,values_from = gene_symbol,values_fn = list(gene_symbol = wrangle_bool_empty_set)) %>%
      mutate_all(funs(replace_na(., FALSE))) %>%
      select(one_of(c("category","gene",pheno_single)))
    out <- ecm_categories %>% map(~occurance %>% 
                                    mutate(category = ifelse(category == .x,.x,"Other matrisome"),category_to_spread = category,gene_to_spread = gene) %>%
                                    pivot_wider(names_from = category_to_spread,values_from = gene_to_spread,values_fn = list(gene_to_spread = wrangle_bool_empty_set)) %>%
                                    select(-`Other matrisome`) %>%
                                    mutate_all(funs(replace_na(., FALSE))) %>%
                                    mutate(category = .x) %>% 
                                    mutate_all(funs(replace_na(., FALSE)))) 
    to_plot <- out %>% 
      reduce(.f = bind_rows) %>%
      mutate_all(funs(replace_na(., FALSE))) %>%
      select(-gene) %>% 
      select(-one_of(pheno_single), everything(),one_of(pheno_single)) %>% # The order of the columns in the dataframe is absolutely crucial
      arrange(category)
    # test enrichement for overlap between phenotype and each matrisome category 
    stats_venn[[pheno_single]] <-  out %>% map_df(~.x %>% wrangle_stats_overlap(df = .,ecm_categories = ecm_categories,pheno_of_interest = pheno_of_interest))
    pval_labs_plot <- stats_venn[[pheno_single]] %>% select(category,pval) %>% mutate(pval = formatC(x = pval,format = "e", digits = 2) %>% paste("\n Pval:",.,sep = " "))
    df_for_venn <- left_join(to_plot,pval_labs_plot,by="category") %>% mutate(category_pval = paste(category,pval)) %>% select(-pval,-category) %>% select(category_pval,everything())
    # Plot Venn diagram for overlap between phenotype and each matrisome category 
    gridfit <- eulerr::euler(df_for_venn, by = category_pval, shape = "circle")
    n_ECM_categories <- ecm_categories %>% n_distinct()
    p_venns[[pheno_single]] <- plot(gridfit,
                                    quantities = list(type = c("counts"),fontsize = venn_fontsize), #c("percent", "counts")
                                    fills = c(palette_categories_venn,palette_phenotypes_single), # works best
                                    edges = list(lty = c(1), col = c("black"),lex = c(rep(3,n_ECM_categories),rep(1,length(palette_phenotypes_single)))),
                                    legend = show_venn_legend)  %>% as.grob()
  }
  stats_venn_table <- stats_venn %>% reduce(.f = bind_rows,.x = .)
  if(do_not_plot) return(stats_venn_table %>% mutate(species = taxon_of_interest))
  ########################################################################
  ############################### Alluvial ###############################
  ########################################################################
  alluvial_label_size <- 5
  sub_species <- df_gene_pheno %>% 
    filter(category %in% names(palette_category_viridis)) %>%
    filter(subject_taxon_label == taxon_of_interest,object_label %in% pheno_of_interest) %>%
    select(gene_symbol,division,category, object_label) %>%
    dplyr::distinct()
  df_alluvial_wide <- sub_species %>% group_by(gene_symbol,division,category) %>% summarise(n_phenotypes = object_label %>% n_distinct(),
                                                                                            phenotype = ifelse(n_phenotypes == 1, object_label,"Involved in \n multiple phenotypes"),
                                                                                            freq = 1) %>%ungroup()
  df_alluvial_long<- df_alluvial_wide %>% pivot_longer(c(-gene_symbol,-n_phenotypes,-freq),names_to = "main_group",values_to = "stratum") %>%
    mutate(main_group = factor(main_group,levels = c("division","category","phenotype")),
           stratum = factor(stratum,levels = names(c(palette_category_viridis,palette_phenotypes_with_multimap_grey))) ) %>%
    filter(main_group != "division") %>%
    mutate(label = case_when(main_group == "phenotype" ~ "right", 
                             main_group == "category" ~ "left", 
                             TRUE ~ "none")) %>% arrange(main_group)
  p_alluvial <- df_alluvial_long %>%
    ggplot(aes(x = main_group, stratum = stratum, alluvium = gene_symbol, y = freq,
               fill = stratum, color = stratum)) +
    scale_x_discrete(expand = c(.4, 0)) +
    geom_flow(width = 1/7, aes.flow = "forward") + #backward
    geom_stratum(width = 1/7,color = "white",size = 1.5) +
    scale_linetype_manual(values = c("blank", "solid")) +
    geom_text_repel(
      aes(label = ifelse(label == "right", as.character(stratum), NA)),
      stat = "stratum", size = alluvial_label_size, direction = "y", nudge_x = .5, color = "black") +
    # geom_text_repel(
    #   aes(label = ifelse(label == "left", as.character(stratum), NA)),
    #   stat = "stratum", size = alluvial_label_size, direction = "y", nudge_x = -.5, color = NA, alpha = 0.3) +
    geom_text_repel(
      aes(label = ifelse(label == "left", as.character(stratum), NA)),
      stat = "stratum", size = alluvial_label_size, direction = "y", nudge_x = -.5, color = "black", alpha = 1) +
    labs(title = "",x = "Matrisome to phenotype relation") +
    scale_fill_manual(values = c(palette_category_viridis,palette_phenotypes_with_multimap_grey)) +
    scale_color_manual(values = c(palette_category_viridis,palette_phenotypes_with_multimap_grey)) +
    theme_void() +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  ########################################################################
  ###############################  UpSetR  ###############################
  ########################################################################
  for_upset <- sub_species %>% 
    mutate(category_to_spread = category,gene_to_spread_cat = gene_symbol,gene_to_spread_pheno = gene_symbol) %>%
    pivot_wider(names_from = category_to_spread,values_from = gene_to_spread_cat,values_fn = list(gene_to_spread_cat = wrangle_bool_empty_set))  %>%
    pivot_wider(names_from = object_label,values_from = gene_to_spread_pheno,values_fn = list(gene_to_spread_pheno = wrangle_bool_empty_set)) %>%
    mutate_all(funs(replace_na(., FALSE))) %>%
    as.data.frame() %>%
    select(-category) %>%
    mutate_if(is.logical,as.numeric) %>%
    column_to_rownames("gene_symbol") %>%
    select(one_of(upsetr_cols <- c(ecm_categories,pheno_of_interest))) 
  # meta table & upsetr plot
  sets <- for_upset %>% colnames() %>% rev()
  meta <- tibble(sets,palette = c(palette_category_viridis,palette_phenotypes_grey), datatype = ifelse(sets %in% names(palette_phenotypes_grey),"Phenotype","Matrisome"))
  palette_upsetr <- c(palette_phenotypes_grey %>% rev(),palette_category_viridis %>% rev())
  bg_color_matrisome_part <- palette_category_viridis[[ceiling(length(palette_category_viridis)/2)]]
  bg_color_pheno_part <- palette_phenotypes_grey[[ceiling(length(palette_phenotypes_grey)/2)]]
  p_upsetr <- upset(data = for_upset,point.size = 5.5,text.scale = 2,mainbar.y.label = "# genes which are shared",sets.x.label = "# genes in each group",
                    set.metadata = list(data = meta, plots = list(list(type = "matrix_rows",column = "datatype",keep.order = TRUE, colors = c("Phenotype" = bg_color_pheno_part, "Matrisome" = "white"),alpha = 0.5))),
                    sets = sets, keep.order = TRUE, mb.ratio = c(0.3, 0.7),sets.bar.color=palette_upsetr, order.by = "freq")
  p_upsetr_cowplot_ready <- cowplot::plot_grid(NULL, 
                                             p_upsetr$Main_bar,
                                             p_upsetr$Sizes, 
                                             p_upsetr$Matrix,
                                             nrow=2, align='hv', rel_heights = c(0.5,1),rel_widths = c(1.4,3))
  ########################################################################
  ###############################  Export  ###############################
  ########################################################################
  output_list <- list(venn = p_venns,
                      alluvial = p_alluvial,
                      upsetr = p_upsetr,
                      p_upsetr_cowplot_ready = p_upsetr_cowplot_ready,
                      stats_table = stats_venn_table)
  return(output_list)
}


plot_ecmcategory_enriched_phenotypes <- function(output_list, path_output, plot_hight = 28,plot_width = 20){
  print("Export function")
  plot_venns <- plot_grid(plotlist = output_list$venn, ncol = 1)
  row1 <- plot_grid(output_list$alluvial,plot_venns, labels = c('A', 'B'),align = "hv", label_size = 30,rel_widths = c(1,2.2))
  row2 <- plot_grid(output_list$p_upsetr_cowplot_ready, labels = c('C'),align = "hv", label_size = 30)
  p_out <- plot_grid(row1,row2, label_size = 12, ncol = 1,rel_heights = c(1.3,1))
  p_out %>% ggexport(filename = path_output,width = plot_width,height = plot_hight,units = "mm")
}

plot_ecmcategory_enriched_phenotypes_separate <- function(output_list, path_output, 
                                                          venn_hight = 13,venn_width = 16,
                                                          alluvial_height = 10, alluvial_width = 6, 
                                                          upset_height = 14, upset_width = 20){
  print("Export function for separate plots")
  plot_venns <- plot_grid(plotlist = output_list$venn, ncol = 1)
  plot_venns %>% ggexport(filename = paste0(path_output,"_venn.pdf"),width = venn_width,height = venn_hight,units = "mm")
  plot_alluvial <- plot_grid(output_list$alluvial)
  plot_alluvial %>% ggexport(filename = paste0(path_output,"_alluvial.pdf"),width = alluvial_width,height = alluvial_height,units = "mm")
  plot_upset <- plot_grid(output_list$p_upsetr_cowplot_ready)
  plot_upset %>% ggexport(filename = paste0(path_output,"_upset.pdf"),width = upset_width,height = upset_height,units = "mm")
}

wrangle_file_path <- function(parent_dir,species,phenotype){
  species <- species %>% str_remove_all(string = .,pattern = " ") 
  phenotype <- phenotype %>% str_remove_all(string = .,pattern = " ") %>% str_remove_all(string = .,pattern = "/")
  filename <- paste0(parent_dir,species, phenotype,"fig.pdf") %>% str_remove_all(string = .,pattern = " ") 
  return(filename)
}



wrangle_decorate_table_ecmcontribution <- function(category,pval,pval_thr = 0.05){
  df <- tibble(category,pval)
  only_matrisome <- df %>% filter(category != "Non-matrisome")
  n_matrisome_categories <- only_matrisome %>% filter(pval < pval_thr) %>% pull(category) %>% n_distinct()
  case_when(n_matrisome_categories == 0 ~ "No matrisome signature",
            n_matrisome_categories == 1 ~ paste0("Matrisome enriched (",n_matrisome_categories," category)"),
            n_matrisome_categories > 1 ~ paste0("Matrisome enriched (",n_matrisome_categories," categories)"),
            TRUE ~ "Other effect")
}


wrangle_overall_phenotype_enrichment <- function(df_gene_pheno,taxon_of_interest,pheno_of_interest,ecm_and_nonecm_categories){
  # Generate overall parameters
  message(paste("processing:",pheno_of_interest))
  covered_pheno <- df_gene_pheno %>% pull(object_label) %>% unique()
  pheno_of_interest <- pheno_of_interest[pheno_of_interest %in% covered_pheno]
  ########################################################################
  ###############################   Venn   ###############################
  ########################################################################
  stats_venn <- list()
  for (pheno_single in pheno_of_interest) {
    message(pheno_single)
    occurance <- df_gene_pheno %>% 
      select(category, object_label,subject_label) %>%
      mutate(gene = subject_label) %>%
      pivot_wider(names_from = object_label,values_from = subject_label,values_fn = list(subject_label = wrangle_bool_empty_set)) %>%
      mutate_all(funs(replace_na(., FALSE))) %>%
      select(one_of(c("category","gene",pheno_single)))
    out <- ecm_and_nonecm_categories %>% map(~occurance %>% 
                                               mutate(category = ifelse(category == .x,.x,"Other"),category_to_spread = category,gene_to_spread = gene) %>%
                                               pivot_wider(names_from = category_to_spread,values_from = gene_to_spread,values_fn = list(gene_to_spread = wrangle_bool_empty_set)) %>%
                                               select(-`Other`) %>%
                                               mutate_all(funs(replace_na(., FALSE))) %>%
                                               mutate(category = .x) %>% 
                                               mutate_all(funs(replace_na(., FALSE)))) 
    to_plot <- out %>% 
      reduce(.f = bind_rows) %>%
      mutate_all(funs(replace_na(., FALSE))) %>%
      select(-gene) %>% 
      select(-one_of(pheno_single), everything(),one_of(pheno_single)) %>% # The order of the columns in the dataframe is absolutely crucial
      arrange(category) 
    # test enrichement for overlap between phenotype and each matrisome category 
    stats_venn[[pheno_single]] <-  out %>% map_df(~.x %>% wrangle_stats_overlap(df = .,ecm_categories = ecm_and_nonecm_categories,pheno_of_interest = pheno_of_interest) %>% mutate(species = taxon_of_interest))
  }
  stats_venn_table <- stats_venn %>% reduce(.f = bind_rows,.x = .)
  return(stats_venn_table)
}