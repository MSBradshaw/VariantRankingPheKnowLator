mln_files <- function(triples_ids, dir, gene2gene = FALSE, geneyobo = FALSE){
  # packages
  require(tidyverse)
  
  # check if directory exists
  if(dir.exists(dir)) stop('Directory exists. Read files or rename directory.')

  # make directory
  dir.create(dir, showWarnings = FALSE)
  
  # read and clean triples
  trip_ids <- triples_ids %>%
    # read data
    read_tsv() %>%
    janitor::clean_names() %>%
    # create gene indicator for subjects and objects 
    mutate(., across(.cols = c(subject, object), .fns = ~str_detect(.x, 'gene/'), .names = '{.col}_gene'),
           across(.cols = c(subject, object), .fns = ~str_detect(.x, 'obo/'), .names = '{.col}_obo'),
           across(.cols = c(subject, predicate, object), .fns = ~basename(.x), .names = '{.col}_abbr'))
  
  # gene to gene?
  if(gene2gene){
    trip_ids <- filter(trip_ids, subject_gene, object_gene)
  }
  
  # genes and obos?
  else if(geneyobo){
    trip_ids <- filter(trip_ids, (subject_gene | subject_obo), (object_gene | object_obo))   
  }
  
  # create files
  layout <- data.frame(node = c(trip_ids$subject_abbr, trip_ids$object_abbr)) %>%
    distinct(., node) %>%
    mutate(., nodeID = row_number()) %>%
    select(., nodeID, nodeLabel = node)
  write_delim(layout, file = str_c(dir, '/layout.txt'))
  
  layers <- data.frame(layer = trip_ids$predicate_abbr) %>%
    distinct(., layer) %>%
    mutate(., layerID = row_number()) %>%
    select(., layerID, layerLabel = layer)
  write_delim(layers, file = str_c(dir, '/layers.txt'))
  
  edges <- trip_ids %>%
    select(., subject_abbr, predicate_abbr, object_abbr) %>%
    # get parent node and layer
    left_join(., select(layout, subject_abbr = nodeLabel, node1 = nodeID)) %>%
    left_join(., select(layers, predicate_abbr = layerLabel, layer1 = layerID)) %>%
    # get child node and layer
    left_join(., select(layout, object_abbr = nodeLabel, node2 = nodeID)) %>%
    left_join(., select(layers, predicate_abbr = layerLabel, layer2 = layerID)) %>%
    # create weight
    mutate(., wgt = 1L) %>%
    select(., node1:wgt) %>%
    set_names(., c('from', 'from_layer', 'to', 'to_layer', 'wgt'))
  write_delim(edges, file = str_c(dir, '/edgelist.edges'))
  
  config <- str_c(c(str_c(dir, '/edgelist.edges'), 
                    str_c(dir, '/layers.txt'), 
                    str_c(dir, '/layout.txt')), 
                  collapse = ';')
  cat(config, file = str_c(dir, '/config.txt'))
}
