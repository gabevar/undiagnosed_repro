library(tidyverse)
library(panelr)
library(kableExtra)
library(mice)
library(broom.mixed)
# library(brms)
library(future)
library(here)

# General use #

mymode <- function(x) {t <- table(x) 
names(t)[ which.max(t) ]}

# Communication events #

comm_event_fn <- function(comm, comm_type = "") {
  
  waves <- comm %>% ungroup() %>% distinct(wave) %>% arrange(wave) %>% pull() 
  res_list <- list()
  for (i in 1:length(waves)) {
    
    wave_lim <- waves[[i]]
    comm_temp <- comm %>% filter(wave == wave_lim)
    
    directed_graph <- graph.data.frame(comm_temp, directed = TRUE)
    undirected_graph <- as.undirected(directed_graph, mode = "collapse", edge.attr.comb = "sum") # Turn to undirected and sum up values.
    undirected_graph <- igraph::as_data_frame(undirected_graph)
    undirected_graph_rev <- undirected_graph %>% select(from = to, to = from, wave, weight, weight_log, weight_scale) # reverse of each edge.
    comm_temp <- rbind(undirected_graph, undirected_graph_rev) # directed edgelist of undirected edges. 
    
    if (wave_lim == 0) {net <- net_df %>% filter(wave == "Wave1" | wave == "Wave2") %>% select(egoid, alterid) %>% mutate(nom = 1)
    } else {net <- net_df %>% filter(wave == paste0("Wave", wave_lim)) %>% select(egoid, alterid) %>% mutate(nom = 1)}
    
    node_vars <- long_var_df %>% 
      filter(wave == wave_lim & (variable == "undiag")) %>% 
      select(-wave) %>% 
      pivot_wider(values_from = value, names_from = variable)
    
    comm_temp <- comm_temp %>%
      mutate(to = as.numeric(to), from = as.numeric(from)) |> 
      left_join(node_vars, by = c("to" = "egoid")) %>% 
      left_join(net, by = c("from" = "egoid", "to" = "alterid")) %>%
      mutate(nom = case_when(is.na(nom) ~ 0, T ~ nom)) %>% 
      mutate(undiag_weighted = case_when(nom == 1 ~ undiag*2, T ~ undiag), # doubling the value of undiagnosed if the person was also nominated.
             undiag_commweighted = undiag*weight_log) # multiply by the logged number of communication events.
             
    comm_temp <- comm_temp %>% 
      group_by(from) %>% 
      summarise(num_comm = n(), # number diff comm partners
                weight_comm = sum(weight), # number of communications
                # weightscale_comm = mean(weight_scale, na.rm = T), # average communication, scaled.
                sum_comm_undiag = ifelse(all(is.na(undiag)), NA, sum(undiag, na.rm=TRUE)),
                sum_comm_undiag_weighted = ifelse(all(is.na(undiag_weighted)), NA, sum(undiag_weighted, na.rm=TRUE)),
                sum_comm_undiag_commweighted = ifelse(all(is.na(undiag_commweighted)), NA, sum(undiag_commweighted, na.rm=TRUE))) %>% 
      mutate(wave = wave_lim)
    
    res_list[[i]] <- comm_temp
    
  }
  
  res_df <- bind_rows(res_list) %>% 
    pivot_longer(c(-from, -wave)) %>% 
    mutate(variable = paste(name, wave, sep = "_"),
           variable = str_remove(variable, "_NA")) %>%
    select(-wave, -name) %>%
    pivot_wider(names_from = "variable", values_from = "value")
  
  res_df <- res_df %>% rename_with(.fn = ~ paste0(comm_type, .x), .cols = -from)
  
  return(res_df)
  
}

# Bipartite networks #

bi_calc <- function(wave_select, input_df, text_input) {
  
  print(paste0("Working through wave ", wave_select, " of the ", text_input, " foci network."))
  
  ## Setup ##
  
  course_size <- input_df %>%  # get foci size and weight by inverse of size
    filter(wave == wave_select) %>% 
    count(eval(parse(text = text_input))) %>% drop_na() %>%
    mutate(n = sqrt(1 + (1/n))) %>% # square root of inverse
    rename(!!text_input := 1)
  
  node_vars <- long_var_df %>% # get full nodelist
    filter(wave == wave_select & (variable == "undiag")) %>% 
    select(-wave) %>% 
    pivot_wider(values_from = value, names_from = variable)
  
  ## Person-person edgelist ##
  
  input_df1 <- input_df %>% # create person-course edgelist
    drop_na() %>%
    filter(wave == wave_select) %>% 
    select(egoid, !!text_input) %>% 
    distinct() %>% 
    left_join(course_size, by = text_input)
  
  countfoci <- input_df1 %>% count(egoid) %>% rename(!!(paste0(text_input, "_num")) := n)
  
  input_df2 <- input_df1 %>% select(c(1,2)) %>% rename(alterid = egoid) # reverse 
  
  person_mode <- input_df1 %>% 
    left_join(input_df2, relationship = "many-to-many") %>% # join both and create indices
    filter(egoid != alterid) %>% # remove self-loops
    mutate(id1 = case_when(egoid <= alterid ~ egoid, TRUE ~ alterid),
           id2 = case_when(egoid <= alterid ~ alterid, TRUE ~ egoid)) %>% 
    group_by(egoid, alterid) %>%
    mutate(directed_id = cur_group_id()) %>% # Directed ids
    ungroup() %>%
    group_by(id1, id2) %>%
    mutate(undirected_id = cur_group_id()) %>% # Undirected ids
    ungroup() %>%
    group_by_at(c(5, 6, 2)) %>% slice(1) %>% # remove duplicates 
    group_by(undirected_id, egoid = id1, alterid = id2) %>% 
    summarise(weight = sum(n), .groups = 'drop') %>% ungroup() %>% select(-undirected_id) # Sum weights of shared foci
  
  ## Alter characteristics ##
  
  res <- person_mode %>% left_join(node_vars, by = c("alterid" = "egoid")) %>% # create weighted characteristics
    mutate(undiag_weighted = undiag*weight) 
  
  res1 <- res %>% group_by(egoid) %>% bi_connect(., text = text_input) # Create summarizations for egoid
  
  res <- person_mode %>% left_join(node_vars, by = c("egoid" = "egoid")) %>% 
    mutate(undiag_weighted = undiag*weight) 
  
  res2 <- res %>% group_by(alterid) %>% bi_connect(., text = text_input) %>% rename(egoid = alterid) # Create summarizations for alterid
  
  res <- res1 %>% rbind(res2) %>% distinct(egoid, .keep_all = T) %>% left_join(countfoci, by = "egoid") # Merge and keep only uniques
  
  person_net <- graph_from_data_frame(person_mode, directed = F)
  
  return(list(person_net, res))
}

bi_connect <- function(edges, text){
  
  site <- text
  
  edges %>% 
    summarise(!!paste0("undiag_", site) := ifelse(all(is.na(undiag)), NA, sum(undiag, na.rm=TRUE)), # Non - imputed
              !!paste0("undiag_weighted_", site) := ifelse(all(is.na(undiag_weighted)), NA, sum(undiag_weighted, na.rm=TRUE)))
}

bi_transform <- function(x, nets = NULL) {
  net <- nets[[x]]
  sums <- res_list[[x]]
  wav <- waves[[x]]
  
  wave_df <- as_data_frame(net, "vertices") %>% 
    rename(egoid = name) %>% as_tibble() %>% 
    mutate(egoid = as.numeric(egoid)) %>% 
    remove_rownames() %>% left_join(sums, by = c("egoid")) %>% 
    mutate(across(where(is.character), as.numeric))
  
  wave_df <- wave_df %>% pivot_longer(c(-egoid)) %>% mutate(wave = wav) %>% 
    select(egoid, wave, variable = name, value)
  
  return(wave_df)
}

# Friendship nomination network #

build_basic_rels <- function(close) {
  
  included <- close %>%
    rename_at(., vars(contains(".x")), list(~ str_replace(.,fixed(".x"), "_ego"))) %>%
    rename_at(., vars(contains(".y")), list(~ str_replace(.,fixed(".y"), "_alter"))) %>%
    filter(in_surv_wave_ego == 1)
  
  # Aggregate of alter characteristics.
  outsurv <- included %>% # number of alters who are also in the survey
    group_by(egoid) %>% 
    filter(in_surv_wave_alter == 1) %>%  
    summarise(num_alter_in = n())
    
  ego_connect <- included %>% 
    group_by(egoid) %>% 
    summarise(num_alter = n(),
              sum_alter_undiag = ifelse(all(is.na(undiag_alter)), NA, sum(undiag_alter, na.rm=TRUE))) %>% 
    left_join(outsurv, by = "egoid") %>% 
    mutate(num_alter_in = replace_na(num_alter_in, 0))
  
  return(ego_connect)
  
}

build_net_vars <- function(x) {
  
  print(x)
  wave_lim <- str_remove(x, "[a-zA-Z]*") %>% as.numeric()
  
  # Build nodelist
  sources <- net_df %>% filter(wave == x) %>% distinct(egoid) %>% rename(label = egoid)
  if (wave_lim == 0) {sources <- net_df %>% filter(wave == "Wave1" | wave == "Wave2") %>% distinct(egoid) %>% rename(label = egoid)}
  destinations <- net_df %>% filter(wave == x) %>% distinct(alterid) %>% rename(label = alterid)
  if (wave_lim == 0) {destinations <- net_df %>% filter(wave == "Wave1" | wave == "Wave2") %>% distinct(alterid) %>% rename(label = alterid)}
  nodes <- full_join(sources, destinations, by = "label") %>% rowid_to_column("id") 
  is_surv <- long_var_df %>% filter(variable == "in_survey_wave" & wave == wave_lim & value == 1) %>% pull(egoid)
  
  # Populate node attributes
  attr_lim <- c("undiag")
  
  in_wave <- sources %>% pull()
  attrs <- long_var_df %>% 
    filter(wave %in% c(wave_lim, NA) & variable %in% attr_lim) %>% select(-wave) %>% 
    pivot_wider(names_from = variable, values_from = value) %>% 
    mutate_at(all_of(attr_lim), as.numeric)
  
  nodes <- nodes %>% 
    left_join(attrs, by = c("label" = "egoid")) %>% 
    mutate(in_net_wave = case_when(label %in% in_wave ~ 1, TRUE ~ 0),
           in_surv_wave = case_when(label %in% is_surv ~ 1, TRUE ~ 0))
  
  # Create edgelist (directed)
  close <- net_df %>% filter(wave == x)
  if (wave_lim == 0){close <- net_df %>% filter(wave == "Wave1" | wave == "Wave2")}
  
  close <- close  %>% select(egoid, alterid) %>%
    left_join(nodes, by = c("egoid" = "label")) %>% rename(from = id) %>% 
    left_join(nodes, by = c("alterid" = "label")) %>% rename(to = id) %>% 
    distinct(from, to, .keep_all = T)
  
  edges <- close %>% select(from, to, egoid, alterid)
  ego_connect <- build_basic_rels(close) # Create basic relational measurements
  
  # Create edgelist (undirected)
  edges1 <- edges %>% 
    mutate(id1 = case_when(egoid <= alterid ~ egoid, TRUE ~ alterid),
           id2 = case_when(egoid <= alterid ~ alterid, TRUE ~ egoid))
  
  edges2 <- edges %>% 
    mutate(id1 = case_when(egoid >= alterid ~ egoid, TRUE ~ alterid),
           id2 = case_when(egoid >= alterid ~ alterid, TRUE ~ egoid))
  
  close2 <- edges1 %>% 
    bind_rows(edges2) %>% 
    select(egoid = id1, alterid = id2) %>% 
    left_join(nodes, by = c("egoid" = "label")) %>% rename(from = id) %>% 
    left_join(nodes, by = c("alterid" = "label")) %>% rename(to = id) %>% 
    distinct(from, to, .keep_all = T)
  
  ego_connect2 <- build_basic_rels(close2)
  ego_connect2 <- ego_connect2 %>% 
    rename_with(~ paste0("undir_", .), -egoid) %>% 
    left_join(ego_connect, by = "egoid")
  
  # Create network and calculate network measurements on full network
  wave_network <- graph_from_data_frame(edges, nodes, directed = T)
  V(wave_network)$wave <- wave_lim # adding wave indicator
  wave_network <- wave_network %>% 
    set_vertex_attr(name = "wave", value = wave_lim)
  
  wave_df <- as_data_frame(wave_network, "vertices") %>% 
    select(label, in_net_wave, in_surv_wave, wave) %>%
    remove_rownames() 
  
  # Cleanup and turn into long dataframe (individuals)
  wave_df <- wave_df %>% left_join(ego_connect2, by = c("label" = "egoid"))
  wave_df <- wave_df %>% 
    filter(in_surv_wave == 1) %>% 
    select(-in_surv_wave) %>% 
    pivot_longer(c(-label, -wave)) %>% 
    rename(variable = name, egoid = label) 
  
  return(wave_df)
  
  # Cleanup and turn into dataframe (dyads)
  # full_dyad_df <<- bind_rows(full_dyad_df, wave_df_dyad)
  
}
