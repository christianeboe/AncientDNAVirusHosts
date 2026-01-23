###########################################
## Author: Christiane Böckel
## Date: 03.06.25
## Damage pattern analysis for DNA virus–host patterns in lake and marine environments over the last glacial cycle
###########################################

library(here)
library(tidyverse)

lineage <- read.delim(here("data/Siberian_rankedlineage.txt"), sep = "\t", header = F) %>%
  select(1,3,5,7,9,11,13,15,17,19,21) %>%
  setNames(c("taxid", "tax_name", "species", "genus","family", "order", "class", "phylum", "kingdom", "clade", "superkingdom")) %>%
  mutate(across(where(is.character), ~na_if(., "")))

meta <- read.delim(here("data/metadata_all.tsv"), sep = "\t", header = T) %>%
  rename(sample_name = sample,
         age = years) %>%
  mutate(epocheII = case_when(core == "PS97" & ka >= 12.7 ~ "Pleistocene",
                              ka >= 14 | age >= 14000 ~ "Pleistocene",
                              ka < 0 | age < 0 ~ "blank",
                              T ~ "Holocene")) 

#========================combined plot========================================
tax_lev <- "Viruses_superkingdom" 
for(tax_lev in c("Bacteria_superkingdom", "Viridiplantae_kingdom", "Metazoa_kingdom", "Viruses_superkingdom")) {   
  
  taxa <- str_split(tax_lev, "_")[[1]][1]
  taxlevel <- str_split(tax_lev, "_")[[1]][2]
  
  print(taxa)
  
  #-------------------------------------------------------------------------
  all_lakes_data <- list()
  
  for (lake in order_lakes) {
    print(lake)
  if (lake == "lama") {
    load(file = here("data/pydamage_all_result_taxon_wide_age_numeric.RData"))
    lake_data <- pydamage_all_numeric %>%
      mutate(age = as.numeric(age)) %>%
      mutate(core = lake) %>%
      filter(get(taxlevel) == taxa) %>%
      filter(!is.na(sample_name.y)) %>%
      rename(sample_name = sample_name.x) %>%
      mutate(reflen = as.integer(reflen)) %>%
      select(sample_name, predicted_accuracy, reflen, superkingdom, kingdom, class, matches("^CtoT-[0-9]$"), nb_reads_aligned) %>%
      left_join(meta, by = "sample_name") %>%
      select(core, ka, age, ecosystem, epocheII, sample_name, predicted_accuracy, reflen, superkingdom, kingdom, class, matches("^CtoT-[0-9]$"), nb_reads_aligned)
    rm(pydamage_all_numeric)
  } else {
    
    pydamage <- read.delim(list.files(here("data", lake), pattern = "_all_name_added_pydamage_result.csv", recursive = F, full.names = T), sep = ",", header = T)
    pydamage <- pydamage %>%
      setNames(sapply(colnames(pydamage), function(x) {gsub("CtoT.","CtoT-", x)})) %>% 
      {if("sample" %in% colnames(pydamage)) {
        rename(., sample_name = sample) 
      } else { . } }
    
    lake_data <- pydamage %>%
      left_join(combined_kraken %>% filter(core == lake) %>% select(-core), by = c("sample_name", "reference")) %>%
      left_join(lineage, by = "taxid") %>% 
      filter(get(taxlevel) == taxa) %>%
      left_join(meta, by = "sample_name") %>%
      select(core, ka, age, ecosystem, epocheII, sample_name, predicted_accuracy, reflen, superkingdom, kingdom, class, matches("^CtoT-[0-9]$"), nb_reads_aligned)
    rm(pydamage)
  }
    
    all_lakes_data[[lake]] <- lake_data
    rm(lake_data)
  }
  
  combined_taxon_data <- bind_rows(all_lakes_data)

  accuracy_filter <- 0.6
  reflen_filter   <- 1000
  
  combined_taxon_data %>%
    filter(!is.na(core)) %>%
    filter(predicted_accuracy >= accuracy_filter & reflen >= reflen_filter) %>%
    pivot_longer(cols = starts_with("CtoT"), names_to = "position", values_to = "freq") %>%
    separate(position, into = c("pos_prefix", "pos_num"), sep = "-", convert = TRUE) %>% 
    { ggplot(., aes(x = pos_num, y = freq)) +
        stat_summary(
          aes(color = factor(core, levels = order_lakes, sapply(names_lakes, function(x) { paste0(x, " (", combined_taxon_data %>% 
                                                                                                      filter(!is.na(core)) %>%
                                                                                                      filter(predicted_accuracy >= accuracy_filter & reflen >= reflen_filter) %>%
                                                                                                      group_by(core) %>%
                                                                                                      summarize(n_reads = sum(nb_reads_aligned, na.rm = T)) %>%
                                                                                                      ungroup %>%
                                                                                                      mutate(core = factor(core, levels = order_lakes, labels = names_lakes)) %>%
                                                                                                      .[.$core == x, "n_reads"], ")")})), group = core),
          fun = mean,
          geom = "line",
          linewidth = 0.8
        ) +
        scale_x_continuous(
          limits = c(0, 9),
          breaks = 0:9,
          labels = 1:10
        ) +
        ylim(0, 0.2) +
        scale_color_manual(values = c('#ffe119', '#4363d8', '#f58231', '#136208', 
                                      '#000075', '#dcbeff', '#a9a9a9', '#e6194B')) +
        labs(
          x = "Position in read",
          y = "C to T substitution frequency"
        ) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank(),
              legend.title = element_blank(),
              legend.position = "bottom")
    }
  ggsave(here("combined_plots", paste0(taxa, "_CtoT_transitions.pdf")), dpi = 900, width = 20, height = 20, unit = "cm")
  rm(all_lakes_data)
}




