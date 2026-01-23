###########################################
## Author: Christiane Böckel
## Date: 06.01.26
## Data analysis for DNA virus–host patterns in lake and marine environments over the last glacial cycle
###########################################

# load libraries
library(here)
library(tidyverse)
library(ggbreak)
library(vegan)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(cowplot)
library(ggVennDiagram)
library(scales)
set.seed(2)

#======read in data=========
order_sites <- c("btoko", "ilirney", "lama", "lele", "ulu", "KL77", "KL12", "PS97")
names_sites <- c("Bolshoe Toko", "Ilirney", "Lama", "Levinson-Lessing", "Ulu", "KL77", "KL12", "PS97")

# read in sample metadata
meta <- read.delim(here("data/metadata.tsv"), sep = "\t", header = T,
                   col.names = c("core", "core_full", "sample", "years", "ka", "total_reads", "mapped_reads", "perc_mapped_reads", "latitude", "longitude", "age-depth_link", "ENA_ProjectID")) %>%
  mutate(core = as.character(factor(site, levels = names_sites, labels = order_sites))) %>% 
  select(core, sample, years, ka, total_reads, mapped_reads, latitude, longitude) %>%
  mutate(years = as.integer(case_when(years %in% c("EB", "LB") ~ "-9999", T ~ years)),
    ka = as.numeric(case_when(ka %in% c("EB", "LB") ~ "-9999", T ~ ka)),     
    epoche = case_when(
      years > 12700 & core == "PS97" ~ "Pleistocene",
      ka > 14 ~ "Pleistocene",
      ka < 0 ~ "blank",
      TRUE ~ "Holocene"),
    epocheII = case_when(
      ka > 14 ~ "Pleistocene",
      ka < 0 ~ "blank",
      TRUE ~ "Holocene"),
    ecosystem = case_when(
      core %in% c("KL12", "KL77", "PS97") ~ "marine",
      T ~ "lake"))
    
  
# Viruses
df_viruses <- read.delim(here("data", "Viruses.tab"), sep = "\t", header = T)

# Correlation host taxa
df_host_species_agg <- read.delim(here("data", "Host_species.tab"), sep = "\t", header = T) %>%
  rename(analysis_group = genus,
         taxon = species) %>%
  mutate(taxlevel = "species") %>%
  left_join(meta)

df_host_genera_agg <- read.delim(here("data", "Host_genera.tab"), sep = "\t", header = T) %>%
  rename(taxon = genus) %>%
  mutate(analysis_group = "host_genera",
         taxlevel = "genus") %>%
  left_join(meta)

df_host_families_agg <- read.delim(here("data", "Host_families.tab"), sep = "\t", header = T) %>%
  rename(taxon = family) %>%
  mutate(analysis_group = "host_families",
         taxlevel = "family") %>%
  left_join(meta)

df_host_phyla_agg <- read.delim(here("data", "Host_phyla.tab"), sep = "\t", header = T) %>%
  rename(taxon = phylum) %>%
  mutate(analysis_group = "host_phyla",
         taxlevel = "phylum") %>%
  left_join(meta)

df_host_superkingdom_agg <- read.delim(here("data", "Host_superkingdom.tab"), sep = "\t", header = T) %>%
  mutate(analysis_group = "host_superkingdoms",
         taxlevel = "superkingdom",
         taxon = superkingdom) %>%
  left_join(meta)

# Bacterial host classes
df_bacterial_hosts_agg <- read.csv(here("data ", "Bacteria_class.csv"), sep = "\t", header = T) %>%
  rename(taxon = class) %>%
  mutate(analysis_group = "bacterial_classes",
         taxlevel = "class") %>%
  left_join(meta)

# read in meta data table from ICTV
ictv_df <- read.csv(here("data/VMR_MSL38_v2.csv"), sep = ";", header = T) %>% 
  dplyr::select(Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species, `Virus.name.s.`) %>%
  rename(Virus = Virus.name.s.)

#=========================define variables===================================
marine <- c("KL12", "PS97", "KL77")

vir_levels <- c("clade", "kingdom", "phylum", "class", "order", "family", "genus", "species")
bac_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

customPalette <- c('#4363d8', '#ffe119', '#800000', '#f58231', '#dcbeff', '#000075',
                   '#a9a9a9', '#e6194B', '#000000', '#ffffff')

bacteriaPalette <- c("Actinomycetes" = '#ffe119', "Alphaproteobacteria" = '#4363d8', "Bacilli" = '#f58231',
                     "Betaproteobacteria" = '#800000', "Cyanophyceae" = '#000075', "Flavobacteriia" = '#dcbeff',
                     "other Bacteria" = '#a9a9a9', "unknown host" = '#9299A0', "Clostridia" = '#000000', "Gammaproteobacteria" = '#e6194B', "Bacteroidia" = '#136208')

host_groupPalette_lakes <- c("Algae" = '#ffe119', "Protists" = '#4363d8', "Plants" = '#f58231',
                             "Cyanobacteria" = '#800000', "Metazoa" = '#000075', "Bacteria/Archaea" = '#dcbeff', "Bacteria" = '#e6194B',
                             "Unknown" = '#a9a9a9', "Invertebrates" = '#000000', "Mimiviridae" = '#136208')

host_groupPalette_marine <- c("Pelagibacter" = '#ffe119', "Protists" = '#4363d8', "Synechococcus" = '#f58231',
                              "Haptophyta" = '#800000', "Prasinophytae" = '#000075', "Bacteria/Archaea" = '#dcbeff', 
                              "Unknown" = '#a9a9a9')

shapes_cores_epoch <- c("H_btoko" = 0, "P_btoko" = 22, "H_ilirney" = 2, "P_ilirney" = 24,
                        "H_lama" = 5, "P_lama" = 23, "H_lele" = 1, "P_lele" = 21, "H_ulu" = 6, "P_ulu" = 25,
                        "H_KL77" = 0, "P_KL77" = 22, "H_KL12" = 2, "P_KL12" = 24, "H_PS97" = 1, "P_PS97" = 21)

envPalette <- c("marine" = "#000075", "lake" = "#FE8939", "both" = "#8B8B8B")

#==========================define functions==================================
# function to fill NA values for a specific taxonomic level
fill_na_with_lower_taxon <- function(df, level) {
  # find the index of the current level in the hierarchy
  level_index <- match(level, c("superkingdom", vir_levels))

  # select higher taxonomic levels based on the current level
  higher_levels <- c("superkingdom", vir_levels)[1:level_index-1]

  # dynamically create "combined" column for the current level
  df <- df %>%
    tidyr::unite("combined_path", all_of(higher_levels), sep = ";", remove = FALSE, na.rm = TRUE) %>%
    rowwise() %>%
    mutate(low_tax = str_split(combined_path, pattern = ";")[[1]] %>% last,
           !!level := if_else(is.na(get(level)), paste("unassigned", low_tax, level), get(level))) %>%
    ungroup() %>%
    dplyr::select(-combined_path, -low_tax)

  return(df)
}


compute_similarity <- function(core1, core2, bin_width, min_year, max_year) {
  
  # bin + average rel_abundance within each bin (per core × species)
  df <- df_rel_abundance %>%
    filter(
      core %in% c(core1, core2),
      years >= min_year, years <= max_year,
      superkingdom == "Viruses",
      taxlevel == "species",
      taxon != "uncultured Caudovirales phage",
      taxon != "Virus NIOZ-UU157"
    ) %>%
    mutate(time_bin = floor(years / bin_width) * bin_width) %>%
    group_by(core, time_bin, taxon) %>%
    summarise(rel_abundance = mean(rel_abundance, na.rm = TRUE)) %>%
    ungroup
  
  # both cores in at least one bin
  shared_bins <- df %>%
    distinct(core, time_bin) %>%
    count(time_bin) %>%
    filter(n == 2) %>%
    pull(time_bin)
  
  if (length(shared_bins) == 0) return(NA_real_)
  
  # Hellinger distance per shared bin
  dists <- map_dbl(shared_bins, function(b) {
    wide <- df %>%
      filter(time_bin == b) %>%
      pivot_wider(names_from = taxon, values_from = rel_abundance, values_fill = 0) %>%
      arrange(core)
    
    wide %>%
      select(-core, -time_bin) %>%
      as.matrix() %>%
      decostand(method = "hellinger") %>%
      vegdist(method = "euclidean") %>%
      as.numeric
    
  })
  
  # mean Hellinger distance per core
  mean(dists)
}



hostPalette <- function(n) {
  colorRampPalette(c("#000075", "#D6D6FF"))(n)
}

virusPalette <- function(n) {
  colorRampPalette(c("#833500", "#FFCEAD"))(n)
}


#=============replace NAs with info about higher taxonomic levels===========
viruses_unassigned <- df_viruses
for (level in rev(vir_levels)) {
  viruses_unassigned <- fill_na_with_lower_taxon(viruses_unassigned, level)
}

#========================normalize data=====================================
# TSS normalization
lca_norm <- viruses_unassigned %>%
  left_join(meta, by = c("sample", "years", "ka", "core")) %>%
  group_by(sample) %>%
  mutate(percentage = taxonReads / mapped_reads * 100) %>% # adjust for differences in sequencing depth
  mutate(total_sample = sum(taxonReads)) %>%
  ungroup()

#========================generate table with cumulative reads==================
# cumulative per virus taxon
df_rel_abundance_vir <- lapply(c("superkingdom", vir_levels), function(taxlevel) {
  lca_norm %>%
    filter(years >= 0,
           superkingdom %in% c("Viruses", "Bacteria", "Eukaryota", "Archaea")) %>%
    group_by(core, sample, years, ka, superkingdom, get(taxlevel)) %>%
    reframe(cum_taxReads = sum(taxonReads), # reads per taxon per sample
              cum_taxpercent = sum(percentage), # percentage of reads per taxon per sample
              taxlevel = taxlevel) %>%
    ungroup() %>%
    dplyr::select(`get(taxlevel)`, superkingdom, taxlevel, cum_taxReads, cum_taxpercent, core, sample, years, ka) %>%
    rename(taxon = `get(taxlevel)`) %>%
    group_by(taxon) %>%
    mutate(analysis_group = "Viruses") %>%
    ungroup()
}) %>% do.call(rbind, .) %>%
  group_by(core, sample, ka, years, superkingdom, taxlevel) %>%
  mutate(taxon = case_when(
    is.na(taxon) ~ "unassigned",
    TRUE ~ taxon
  )) %>%
  ungroup()
  
# combine all taxa in one table
combined_taxa <- df_bacterial_hosts_agg %>%
  bind_rows(df_host_families_agg, df_host_genera_agg, df_host_phyla_agg, df_host_species_agg, df_host_superkingdom_agg) %>% 
  rename(cum_taxReads = taxonReads) %>% 
  mutate(cum_taxpercent = cum_taxReads / mapped_reads * 100) %>%
  group_by(core, sample, ka, years, superkingdom, taxlevel) %>%
  mutate(taxon = case_when(
    is.na(taxon) ~ "unassigned",
    TRUE ~ taxon
  )) %>%
  ungroup() %>%
  select(taxon, superkingdom, taxlevel, cum_taxReads, cum_taxpercent, core, sample, years, ka, analysis_group)

df_superkingdoms <- combined_taxa %>%
  filter(taxlevel == "superkingdom") %>%
  select(sample, superkingdom, cum_taxReads) %>%
  rename(superkingdom_taxReads = cum_taxReads)

df_rel_abundance <- combined_taxa %>%
  bind_rows(df_rel_abundance_vir) %>%
  left_join(df_superkingdoms, by = c("superkingdom", "sample")) %>%
  mutate(rel_abundance = (cum_taxReads / superkingdom_taxReads) * 100) %>%
  select(-superkingdom_taxReads) %>%
  filter(years >= 0 | ka >= 0) %>%
  mutate(core = factor(core, levels = order_sites)) %>%
  arrange(taxlevel, superkingdom, core, sample) %>%
  mutate(core = as.character(core))

# additionally filtered
df_filtered <- df_rel_abundance %>%
  group_by(taxon, core) %>%
  mutate(max_reads = max(cum_taxReads, rm.na = T)) %>%
  ungroup %>%
  filter(max_reads >= 5)

# blank table
df_blanks <- lapply(vir_levels, function(taxlevel) {
  lca_norm %>%
    filter(years < 0 | ka < 0 | sample == "200313_A00902_A_L001_APMG-8-7") %>%
    group_by(core, sample, years, ka, superkingdom, get(taxlevel)) %>%
    reframe(cum_taxReads = sum(taxonReads), # reads per taxon per sample
              cum_taxpercent = sum(percentage), # percentage of reads per taxon per sample
              taxlevel = taxlevel) %>%
    ungroup() %>%
    dplyr::select(`get(taxlevel)`, superkingdom, taxlevel, cum_taxReads, cum_taxpercent, core, sample, years, ka) %>%
    rename(taxon = `get(taxlevel)`) %>%
    group_by(taxon) %>%
    mutate(all_taxReads = sum(cum_taxReads)) %>% # reads per taxon across all samples
    ungroup()
}) %>% do.call(rbind, .) %>%
  group_by(core, sample, ka, years, superkingdom, taxlevel) %>%
  mutate(rel_abundance = (cum_taxpercent / sum(cum_taxpercent)) * 100) %>%
  mutate(taxon = case_when(
    is.na(taxon) ~ "unassigned",
    TRUE ~ taxon
  )) %>%
  ungroup()

#=================stats for data overview======
# number of read in dataset
meta %>% 
  pull(total_reads) %>%
  sum

# number of samples
meta %>%
  filter(years >= 0 | ka >= 0) %>%
  pull(sample) %>%
  unique %>%
  length

# number of blanks
meta %>%
  filter(!(years >= 0 | ka >= 0)) %>%
  pull(sample) %>%
  unique %>%
  length

# sum of mapped reads in dataset
meta %>%
  filter(years >= 0 | ka >= 0) %>%
  pull(mapped_reads) %>%
  sum

# sum of viral reads in dataset
df_viruses %>%
  filter(years >= 0 | ka >= 0) %>%
  pull(taxonReads) %>%
  sum

# percentage of viral reads in dataset
100 * (df_viruses %>% 
         filter(years >=0 | ka >= 0) %>% 
         pull(taxonReads) %>% 
         sum) / (meta %>%
                   filter(years >= 0 | ka >= 0) %>%
                   pull(mapped_reads) %>%
                   sum)


# number viral of taxa in dataset
df_viruses %>%
  filter(years >= 0 | ka >= 0) %>%
  pull(taxon) %>%
  unique %>%
  length

# median viral taxa per sample
df_viruses %>%
  filter(years >= 0 | ka >= 0) %>%
  group_by(sample) %>%
  reframe(n_taxa = n_distinct(taxon)) %>%
  ungroup %>%
  pull(n_taxa) %>%
  median

# median viral reads per sample
df_viruses %>%
  filter(years >= 0 | ka >= 0) %>%
  group_by(sample) %>%
  reframe(cum_reads = sum(taxonReads)) %>%
  ungroup %>%
  pull(cum_reads) %>%
  median

# rel read contribution per taxonomic level
df_viruses %>%
  filter(years >= 0 | ka >= 0) %>%
  group_by(rank) %>%
  reframe(cum_reads = sum(taxonReads)) %>%
  ungroup %>%
  mutate(rel_reads = 100 * cum_reads / sum(cum_reads)) %>%
  arrange(desc(rel_reads))

# taxonomic level of identified viral taxa
df_viruses %>%
  filter(years >= 0 | ka >= 0) %>%
  group_by(rank) %>%
  reframe(n_taxa = n_distinct(taxon)) %>%
  ungroup %>%
  mutate(rel_taxa = 100 * n_taxa / sum(n_taxa)) %>%
  arrange(desc(rel_taxa))

# percent reads per group in the Venn diagram
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(ecosystem = case_when(n() == 2 ~ "both",
                                T ~ ecosystem)) %>%
  ungroup %>% 
  unique %>%
  group_by(ecosystem) %>%
  mutate(n = n()) %>%
  ungroup %>%
  left_join(df_rel_abundance) %>%
  group_by(ecosystem) %>%
  reframe(n = mean(n),
          sum_reads = sum(cum_taxReads)) %>%
  ungroup %>%
  mutate(cum_percent = 100 * sum_reads / sum(sum_reads),
         n_percent = 100 * n / sum(n)) 

# Cirlivirales difference between the ecosystems
df_rel_abundance %>%
  filter(taxon == "Cirlivirales") %>%
  select(sample, rel_abundance) %>%
  left_join(meta) %>%
  group_by(ecosystem) %>%
  reframe(median(rel_abundance))

# number samples with more than 50% uncultured Caudovirales phage
df_rel_abundance %>%
  filter(taxon == "uncultured Caudovirales phage") %>%
  select(sample, rel_abundance) %>%
  merge(meta, by = "sample", all = T) %>%
  replace(is.na(.), 0) %>%
  filter(years >= 0 | ka >= 0,
         ecosystem == "lake") %>%
  mutate(phage50 = case_when(rel_abundance > 50 ~ "more than 50%",
                             T ~ "equal or less to 50%")) %>%
  group_by(phage50) %>%
  summarize(n_samples = n()) %>%
  ungroup %>%
  pivot_wider(values_from = n_samples, names_from = phage50) %>%
  mutate(perc_samples = 100 * `more than 50%` / (`more than 50%` + `equal or less to 50%`))

# number samples with more than 90% uncultured Caudovirales phage
df_rel_abundance %>%
  filter(taxon == "uncultured Caudovirales phage") %>%
  select(sample, rel_abundance) %>%
  merge(meta, by = "sample", all = T) %>%
  replace(is.na(.), 0) %>%
  filter(years >= 0 | ka >= 0) %>%
  mutate(phage90 = case_when(rel_abundance > 90 ~ "more than 90%",
                             T ~ "equal or less to 90%")) %>%
  group_by(ecosystem, phage90) %>%
  summarize(n_samples = n()) %>%
  ungroup %>%
  pivot_wider(id_cols = "ecosystem", values_from = n_samples, names_from = phage90) %>%
  mutate(perc_samples = 100 * `more than 90%` / (`more than 90%` + `equal or less to 90%`))

# max relative abundance of uncultured Caudovirales phage
df_rel_abundance %>%
  filter(taxon == "uncultured Caudovirales phage") %>%
  select(sample, rel_abundance) %>%
  merge(meta, by = "sample", all = T) %>%
  replace(is.na(.), 0) %>%
  filter(years >= 0 | ka >= 0) %>%
  group_by(ecosystem) %>%
  reframe(minimum = min(rel_abundance),
            maximum = max(rel_abundance, na.rm = T),
            average = mean(rel_abundance, na.rm = T))

# median relative abundance Pithovirus LCPA304 in PS97
df_rel_abundance %>%
  filter(taxon == "Pithovirus LCPAC304",
         core == "PS97") %>%
  select(sample, rel_abundance) %>%
  merge(meta, by = "sample") %>%
  filter(years >= 0 | ka >= 0) %>%
  group_by(epoche) %>%
  reframe(mean(rel_abundance))

# abundance of Mimiviridae
df_rel_abundance %>%
  filter(taxon == "Mimiviridae") %>%
  select(sample, rel_abundance) %>%
  merge(meta, by = "sample") %>%
  filter(years >= 0 | ka >= 0,
         ecosystem == "lake") %>%
  group_by(core) %>%
  mutate(max_alltime = max(rel_abundance)) %>%
  ungroup %>%
  group_by(core, epoche, max_alltime) %>%
  reframe(mean = mean(rel_abundance)) %>% 
  pivot_wider(id_cols = c(core, max_alltime), names_from = epoche, values_from = mean) %>%
  mutate(fold_change = Holocene / Pleistocene)

# abundance of Casjensviridae
df_rel_abundance %>%
  filter(taxon == "Casjensviridae") %>%
  select(sample, rel_abundance) %>%
  merge(meta, by = "sample") %>%
  filter(years >= 0 | ka >= 0,
         ecosystem == "lake") %>%
  group_by(core) %>%
  reframe(round(max(rel_abundance), 1))

# max abundance of Synechococcus strains in KL77
df_rel_abundance %>%
  filter(taxon %in% c("Synechococcus sp. Ace-Pa", "Synechococcus sp. UW69"),
         core == "KL77") %>%
  group_by(taxon) %>%
  reframe(max = max(rel_abundance))

# abundance of Cyanophyceae
df_rel_abundance %>% 
  filter(taxon == "Cyanophyceae") %>%
  left_join(meta) %>%
  group_by(ecosystem) %>%
  reframe(median(rel_abundance))

#==================ecosystem difference=======================================
# percentage of unique taxa
df_rel_abundance %>% 
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>% 
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(ecosystem = case_when(n() == 2 ~ "both",
                                  T ~ ecosystem)) %>%
  ungroup %>% 
  unique %>%
  group_by(ecosystem) %>%
  mutate(n = n()) %>%
  ungroup %>%
  left_join(df_rel_abundance) %>%
  group_by(ecosystem) %>%
  reframe(n = mean(n),
            sum_reads = sum(cum_taxReads)) %>%
  ungroup %>%
  mutate(cum_percent = 100 * sum_reads / sum(sum_reads),
         n_percent = 100 * n / sum(n)) %>% 
  ggplot(aes(x = "", y = n, fill = ecosystem)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = envPalette)
ggsave(here("Figures", "all_cores", "shared_taxa", "Piechart_unique_taxa.pdf"), dpi = 900) # Figure 1C

#top 5 species - unique to environment
# marine
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(presence = case_when(all(ecosystem == "lake") ~ "lake",
                               all(ecosystem == "marine") ~ "marine",
                               TRUE ~ "both")) %>%
  filter(presence == "marine") %>%
  left_join(df_rel_abundance %>% filter(taxlevel == "species")) %>%
  left_join(meta) %>%
  dplyr::select(presence, taxon, sample, ecosystem, rel_abundance) %>%
  pivot_wider(id_cols = c("presence", "taxon", "ecosystem"), names_from = "sample", values_from = "rel_abundance", values_fill = 0) %>% 
  pivot_longer(cols = -c("presence", "taxon", "ecosystem"), names_to = "sample", values_to = "rel_abundance") %>%
  group_by(presence, taxon, ecosystem) %>%
  reframe(mean_abundance = mean(rel_abundance)) %>%
  ungroup %>%
  group_by(ecosystem) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:5) %>%
  ungroup %>% View

# lake
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(presence = case_when(all(ecosystem == "lake") ~ "lake",
                               all(ecosystem == "marine") ~ "marine",
                               TRUE ~ "both")) %>%
  filter(presence == "lake") %>%
  left_join(df_rel_abundance %>% filter(taxlevel == "species")) %>%
  left_join(meta) %>%
  filter(ecosystem == "lake") %>%
  dplyr::select(presence, taxon, sample, ecosystem, rel_abundance) %>%
  pivot_wider(id_cols = c("presence", "taxon", "ecosystem"), names_from = "sample", values_from = "rel_abundance", values_fill = 0) %>% 
  pivot_longer(cols = -c("presence", "taxon", "ecosystem"), names_to = "sample", values_to = "rel_abundance") %>%
  group_by(presence, taxon, ecosystem) %>%
  reframe(mean_abundance = mean(rel_abundance)) %>%
  ungroup %>%
  group_by(ecosystem) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:5) %>%
  ungroup %>% View

# abundance of unique taxa
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(ecosystem = case_when(n() == 2 ~ "both",
                                  T ~ ecosystem)) %>%
  ungroup %>% 
  unique %>%
  group_by(ecosystem) %>%
  mutate(n = n()) %>%
  ungroup %>%
  left_join(df_rel_abundance) %>%
  group_by(ecosystem) %>%
  reframe(n = mean(n),
            sum_reads = sum(cum_taxReads)) %>%
  ungroup %>%
  mutate(cum_percent = 100 * sum_reads / sum(sum_reads)) %>%
  ggplot(aes(x = "", y = cum_percent, fill = ecosystem)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = envPalette)
ggsave(here("Figures", "all_cores", "shared_taxa", "Piechart_abundance_unique_taxa.pdf"), dpi = 900) # Figure 1D


# differential abundance of top shared taxa
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(presence = case_when(all(ecosystem == "lake") ~ "lake",
                                 all(ecosystem == "marine") ~ "marine",
                                 TRUE ~ "both")) %>%
  filter(presence == "both") %>%
  left_join(df_rel_abundance %>% filter(taxlevel == "species")) %>%
  left_join(meta) %>%
  dplyr::select(presence, taxon, sample, ecosystem, rel_abundance) %>%
  pivot_wider(id_cols = c("presence", "taxon", "ecosystem"), names_from = "sample", values_from = "rel_abundance", values_fill = 0) %>% 
  pivot_longer(cols = -c("presence", "taxon", "ecosystem"), names_to = "sample", values_to = "rel_abundance") %>%
  group_by(presence, taxon, ecosystem) %>%
  reframe(mean_abundance = mean(rel_abundance)) %>%
  ungroup %>%
  group_by(ecosystem) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:7) %>%
  ungroup %>%
  dplyr::select(presence, taxon) %>%
  unique %>%
  left_join(df_rel_abundance %>% filter(taxlevel == "species")) %>% 
  dplyr::select(presence, taxon, sample, rel_abundance) %>%
  pivot_wider(id_cols = c("presence", "taxon",), names_from = "sample", values_from = "rel_abundance", values_fill = 0) %>% 
  pivot_longer(cols = -c("presence", "taxon"), names_to = "sample", values_to = "rel_abundance") %>%
  group_by(presence, taxon) %>%
  reframe(mean_abundance = mean(rel_abundance)) %>%
  ungroup %>% 
  group_by(presence) %>%
  mutate(norm_abundance = mean_abundance / max(mean_abundance)) %>%
  ungroup %>% 
  filter(presence == "both") %>%
  left_join(df_rel_abundance, by = "taxon") %>%
  left_join(meta, by = "sample") %>%
  group_by(ecosystem, taxon) %>%
  reframe(mean = mean(rel_abundance),
            sd = sd(rel_abundance)) %>%
  ungroup %>%
  arrange(rev(ecosystem), mean) %>%
  filter(ecosystem == "marine") %>%
  {. ->> df_marine_shared} %>% 
  ggplot(aes(y = factor(taxon, levels = unique(taxon)))) +
  geom_bar(aes(x = mean),
           stat = "identity") +
  geom_errorbar(aes(xmin = mean-0.05, xmax = mean + sd), width = 0.4) +
  theme_minimal() +
  #facet_wrap(~ecosystem, ncol = 2, scales = "free_x") +
  xlab("mean relative abundance [%]") +
  xlim(0,30) +
  scale_x_reverse() +
  ylab(NULL)
ggsave(here("Figures", "all_cores", "shared_taxa", "Abundance_shared_taxa_marine.pdf"), dpi = 900) # Figure 1E


df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(ecosystem, taxon) %>%
  unique %>%
  group_by(taxon) %>%
  reframe(presence = case_when(all(ecosystem == "lake") ~ "lake",
                                 all(ecosystem == "marine") ~ "marine",
                                 TRUE ~ "both")) %>%
  filter(presence == "both") %>%
  left_join(df_rel_abundance %>% filter(taxlevel == "species")) %>%
  left_join(meta) %>%
  dplyr::select(presence, taxon, sample, ecosystem, rel_abundance) %>%
  pivot_wider(id_cols = c("presence", "taxon", "ecosystem"), names_from = "sample", values_from = "rel_abundance", values_fill = 0) %>% 
  pivot_longer(cols = -c("presence", "taxon", "ecosystem"), names_to = "sample", values_to = "rel_abundance") %>%
  group_by(presence, taxon, ecosystem) %>%
  reframe(mean_abundance = mean(rel_abundance)) %>%
  ungroup %>%
  group_by(ecosystem) %>%
  arrange(desc(mean_abundance)) %>%
  slice(1:7) %>%
  ungroup %>%
  dplyr::select(presence, taxon) %>%
  unique %>%
  left_join(df_rel_abundance %>% filter(taxlevel == "species")) %>% 
  dplyr::select(presence, taxon, sample, rel_abundance) %>%
  pivot_wider(id_cols = c("presence", "taxon",), names_from = "sample", values_from = "rel_abundance", values_fill = 0) %>% 
  pivot_longer(cols = -c("presence", "taxon"), names_to = "sample", values_to = "rel_abundance") %>%
  group_by(presence, taxon) %>%
  reframe(mean_abundance = mean(rel_abundance)) %>%
  ungroup %>% 
  group_by(presence) %>%
  mutate(norm_abundance = mean_abundance / max(mean_abundance)) %>%
  ungroup %>% 
  filter(presence == "both") %>%
  left_join(df_rel_abundance, by = "taxon") %>%
  left_join(meta, by = "sample") %>%
  group_by(ecosystem, taxon) %>%
  reframe(mean = mean(rel_abundance),
            sd = sd(rel_abundance)) %>%
  ungroup %>%
  filter(ecosystem == "lake") %>%
  ggplot(aes(y = factor(taxon, levels = unique(df_marine_shared$taxon)))) +
  geom_bar(aes(x = mean),
           stat = "identity") +
  geom_errorbar(aes(xmin = mean-0.05, xmax = mean + sd), width = 0.4) +
  theme_minimal() +
  #facet_wrap(~ecosystem, ncol = 2, scales = "free_x") +
  xlab("mean relative abundance [%]") +
  scale_x_break(c(15, 70)) +
  ylab(NULL)
ggsave(here("Figures", "all_cores", "shared_taxa", "Abundance_shared_taxa_lake.pdf"), dpi = 900) # Figure 1E


#=======================NMDS plot all viruses===================================
nmds_meta <- meta %>%
  filter(years >= 0 | ka >= 0) %>%
  dplyr::select(sample, years, ecosystem, epocheII, epoche, core) %>%
  mutate(ecosystem = factor(ecosystem, levels = c("marine", "lake"))) %>%
  mutate(epoche_core = case_when(epocheII == "Holocene" ~ paste0("H_", core),
                                 epocheII == "Pleistocene" ~ paste0("P_", core))) %>%
  mutate(epoche_env = case_when(epocheII == "Holocene" ~ paste0("H_", ecosystem),
                                epocheII == "Pleistocene" ~ paste0("P_", ecosystem)))

# multiple lakes
for (eco in c("all", "lake", "marine")) {
  for (caudovirales in c(T, F)) {
        tmp_path <- here("Figures", "all_cores", "NMDS", eco)
        dir.create(tmp_path, recursive = T)
        set.seed(2)
        
        NMDS <- df_rel_abundance

        if(eco == "lake") {
          NMDS <- NMDS %>%
            filter(!core %in% marine)

        } else if (eco == "marine") {
          NMDS <- NMDS %>%
            filter(core %in% marine)
        }

        if(!caudovirales) {
          NMDS <- NMDS %>%
            filter(taxon != "uncultured Caudovirales phage")
        }

        nmds_meta_filtered <- nmds_meta %>%
          filter(sample %in% unique(NMDS$sample)) %>%
          column_to_rownames("sample")

        NMDS <- NMDS %>%
          filter(superkingdom == "Viruses") %>%
          filter(taxlevel == "species") %>%
          filter(ka >= 0) %>%
          dplyr::select(taxon, sample, rel_abundance) %>%
          pivot_wider(names_from = sample, values_from = rel_abundance) %>%
          column_to_rownames("taxon") %>%
          as.matrix %>%
          t %>%
          replace(is.na(.), 0) %>%
          decostand(method = "hellinger") %>%
          metaMDS(k = 2, trymax = 100, trace = F, autotransform = F, distance = "bray")

        nmds_samples <- NMDS %>%
          scores %>%
          .[["sites"]] %>%
          as.data.frame
        
        top_taxa <- df_filtered %>%
          left_join(meta) %>%
          filter(superkingdom == "Viruses") %>%
          {if (eco == "all") {
            .
          } else {
            filter(., ecosystem == eco)
          }} %>%
          {if (caudovirales) {
            .
          } else {
            filter(., taxon != "uncultured Caudovirales phage")
          }} %>%
          group_by(taxon, core, epocheII) %>%
          mutate(n_samples_present = n(),
                 mean_abund = mean(rel_abundance),
                 max_rel_abund = max(rel_abundance),
                 n_abundant = sum(rel_abundance >= (0.3 * max_rel_abund))) %>%
          ungroup() %>%
          filter(n_samples_present > 3, 
                 n_abundant > 3) %>%
          filter(taxlevel == level,
                 #!taxon %in% core_virome$taxon,
                 !grepl("unassigned", taxon)) %>%
          group_by(core, epocheII) %>%
          arrange(taxlevel, desc(mean_abund)) %>%
          slice(1:20) %>%
          ungroup %>%
          pull(taxon) %>%
          unique %>%
          c(., "Timirovirus borboradaptatum",
            "Deltanudivirus tipoleraceae",
            "Sphingomonas phage vB_StuS_MMDA13",
            "Pseudomonas phage SC_9_H1H8_2017",
            "Microviridae sp.",
            "Pandoravirus dulcis")

        nmds_taxa <- NMDS %>%
          scores %>%
          .[["species"]] %>%
          as.data.frame %>%
          rownames_to_column("taxon") %>%
          filter(taxon %in% top_taxa)


        plot_epoch <- nmds_samples %>%
          rownames_to_column("sample") %>%
          left_join(nmds_meta, by = "sample") %>%
          {. ->> NMDS_plot_data } %>%
          ggplot() +
          stat_ellipse(aes(x = NMDS1, y = NMDS2, linetype = factor(epocheII, levels = c("Pleistocene", "Holocene")), color = ecosystem), geom="polygon", level=0.95, alpha=0) +
          geom_point(aes(x = NMDS1, y = NMDS2, color = ecosystem, fill = ecosystem, shape = epoche_core), size=2, alpha = 0.5) +
          theme_minimal() +
          labs(color = "Site",
               shape = "Geological epoch") +
          scale_fill_manual(values = envPalette) +
          scale_color_manual(values = envPalette) +
          scale_shape_manual(values = shapes_cores_epoch) +
          geom_point(
            data = nmds_taxa,
            aes(NMDS1, NMDS2),
            inherit.aes = FALSE,
            shape = 4,
            size = 2.5,
            stroke = 0.9,
            color = "red"
          ) +
          geom_text_repel(data = nmds_taxa, aes(x = NMDS1, y = NMDS2, label = taxon),
                          size = 3,
                          min.segment.length = unit(0, 'lines'),
                          seed = 2) +
          theme(legend.position = "none")

        ggsave(here(tmp_path, paste0("NMDS_epochII_", eco, "_caudovirales_", caudovirales, "_top_taxa_repel.pdf")), plot_epoch, width = 10, height = 6, dpi = 900) # Figure 2A 
        
  }
}


#==============RDA================
rda_data <- df_rel_abundance %>%
  filter(superkingdom == "Viruses") %>%
  filter(taxlevel == "species") %>%
  filter(ka >= 0 | years >= 0) %>%
  dplyr::select(taxon, sample, rel_abundance) %>%
  pivot_wider(id_cols = sample, names_from = taxon, values_from = rel_abundance) %>%
  column_to_rownames("sample") %>%
  as.matrix %>%
  replace(is.na(.), 0) %>% 
  decostand(method = "hellinger")

rda_env <- meta %>%
  filter(sample %in% (df_rel_abundance %>%   filter(superkingdom == "Viruses",
                                                    taxlevel == "species") %>% pull(sample) %>% unique)) %>%
  dplyr::select(sample, ecosystem, epocheII) %>%
  column_to_rownames("sample")

rda <- rda(rda_data ~ ., data = rda_env)
summary(rda)
RsquareAdj(rda)
anova.cca(rda, permutations = 1000)

rda_unique_env <- rda(rda_data ~ ecosystem + Condition(epocheII), data = rda_env)
RsquareAdj(rda_unique_env)
anova.cca(rda_unique_env, permutations = 1000)

rda_unique_epoch <- rda(rda_data ~ epocheII + Condition(ecosystem), data = rda_env)
RsquareAdj(rda_unique_epoch)
anova.cca(rda_unique_epoch, permutations = 1000)


#========================virus composition==================================
df_rel_abundance %>%
  filter(taxlevel == "order" | taxon == "unassigned Viruses species",
         superkingdom == "Viruses") %>%
  dplyr::select(sample, taxon, rel_abundance) %>%
  pivot_wider(id_cols = "sample", names_from = "taxon", values_from = "rel_abundance") %>%
  replace(is.na(.), 0) %>% 
  mutate(`unassigned Viruses order` = `unassigned Viruses order` - `unassigned Viruses species`) %>%
  pivot_longer(cols = -"sample", names_to = "taxon", values_to = "rel_abundance") %>%
  left_join(meta) %>%
  dplyr::select(core, ecosystem, taxon, rel_abundance) %>%
  group_by(core, ecosystem, taxon) %>% 
  reframe(mean_percent = mean(rel_abundance)) %>%
  ungroup %>% 
  group_by(core, ecosystem) %>%
  mutate(rel_abundance = (mean_percent / sum(mean_percent)) * 100) %>%
  ungroup %>%
  group_by(taxon) %>%
  mutate(n_abundant = sum(rel_abundance > 3)) %>%
  ungroup %>%
  mutate(taxon = case_when(taxon %in% (df_viruses %>% filter(class == "Caudoviricetes") %>% pull(order) %>% unique) ~ "Caudoviricetes (class)",
                           taxon == "unassigned Caudoviricetes order" ~ "Caudoviricetes (class)",
                           taxon == "unassigned Viruses order" ~ "assigned virus sp. without order",
                           taxon == "unassigned Viruses species" ~ "unassigned Viruses",
                           n_abundant >= 1 ~ paste0(taxon, " (order)"),
                           TRUE ~ "other Viruses order")) %>%
  group_by(core, ecosystem, taxon) %>%
  reframe(rel_abundance = sum(rel_abundance)) %>%
  ungroup %>% 
  mutate(taxon = factor(taxon, levels = c("Algavirales (order)", "Cirlivirales (order)", "Imitervirales (order)", "other Viruses order", "Caudoviricetes (class)", "assigned virus sp. without order", "unassigned Viruses")),
         core = factor(core, levels = rev(order_sites), labels = rev(names_sites))) %>%
  ggplot(aes(y = core, x = rel_abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  xlab("Relative Abundance [%]") +
  ylab(NULL) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = customPalette)

ggsave(here("Figures", "all_cores", "overview_stats", "Virus_composition_order.pdf"), dpi = 900) # Figure 2B

#========================core map===============================
# Northern Hemisphere map
window_crs <- "+proj=robin +datum=WGS84"

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = window_crs)

siberia_crs <- "+proj=aeqd +lat_0=60 +lon_0=135 +x_0=0 +y"

siberia_crop_box <- st_as_sfc(st_bbox(c(
  xmin = -2500000,
  xmax =  3500000,
  ymin = -3000000,
  ymax =  3000000
), crs = siberia_crs))

world_siberia <- world %>%
  st_transform(crs = siberia_crs) %>%
  st_crop(siberia_crop_box)

core_siberia <- meta %>%
  filter(!is.na(core)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = siberia_crs) %>%
  st_crop(siberia_crop_box)

ggplot() +
  geom_sf(data = world_siberia, fill = "gray95", color = "gray60", size = 0.3) +
  geom_sf(data = core_siberia, aes(color = ecosystem, fill = ecosystem), size = 4, shape = 23) +
  coord_sf(
    crs = siberia_crs,
    xlim = c(-2500000, 3500000),
    ylim = c(-3000000, 3000000),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    legend.position = "none"
  ) +
  scale_fill_manual(values = envPalette) +
  scale_color_manual(values = envPalette)
ggsave(here("Figures", "all_cores", "map", paste0("map_Siberia.pdf")), dpi = 900) # Figure 1A


# Southern Hemisphere map
window_crs <- "+proj=robin +datum=WGS84"

world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = window_crs)

tilted_crs <- "+proj=aeqd +lat_0=-65 +lon_0=-55 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

crop_box <- st_as_sfc(st_bbox(c(
  xmin = -2000000,
  xmax =  2000000,
  ymin = -2000000,
  ymax =  2000000
), crs = tilted_crs))

core_points <- meta %>%
  filter(!is.na(core)) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = tilted_crs) %>%
  st_crop(crop_box)

world_tilted_clipped <- world %>%
  st_transform(crs = tilted_crs) %>%
  st_intersection(crop_box)

ggplot() +
  geom_sf(data = world_tilted_clipped, fill = "gray95" , color = "gray60", size = 0.3, alpha = 1) +
  geom_sf(data = core_points, color = "#000075", size = 4, shape = 23, fill = "#000075") +
  coord_sf(
    crs = tilted_crs,
    xlim = c(-2000000, 2000000),
    ylim = c(-2000000, 2000000),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),  # Remove map background
    plot.background = element_rect(fill = NA, color = NA)    # Remove outer background
  )
ggsave(here("Figures", "all_cores", "map", paste0("map_Antarctica.pdf")), dpi = 900) # Figure 1B


#================heatmap shared taxa between cores========================
# make core pairing table
df_dissimilarity <- crossing(core1 = order_sites, core2 = order_sites) %>% 
  filter(core1 != core2) %>%
  rowwise() %>%
  mutate(
    # Pull data for each core
    c1_years = list(meta %>% filter(core == core1, years >= 0) %>% pull(years)),
    c2_years = list(meta %>% filter(core == core2, years >= 0) %>% pull(years)),
    
    # Overlap range
    min_years = max(min(c1_years), min(c2_years)),
    max_years = min(max(c1_years), max(c2_years)),
    
    # Median spacing
    med1 = median(diff(sort(c1_years))),
    med2 = median(diff(sort(c2_years))),
    
    # Choose bin width per pair
    bin_width = max(med1, med2, na.rm = TRUE),
    filter = T) %>%
  ungroup %>% 
  select(core1, core2, bin_width, min_years, max_years) %>%
  mutate(similarity = pmap_dbl(
    list(core1, core2, bin_width, min_years, max_years),
    ~ compute_similarity(..1, ..2, ..3, ..4, ..5)
  )) %>%
  mutate(core1 = factor(core1, levels = order_sites, labels = names_sites),
         core2 = factor(core2, levels = rev(order_sites), labels = rev(names_sites)))

df_dissimilarity %>%
  ggplot(aes(core1, core2, fill = similarity)) +
  geom_tile() +
  theme_minimal() +
  labs(
    x = element_blank(),
    y = element_blank(),
    fill = 'Similarity'
  ) +
  scale_fill_gradient(
    low = "#000075",
    high = 'white'
  ) +
  theme(legend.position = 'top')
ggsave(here("Figures", "all_cores", "similarity", "Similarity_heatmap_hellinger.pdf"), dpi = 300, width = 10, height = 10) # Figure 2C and Supplementary Figure 6



#===============make virus host table===================
virus_host_df <- df_rel_abundance %>%
  left_join(meta) %>%
  filter(superkingdom == "Viruses") %>%
  group_by(taxon, core, epoche) %>%
  mutate(n_samples_present = n(),
         mean_abund = mean(rel_abundance),
         max_rel_abund = max(rel_abundance),
         n_abundant = sum(rel_abundance >= (0.3 * max_rel_abund))) %>%
  ungroup() %>%
  filter(n_samples_present > 3, 
         n_abundant > 3) %>%
  arrange(taxlevel, desc(mean_abund)) %>%
  dplyr::select(taxon, taxlevel) %>%
  unique
write.table(virus_host_df, file = here("results/Top_virus-host_pairs_prefiltered.tsv"), sep = "\t", quote = F, row.names = F) # preselection for literature search

#### this table was filled with corresponding hosts based on literature search
df_vir_hosts <- read.delim(here("data/SupplementaryTable3.tsv"), sep = "\t", header = T) %>%
  mutate(virhosts = paste0(taxon, ";", host_taxon_correlation)) %>%
  arrange(host_taxon)

# in which cores abundant?
df_vir_hosts_cores <- df_rel_abundance %>%
  left_join(meta, by = c("sample", "core")) %>%
  filter(superkingdom == "Viruses") %>%
  filter(taxlevel == "species") %>%
  dplyr::select(taxon, core, epoche, rel_abundance, sample) %>%
  pivot_wider(id_cols = c("core", "epoche", "sample"), names_from = "taxon", values_from = "rel_abundance") %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols = -c("core", "epoche", "sample"), names_to = "taxon", values_to = "rel_abundance") %>%
  replace(is.na(.), 0) %>%
  group_by(taxon, core, epoche) %>%
  mutate(n_samples_present = sum(rel_abundance > 0),
         mean_abund = mean(rel_abundance),
         max_rel_abund = max(rel_abundance),
         n_abundant = sum(rel_abundance >= (0.1 * max_rel_abund))) %>%
  ungroup() %>%
  filter(n_samples_present > 3,
         n_abundant > 3) %>%
  arrange(desc(mean_abund)) %>%
  dplyr::select(taxon, core) %>%
  unique %>%
  mutate(abundant = TRUE) %>%
  pivot_wider(id_cols = c("taxon"), names_from = "core", values_from = "abundant") %>%
  replace(is.na(.), FALSE) %>%
  rowwise %>%
  mutate(num_cores = sum(c(ulu, btoko, ilirney, lama, lele, PS97, KL77, KL12), na.rm = T)) %>%
  mutate(ecosystem_presence = case_when(!ulu & !btoko & !ilirney & !lama & !lele ~ "marine",
                                        !PS97 & !KL77 & !`KL12` ~ "lake",
                                        TRUE ~ "both")) %>%
  merge(df_vir_hosts, by = "taxon", all.x = F, all.y = T)

#=========================virus-host co-occurrence patterns====================
host_tax <- unique(df_vir_hosts_cores$host_taxon)[1]
df_cor_spearman_whole_core_combined_taxlevels <- lapply(unique(df_vir_hosts_cores$host_taxon), function(host_tax) {
    print(host_tax)
    host_correlation <- df_vir_hosts_cores$host_taxon_correlation[df_vir_hosts_cores$host_taxon == host_tax] %>% unique
    print(host_correlation)
    viruses <- df_vir_hosts_cores %>% filter(host_taxon == host_tax) %>% pull(taxon)
    
    if (!is.na(host_correlation) & any(sapply(viruses, function(x) { x %in% df_rel_abundance$taxon}))) {
      df_vir <- df_rel_abundance %>%
        filter((taxon %in% viruses & taxlevel == "species") | taxon == host_correlation) %>%
        dplyr::select(taxon, core, years, cum_taxpercent, sample) %>%
        mutate(vir_host = case_when(taxon == host_correlation ~ "host",
                                    taxon %in% viruses ~ "virus")) %>%
        group_by(sample, vir_host) %>%
        reframe(cum_taxpercent = sum(cum_taxpercent, na.rm = T)) %>%
        ungroup %>% 
        complete(sample, vir_host, fill = list(cum_taxpercent = 0)) %>%
        left_join(meta) %>%
        group_by(vir_host, core, years) %>%
        reframe(cum_taxpercent = mean(cum_taxpercent)) %>% 
        ungroup() %>% 
        pivot_wider(id_cols = c("core", "years"), names_from = "vir_host", values_from = "cum_taxpercent") %>%
        replace(is.na(.), 0) %>% 
        group_by(core) %>% 
        reframe(
          # Pearson
          Pearson = if (dplyr::n() >= 3 && sd(virus) > 0 && sd(host) > 0)
            cor(virus, host, method = "pearson") else NA_real_,
          Pearson_p = if (dplyr::n() >= 3 && sd(virus) > 0 && sd(host) > 0)
            cor.test(virus, host, method = "pearson")$p.value else NA_real_,
          
          # Spearman
          Spearman = if (dplyr::n() >= 3 && sd(virus) > 0 && sd(host) > 0)
            cor(virus, host, method = "spearman") else NA_real_,
          Spearman_p = if (dplyr::n() >= 3 && sd(virus) > 0 && sd(host) > 0)
            cor.test(virus, host, method = "spearman", exact = FALSE)$p.value else NA_real_
        ) %>%
        ungroup() %>%
        dplyr::select(core, Pearson, Pearson_p, Spearman, Spearman_p) %>%
        mutate(viruses = list(viruses),
               host_taxon = host_tax)
    } else { NULL }
  }) %>% do.call(rbind, .) %>%
    mutate(ecosystem = case_when(
      core %in% marine ~ "marine",
      TRUE ~ "lake"
    )) %>%
    merge(df_vir_hosts_cores, by = "host_taxon") %>%
    unique %>%
    pivot_longer(cols = c("btoko", "ilirney", "KL77", "lama", "lele", "KL12", "PS97", "ulu"), names_to = "core_presence", values_to = "core_abspres") %>%
    filter(core_abspres,
           core == core_presence)

plot_correlations <- 
  df_cor_spearman_whole_core_combined_taxlevels %>%
  filter(!grepl("unassigned", host_taxon_correlation)) %>%
  dplyr::select(host_superkingdom, ecosystem, core, host_group, Spearman, host_taxon) %>%
  unique %>% 
  group_by(host_superkingdom, ecosystem, host_group) %>%
  mutate(nsize = n(),
         host_group = factor(host_group, 
                               levels = rev(c("Bacteria", "cyanobacteria", "pseudomonadota", "Psychrobacillus", "Flavobacterium", "Eukaryota", "algae", "amoeba", "metazoa", "archaea", "viruses")),
                               labels = rev(c("Bacteria", "Cyanobacteria", "Pseudomonadota", "Psychrobacillus", "Flavobacterium", "Eukaryota", "Algae", "Amoeba", "Metazoa", "Archaea", "Viruses")))) %>%
  ungroup() %>%
  arrange(host_group) %>%
  mutate(host_group = factor(host_group, levels = unique(host_group), 
                               labels = (dplyr::select(., host_group, nsize, ecosystem) %>%
                                           unique %>%
                                           pivot_wider(id_cols = "host_group", names_from = "ecosystem", values_from = "nsize") %>%
                                           replace(is.na(.), 0) %>%
                                           mutate(host_group = paste0(host_group, " (n = ", marine, "; n = ", lake, ")")) %>%
                                           pull(host_group)))) %>% #mutate(core = factor(core, levels = order_sites)) %>%arrange(host_superkingdom, core, host_taxon) %>%
  ggplot(aes(y = host_group, 
             x = Spearman)) +
  geom_point(aes(color = host_superkingdom), alpha = 0.7, size = 3) +
  geom_boxplot(data = (. %>% filter(nsize > 4)), fill = alpha("white", 0), outlier.size = 1) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted") +
  geom_point(data = (. %>% filter(nsize <= 4 & nsize > 2) %>% 
                       group_by(host_superkingdom, ecosystem, host_group) %>%
                       reframe(median_val = median(Spearman),
                                 Q1 = quantile(Spearman, 0.25),
                                 Q3 = quantile(Spearman, 0.75))),
             aes(x = median_val, y = host_group), shape = 23, size = 4, fill = alpha("white", 0)) +
  geom_errorbarh(data = (. %>% filter(nsize <= 4 & nsize > 2) %>% 
                           group_by(host_superkingdom, ecosystem, host_group) %>%
                           reframe(median_val = median(Spearman),
                                     Q1 = quantile(Spearman, 0.25),
                                     Q3 = quantile(Spearman, 0.75))),
                 aes(xmin = Q1, xmax = Q3, y = host_group, x = NULL), height = 0, linewidth = 0.7) +
  theme_minimal() +
  facet_wrap(~factor(ecosystem, levels = c("marine", "lake"), labels = c("Marine", "Lake"))) +
  ylab(NULL) +
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values = customPalette[c(2,3,6,4)]) +
  xlim(-1,1) 

plot_correlations_combined <-
  df_cor_spearman_whole_core_combined_taxlevels %>%
  filter(!grepl("unassigned", host_taxon_correlation)) %>%
  dplyr::select(host_superkingdom, ecosystem, core, host_group, Spearman, host_taxon) %>%
  unique %>%
  ggplot(aes(y = factor("All"), x = Spearman)) +
  geom_point(alpha = 0.2, size = 3) +
  geom_boxplot(fill = alpha("white", 0), outlier.size = 1) +
  geom_vline(xintercept = 0, col = "black", linetype = "dotted") +
  theme_minimal() +
  facet_wrap(~factor(ecosystem, levels = c("marine", "lake"), labels = c("Marine", "Lake"))) +
  theme(panel.spacing = unit(1.5, "lines"),
        legend.position="none",
        strip.text = element_blank()) +
  xlim(-1,1) +
  ylab(NULL) +
  xlab("Spearman coefficient")

# combine plots
(correlations_plot <- plot_grid(
  plot_correlations, plot_correlations_combined, 
  ncol = 1,
  rel_heights = c(5,1), 
  align = "v"))

ggsave(here("Figures", "all_cores", "correlations", "Taxonomic_levels_correlations_combined_final_selection.pdf"), dpi = 300, width = 10, height = 4) # Figure 3A


wilcox_results <- df_cor_spearman_whole_core_combined_taxlevels %>%
  filter(!grepl("unassigned", host_taxon_correlation)) %>%
  dplyr::select(host_superkingdom, ecosystem, core, host_group, Spearman, host_taxon) %>%
  unique %>%
  group_by(host_group, ecosystem) %>%
  summarise(
    n = sum(!is.na(Spearman)),
    wilcox_p = if (n >= 1 && any(Spearman != 0, na.rm = TRUE)) {
      wilcox.test(Spearman, mu = 0, alternative = "greater")$p.value
    } else {
      NA_real_
    }) %>%
  ungroup %>%
  mutate(wilcox_p = case_when(wilcox_p > 0.05 ~ "ns",
                              wilcox_p < 0.001 ~ "***",
                              wilcox_p < 0.01 ~ "**",
                              wilcox_p < 0.05 ~ "*")) %>%
  pivot_wider(id_cols = host_group, names_from = ecosystem, values_from = wilcox_p)

df_cor_spearman_whole_core_combined_taxlevels %>%
  filter(ecosystem == "marine") %>%
  filter(!grepl("unassigned", host_taxon_correlation)) %>%
  dplyr::select(host_superkingdom, ecosystem, core, host_group, Spearman, host_taxon) %>%
  unique %>%
  pull(Spearman) %>%
  wilcox.test(., alternative = "greater")

df_cor_spearman_whole_core_combined_taxlevels %>%
  filter(ecosystem == "lake") %>%
  filter(!grepl("unassigned", host_taxon_correlation)) %>%
  dplyr::select(host_superkingdom, ecosystem, core, host_group, Spearman, host_taxon) %>%
  unique %>%
  pull(Spearman) %>% 
  wilcox.test(., alternative = "greater")

#=======================Micromonas=======================
for (lake in marine) {
  tmp_path <- here("Figures", lake, "vir_hosts", "micromonas")
  dir.create(tmp_path, recursive = T)
  
  plot_viruses <- df_rel_abundance %>%
    filter(core == lake,
           grepl("Micromonas", taxon),
           superkingdom == "Viruses") %>%
    mutate(taxon = "Micromonas-infecting viruses") %>%
    dplyr::select(sample, taxon, cum_taxpercent) %>%
    group_by(sample, taxon) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>% 
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>% 
    mutate(taxon = case_when(is.na(taxon) ~ "Micromonas-infecting viruses",
                             T ~ taxon),
           cum_taxpercent = case_when(is.na(cum_taxpercent) ~ 0,
                                      T ~ cum_taxpercent)) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    scale_fill_manual(values = virusPalette(tmp %>% pull(taxon) %>% unique %>% length)) +
    ylab("cumulative virus abundance [%]") +
    labs(fill = NULL) +
    theme_minimal() +
    xlab("cal ka BP") +
    theme(axis.text = element_text(size = 12))
  
  
  plot_host <- df_rel_abundance %>%
    filter(core == lake,
           taxon == "Micromonas") %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>% 
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    mutate(taxon = case_when(is.na(taxon) ~ "Micromonas",
                             T ~ taxon),
           cum_taxpercent = case_when(is.na(cum_taxpercent) ~ 0,
                                      T ~ cum_taxpercent)) %>%
    group_by(taxon, ka, core) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative host abundance [%]") +
    labs(fill = NULL) +
    theme_minimal() +
    xlab(NULL) +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank())
  
  
  
  overview_plot <- plot_grid(
    plot_host, plot_viruses, ncol = 1,
    rel_heights = c(1,1), 
    align = "v")
  
  ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_virus_all_vs_host_micromonas_grey.pdf")), overview_plot, dpi = 900, width = 15, height = 10) # Figure 3B
  
}

#=======================Chrysochromulina=================
for (lake in order_sites) {
  tmp_path <- here("Figures", lake, "vir_hosts", "algae")
  dir.create(tmp_path, recursive = T)
  
  plot_viruses <- df_rel_abundance %>%
    filter(core == lake,
           (taxon %in% df_vir_hosts_cores$taxon[df_vir_hosts_cores$host_taxon == "Chrysochromulina"]) | grepl("Chrysochromulina", taxon),
           superkingdom == "Viruses") %>%
    mutate(taxon = "Chrysochromulina-infecting viruses") %>%
    dplyr::select(sample, taxon, cum_taxpercent) %>%
    group_by(taxon, sample) %>% 
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>% 
    mutate(taxon = case_when(is.na(taxon) ~ "Chrysochromulina-infecting viruses",
                             T ~ taxon),
           cum_taxpercent = case_when(is.na(cum_taxpercent) ~ 0,
                                      T ~ cum_taxpercent)) %>%
    mutate(taxon = "Chrysochromulina-infecting viruses") %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp } %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative virus abundance [%]") +
    labs(fill = NULL) +
    xlab("cal ka BP") +
    theme_minimal() +
    scale_fill_manual(values = virusPalette(tmp %>% pull(taxon) %>% unique %>% length))
  
  
  plot_host <- df_rel_abundance %>%
    filter(core == lake,
           taxon == "Haptophyta") %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>% 
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp } %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative host abundance [%]") +
    labs(fill = NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank())
  
  
  overview_plot <- plot_grid(
    plot_host, plot_viruses, ncol = 1,
    rel_heights = c(1,1), 
    align = "v")
  
  ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_Chrysochromulina_virus_all_vs_host_haptophyta_grey.pdf")), overview_plot, dpi = 900, width = 15, height = 10) # Figure 3C
  
}

df_cor_spearman_whole_core_combined_taxlevels %>%
  filter(host_taxon == "Chrysochromulina") %>%
  select(host_taxon, core, Spearman, Spearman_p, ecosystem, host_taxon_correlation, host_group, host_superkingdom) %>%
  unique

#=======================Pelagibacter=====================
for (lake in marine) {
  tmp_path <- here("Figures", lake, "vir_hosts", "pelagibacter")
  dir.create(tmp_path, recursive = T)
  
  pelagibacterviruses <- c(df_viruses %>% filter(grepl("Pelagibacter phage", taxon), !grepl("virophage", taxon), superkingdom == "Viruses", core == lake) %>% pull(species),
                           df_vir_hosts_cores$taxon[df_vir_hosts_cores$host_taxon == "Pelagibacter"]) %>% unique
  pelagibacterviruses <- df_rel_abundance %>% 
    filter(taxon %in% pelagibacterviruses, core == lake) %>%
    group_by(taxon) %>%
    reframe(mean_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(desc(mean_taxpercent)) %>% 
    pull(taxon) %>% 
    unique
  
  plot_viruses <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% pelagibacterviruses ,
           superkingdom == "Viruses") %>%
    mutate(taxon = "Pelagibacter-infecting viruses") %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative virus abundance [%]") +
    scale_fill_manual(values = virusPalette(tmp %>% pull(taxon) %>% unique %>% length)) +
    labs(fill = NULL) +
    theme_minimal() +
    xlab("cal ka BP") +
    geom_vline(xintercept = 13, color = "black", linetype = "dotted")
  
  pelagibacter <- df_rel_abundance %>% filter(analysis_group == "Pelagibacter") %>% pull(taxon) %>% unique
  pelagibacter <- df_rel_abundance %>% 
    filter(taxon %in% pelagibacter, core == lake) %>%
    group_by(taxon) %>%
    reframe(mean_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(desc(mean_taxpercent)) %>% 
    pull(taxon) %>% 
    unique
  
  plot_host1 <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% pelagibacter) %>%
    mutate(taxon = case_when(taxon == "Candidatus Pelagibacter sp. IMCC9063" ~ taxon,
                             TRUE ~ "other Pelagibacter")) %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>% 
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent, fill = taxon)) +
    geom_area(position = "stack", stat = "identity") +
    ylab("cumulative host abundance [%]") +
    scale_fill_manual(values = c("Candidatus Pelagibacter sp. IMCC9063" = "#800000", "other Pelagibacter" = "#707070")) +
    labs(fill = NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank()) +
    geom_vline(xintercept = 13, color = "black", linetype = "dotted")
  
  plot_host2 <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% pelagibacter,
           taxon != "Candidatus Pelagibacter sp. IMCC9063") %>%
    mutate(taxon = "other Pelagibacter") %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>% 
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative host abundance [%]") +
    labs(fill = NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank()) +
    geom_vline(xintercept = 13, color = "black", linetype = "dotted")
  
  
  
  overview_plot <- plot_grid(
    plot_host1, plot_host2, plot_viruses, ncol = 1,
    rel_heights = c(1,1,1), 
    align = "v") 
  
  ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_Pelagibacter_virus_all_vs_host_woIMCC9063_grey.pdf")), overview_plot, dpi = 900, width = 10, height = 8) # Figure 4
  
}

df_rel_abundance %>%
  filter(core == "KL77",
         taxon %in% { c(df_viruses %>% filter(grepl("Pelagibacter phage", taxon), !grepl("virophage", taxon), superkingdom == "Viruses", core == lake) %>% pull(species),
                      df_vir_hosts_cores$taxon[df_vir_hosts_cores$host_taxon == "Pelagibacter"]) %>% unique },
         taxlevel == "species") %>%
  merge(meta[meta$core == "KL77", ], all.y = T) %>%
  filter(years >= 0) %>%
  mutate(taxon = "virus") %>%
  group_by(sample, years, ka, core, epocheII, taxon) %>%
  summarize(cum_taxpercent = sum(cum_taxpercent)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("core", "sample", "years", "ka", "epocheII"), names_from = "taxon", values_from = "cum_taxpercent") %>%
  replace(is.na(.), 0) %>% 
  #pivot_longer(cols = -c("core", "sample", "years", "ka", "epocheII"), names_to = "virus", values_to = "cum_taxpercent") %>% 
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon %in% (df_rel_abundance %>% filter(analysis_group == "Pelagibacter", taxon != "Candidatus Pelagibacter sp. IMCC9063") %>% pull(taxon) %>% unique)) %>% group_by(sample) %>% summarize(woIMCC9063 = sum(cum_taxpercent)) %>% ungroup(), by = "sample") %>%
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon == "Pelagibacter") %>% rename(all = cum_taxpercent) %>% select(sample, all), by = "sample") %>% 
  pivot_longer(cols = c("woIMCC9063", "all"), values_to = "host", names_to = "host_group") %>%
  group_by(epocheII, host_group) %>%
  summarize(cor_test = list(cor.test(virus, host, method = "spearman"))) %>% 
  ungroup() %>%
  mutate(
    cor = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  arrange(epocheII, desc(cor)) %>%
  mutate(significant = case_when(p_value < 0.05 ~ T,
                                 T ~ F)) %>% print


#=======================Bacteriophages correlations======
df_bac_paths <- read.delim(here('data/Bacteria_taxPath.tab'), sep = "\t", header = T)
viruses <- df_rel_abundance %>% filter(superkingdom == "Viruses") %>% pull(taxon) %>% unique %>% .[!is.na(.)]

# filter table for present viruses
df_viruses_ictv <- ictv_df %>% 
  filter(Family %in% viruses) %>% 
  dplyr::select(Family, Virus) %>% 
  unique()

df_vir_genus_host_ictv <- ictv_df %>%
  select(Family, Genus, Virus) %>%
  unique %>%
  mutate(Host = sapply(Virus, function(x) {strsplit(x, split = " ")[[1]][1]})) %>%
  merge(df_bac_paths, by.x = "Host", by.y = "genus", all.x = T, all.y = F) %>%
  dplyr::select(Family, Genus, class) %>%
  unique %>%
  filter(!is.na(class)) %>%
  group_by(Family, Genus) %>%
  mutate(n_hosts = n_distinct(class)) %>%
  ungroup %>%
  mutate(class = case_when(n_hosts == 1 ~ class, 
                           T ~ "Uncertain")) %>%
  filter(!is.na(Genus)) %>%
  rename(Host = class)


df_vir_ncbi_host <- df_viruses %>%
  select(taxon, ncbi_name, genus, family) %>%
  rename(Genus = genus,
         Family = family) %>%
  unique %>%
  mutate(Host = sapply(ncbi_name, function(x) {strsplit(x, split = " ")[[1]][1]})) %>%
  merge(df_bac_paths, by.x = "Host", by.y = "genus", all.x = T, all.y = F) %>%
  dplyr::select(Genus, Family, class) %>%
  unique %>%
  filter(!is.na(class)) %>%
  group_by(Family, Genus) %>%
  mutate(n_hosts = n_distinct(class)) %>%
  ungroup %>%
  mutate(class = case_when(n_hosts == 1 ~ class, 
                           T ~ "Uncertain")) %>%
  filter(!is.na(Genus)) %>%
  rename(Host = class)
  
lake <- "lama"
{
  bacterial_host_df <- df_vir_genus_host_ictv %>% select(Genus, Host) %>%
    rbind(df_vir_ncbi_host %>% select(Genus, Host))
  
  bacteriophages <- c(viruses_unassigned %>% 
                 filter(family %in% c("Autographiviridae", "Mesyanzhinovviridae", "Shitoviridae", "Zobellviridae", "Casjensviridae", "Peduoviridae"), 
                        family != "unassigned Caudoviricetes genus") %>% pull(genus),
               "Ackermannviridae", "Straboviridae", "Guelinviridae", "Steigviridae", 
               "Grimontviridae", "Herelleviridae", "Rountreeviridae", "Suoliviridae",
               "Aggregaviridae", "Demerecviridae",
               viruses_unassigned %>% filter(class == "Caudoviricetes", family != "unassigned Caudoviricetes family",
                                                       !family %in% c("Autographiviridae", "Mesyanzhinovviridae", "Shitoviridae", "Zobellviridae", "Casjensviridae", "Peduoviridae")) %>% pull(family) %>% unique)
  
  for (lake in order_sites) {
    tmp_path <- here("Figures", lake, "vir_hosts", "bacteria")
    dir.create(tmp_path, recursive = T)
    
    overview_plot <- df_rel_abundance %>%
        filter(core == lake,
               (taxlevel == "class" & superkingdom == "Bacteria") | taxon %in% bacteriophages | taxon %in% bacterial_host_df$Genus ) %>% 
        merge(bacterial_host_df, by.x = "taxon", by.y = "Genus", all.x = T, all.y = F) %>% 
        mutate(virhost = case_when(superkingdom == "Viruses" ~ "Inferred Hosts of Viral Community",
                                   superkingdom == "Bacteria" ~ "Bacteria")) %>% 
        mutate(taxon = case_when(superkingdom == "Viruses" & is.na(Host) ~ "unknown host",
                                 superkingdom == "Viruses" & Host == "Uncertain" ~ "unknown host",
                                 superkingdom == "Viruses" ~ Host,
                                 superkingdom == "Bacteria" & taxon %in% (c("Actinomycetes", "Alphaproteobacteria", "Bacilli",
                                                                            "Betaproteobacteria", "Cyanophyceae", "Flavobacteriia",
                                                                            "Clostridia", "Gammaproteobacteria", "Bacteroidia") %>% unique) ~ taxon,
                                 superkingdom == "Bacteria" ~ "other Bacteria")) %>% 
        mutate(taxon = case_when(superkingdom == "Viruses" & !(taxon %in% c("Actinomycetes", "Alphaproteobacteria", "Bacilli",
                                                                           "Betaproteobacteria", "Cyanophyceae", "Flavobacteriia",
                                                                           "Clostridia", "Gammaproteobacteria", "Bacteroidia", "unknown host")) ~ "other Bacteria",
                                 T ~ taxon)) %>%
        group_by(taxon, sample, core, virhost) %>%
        reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
        merge(meta[meta$core == lake, ], all.y = T) %>%
        filter(years >= 0 | ka >= 0) %>% 
        dplyr::select(taxon, sample, years, ka, core, virhost, cum_taxpercent) %>%
        pivot_wider(id_cols = c("core", "sample", "years", "ka", "virhost"), names_from = "taxon", values_from = "cum_taxpercent") %>%
        replace(is.na(.), 0) %>%
        pivot_longer(cols = -c("core", "sample", "years", "ka", "virhost"), names_to = "taxon", values_to = "cum_taxpercent") %>%
        filter(taxon != "NA") %>%
        group_by(ka, taxon, virhost) %>%
        reframe(cum_taxpercent = mean(cum_taxpercent, na.rm = T)) %>%
        ungroup %>% 
        complete(ka, virhost, taxon, fill = list(cum_taxpercent = 0)) %>%
        ggplot(aes(x = ka, y = cum_taxpercent, fill = taxon)) +
        geom_area(position = "stack", stat = "identity") +
        ylab("cumulative abundance [%]") +
        xlab("cal ka BP") +
        labs(fill = NULL) +
        scale_fill_manual(values = bacteriaPalette) +
        facet_wrap(~virhost, ncol = 1, scales = "free_y") +
        theme_minimal() +
        theme(strip.text = element_text(face = "bold", size = 13),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 13),
              legend.position = "bottom")
    
    
    ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_virus_host_classes_vs_host_Bacteria_classes_selected.pdf")), overview_plot, dpi = 900) # Figure 5A and Figure 5B and Supplementary Figures 18, 21 and 22
    
  }
}

{
  df_phage_correlations <- df_rel_abundance %>%
    filter(taxon %in% bacterial_host_df$Host | taxon %in% bacterial_host_df$Genus) %>%         
    merge(bacterial_host_df, by.x = "taxon", by.y = "Genus", all.x = T, all.y = F) %>%
    mutate(virhost = case_when(superkingdom == "Viruses" ~ "virus",
                               superkingdom == "Bacteria" ~ "host")) %>%
    mutate(taxon = case_when(superkingdom == "Viruses" ~ Host,
                             superkingdom == "Bacteria" ~ taxon)) %>%
    group_by(taxon, sample, core, virhost) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>% 
    merge(meta, all.y = T) %>%
    filter(years >= 0) %>% 
    dplyr::select(sample, years, ka, core, taxon, virhost, cum_taxpercent, ecosystem) %>%
    group_by(ka, virhost, ecosystem, taxon, core) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent, na.rm = T)) %>%
    ungroup %>%
    filter(virhost != "NA", taxon != "Uncertain") %>%
    pivot_wider(id_cols = c("core", "ecosystem", "ka", "taxon"), names_from = "virhost", values_from = "cum_taxpercent", values_fill = 0) %>%
    group_by(core, ecosystem, taxon) %>%
    reframe(Pearson = ifelse(sd(virus) == 0 | sd(host) == 0, NA, cor(virus, host)),
              Spearman = ifelse(sd(virus) == 0 | sd(host) == 0, NA, cor(virus, host, method = "spearman"))) %>%
    ungroup() %>%
    dplyr::select(core, ecosystem, taxon, Pearson, Spearman) %>%
    filter(!is.na(Spearman)) %>%
    group_by(taxon, ecosystem) %>%
    mutate(nsize = n()) %>%
    ungroup() %>%
    arrange(taxon) %>%
    filter(taxon %in% c("Gammaproteobacteria", "Flavobacteriia", "Cyanophyceae", "Clostridia", "Betaproteobacteria", "Bacteroidia", "Bacilli", "Alphaproteobacteria", "Actinomycetes"))
  
  plot <- df_phage_correlations %>%
    mutate(taxon_col = taxon,
           taxon = factor(taxon, levels = unique(taxon), 
                          labels = (dplyr::select(., taxon, nsize, ecosystem) %>%
                                      unique %>%
                                      pivot_wider(id_cols = "taxon", names_from = "ecosystem", values_from = "nsize") %>%
                                      replace(is.na(.), 0) %>%
                                      mutate(taxon = paste0(taxon, " (n = ", marine, "; n = ", lake, ")")) %>%
                                      pull(taxon)))) %>%
    ggplot(aes(y = taxon, 
               x = Spearman)) +
    geom_point(aes(col = taxon_col), alpha = 0.7, size = 3) +
    geom_boxplot(data = (. %>% filter(nsize > 4)), alpha = 0, outlier.shape = 19, outlier.size = 1) +
    geom_vline(xintercept = 0, col = "black", linetype = "dotted") +
    geom_point(data = (. %>% filter(nsize <= 4 & nsize > 2) %>% 
                         group_by(ecosystem, taxon, taxon_col) %>%
                         reframe(median_val = median(Spearman),
                                   Q1 = quantile(Spearman, 0.25),
                                   Q3 = quantile(Spearman, 0.75))),
               aes(x = median_val, y = taxon), shape = 23, size = 4, fill = alpha("white", 0)) +
    geom_errorbarh(data = (. %>% filter(nsize <= 4 & nsize > 2) %>% 
                             group_by(ecosystem, taxon, taxon_col) %>%
                             reframe(median_val = median(Spearman),
                                       Q1 = quantile(Spearman, 0.25),
                                       Q3 = quantile(Spearman, 0.75))),
                   aes(xmin = Q1, xmax = Q3, y = taxon, x = NULL), height = 0, linewidth = 0.7) +
    theme_minimal() +
    facet_wrap(~factor(ecosystem, levels = c("marine", "lake"))) +
    ylab(NULL) +
    scale_color_manual(values = bacteriaPalette) +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlim(-1,1)
  
  plot_combined <- df_phage_correlations %>% 
    ggplot(aes(y = factor(paste0("All (n = ", df_phage_correlations %>% filter(ecosystem == "marine") %>% pull(core) %>% length, "; n = ", df_phage_correlations %>% filter(ecosystem == "lake") %>% pull(core) %>% length, ")")), x = Spearman)) +
    geom_point(alpha = 0.2, shape = 19, size = 3) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1, fill = alpha("white", 0)) +
    geom_vline(xintercept = 0, col = "black", linetype = "dotted") +
    theme_minimal() +
    facet_wrap(~factor(ecosystem, levels = c("marine", "lake"))) +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position="none",
          strip.text = element_blank()) +
    xlim(-1,1) +
    ylab(NULL) +
    xlab("Spearman correlation coefficient")
  
  # combine plots
  (plot_correlations <- plot_grid(
    plot, plot_combined, 
    ncol = 1,
    rel_heights = c(5,1), 
    align = "v"))
  
  ggsave(here("Figures", "all_cores", "correlations", "Correlations_virus_host_classes_vs_host_Bacteria_classes.pdf"), plot_correlations, dpi = 900, width = 10, height = 4) # Figure 5C
}

wilcox_results_bacteria <- df_phage_correlations %>%
  dplyr::select(taxon, ecosystem, core, Spearman) %>%
  unique %>%
  group_by(taxon, ecosystem) %>%
  summarise(
    n = sum(!is.na(Spearman)),
    wilcox_p = if (n >= 1 && any(Spearman != 0, na.rm = TRUE)) {
      wilcox.test(Spearman, mu = 0)$p.value #, alternative = "greater"
    } else {
      NA
    }) %>%
  ungroup %>%
  mutate(wilcox_p = case_when(wilcox_p > 0.05 ~ "ns",
                              wilcox_p < 0.001 ~ "***",
                              wilcox_p < 0.01 ~ "**",
                              wilcox_p < 0.05 ~ "*")) %>%
  select(-n) %>%
  pivot_wider(names_from = ecosystem, values_from = wilcox_p)

df_phage_correlations %>%
  filter(ecosystem == "marine") %>%
  dplyr::select(taxon, ecosystem, core, Spearman) %>%
  unique %>%
  pull(Spearman) %>%
  wilcox.test(., alternative = "greater")

df_phage_correlations %>%
  filter(ecosystem == "lake") %>%
  dplyr::select(taxon, ecosystem, core, Spearman) %>%
  unique %>%
  pull(Spearman) %>% 
  wilcox.test(., alternative = "greater")

#===========================Supplementary Figures=============================
#========================Rarefaction==========================================
df_viruses_rarefied <- read.delim(here("data", "Viruses_rarefied.tab"), sep = "\t", header = T)

(rarefaction_plot <- df_viruses_rarefied %>%
    group_by(depth) %>%
    mutate(perc = n_species / sum(n_species)) %>%
    ungroup %>%
    { . ->> tmp } %>%
    ggplot() +
    geom_bar(aes(x = depth, y = perc, fill = ecosystem),
             stat = "identity", alpha = 0.85) +
    scale_x_log10() +
    geom_line(data = tmp %>% 
                group_by(depth) %>%
                summarize(n_species = sum(n_species)) %>%
                ungroup() %>%
                { . ->> tmp2 } %>%
                mutate(n_species = n_species / 2000),
              aes(x = depth, y = n_species), col = "black", linewidth = 0.75) +
    scale_y_continuous(
      name = "proportion of distinct viral species per ecosystem",
      sec.axis = sec_axis( trans= ~ . * 2000, 
                           name="number of distinct viral species")) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    xlab("rarefaction read depth") +
    scale_fill_manual(values = envPalette))
ggsave(here("Figures", "all_cores", "shared_taxa", "Rarefaction_analysis.pdf"), rarefaction_plot, dpi = 900, width = 15, height = 10) # Supplementary Figure 3

#===========================Caudovirales======================================
df_rel_abundance %>%
  filter(taxon == "uncultured Caudovirales phage") %>%
  group_by(ka, core) %>%
  summarize(rel_abundance = mean(rel_abundance)) %>%
  ungroup %>%
  ggplot(aes(x = ka, y = rel_abundance)) +
  geom_area() +
  xlab("cal ka BP") +
  ylab("percent viral reads [%]") +
  theme_minimal() +
  facet_wrap(~factor(core, levels = order_sites, labels = names_sites), ncol = 1, scales = "free_y") +
  theme(strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13))
ggsave(here("Figures", "all_cores", "taxa", "Stratigraphic_plot_Uncultured_Caudovirales_phage.pdf"), dpi = 900, height = 20, width = 15.93, units = "cm") # Supplementary Figure 4

#========================Shared taxa==========================================
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species") %>%
  left_join(meta, by = c("sample", "ka", "years", "core")) %>%
  dplyr::select(core, taxon) %>%
  unique %>%
  group_by(core) %>%
  summarize(unique_taxa = list(taxon)) %>%
  mutate(core = factor(core, levels = order_sites, labels = names_sites)) %>%
  deframe() %>%
  ggVennDiagram(force_upset = T)
ggsave(here("Figures", "all_cores", "venn", "UpsetPlot_lakes.pdf"), dpi = 900, height = 15, width = 15.93, units = "cm") # Supplementary Figure 5

#=======================Stratigraphic plots===================================
for(lake in marine) {
  df_rel_abundance %>%
    filter(core == lake,
           taxlevel == "genus" | taxon == "Pithovirus LCPAC304",
           superkingdom == "Viruses") %>%
    dplyr::select(sample, taxon, rel_abundance) %>% 
    left_join(meta %>% filter(core == lake, ka >= 0) %>% dplyr::select(sample)) %>% 
    complete(sample, taxon, fill = list(rel_abundance = 0)) %>%
    left_join(meta %>% filter(core == lake, ka >= 0) %>% dplyr::select(sample, ka)) %>% 
    group_by(taxon, ka) %>% 
    summarize(rel_abundance = mean(rel_abundance)) %>%
    ungroup %>% 
    mutate(taxon = case_when(lake != "PS97" & taxon %in% c("Prymnesiovirus", "Prasinovirus",
                                                           "Pelagivirus", "Stopavirus", "Llyrvirus", "Mazuvirus") ~ taxon,
                             lake == "PS97" & taxon %in% c("Prymnesiovirus", "Prasinovirus", "Theiavirus", "Pithovirus LCPAC304") ~ taxon,
                             taxon == "unassigned Caudoviricetes genus" ~ "unassigned Caudoviricetes", 
                             taxon == "unassigned Viruses genus" ~ "unassigned Viruses",
                             TRUE ~ paste0("other virus genera"))) %>%
    group_by(taxon, ka) %>%
    summarize(rel_abundance = sum(rel_abundance)) %>%
    ungroup %>%
    {if (lake == "PS97") {
    pivot_wider(., id_cols = "ka", names_from = "taxon", values_from = "rel_abundance") %>%
    mutate(`unassigned Viruses` = `unassigned Viruses` - `Pithovirus LCPAC304`) %>%
    pivot_longer(cols = -c("ka"), names_to = "taxon", values_to = "rel_abundance")
    } else {
      .
    }} %>%
    mutate(host = case_when(taxon %in% c("Prymnesiovirus") ~ "Haptophyta",
                              taxon %in% c("Prasinovirus") ~ "Prasinophytae",
                              taxon %in% c("Theiavirus") ~ "Protists",
                              taxon %in% c("Pelagivirus", "Stopavirus") ~ "Pelagibacter",
                              taxon %in% c("Llyrvirus", "Mazuvirus") ~ "Synechococcus",
                              taxon == "unassigned Caudoviricetes" ~ "Bacteria/Archaea",
                              TRUE ~ "Unknown")) %>%
    ggplot(aes(x = ka, y = rel_abundance, fill = host)) +
    geom_area(stat = "identity") +
    {if (lake != "PS97") {
      geom_vline(xintercept = 14, linetype = "dotted")
    } else {
      geom_vline(xintercept = 12.7, linetype = "dotted")
    }} +
    theme_minimal() +
    ylab("relative virus abundance [%]") +
    xlab("cal ka BP") +
    {if (lake != "PS97") {
      facet_wrap(~factor(taxon, 
                         levels = c("unassigned Viruses", "other virus genera", "Prasinovirus", "unassigned Caudoviricetes",
                                    "Pelagivirus", "Stopavirus", "Llyrvirus", "Mazuvirus", "Prymnesiovirus")), 
                 ncol = 1, scales = "free_y", strip.position = "right")
    } else {
      facet_wrap(~factor(taxon, 
                         levels = c("unassigned Viruses", "other virus genera", "Prymnesiovirus", "unassigned Caudoviricetes", "Prasinovirus", "Theiavirus", "Pithovirus LCPAC304")), 
                 ncol = 1, scales = "free_y", strip.position = "right")
    }} +
    theme(axis.text.y = element_text(size = 8),
          strip.text.y = element_text(angle = 0),
          legend.position = "bottom") +
    scale_fill_manual(values = host_groupPalette_marine) +
    scale_y_continuous(labels = label_number(sci = FALSE))
  
  ggsave(here("Figures", lake, "stratigraphic_plots", paste0("Virus_composition_", lake, "_genus_selected.pdf")), dpi = 900, height = 15, width = 20, units = "cm") # Supplementary Figures 7-9
  
}

for(lake in order_sites[!order_sites %in% marine]) {
  df_rel_abundance %>%
    filter(core == lake,
           taxlevel == "family",
           superkingdom == "Viruses") %>%
    dplyr::select(sample, taxon, rel_abundance) %>% 
    left_join(meta %>% filter(core == lake, ka >= 0) %>% dplyr::select(sample)) %>% 
    complete(sample, taxon, fill = list(rel_abundance = 0)) %>%
    left_join(meta %>% filter(core == lake, ka >= 0) %>% dplyr::select(sample, ka)) %>% 
    group_by(taxon, ka) %>% 
    summarize(rel_abundance = mean(rel_abundance)) %>%
    ungroup %>% 
    mutate(taxon = case_when(lake == "ilirney" & taxon %in% c("Mimiviridae", "Phycodnaviridae", "Mesyanzhinovviridae", "Casjensviridae") ~ taxon,
                             lake == "lele" & taxon %in% c("Lavidaviridae", "Mimiviridae", "Phycodnaviridae",
                                                           "Geminiviridae", "Iridoviridae",
                                                           "Kyanoviridae", "Marseilleviridae") ~ taxon,
                             lake == "lama" & taxon %in% c("Mimiviridae", "Phycodnaviridae", "Lavidaviridae",
                                                           "Theiavirus", "Geminiviridae",
                                                           "Asfarviridae", "Herelleviridae", "Casjensviridae") ~ taxon,
                             lake == "btoko" & taxon %in% c("Mimiviridae", "Phycodnaviridae", "Casjensviridae") ~ taxon,
                             lake == "ulu" & taxon %in% c("Mimiviridae", "Phycodnaviridae", 
                                                          "Casjensviridae", "Baculoviridae") ~ taxon,
                             taxon == "unassigned Caudoviricetes family" ~ "unassigned Caudoviricetes", 
                             taxon == "unassigned Viruses family" ~ "unassigned Viruses",
                             !grepl("viridae", taxon) ~ paste0("other virus families"))) %>%
    filter(!is.na(taxon)) %>%
    mutate(host = case_when(taxon %in% c("Mimiviridae", "Marseilleviridae") ~ "Protists",
                            taxon %in% c("Lavidaviridae") ~ "Mimiviridae",
                            taxon %in% c("Phycodnaviridae") ~ "Algae",
                            taxon %in% c("Casjensviridae", "Mesyanzhinovviridae",
                                         "Herelleviridae") ~ "Bacteria",
                            taxon %in% c("Baculoviridae") ~ "Invertebrates",
                            taxon %in% c("Geminiviridae") ~ "Plants",
                            taxon %in% c("Asfarviridae", "Iridoviridae") ~ "Metazoa",
                            taxon %in% c("Kyanoviridae") ~ "Cyanobacteria",
                            taxon == "unassigned Caudoviricetes" ~ "Bacteria/Archaea",
                            TRUE ~ "Unknown")) %>%
    group_by(taxon, ka, host) %>%
    summarize(rel_abundance = sum(rel_abundance)) %>%
    ungroup %>%
    ggplot(aes(x = ka, y = rel_abundance, fill = host)) +
    geom_area(stat = "identity") +
    geom_vline(xintercept = 14, linetype = "dotted") +
    theme_minimal() +
    ylab("relative virus abundance [%]") +
    xlab("cal ka BP") +
    {if (lake == "btoko") {
      facet_wrap(~factor(taxon,
                         levels = c("unassigned Viruses", "other virus families", "unassigned Caudoviricetes", "Mimiviridae",
                                    "Phycodnaviridae", "Casjensviridae")),
                 ncol = 1, scales = "free_y", strip.position = "right")
    } else if (lake == "ilirney") {
      facet_wrap(~factor(taxon,
                         levels = c("unassigned Viruses", "other virus families", "unassigned Caudoviricetes",
                                    "Mimiviridae", "Phycodnaviridae", "Casjensviridae", "Mesyanzhinovviridae")),
                 ncol = 1, scales = "free_y", strip.position = "right")
    } else if (lake == "lama") {
      facet_wrap(~factor(taxon,
                         levels = c("unassigned Viruses", "other virus families", "unassigned Caudoviricetes", "Geminiviridae",
                                    "Mimiviridae", "Asfarviridae",
                                    "Phycodnaviridae", "Herelleviridae", "Casjensviridae", "Lavidaviridae")),
                 ncol = 1, scales = "free_y", strip.position = "right")
    } else if (lake == "lele") {
      facet_wrap(~factor(taxon,
                         levels = c("unassigned Viruses", "other virus families", "unassigned Caudoviricetes", "Mimiviridae",
                                    "Marseilleviridae", "Phycodnaviridae", "Geminiviridae", "Kyanoviridae", "Iridoviridae",
                                    "Lavidaviridae")),
                 ncol = 1, scales = "free_y", strip.position = "right")
    } else if (lake == "ulu") {
      facet_wrap(~factor(taxon,
                         levels = c("unassigned Viruses", "other virus families", "unassigned Caudoviricetes", "Mimiviridae",
                                    "Phycodnaviridae", "Baculoviridae", "Casjensviridae")),
                 ncol = 1, scales = "free_y", strip.position = "right")
    }} +
    theme(panel.spacing = unit(1, "lines"),
          axis.text.y = element_text(size = 8),
          strip.text.y = element_text(angle = 0),
          legend.position = "bottom") +
    scale_fill_manual(values = host_groupPalette_lakes) +
    scale_y_continuous(labels = label_number(sci = FALSE))
  
  ggsave(here("Figures", lake, "stratigraphic_plots", paste0("Virus_abundance_", lake, "_family_selected.pdf")), dpi = 900, height = 20, width = 20, units = "cm") # Supplementary Figures 10-14
  
}

#===========================marine taxa Levinson-Lessing======================
df_rel_abundance %>%
  filter(superkingdom == "Viruses",
         taxlevel == "species") %>%
  group_by(taxon) %>%
  summarize(presence = case_when(any(core == "KL77") & any(core == "KL12") & any(core == "PS97") & any(core == "lele") ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  ungroup %>% 
  filter(presence,
         taxon != "uncultured Caudovirales phage") %>%
  left_join(df_rel_abundance %>% filter(core == "lele")) %>%
  group_by(ka, years) %>%
  summarize(cum_abundance = sum(rel_abundance),
            cum_reads = sum(cum_taxReads)) %>%
  ungroup %>%
  ggplot(aes(x = ka, y = cum_abundance)) +
  geom_area() + 
  ylab("cumulative relative abundance [%]") +
  xlab("cal ka BP") +
  theme_minimal()
ggsave(here("Figures", "lele", "stratigraphic_plots", "Marine_taxa_anylake_lele.pdf"), dpi = 900, width = 10, height = 5) # Supplementary Figure 15

#========================correlations Pelagibacter viruses====================
df_filtered %>%
  filter(core == "KL77",
         grepl("Pelagivirus", taxon) | grepl("Pelagibacter phage", taxon) | taxon == "Stopavirus HTVC011P",
         taxlevel == "species",
         !grepl("unassigned", taxon)) %>% 
  merge(meta[meta$core == "KL77", ], all.y = T) %>%
  filter(years >= 0) %>% 
  dplyr::select(sample, years, ka, core, cum_taxpercent, taxon, epocheII) %>% 
  pivot_wider(id_cols = c("core", "sample", "years", "ka", "epocheII"), names_from = "taxon", values_from = "cum_taxpercent") %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols = -c("core", "sample", "years", "ka", "epocheII"), names_to = "virus", values_to = "cum_taxpercent") %>% 
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon == "Candidatus Pelagibacter sp. IMCC9063") %>% rename(host = cum_taxpercent) %>% select(sample, host), by = "sample") %>%
  group_by(virus, epocheII) %>%
  summarize(cor_test = list(cor.test(cum_taxpercent, host, method = "spearman"))) %>%
  ungroup() %>%
  mutate(
    cor = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  arrange(epocheII, desc(cor)) %>%
  mutate(significant = case_when(p_value < 0.05 ~ T,
                             T ~ F)) %>%
  rename(`geological epoch` = epocheII) %>%
  ggplot(aes(y = virus, x = cor, col = `geological epoch`, shape = `geological epoch`, alpha = significant)) +
  geom_point(size = 3) +
  xlim(-1, 1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlab("spearman correlation coefficient") +
  ylab(NULL) +
  theme_minimal() +
  scale_color_manual(values = c("#800000", "#000075")) +
  scale_alpha_manual(values = c(0.4, 1)) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13))
ggsave(here(here("Figures", "KL77", "vir_hosts", "pelagibacter"), "Correlations_vs_Candidatus Pelagibacter_sp._IMCC9063.pdf"), dpi = 900) # Supplementary Figure 16


df_rel_abundance %>%
  filter(core == "KL77",
         grepl("Pelagivirus", taxon) | grepl("Pelagibacter phage", taxon) | taxon == "Stopavirus HTVC011P",
         taxlevel == "species") %>% 
  merge(meta[meta$core == "KL77", ], all.y = T) %>%
  filter(years >= 0) %>%
  mutate(taxon = "virus") %>%
  group_by(sample, years, ka, core, epocheII, taxon) %>%
  summarize(cum_taxpercent = sum(cum_taxpercent)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("core", "sample", "years", "ka", "epocheII"), names_from = "taxon", values_from = "cum_taxpercent") %>%
  replace(is.na(.), 0) %>% 
  #pivot_longer(cols = -c("core", "sample", "years", "ka", "epocheII"), names_to = "virus", values_to = "cum_taxpercent") %>% 
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon %in% (df_rel_abundance %>% filter(analysis_group == "Pelagibacter", taxon != "Candidatus Pelagibacter sp. IMCC9063") %>% pull(taxon) %>% unique)) %>% group_by(sample) %>% summarize(woIMCC9063 = sum(cum_taxpercent)) %>% ungroup(), by = "sample") %>% 
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon == "Pelagibacter") %>% rename(all = cum_taxpercent) %>% select(sample, all), by = "sample") %>% 
  pivot_longer(cols = c("woIMCC9063", "all"), values_to = "host", names_to = "host_group") %>%
  group_by(epocheII, host_group) %>%
  summarize(cor_test = list(cor.test(virus, host, method = "spearman"))) %>% 
  ungroup() %>%
  mutate(
    cor = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  arrange(epocheII, desc(cor)) %>%
  mutate(significant = case_when(p_value < 0.05 ~ T,
                                 T ~ F)) %>% print
  
#==========================Synechococcus===================================
for (lake in marine) {
  tmp_path <- here("Figures", lake, "vir_hosts", "synechococcus")
  dir.create(tmp_path, recursive = T)
  
  virus <- c(df_vir_hosts$taxon[df_vir_hosts$host_taxon == "Synechococcus"], df_rel_abundance %>% filter(grepl("Synechococcus", taxon), superkingdom == "Viruses")) %>% unique
  
  plot_viruses <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% virus) %>%
    mutate(taxon = "Synechococcus-infecting viruses") %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative virus abundance [%]") +
    scale_fill_manual(values = virusPalette(tmp %>% pull(taxon) %>% unique %>% length)) +
    labs(fill = NULL) +
    theme_minimal() +
    xlab("cal ka BP") +
    geom_vline(xintercept = 13, color = "black", linetype = "dotted")
  
  synechococcus <- df_rel_abundance %>% filter(analysis_group == "Synechococcus") %>% pull(taxon) %>% unique
  
  synechococcus <- df_rel_abundance %>% 
    filter(taxon %in% synechococcus, 
           core == lake) %>%
    group_by(taxon) %>%
    summarize(mean_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(desc(mean_taxpercent)) %>% 
    pull(taxon) %>% 
    unique
  
  plot_host1 <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% synechococcus) %>%
    mutate(taxon = case_when(taxon %in% c("Synechococcus sp. Ace-Pa", "Synechococcus sp. UW69") ~ taxon,
                             TRUE ~ "other Synechococcus")) %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>% 
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent, fill = taxon)) +
    geom_area(position = "stack", stat = "identity") +
    ylab("cumulative host abundance [%]") +
    scale_fill_manual(values = c("Synechococcus sp. Ace-Pa" = "#800000", "Synechococcus sp. UW69" = "#000075", "other Synechococcus" = "#707070")) +
    labs(fill = NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank()) +
    geom_vline(xintercept = 13, color = "black", linetype = "dotted")
  
  plot_host2 <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% synechococcus,
           taxon != "Synechococcus sp. Ace-Pa") %>%
    mutate(taxon = case_when(taxon == "Synechococcus sp. UW69" ~ taxon,
                             TRUE ~ "other Synechococcus")) %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>% 
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp} %>%
    ggplot(aes(x = ka, y = cum_taxpercent, fill = taxon)) +
    geom_area(position = "stack", stat = "identity") +
    ylab("cumulative host abundance [%]") +
    scale_fill_manual(values = c("Synechococcus sp. UW69" = "#000075", "other Synechococcus" = "#707070")) +
    labs(fill = NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank()) +
    geom_vline(xintercept = 13, color = "black", linetype = "dotted")
  
  overview_plot <- plot_grid(
    plot_host1, plot_host2, plot_viruses, ncol = 1,
    rel_heights = c(1,1,1), 
    align = "v") 
  
  ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_virus_vs_host_Synechococcus_wo_Ace-Pa_grey.pdf")), overview_plot, dpi = 900, width = 10, height = 8) # Supplementary Figure 17
  
}

df_rel_abundance %>%
  filter(core == "KL77",
       taxon %in% c(df_vir_hosts$taxon[df_vir_hosts$host_taxon == "Synechococcus"], df_rel_abundance %>% filter(grepl("Synechococcus", taxon), superkingdom == "Viruses") %>% unique),
         taxlevel == "species") %>%
  merge(meta[meta$core == "KL77", ], all.y = T) %>%
  filter(years >= 0) %>%
  mutate(taxon = "virus") %>%
  group_by(sample, years, ka, core, epocheII, taxon) %>%
  summarize(cum_taxpercent = sum(cum_taxpercent)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c("core", "sample", "years", "ka", "epocheII"), names_from = "taxon", values_from = "cum_taxpercent") %>%
  replace(is.na(.), 0) %>% 
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon %in% (df_rel_abundance %>% filter(analysis_group == "Synechococcus", taxon != "Synechococcus sp. Ace-Pa") %>% pull(taxon) %>% unique)) %>% group_by(sample) %>% summarize(woAcePa = sum(cum_taxpercent)) %>% ungroup(), by = "sample") %>%
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon %in% (df_rel_abundance %>% filter(analysis_group == "Synechococcus", taxon != "Synechococcus sp. Ace-Pa", taxon != "Synechococcus sp. UW69") %>% pull(taxon) %>% unique)) %>% group_by(sample) %>% summarize(woAcePaUW69 = sum(cum_taxpercent)) %>% ungroup(), by = "sample") %>%
  left_join(df_rel_abundance %>% filter(core == "KL77", taxon == "Synechococcus") %>% rename(all = cum_taxpercent) %>% select(sample, all), by = "sample") %>% 
  pivot_longer(cols = c("woAcePa", "woAcePaUW69", "all"), values_to = "host", names_to = "host_group") %>%
  group_by(epocheII, host_group) %>%
  summarize(cor_test = list(cor.test(virus, host, method = "spearman"))) %>% 
  ungroup() %>%
  mutate(
    cor = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  arrange(epocheII, desc(cor)) %>%
  mutate(significant = case_when(p_value < 0.05 ~ T,
                                 T ~ F)) %>% print

#========================Bacterial classes===================================
{
  bacterial_host_df <- df_vir_genus_host_ictv %>% select(Genus, Host) %>%
    rbind(df_vir_ncbi_host %>% select(Genus, Host))
  
  for (lake in order_sites) {
    for(host in c("Actinomycetes", "Alphaproteobacteria", "Bacilli",
                 "Betaproteobacteria", "Cyanophyceae", "Flavobacteriia",
                 "Clostridia", "Gammaproteobacteria", "Bacteroidia")) { 
      tmp_path <- here("Figures", lake, "vir_hosts", "bacteria_classes")
      dir.create(tmp_path, recursive = T)
      
      overview_plot <- df_rel_abundance %>%
        filter(core == lake,
               taxon == host | taxon %in% bacterial_host_df$Genus[bacterial_host_df$Host == host]) %>%        
        merge(bacterial_host_df, by.x = "taxon", by.y = "Genus", all.x = T, all.y = F) %>%
        mutate(virhost = case_when(superkingdom == "Viruses" ~ paste("Viruses infecting", host),
                                   superkingdom == "Bacteria" ~ host)) %>%
        mutate(taxon = case_when(superkingdom == "Viruses" ~ Host,
                                 superkingdom == "Bacteria" ~ taxon)) %>%
        group_by(taxon, sample, core, virhost) %>%
        summarize(cum_taxpercent = sum(cum_taxpercent)) %>%
        ungroup %>% 
        merge(meta[meta$core == lake, ], all.y = T) %>%
        filter(years >= 0) %>% 
        dplyr::select(sample, years, ka, core, virhost, cum_taxpercent) %>%
        group_by(ka, virhost) %>%
        summarize(cum_taxpercent = mean(cum_taxpercent, na.rm = T)) %>% 
        ungroup %>%
        complete(ka, virhost, fill = list(cum_taxpercent = 0)) %>%
        filter(virhost != "NA") %>%
        ggplot(aes(x = ka, y = cum_taxpercent)) +
        geom_area(position = "stack", stat = "identity") +
        ylab("cumulative abundance [%]") +
        xlab("cal ka BP") +
        labs(fill = NULL) +
        scale_fill_manual(values = bacteriaPalette) +
        facet_wrap(~virhost, ncol = 1, scales = "free_y") +
        theme_minimal() +
        theme(strip.text = element_text(face = "bold", size = 13),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 13))
      
      
      ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_virus_host_classes_vs_host_Bacteria_classes_", host, ".pdf")), overview_plot, dpi = 300) # Supplementary Figures 19, 20, and 24-30
    }
  }
}

#===========================Pseudomonas=====================================
for (lake in order_sites) {
  tmp_path <- here("Figures", lake, "vir_hosts", "pseudomonas")
  dir.create(tmp_path, recursive = T)
  
  pseudomonas_viruses <- unique(c(df_vir_hosts_cores$taxon[df_vir_hosts_cores$host_taxon == "Pseudomonas"],
                                  df_rel_abundance %>% filter(grepl("Pseudomonas", taxon), !grepl("virophage", taxon), superkingdom == "Viruses", taxlevel == "species") %>% pull(taxon)))
  
  plot_viruses <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% pseudomonas_viruses,
           superkingdom == "Viruses") %>%
    mutate(taxon = "Pseudomonas-infecting viruses") %>%
    dplyr::select(sample, taxon, cum_taxpercent) %>%
    group_by(taxon, sample) %>% 
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>% 
    mutate(taxon = case_when(is.na(taxon) ~ "Pseudomonas-infecting viruses",
                             T ~ taxon),
           cum_taxpercent = case_when(is.na(cum_taxpercent) ~ 0,
                                      T ~ cum_taxpercent)) %>%
    mutate(taxon = "Pseudomonas-infecting viruses") %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp } %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative virus abundance [%]") +
    labs(fill = NULL) +
    xlab("cal ka BP") +
    theme_minimal() +
    scale_fill_manual(values = virusPalette(tmp %>% pull(taxon) %>% unique %>% length))
  
  
  plot_host <- df_rel_abundance %>%
    filter(core == lake,
           taxon == "Pseudomonas") %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>% 
    arrange(core, ka, taxon) %>%
    filter(taxon != "NA") %>%
    {. ->> tmp } %>%
    ggplot(aes(x = ka, y = cum_taxpercent)) +
    geom_area(position = "stack", stat = "identity", fill = "#707070") +
    ylab("cumulative host abundance [%]") +
    labs(fill = NULL) +
    xlab(NULL) +
    theme_minimal() +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank())
  
  
  overview_plot <- plot_grid(
    plot_host, plot_viruses, ncol = 1,
    rel_heights = c(1,1), 
    align = "v")
  
  ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_Pseudomonas_virus_all_vs_host_Pseudomonas_grey.pdf")), overview_plot, dpi = 900, width = 15, height = 10) # Supplementary Figure 23
  
}

lake <- "ilirney"
for (lake in order_sites) {
  tmp_path <- here("Figures", lake, "vir_hosts", "pseudomonas")
  dir.create(tmp_path, recursive = T)
  
  pseudomonas_viruses <- unique(c(df_vir_hosts_cores$taxon[df_vir_hosts_cores$host_taxon == "Pseudomonas"],
                    df_rel_abundance %>% filter(grepl("Pseudomonas", taxon), !grepl("virophage", taxon), superkingdom == "Viruses", taxlevel == "species") %>% pull(taxon)))
  
  plot_viruses <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% pseudomonas_viruses,
           superkingdom == "Viruses") %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    {. ->> tmp} %>% 
    filter(!taxon == "NA") %>%
    ggplot(aes(x = ka, y = cum_taxpercent, fill = taxon)) +
    geom_area(position = "stack", stat = "identity") +
    scale_fill_manual(values = virusPalette(tmp %>% pull(taxon) %>% unique %>% length)) +
    ylab("cumulative virus abundance [%]") +
    labs(fill = NULL) +
    theme_minimal() +
    xlab("cal ka BP") +
    theme(axis.text = element_text(size = 12))
  
  pseudomonas <- df_rel_abundance %>% filter(analysis_group == "Pseudomonas", !grepl("unassigned", taxon)) %>% pull(taxon) %>% unique
  pseudomonas <- df_rel_abundance %>% 
    filter(taxon %in% pseudomonas, core == lake) %>%
    group_by(taxon) %>%
    reframe(mean_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(desc(mean_taxpercent)) %>% 
    pull(taxon) %>% 
    unique
  
  plot_host <- df_rel_abundance %>%
    filter(core == lake,
           taxon %in% pseudomonas) %>%
    mutate(taxon = case_when(taxon %in% pseudomonas[1:7] ~ taxon,
                             TRUE ~ "other Pseudomonas")) %>%
    group_by(taxon, sample, core) %>%
    reframe(cum_taxpercent = sum(cum_taxpercent)) %>%
    ungroup %>%
    merge(meta[meta$core == lake, ], all.y = T) %>%
    filter(years >= 0) %>%
    dplyr::select(taxon, sample, years, ka, core, cum_taxpercent) %>%
    pivot_wider(id_cols = c("core", "sample", "years", "ka"), names_from = "taxon", values_from = "cum_taxpercent") %>%
    replace(is.na(.), 0) %>%
    pivot_longer(cols = -c("core", "sample", "years", "ka"), names_to = "taxon", values_to = "cum_taxpercent") %>%
    group_by(core, sample, years, ka) %>%
    mutate(rel_abundance = cum_taxpercent / sum(cum_taxpercent) * 100) %>%
    ungroup() %>%
    group_by(core, ka, taxon) %>%
    reframe(cum_taxpercent = mean(cum_taxpercent)) %>%
    ungroup() %>%
    arrange(core, ka, taxon) %>%
    {. ->> tmp} %>% 
    ggplot(aes(x = ka, y = cum_taxpercent, fill = taxon)) +
    geom_area(position = "stack", stat = "identity") +
    ylab("cumulative host abundance [%]") +
    labs(fill = NULL) +
    theme_minimal() +
    xlab(NULL) +
    scale_fill_manual(values = hostPalette(tmp %>% pull(taxon) %>% unique %>% length)) +
    theme(axis.text = element_text(size = 10),
          axis.text.x = element_blank())
  
  
  overview_plot <- plot_grid(
    plot_host, plot_viruses, ncol = 1,
    rel_heights = c(1,1), 
    align = "v")
  
  ggsave(here(tmp_path, paste0("Relative_abundance_", lake, "_virus_all_vs_host_pseudomonas.pdf")), overview_plot, dpi = 900, width = 15, height = 10) # Supplementary Figure 22
  
}

