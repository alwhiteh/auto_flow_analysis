######################################

# Author: Alex Whitehead, Ph.D.
# Contact: ajwhitehead@wisc.edu

#!/usr/bin/env Rscript

######################################

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
 stop("Usage of the script is as follows, including these arguments:
      Rscript get_flow_stats_with_ICS.R <path to FlowJo wsp file> 
      <path/to/fcs_files> <output/path/folder> <surface dilution> <count beads? 1/0> <Bead Gate Name> <path/to/ICS.wsp>
      <ICS dilution>
      
      Directories must not contain spaces
      FlowJo wsp needs to have samples arranged in groups - statistical testing will compare groups
      Count beads assumes 50k Bead were added to each sample in 50uL
      Dilution is an number (i.e. 1/4 lung would be 4), undiluted = 1
      
      Note: if any gates are off-scale it may generate problems - !check if FSC3.1 fixes this!
      TODO: add MFI
      
      
      Example:
      Rscript get_flow_stats_with_ICS.R Surface.wsp /CEng2_Resting_Lung_Tet_NR4A1KO_AW_1_5_23 /out 4 1 Beads ICS.wsp 2.666 "
      , call.=FALSE)
}

# Read in arguments from command line
ws_path <- args[1]
fcs_path <- args[2]
output_path <- as.character(args[3])
dilution <- as.numeric(args[4])
count_beads <- as.integer(args[5])
bead_gate <- as.character(args[6])
ics_path <- as.character(args[7])
ics_dilution <- as.numeric(args[8])

# Load Libraries
# List Packages
packages <- c("dplyr","stringr","ggplot2","CytoML","flowWorkspace",
              "parallel","dunn.test","tidyr","purrr","ggpubr","broom",
              "gdata","shiny","plotly","gridlayout","bslib","DT","fuzzyjoin")
# If package is not installed, then install
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Quietly load all packages
invisible(lapply(packages, library, character.only = TRUE))
# suppressMessages(library("dplyr"))
# suppressMessages(library("stringr"))
# suppressMessages(library("ggplot2"))
# suppressMessages(library("CytoML"))
# suppressMessages(library("flowWorkspace"))
# library("parallel")
# library("dplyr")
# library("dunn.test") 
# library("tidyr")
# library("purrr")
# library("ggplot2")
# library("ggpubr")
# library("broom")
#suppressMessages(library("gdata"))
# For ICS 
# library("shiny")
#devtools::install_github("rstudio/gridlayout", force = TRUE)
#options(timeout=4000) # use this if your connection is slow and times out
#devtools::install_github("rstudio/shinyuieditor")
#suppressMessages(library("shinyuieditor"))
# suppressMessages(library("plotly"))
# suppressMessages(library("gridlayout"))
# suppressMessages(library("bslib"))
# suppressMessages(library("DT"))
# suppressMessages(library("fuzzyjoin"))
cores <- as.integer(detectCores())

#Create Cleanup Function to set Theme for Clean Graphs
cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = "black"))


# Read in the data files
ws <- open_flowjo_xml(ws_path)

# Get sample groups and generate gatingset from the worksheet
fj_ws_get_sample_groups(ws) %>% filter(groupName != "All Samples") -> groups
samples <- fj_ws_get_samples(ws)
sample_groups <- inner_join(samples, groups, by = "sampleID")

# Get the gating set without loading FCS files first
gs <- flowjo_to_gatingset(ws,
                           mc.cores = cores,
                           name = 1,
                           path = fcs_path,
                           execute = F) # name = 1 means all samples, execute = F means don't load .fcs files

# Extract the comp matrix
gh <- gs[[1]]
comp_mat <- gh_get_compensations(gh) # pick from the first sample

# Now repeat generating the gatingset with the FCS files 
print("Now Reading in Surface FCS Files, this may take a minute")
# We need to edit this in case someone changed the axes
gs <- flowjo_to_gatingset(ws, 
                          mc.cores = cores, 
                          name = 1, 
                          path = fcs_path,
                          compensation = comp_mat,
                          channel.ignore.case = T,
                          execute = T) # name = 1 means all samples, execute = F means don't load .fcs files
print("Finished reading in Surface FCS Files")

# get the counts from each gate from each sample
populations <- gs_pop_get_count_fast(gs, path = "full", XML = T)
populations$pct <- populations$Count/populations$ParentCount *100
populations$pct[is.nan(populations$pct)] <- 0 # Replace 0/0 NaNs with Zero
populations$name <- gsub("_\\w+$", "", populations$name)

# append the sample group data to each line 
populations <- left_join(populations, sample_groups)
# Remove any samples that are not in a group
populations <- populations[populations$groupName != 'NA']

# add the MFI?
# use the pop.MFI() function


# Calculate the total number of beads
# Extract beads measured
if (count_beads == 1) {
  nbeads <- populations[grepl(paste0("\\", bead_gate, "$"),populations$Population), c("name","Count")]
  nbeads %>% group_by(name) %>% slice(which.min(Count)) -> nbeads
  colnames(nbeads)[2] <-"nbeads"
  populations <- left_join(populations, nbeads)
  # add calculated cells to each row using the count beads and sample_volume
  populations$calc_cells <- ((1/(populations$nbeads/50000))*dilution) * populations$Count
} else {
  populations$calc_cells <- dilution * populations$Count 
}

# Perform stats and ID what populations are different 
if (max(sample_groups$groupID >2)){ #Do Kruskal-Wallis with Dunn post-test and bonferroni correction if >2 groups
  results_tbl_ncells_sig <- populations %>%
    group_by(Population) %>%
    filter(var(Count) != 0) %>% # Remove populations with zero variance
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$calc_cells, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted <= 0.05) %>% arrange(P.adjusted)
  results_tbl_ncells_not_sig <- populations %>% # Get non-significant populations, too
    group_by(Population) %>%
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$calc_cells, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted > 0.05) %>% arrange(P.adjusted)
  
  # Repeat for pct comparisons
  results_tbl_pct_sig <- populations %>% # Get significant comparisons
    group_by(Population) %>%
    filter(var(Count) != 0) %>%
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$pct, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted <= 0.05) %>% arrange(P.adjusted)
  results_tbl_pct_not_sig <- populations %>% # Get non-significant comparisons
    group_by(Population) %>%
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$pct, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted > 0.05) %>% arrange(P.adjusted)
} else{ # Otherwise do Whitney-Mann U with wilcox.test and bonferroni correction
  results_tbl_ncells_sig <- populations %>%
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$calc_cells~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value <= 0.05) %>% arrange(p.value)
  results_tbl_ncells_not_sig <- populations %>% # get non-significant comparisons
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$calc_cells~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value > 0.05) %>% arrange(p.value)

  # Repeat for pct comparisons
  results_tbl_pct_sig <- populations %>%
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$pct~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value <= 0.05) %>% arrange(p.value)
  results_tbl_pct_not_sig <- populations %>% # Get non-significant comparisons
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$pct~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value > 0.05) %>% arrange(p.value)
}


# Generate plots for significant comparisons
filtered_pct_sig <- results_tbl_pct_sig %>%
  distinct(Population, .keep_all = T) # remove duplicate populations
filtered_pct_ns <- results_tbl_pct_not_sig %>%
  distinct(Population, .keep_all = T) # remove duplicate populations

pct_plots <- list() # Generate plots - later would be good to add sig to plot
# for (f in 1:nrow(filtered_pct_sig)){
#   as.data.frame(filtered_pct_sig$data[f]) %>% select(c("groupName","pct")) %>%
#     ggplot(aes(groupName, pct)) +
#     geom_boxplot(outlier.shape = NA) + 
#     geom_jitter() +
#     ggtitle(paste0(sub(".*/", "",filtered_pct_sig[f,1])),
#             subtitle = str_wrap(paste0(filtered_pct_sig$data[[f]]$Parent), 
#                                 width = 80, whitespace_only = F)) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(plot.subtitle=element_text(size=12, hjust=0.5)) +
#     scale_x_discrete(guide = guide_axis(n.dodge=2)) + # prevent overlapping labs
#     ylab("Percentage of Parent") +
#     xlab("Group") +
#     cleanup-> p1
#     pct_plots[[f]] <- p1
# }

# make a list of all possible comparisons
#all_comparisons <- combn(unique(groups$groupName), 2, simplify = F)

for (f in 1:nrow(filtered_pct_sig)){ # For each population, let's make a plot
  results_tbl_pct_sig %>% 
    filter(Population == filtered_pct_sig$Population[f]) %>% 
    ungroup() %>%
    select(comparisons)-> comparisons #Find the sig comparisons for each pop
  plot_comparisons <- list()
  for (i in 1:nrow(comparisons)){ #format the comparisons for stat_compare_means
    plot_comparisons[i] <- str_split(comparisons$comparisons[i]," - ")
  }
  if (max(sample_groups$groupID >2)){
    post_test <- "wilcox.test" 
  } else {post_test <- "t.test"
    }

  as.data.frame(filtered_pct_sig$data[f]) %>% select(c("groupName","pct")) %>%
    ggplot(aes(groupName, pct)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    ggtitle(paste0(sub(".*/", "",filtered_pct_sig[f,1])),
            subtitle = str_wrap(paste0(filtered_pct_sig$data[[f]]$Parent), 
                                width = 80, whitespace_only = F)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5)) +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) + # prevent overlapping labs
    # We need to only feed the comparisons we already know are significant
    stat_compare_means(comparisons = plot_comparisons, hide.ns = F,
                       label = "p.signif",p.adjust.methods = "bonferonni", #Perform bonferroni correction
                       method = post_test) + 
    # FYI, we the p values from wilcox test here are likely to be smaller than from Dunn. 
    # When in doubt, use Dunn test results
    ylab("Percentage of Parent") +
    xlab("Group") +
    cleanup-> p1
  pct_plots[[f]] <- p1
}
pdf(paste0(output_path,"/Surface_Percentage_of_Cells_Sig.pdf"), width = 20, height = 20)
ggarrange(plotlist = pct_plots, ncol = 3, nrow = 4)
dev.off()

# repeat for calculated number of cells comparisons
filtered_ncell_sig <- results_tbl_ncells_sig %>%
  distinct(Population, .keep_all = T) # remove duplicate populations
filtered_ncell_ns <- results_tbl_ncells_not_sig %>%
  distinct(Population, .keep_all = T) # remove duplicate populations

ncell_plots <- list() # Generate plots
for (f in 1:nrow(filtered_ncell_sig)){
  results_tbl_ncells_sig %>% 
    filter(Population == filtered_ncell_sig$Population[f]) %>% 
    ungroup() %>%
    select(comparisons)-> comparisons #Find the sig comparisons for each pop
  plot_comparisons <- list()
  for (i in 1:nrow(comparisons)){ #format the comparisons for stat_compare_means
    plot_comparisons[i] <- str_split(comparisons$comparisons[i]," - ")
  }
  as.data.frame(filtered_ncell_sig$data[f]) %>% select(c("groupName","calc_cells")) %>%
    ggplot(aes(groupName, calc_cells)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    ggtitle(paste0(sub(".*/", "",filtered_ncell_sig[f,1])),
            subtitle = str_wrap(paste0(filtered_ncell_sig$data[[f]]$Parent), 
                                width = 80, whitespace_only = F)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle=element_text(size=10, hjust=0.5)) +
    scale_y_log10(expand = c(0, 0.5), na.value=0) + # adds half a log of buffer for Y axis, makes 0's into 0's
    scale_x_discrete(guide = guide_axis(n.dodge=2)) + # prevent overlapping labs
    stat_compare_means(comparisons = plot_comparisons, hide.ns = F,
                       label = "p.signif",p.adjust.methods = "bonferonni", #Perform bonferroni correction
                       method = post_test) + 
    ylab("Number of Cells") +
    xlab("Group") +
    cleanup -> p1
  ncell_plots[[f]] <- p1
}
pdf(paste0(output_path,"/Surface_Calculated_Cells_Sig.pdf"),width = 20, height = 20)
suppressWarnings(ggarrange(plotlist = ncell_plots, ncol = 3, nrow = 4))
dev.off()

# Send the data to the output folder
# Plot the gating strategy
pdf(paste0(output_path,"/Surface_Gating_Tree.pdf"),width = 20, height = 20)
plot(gs, fontsize = 15)
dev.off()

# Show hierarchy - only works with parsed .fcs files
gh <- gs[[1]]
pdf.options(reset = TRUE, onefile = FALSE) # removes blank page from autoplot
pdf(paste0(output_path,"/Surface_Gating_Strategy.pdf"),width = 20, height = 20)
suppressMessages(autoplot(gh, bins = 128))
dev.off()

# Write output data tables
if (max(sample_groups$groupID >2)){ # Use these columns if you did an ANOVA
write.table(as.data.frame(results_tbl_ncells_sig[,c(1,3:7)]), paste0(output_path,"/Surface_Calculated_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
write.table(as.data.frame(results_tbl_pct_sig[,c(1,3:7)]), paste0(output_path,"/Surface_Pct_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

# Print the populations that are significant for both numbers and percentages
pct_and_ncell_sig <- merge(results_tbl_ncells_sig,results_tbl_pct_sig, by = c("Population", "comparisons"), suffixes = c("_ncells","_pct"))
write.table(as.data.frame(pct_and_ncell_sig[,c(1,2,4:7,9:12)]), paste0(output_path,"/Surface_Pct_and_nCells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

} else{ # Use these columns if you used a Whitney-Mann U with wilcox.test
  write.table(as.data.frame(results_tbl_ncells_sig[,c(1,3:5)]), paste0(output_path,"/Surface_Calculated_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  write.table(as.data.frame(results_tbl_pct_sig[,c(1,3:5)]), paste0(output_path,"/Surface_Pct_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  pct_and_ncell_sig <- merge(results_tbl_ncells_sig,results_tbl_pct_sig, by = c("Population"), suffixes = c("_ncells","_pct"))
  write.table(as.data.frame(pct_and_ncell_sig[,c(1,3:6,8:11)]), paste0(output_path,"/Surface_Pct_and_nCells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

# Now write the raw data for the significant plots into a table in case people want to change plotting settings in Prism
# Do this for pct
pct_raw_data_sig <- list() # Initialize a list that will hold all of our dataframes
for (f in 1:nrow(filtered_pct_sig)){ # For each significant comparison
  as.data.frame(filtered_pct_sig$data[f]) %>% #Extract the raw values
    select(c("groupName","pct")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "pct") %>%
    select(-id) %>% 
    select(-'"groupName"') %>%
    as.data.frame() -> hold
  pct_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    pct_df[[i]] <- na.omit(hold[i])
  }
  pct_df <- do.call("cbindX",pct_df) # unlist into a data frame
  #pct_df <- as.data.frame(do.call("paste",pct_df))
  # Store the dataframe in the list
  pct_raw_data_sig[[f]] <- pct_df
  # Name each item in the list so we know what gate the raw data came from
  names(pct_raw_data_sig)[f] <- paste0(filtered_pct_sig$Population[f]) 
}
# Write the raw data to an excel file
pct_raw_data_sig <- do.call("rbind", pct_raw_data_sig)
write.csv(pct_raw_data_sig, paste0(output_path,"/Surface_Pct_Sig_Raw_Data.csv"), row.names = T)

# Repeat this for non-significant data
pct_raw_data_ns <- list()
for (f in 1:nrow(filtered_pct_ns)){ # For each significant comparison
  as.data.frame(filtered_pct_ns$data[f]) %>% #Extract the raw values
    select(c("groupName","pct")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "pct") %>%
    select(-id) %>% 
    select(-'"groupName"') %>%
    as.data.frame() -> hold
  pct_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    pct_df[[i]] <- na.omit(hold[i])
  }
  pct_df <- do.call("cbindX",pct_df) # unlist into a data frame
  #pct_df <- as.data.frame(do.call("paste",pct_df))
  # Store the dataframe in the list
  pct_raw_data_ns[[f]] <- pct_df
  # Name each item in the list so we know what gate the raw data came from
  names(pct_raw_data_ns)[f] <- paste0(filtered_pct_ns$Population[f]) 
}
# Write the raw data to an excel file
pct_raw_data_ns <- do.call("rbind", pct_raw_data_ns)
write.csv(pct_raw_data_ns, paste0(output_path,"/Surface_Pct_Not_Sig_Raw_Data.csv"), row.names = T)

# Then do this for ncell
ncell_raw_data_sig <- list() 
for (f in 1:nrow(filtered_ncell_sig)){
  as.data.frame(filtered_ncell_sig$data[f]) %>% 
    select(c("groupName","calc_cells")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "calc_cells") %>%
    select(-id) %>% 
    as.data.frame() -> hold
  hold <- hold[,2:ncol(hold)] # remove the "groupName" dummy column
  ncell_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    ncell_df[[i]] <- na.omit(hold[i])
  }
  ncell_df <- do.call("cbindX",ncell_df) # unlist into a data frame
  # Store the dataframe in the list
  ncell_raw_data_sig[[f]] <- ncell_df
  # Name each item in the list so we know what gate the raw data came from
  names(ncell_raw_data_sig)[f] <- paste0(filtered_ncell_sig$Population[f]) 
}
rm(hold) # Remove dummy variables
# Write the raw data to an excel file 
ncell_raw_data_sig <- do.call("rbind", ncell_raw_data_sig)
ncell_raw_data_sig <- round(ncell_raw_data_sig)
write.csv(ncell_raw_data_sig, paste0(output_path,"/Surface_Calculated_Cells_Sig_Raw_Data.csv"), row.names = T)

# Repeat this for non-significant populations
ncell_raw_data_ns <- list() 
for (f in 1:nrow(filtered_ncell_ns)){
  as.data.frame(filtered_ncell_ns$data[f]) %>% 
    select(c("groupName","calc_cells")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "calc_cells") %>%
    select(-id) %>% 
    as.data.frame() -> hold
  hold <- hold[,2:ncol(hold)] # remove the "groupName" dummy column
  ncell_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    ncell_df[[i]] <- na.omit(hold[i])
  }
  ncell_df <- do.call("cbindX",ncell_df) # unlist into a data frame
  # Store the dataframe in the list
  ncell_raw_data_ns[[f]] <- ncell_df
  # Name each item in the list so we know what gate the raw data came from
  names(ncell_raw_data_ns)[f] <- paste0(filtered_ncell_ns$Population[f]) 
}
rm(hold) # Remove dummy variables
# Write the raw data to an excel file 
ncell_raw_data_ns <- do.call("rbind", ncell_raw_data_ns)
ncell_raw_data_ns <- round(ncell_raw_data_ns)
write.csv(ncell_raw_data_ns, paste0(output_path,"/Surface_Calculated_Cells_Not_Sig_Raw_Data.csv"), row.names = T)



############# ICS Analysis ##############

# IF there was no ICS, stop the script here
if(is.na(ics_path == TRUE)){
  print("No ICS worspace provided, exiting script")
  stop()
}
# Read in the ICS files 
ics_ws <- open_flowjo_xml(ics_path)

# Get sample groups and generate gatingset from the worksheet
fj_ws_get_sample_groups(ics_ws) %>% filter(groupName != "All Samples") -> ics_groups
ics_samples <- fj_ws_get_samples(ics_ws)
ics_sample_groups <- inner_join(ics_samples, ics_groups, by = "sampleID")

# Get the gating set without loading FCS files first
ics_gs <- flowjo_to_gatingset(ics_ws,
                          mc.cores = cores,
                          name = 1,
                          path = fcs_path,
                          execute = F) # name = 1 means all samples, execute = F means don't load .fcs files

# Extract the comp matrix
ics_gh <- ics_gs[[1]]
ics_comp_mat <- gh_get_compensations(ics_gh) # pick from the first sample

# Now repeat generating the gatingset with the FCS files 
print("Now Reading in ICS FCS Files, this may take a minute")
# We need to edit this in case someone changed the axes
ics_gs <- flowjo_to_gatingset(ics_ws, 
                          mc.cores = cores, 
                          name = 1, 
                          path = fcs_path,
                          compensation = ics_comp_mat,
                          channel.ignore.case = T,
                          execute = T) # name = 1 means all samples, execute = F means don't load .fcs files
print("Finished reading in ICS FCS Files")

# get the counts from each gate from each sample
ics_populations <- gs_pop_get_count_fast(ics_gs, path = "full", XML = T)
ics_populations$pct <- ics_populations$Count/ics_populations$ParentCount *100
ics_populations$pct[is.nan(ics_populations$pct)] <- 0 # Replace 0/0 NaNs with Zero
ics_populations$name <- gsub("_\\w+$", "", ics_populations$name)

# append the sample group data to each line 
ics_populations <- left_join(ics_populations, ics_sample_groups)
# Remove any samples that are not in a group
ics_populations <- ics_populations[ics_populations$groupName != 'NA']

# Create the variables we will use in the GUI 
ics_file_list <- as.list(unique(ics_populations$name))
surface_file_list <- as.list(unique(populations$name))

if(count_beads == 1){
  # Here we need to import the number of count beads for each sample from the surface samples
  # nbeads contains the samples and number of beads - we need to link each one to an ICS sample
  ui <- grid_page(
    layout = c(
      "header header header",
      "table  table  plotly",
      "table  table  plotly",
      "table  table  plotly"
    ),
    row_sizes = c(
      "100px",
      "1fr",
      "1fr",
      "1fr"
    ),
    col_sizes = c(
      "250px",
      "0.27fr",
      "1.73fr"
    ),
    gap_size = "1rem",
    grid_card_text(
      area = "header",
      content = "Please pair the ICS to Surface Files, Close Window When Done",
      alignment = "center",
      is_title = FALSE
    ),
    grid_card(
      area = "table",
      card_header("ICS Files"),
      card_body(
        radioButtons(
          inputId = "ICS_Radio_Button",
          label = "ICS File",
          choices = ics_file_list,
          width = "100%"
        )
      )
    ),
    grid_card(
      area = "plotly",
      card_header(
        actionButton("Link","Link Files"),
        actionButton(
          inputId = "Delete",
          label = "Delete Selected Link"
        )
      ),
      card_body(
        "Surface Files",
        selectInput(
          inputId = "Surface_Dropdown",
          label = "Select Surface File ",
          choices = surface_file_list
        ),
        DTOutput(outputId = "Linked_Files", width = "100%")
      )
    )
  )
  
  
  server <- function(input, output,session){ # session allows us to pull from the Shiny envt
    # We want to create a "guess" for the pairings of the files to save time
    guess <- stringdist_join(unique(ics_populations$name) %>%
                               as.data.frame %>%
                               set_names("files"),
                             unique(populations$name) %>%
                               as.data.frame %>%
                               set_names("files"),
                             mode = "left",
                             method = "osa", #osa seems to work well
                             max_dist = 99, 
                             distance_col = "dist") %>%
      group_by(files.x) %>%
      slice_min(order_by = dist, n = 1) %>%
      select(-dist)
    
    #guess <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(guess) <- c("ICS_Files","Surface_Files")
    
    values <- reactiveValues() 
    values$data <- guess
    
    output$Linked_Files <- renderDT(
      values$data, options = list(pageLength = 100)
    )
    
    observeEvent(input$Link,{
      new_row=data.frame(
        ICS_Files=as.character(input$ICS_Radio_Button),
        Surface_Files=as.character(input$Surface_Dropdown))
      output$text1 <- renderText({input$ICS_Radio_Button})
      values$data <- rbind(values$data, new_row)
    })
    
    observeEvent(input$Delete,{
      if (!is.null(input$Linked_Files_rows_selected)) {
        values$data <- values$data[-as.numeric(input$Linked_Files_rows_selected),]
      }
    })
    
    # When the session is ended, the table gets sent to our working environment
    session$onSessionEnded(function() {
      linked_files_table <<- isolate(values$data) # The double arrow sends to global env
      stopApp() # This is needed to keep running the code after closing the app
    })
  }
  
  # Run the app
  app <- shinyApp(ui, server)
  print("Running App")
  runApp(app)
  
  # wait for UI to close and define var to continue
  while(exists("linked_files_table") == FALSE){ 
    print("waiting on input")
    Sys.sleep(1)
  }
  #Use the linked_files_table to join the nbeads and ics_populations
  linked_files_table <- merge(linked_files_table, nbeads, by.x ="Surface_Files",by.y ="name", all = T)
  
  # now add a beads column to the populations based on the linked files table
  ics_populations <- left_join(ics_populations, linked_files_table,
                               join_by("name" == "ICS_Files"))
  
  # add calculated cells to each row using the count beads and sample_volume
  # Check the math here
  ics_populations$calc_cells <- (ics_populations$Count/ics_populations$nbeads)*50000*ics_dilution
  print("Cell Counts were Calculated from Count Beads and Dilution")
} else{
  ics_populations$calc_cells <- ics_dilution * ics_populations$Count 
  print("Cell Counts were Calculated from Dilution")
}

# Now we can calculate stats for the ICS files based on pct and ncell
# Perform stats and ID what populations are different 
if (max(ics_sample_groups$groupID >2)){ #Do ANOVA with Kruskal-Wallis if >2 groups
  ics_results_tbl_ncells <- ics_populations %>%
    group_by(Population) %>%
    filter(var(Count) != 0) %>% # Remove populations with zero variance
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$calc_cells, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted <= 0.05) %>% arrange(P.adjusted)
  ics_results_tbl_ncells_ns <- ics_populations %>% # pull the non-significant pops too
    group_by(Population) %>%
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$calc_cells, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted > 0.05) %>% arrange(P.adjusted)
  
  # Repeat for pct comparisons
  ics_results_tbl_pct <- ics_populations %>%
    group_by(Population) %>%
    filter(var(Count) != 0) %>%
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$pct, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted <= 0.05) %>% arrange(P.adjusted)
  ics_results_tbl_pct_ns <- ics_populations %>% # Get non-significant pops too
    group_by(Population) %>%
    nest() %>% 
    mutate(model = map(data, ~dunn.test(x = .x$pct, g = .x$groupName, method = "bonferroni") %>% as_tibble())) %>% 
    unnest(model) %>% filter(P.adjusted > 0.05) %>% arrange(P.adjusted)
  
} else{ # Otherwise do Whitney-Mann U with wilcox.test
  ics_results_tbl_ncells <- ics_populations %>%
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$calc_cells~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value <= 0.05) %>% arrange(p.value)
  ics_results_tbl_ncells_ns <- ics_populations %>% # Get non-significant pops too
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$calc_cells~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value > 0.05) %>% arrange(p.value)
  
  # Repeat for pct comparisons
  ics_results_tbl_pct <- ics_populations %>%
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$pct~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value <= 0.05) %>% arrange(p.value)
  ics_results_tbl_pct_ns <- ics_populations %>% # Get non-significant pops too
    group_by(Population) %>%
    nest() %>%
    mutate(model = map(data, ~wilcox.test(.x$pct~.x$groupName, var.equal = F) %>% tidy())) %>%
    unnest(model) %>% filter(p.value > 0.05) %>% arrange(p.value)
}

# Generate plots for significant ICS comparisons
ics_filtered_pct <- ics_results_tbl_pct %>%
  distinct(Population, .keep_all = T) # remove duplicate populations
ics_filtered_pct_ns <- ics_results_tbl_pct_ns %>%
  distinct(Population, .keep_all = T) # remove duplicate populations

ics_pct_plots <- list() # Generate plots 
for (f in 1:nrow(ics_filtered_pct)){
  ics_results_tbl_pct %>% 
    filter(Population == ics_filtered_pct$Population[f]) %>% 
    ungroup() %>%
    select(comparisons)-> comparisons #Find the sig comparisons for each pop
  plot_comparisons <- list()
  for (i in 1:nrow(comparisons)){ #format the comparisons for stat_compare_means
    plot_comparisons[i] <- str_split(comparisons$comparisons[i]," - ")
  }
  as.data.frame(ics_filtered_pct$data[f]) %>% select(c("groupName","pct")) %>%
    ggplot(aes(groupName, pct)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    ggtitle(paste0(sub(".*/", "",ics_filtered_pct[f,1])),
            subtitle = str_wrap(paste0(ics_filtered_pct$data[[f]]$Parent),
            width = 80, whitespace_only = F)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5)) +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) + # prevent overlapping labs
    stat_compare_means(comparisons = plot_comparisons, hide.ns = F,
                       label = "p.signif",p.adjust.methods = "bonferonni", #Perform bonferroni correction
                       method = post_test) + 
    ylab("Percentage of Parent") +
    xlab("Group") +
    cleanup-> p1
  ics_pct_plots[[f]] <- p1
}
pdf(paste0(output_path,"/ICS_Percentage_of_Cells_Sig.pdf"), width = 20, height = 20)
ggarrange(plotlist = ics_pct_plots, ncol = 3, nrow = 4)
dev.off()

# repeat for calculated number of cells comparisons
ics_filtered_ncell <- ics_results_tbl_ncells %>%
  distinct(Population, .keep_all = T) # remove duplicate populations
ics_filtered_ncell_ns <- ics_results_tbl_ncells_ns %>%
  distinct(Population, .keep_all = T) # remove duplicate populations

ics_ncell_plots <- list() # Generate plots
for (f in 1:nrow(ics_filtered_ncell)){
  ics_results_tbl_ncells %>% 
    filter(Population == ics_filtered_ncell$Population[f]) %>% 
    ungroup() %>%
    select(comparisons)-> comparisons #Find the sig comparisons for each pop
  plot_comparisons <- list()
  for (i in 1:nrow(comparisons)){ #format the comparisons for stat_compare_means
    plot_comparisons[i] <- str_split(comparisons$comparisons[i]," - ")
  }
  as.data.frame(ics_filtered_ncell$data[f]) %>% select(c("groupName","calc_cells")) %>%
    ggplot(aes(groupName, calc_cells)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    ggtitle(paste0(sub(".*/", "",ics_filtered_ncell[f,1])),
            subtitle = str_wrap(paste0(ics_filtered_ncell$data[[f]]$Parent),
                                width = 80, whitespace_only = F)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5)) +
    scale_y_log10(expand = c(0, 0.5), na.value=0) + # adds half a log of buffer for Y axis, makes 0's into 0's
    scale_x_discrete(guide = guide_axis(n.dodge=2)) + # prevent overlapping labs
    stat_compare_means(comparisons = plot_comparisons, hide.ns = F,
                       label = "p.signif", p.adjust.methods = "bonferonni", #Perform bonferroni correction
                       method = post_test) + 
    ylab("Number of Cells") +
    xlab("Group") +
    cleanup -> p1
  ics_ncell_plots[[f]] <- p1
}
pdf(paste0(output_path,"/ICS_Calculated_Cells_Sig.pdf"),width = 20, height = 20)
ggarrange(plotlist = ics_ncell_plots, ncol = 3, nrow = 4)
dev.off()

# Plot the gating strategy
pdf(paste0(output_path,"/ICS_Gating_Tree.pdf"),width = 20, height = 20)
plot(ics_gs, fontsize = 15)
dev.off()

# Show hierarchy - only works with parsed .fcs files
ics_gh <- ics_gs[[1]]
pdf.options(reset = TRUE, onefile = FALSE) # removes blank page from autoplot
pdf(paste0(output_path,"/ICS_Gating_Strategy.pdf"),width = 20, height = 20)
autoplot(ics_gh, bins = 128)
dev.off()

# Write output data tables
if (max(ics_sample_groups$groupID >2)){ # Use these columns if you did an ANOVA
  write.table(as.data.frame(ics_results_tbl_ncells[,c(1,3:7)]), paste0(output_path,"/ICS_Calculated_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  write.table(as.data.frame(ics_results_tbl_pct[,c(1,3:7)]), paste0(output_path,"/ICS_Pct_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  # Print the populations that are significant for both numbers and percentages
  ics_pct_and_ncell_sig <- merge(ics_results_tbl_ncells,ics_results_tbl_pct, by = c("Population", "comparisons"), suffixes = c("_ncells","_pct"))
  write.table(as.data.frame(ics_pct_and_ncell_sig[,c(1,2,4:7,9:12)]), paste0(output_path,"/ICS_Pct_and_nCells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
} else{ # Use these columns if you used a Whitney-Mann U with wilcox.test
  write.table(as.data.frame(ics_results_tbl_ncells[,c(1,3:5)]), paste0(output_path,"/ICS_Calculated_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  write.table(as.data.frame(ics_results_tbl_pct[,c(1,3:5)]), paste0(output_path,"/ICS_Pct_Cells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  ics_pct_and_ncell_sig <- merge(ics_results_tbl_ncells,ics_results_tbl_pct, by = c("Population"), suffixes = c("_ncells","_pct"))
  write.table(as.data.frame(ics_pct_and_ncell_sig[,c(1,3:6,8:11)]), paste0(output_path,"/ICS_Pct_and_nCells_Stats.txt"),  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

# Now write the raw data for the significant plots into a table in case people want to change plotting settings in Prism
# Do this for pct
ics_pct_raw_data <- list() # Initialize a list that will hold all of our dataframes
for (f in 1:nrow(ics_filtered_pct)){ # For each significant comparison
  as.data.frame(ics_filtered_pct$data[f]) %>% #Extract the raw values
    select(c("groupName","pct")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "pct") %>%
    select(-id) %>% 
    as.data.frame() -> hold
  hold <- hold[,2:ncol(hold)] # remove the "groupName" dummy column
  pct_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    pct_df[[i]] <- na.omit(hold[i])
  }
  ics_pct_df <- do.call("cbindX",pct_df) # unlist into a data frame
  # Store the dataframe in the list
  ics_pct_raw_data[[f]] <- ics_pct_df
  # Name each item in the list so we know what gate the raw data came from
  names(ics_pct_raw_data)[f] <- paste0(ics_filtered_pct$Population[f]) 
}
# Write the raw data to an excel file
ics_pct_raw_data <- do.call("rbind", ics_pct_raw_data)
write.csv(ics_pct_raw_data, paste0(output_path,"/ICS_Pct_Sig_Raw_Data.csv"), row.names = T)

# Repeat for non-significant populations
ics_pct_raw_data_ns <- list() # Initialize a list that will hold all of our dataframes
for (f in 1:nrow(ics_filtered_pct_ns)){ # For each significant comparison
  as.data.frame(ics_filtered_pct_ns$data[f]) %>% #Extract the raw values
    select(c("groupName","pct")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "pct") %>%
    select(-id) %>% 
    as.data.frame() -> hold
  hold <- hold[,2:ncol(hold)] # remove the "groupName" dummy column
  pct_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    pct_df[[i]] <- na.omit(hold[i])
  }
  ics_pct_df <- do.call("cbindX",pct_df) # unlist into a data frame
  # Store the dataframe in the list
  ics_pct_raw_data_ns[[f]] <- ics_pct_df
  # Name each item in the list so we know what gate the raw data came from
  names(ics_pct_raw_data_ns)[f] <- paste0(ics_filtered_pct_ns$Population[f]) 
}
# Write the raw data to an excel file
ics_pct_raw_data_ns <- do.call("rbind", ics_pct_raw_data_ns)
write.csv(ics_pct_raw_data_ns, paste0(output_path,"/ICS_Pct_Not_Sig_Raw_Data.csv"), row.names = T)

# Then do this for ncell
ics_ncell_raw_data <- list() 
for (f in 1:nrow(ics_filtered_ncell)){
  as.data.frame(ics_filtered_ncell$data[f]) %>% 
    select(c("groupName","calc_cells")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "calc_cells") %>%
    select(-id) %>% 
    as.data.frame() -> hold
  hold <- hold[,2:ncol(hold)] # remove the "groupName" dummy column
  ics_ncell_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    ics_ncell_df[[i]] <- na.omit(hold[i])
  }
  ics_ncell_df <- do.call("cbindX",ics_ncell_df) # unlist into a data frame
  # Store the dataframe in the list
  ics_ncell_raw_data[[f]] <- ics_ncell_df
  # Name each item in the list so we know what gate the raw data came from
  names(ics_ncell_raw_data)[f] <- paste0(ics_filtered_ncell$Population[f]) 
}
rm(hold) # Remove dummy variables
# Write the raw data to an excel file
ics_ncell_raw_data <- do.call("rbind", ics_ncell_raw_data)
ics_ncell_raw_data <- round(ics_ncell_raw_data)
write.csv(ics_ncell_raw_data, paste0(output_path,"/ICS_Calculated_Cells_Sig_Raw_Data.csv"), row.names = T)

# Repeat this for negative populations
ics_ncell_raw_data_ns <- list() 
for (f in 1:nrow(ics_filtered_ncell_ns)){
  as.data.frame(ics_filtered_ncell_ns$data[f]) %>% 
    select(c("groupName","calc_cells")) %>%
    group_by("groupName") %>%
    mutate(id = row_number()) %>% 
    ungroup() %>%
    pivot_wider(names_from = "groupName", 
                values_from = "calc_cells") %>%
    select(-id) %>% 
    as.data.frame() -> hold
  hold <- hold[,2:ncol(hold)] # remove the "groupName" dummy column
  ics_ncell_df <- list() # initiate a list
  for (i in 1: ncol(hold)){ # Remove the NAs in each column individually
    ics_ncell_df[[i]] <- na.omit(hold[i])
  }
  ics_ncell_df <- do.call("cbindX",ics_ncell_df) # unlist into a data frame
  # Store the dataframe in the list
  ics_ncell_raw_data_ns[[f]] <- ics_ncell_df
  # Name each item in the list so we know what gate the raw data came from
  names(ics_ncell_raw_data_ns)[f] <- paste0(ics_filtered_ncell_ns$Population[f]) 
}
rm(hold) # Remove dummy variables
# Write the raw data to an excel file
ics_ncell_raw_data_ns <- do.call("rbind", ics_ncell_raw_data_ns)
ics_ncell_raw_data_ns <- round(ics_ncell_raw_data_ns)
write.csv(ics_ncell_raw_data_ns, paste0(output_path,"/ICS_Calculated_Cells_Not_Sig_Raw_Data.csv"), row.names = T)


# Finished!
print("Finished!")