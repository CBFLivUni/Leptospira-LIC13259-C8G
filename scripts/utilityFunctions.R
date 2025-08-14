# Function to process log files
extractQCData <-function(input_files){
  
  # Read in and process data
  log_data <- input_files %>% 
    map(read_delim, col_names = c("date", "time", "model", "pLDDT", "pTM", "blank1", "blank2", "sequencing_depth"))
  
  names(log_data) <- basename(input_files)
  
  # Remove rows of jobs that have failed and
  # create vector of failed jobs
  failed_jobs <- names(log_data[sapply(log_data, ncol) < 8])
  failed_jobs <- gsub(".txt", "", failed_jobs)
  log_data <- log_data[sapply(log_data, ncol)>=8] 
  
  # Bind data
  log_data <- log_data %>% rbindlist(idcol = TRUE) 
  
  # Format columns
  log_data <- log_data %>% select(!(blank1:blank2))
  log_data$ID <- gsub(".txt", "", log_data$.id)
  log_data$pLDDT <- as.numeric(gsub("pLDDT=", "", log_data$pLDDT))
  log_data$pTM <- as.numeric(gsub("pTM=", "", log_data$pTM))
  log_data$sequencing_depth <- as.numeric(log_data$sequencing_depth)
  
  # Rename old ID column to be file name instead
  log_data <- log_data %>% dplyr::rename("file_name" = ".id")
  
  return(list("logData" = log_data, "failedJobs" = failed_jobs))
  
}


# Function to process structural alignment files
extractAlignmentData <-function(input_files, col_names = c("query", "target", "alntmscore", "lddt", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits")){
  
  # Read in and process data
  sa_data <- input_files %>% 
    map(read_delim, col_names = col_names)
  
  names(sa_data) <- basename(input_files)
  
  # Remove (failed) items that have more columns than specified
  sa_data <- sa_data[sapply(sa_data, ncol)==length(col_names)] 
  
  # Bind data
  sa_data <- sa_data %>% rbindlist(idcol = FALSE) 
  
  
  # Remove '.pdb' from query and target column
  sa_data$query <- gsub(".pdb", "", sa_data$query)
  sa_data$target <- gsub(".pdb", "", sa_data$target)
  
  return(sa_data)
  
}



# Function to process database queries 
processDBQuery <-function(db_hits_file_path = "", protein_anno=protein_anno, col_names = c("query", "target", "alntmscore", "lddt", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "taxid", "taxname", "taxlineage"), e.value = 0.001){
  
  # Get list of files from file path
  db_hits_input_files <- list.files(db_hits_file_path, pattern = ".m8", full.names = TRUE) 
  
  # Run the function to extract the alignment data from the file list
  db_hits_data <- extractAlignmentData(db_hits_input_files, col_names = col_names)
  
  # Filter hits to only include fungal species
  db_hits_data <- db_hits_data %>% filter(grepl('Fungi', taxlineage))
  
  # Join to protein annotation to include query effector type, effector score and effector pLDDT
  db_hits_data_anno <- db_hits_data %>% 
    dplyr::inner_join(protein_anno[,c("ID", "effector_type_score", "prediction")], join_by(query == ID)) %>%
    dplyr::rename("query_prediction" = "prediction", "query_effector_type_score" = "effector_type_score") %>% 
    arrange(desc(alntmscore))
  
  # Create filtered version of dataframe that only has
  # TM score >0.5, evalue < 0.001
  db_hits_data_filter <- db_hits_data_anno %>%
    filter(alntmscore >= 0.5 & evalue < e.value) 
  
  return(list(hitsDF = db_hits_data_anno, hitsDF.filtered = db_hits_data_filter))
  
}



# Function to process database queries for effectors
globalSequenceAlignment <-function(seq1 = df$x, seq2 = df$y, local=F){
  
  data(BLOSUM62)
  
  if(local){
    
    global_align <- pairwiseAlignment(as.character(seq1),
                                      as.character(seq2),
                                      substitutionMatrix = BLOSUM62, 
                                      gapOpening = -2,
                                      gapExtension = -8, 
                                      scoreOnly = FALSE,
                                      type="local")
    
  } else {
    
    global_align <- pairwiseAlignment(as.character(seq1),
                                      as.character(seq2),
                                      substitutionMatrix = BLOSUM62, 
                                      gapOpening = -2,
                                      gapExtension = -8, 
                                      scoreOnly = FALSE)
    
    
  }
  
  
  sequence_identity <- pid(global_align)/100
  return(sequence_identity)
  
}

# Function to process log files
extractAlphaPulldown <-function(input_files){
  #, receptor="plasminogen"
  
  # Read in and process data
  pulldown_data <- input_files %>% 
    map(read_delim, col_names = c(1:8, "model", "ptm_iptm"))
  
  names(pulldown_data) <- basename(input_files)
  
  # Bind data
  pulldown_data <- pulldown_data %>% rbindlist(idcol = TRUE) 
  
  # Format columns
  pulldown_data <- pulldown_data %>% select(c(.id, model, ptm_iptm))
  pulldown_data$ID <- gsub(".txt", "", pulldown_data$.id)
  
  # Rename old ID column to be file name instead
  pulldown_data <- pulldown_data %>% dplyr::rename("file_name" = ".id")
  
  # Add columns for plasminogen species and leptospira species
  pulldown_data$lepto_species <- str_extract(pulldown_data$ID, "Leptospira_biflexa|Leptospira_borgpetersenii|Leptospira_interrogans_s.Copenhageni|Leptospira_interrogans|Leptospira_kmetyi|Leptospira_wolffii|Leptospira_ilyithenensis")
  pulldown_data$lepto_species <- gsub("Leptospira_", "L. ", pulldown_data$lepto_species)
  pulldown_data$lepto_species <- factor(pulldown_data$lepto_species, ordered = TRUE, levels=c("L. borgpetersenii", "L. interrogans", "L. biflexa", "L. ilyithenensis", "L. kmetyi", "L. wolffii", "L. interrogans_s.Copenhageni"))
  
  pulldown_data$receptor_species <- str_extract(pulldown_data$ID, "Bos_taurus|Canis_lupus|Homo_sapiens|Macaca_mulatta|Mesocricetus_auratus|Mus_musculus|Rattus_norvegicus|Sus_scrofa")
  pulldown_data$receptor_species <- gsub("_", " ", pulldown_data$receptor_species)
  
  #receptor_column <- paste0(receptor, "_species")
  #pulldown_data <- pulldown_data %>% dplyr::rename(!!receptor_column := "receptor_species")
  
  return(pulldown_data)
  
}



# Function to extract scores from each line
extractColabFold <- function(file) {
  line <- read_lines(file, n_max = 1)
  actifptm <- as.numeric(str_extract(line, '(?<="actifptm": )\\d+\\.\\d+'))
  ptm      <- as.numeric(str_extract(line, '(?<="ptm": )\\d+\\.\\d+'))
  iptm     <- as.numeric(str_extract(line, '(?<="iptm": )\\d+\\.\\d+'))
  
  tibble(
    file = file,
    actifptm = actifptm,
    ptm = ptm,
    iptm = iptm,
    ptm_iptm = 0.2 * ptm + 0.8 * iptm
  )
}



# Function to annotate dataframe
annotate_df <- function(df, id_col = "ID") {
  df %>%
    mutate(
      Subunit = case_when(
        grepl("C8A", .data[[id_col]]) ~ "C8A",
        grepl("C8B", .data[[id_col]]) ~ "C8B",
        grepl("C8G", .data[[id_col]]) ~ "C8G",
        TRUE ~ "Unknown"
      ),
      
      # Extract and clean Leptospira species (case-insensitive)
      lepto_species = str_extract(
        str_to_lower(.data[[id_col]]),
        "leptospira_biflexa|leptospira_borgpetersenii|leptospira_interrogans_s.copenhageni|leptospira_interrogans|leptospira_kmetyi|leptospira_wolfii|leptospira_ilyithenensis"
      ) %>%
        str_replace_all("leptospira_", "L. ") %>%
        str_replace("l\\. interrogans_s\\.copenhageni", "L. interrogans_s.Copenhageni") %>%
        str_to_title() %>%
        factor(
          levels = c("L. Borgpetersenii", "L. Interrogans", "L. Biflexa", 
                     "L. Ilyithenensis", "L. Kmetyi", "L. Wolfii", 
                     "L. Interrogans_s.Copenhageni"),
          ordered = TRUE
        ),
      
      c8g_type = case_when(
        grepl("trunucated", .data[[id_col]], ignore.case = TRUE) ~ "Truncated C8G",
        TRUE ~ "Full length C8G"
        )
      
      )
}


# Function to extract and format interface data
extractInterfaces <-function(interface_input_files, quantile=0.9){
  #, receptor="plasminogen"
  
  # Read in and process data
  interface_data <- interface_input_files %>%
    map(read.csv) %>%
    map(~ filter(.x, Total.Count >= quantile(.x$Total.Count, quantile)))
  
  names(interface_data) <- basename(interface_input_files)
  
  # Bind data
  interface_data <- interface_data %>% rbindlist(idcol = TRUE) 
  
  # Format columns
  interface_data$ID <- gsub(".csv", "", interface_data$.id)
  
  # Rename old ID column to be file name instead
  interface_data <- interface_data %>% dplyr::rename("file_name" = ".id", "Residue" = "Residue.1", "Count" = "Total.Count")
  
  # Add columns for plasminogen species and leptospira species
  interface_data$lepto_species <- str_extract(interface_data$ID, "Leptospira_biflexa|Leptospira_borgpetersenii|Leptospira_interrogans|Leptospira_kmetyi|Leptospira_wolffii|Leptospira_ilyithenensis")
  interface_data$lepto_species <- gsub("_", " ", interface_data$lepto_species)
  interface_data$lepto_species <- factor(interface_data$lepto_species, ordered = TRUE, levels=c("Leptospira biflexa", "Leptospira ilyithenensis", "Leptospira kmetyi", "Leptospira wolffii", "Leptospira borgpetersenii", "Leptospira interrogans"))
  
  interface_data$receptor_species <- str_extract(interface_data$ID, "Bos_taurus|Canis_lupus|Homo_sapiens|Macaca_mulatta|Mesocricetus_auratus|Mus_musculus|Rattus_norvegicus|Sus_scrofa")
  interface_data$receptor_species <- gsub("_", " ", interface_data$receptor_species)
  
  #receptor_column <- paste0(receptor, "_species")
  #interface_data <- interface_data %>% dplyr::rename(!!receptor_column := "receptor_species")
  
  interface_data$Residue <- gsub(" A", "", interface_data$Residue)
  interface_data$Residue_type <- gsub(" .*", "", interface_data$Residue)
  
  return(interface_data)
  
}



# Function to extract interface binding-pair data
extractResiduePairs <-function(residue_pair_files){
  
  # Read in and process data
  residue_pair_data <- residue_pair_files %>%
    map(read.csv) %>%
    map(~ {
      # Calculate the total count for each residue
      total_counts <- .x %>%
        group_by(Residue.2) %>%
        summarise(Candidate_residue_count = sum(Count)) %>%
        ungroup()
      
      # Merge the total counts back to the original data
      .x <- .x %>%
        left_join(total_counts, by = "Residue.2")
      
      # Filter to keep only rows with Count in the upper 0.9 quantile
      #.x <- .x %>%
      #  filter(Count >= quantile(Count, 0.9))
      
      return(.x)
    })
  
  names(residue_pair_data) <- basename(residue_pair_files)
  
  # Bind data
  residue_pair_data <- rbindlist(residue_pair_data, idcol = TRUE) 
  
  # Format columns
  residue_pair_data$ID <- gsub(".csv", "", residue_pair_data$.id)
  
  # Rename old ID column to be file name instead
  residue_pair_data <- residue_pair_data %>% dplyr::rename("file_name" = ".id", "Bait_residue" = "Residue.1", "Candidate_residue" = "Residue.2")
  
  # Add columns for plasminogen species and leptospira species
  residue_pair_data$lepto_species <- str_extract(residue_pair_data$ID, "Leptospira_biflexa|Leptospira_borgpetersenii|Leptospira_interrogans|Leptospira_kmetyi|Leptospira_wolffii|Leptospira_ilyithenensis")
  residue_pair_data$lepto_species <- gsub("_", " ", residue_pair_data$lepto_species)
  residue_pair_data$lepto_species <- factor(residue_pair_data$lepto_species, ordered = TRUE, levels=c("Leptospira borgpetersenii", "Leptospira interrogans", "Leptospira biflexa", "Leptospira ilyithenensis", "Leptospira kmetyi", "Leptospira wolffii"))
  
  residue_pair_data$receptor_species <- str_extract(residue_pair_data$ID, "Bos_taurus|Canis_lupus|Homo_sapiens|Macaca_mulatta|Mesocricetus_auratus|Mus_musculus|Rattus_norvegicus|Sus_scrofa")
  residue_pair_data$receptor_species <- gsub("_", " ", residue_pair_data$receptor_species)
  
  #receptor_column <- paste0(receptor, "_species")
  #residue_pair_data <- residue_pair_data %>% dplyr::rename(!!receptor_column := "receptor_species")
  
  residue_pair_data$Candidate_residue_sequence_position <- stringr::str_extract(residue_pair_data$Candidate_residue, "\\d+")
  residue_pair_data$Candidate_residue_sequence_position <- as.numeric(residue_pair_data$Candidate_residue_sequence_position)
  
  return(residue_pair_data)
  
}

# Plot density plots
candidateDensityPlot <- function(residue_pair_data_format, lepto.species="", receptor.species="", residue.count=1250){
  
  residue_pair_data_species <- residue_pair_data_format %>%
    filter(grepl(lepto.species, lepto_species) & grepl(receptor.species, receptor_species))
  
  candidate_binding_plot <- ggplot(residue_pair_data_species, aes(x=Candidate_residue_sequence_position, fill=lepto_species)) +
    geom_histogram(aes(y = ..density..), binwidth = 1, alpha = 0.8) +
    geom_density(color = "black", size = 0.02, alpha = 0.2) +
    geom_rug(sides = "b", length = unit(0.02, "npc"), alpha = 0.2) +
    theme_minimal() + 
    scale_fill_manual(values = project_palettes$nick_species) + 
    theme(legend.position="none") + 
    ylab("Density") + 
    xlab("Residue #") +
    xlim(0, residue.count)
  
  candidate_binding_plot
  
  return(candidate_binding_plot)
  
} 
  
  

# Code adapted from multinichenetR circos plot
# https://github.com/saeyslab/multinichenetr/blob/main/R/plotting.R
produceCircosPlot = function(circos_df, colors_sender, colors_receiver, title=""){
  
  require("dplyr")
  require("ggplot2")
  require("circlize")
  
  
  # Check if grouping column is in pheno
  required_cols <- c("ligand", "receptor", "weight", "sender", "receiver")
  
  if (!(all(required_cols %in% colnames(circos_df)))) {
    stop("The provided dataframe is missing some or all of the required columns. It must contain the columns 'ligand', 'receptor', 'weight', 'sender', receiver'.")
  }
  
  # Check weight column is numeric 
  if (!is.numeric(circos_df$weight)) {
    stop("The provided dataframe is missing some or all of the required columns. It must contain the columns 'ligand', 'receptor', 'weight', 'sender', receiver'.")
  }
  
  # Link each cell type to a color
  grid_col_tbl_ligand = tibble::tibble(sender = colors_sender %>% names(), color_ligand_type = colors_sender)
  grid_col_tbl_receptor = tibble::tibble(receiver = colors_receiver %>% names(), color_receptor_type = colors_receiver)
  
  
  # Make the plot for condition of interest - title of the plot
  title = title
  
  # Create temp df
  df = circos_df
  
  
  # Each pair of ligand-receptors needs to be unique for circos plot to work
  # Code to make each pair unique by adding spaces after name
  ligand.uni = unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i = df[df$ligand == ligand.uni[i], ]
    sender.uni = unique(df.i$sender)
    for (j in 1:length(sender.uni)) {
      df.i.j = df.i[df.i$sender == sender.uni[j], ]
      df.i.j$ligand = paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] = df.i.j$ligand
    }
  }
  receptor.uni = unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i = df[df$receptor == receptor.uni[i], ]
    receiver.uni = unique(df.i$receiver)
    for (j in 1:length(receiver.uni)) {
      df.i.j = df.i[df.i$receiver == receiver.uni[j], ]
      df.i.j$receptor = paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] = df.i.j$receptor
    }
  }
  
  intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  
  while(length(intersecting_ligands_receptors) > 0){
    df_unique = df %>% dplyr::filter(!receptor %in% intersecting_ligands_receptors)
    df_duplicated = df %>% dplyr::filter(receptor %in% intersecting_ligands_receptors)
    df_duplicated = df_duplicated %>% dplyr::mutate(receptor = paste(" ",receptor, sep = ""))
    df = dplyr::bind_rows(df_unique, df_duplicated)
    intersecting_ligands_receptors = generics::intersect(unique(df$ligand),unique(df$receptor))
  }
  
  circos_df = df
  
  # Link ligands/Receptors to the colors of senders/receivers
  circos_df = circos_df %>% dplyr::inner_join(grid_col_tbl_ligand) %>% dplyr::inner_join(grid_col_tbl_receptor)
  links_circle = circos_df %>% dplyr::distinct(ligand,receptor, weight)
  ligand_color = circos_df %>% dplyr::distinct(ligand,color_ligand_type)
  grid_ligand_color = ligand_color$color_ligand_type %>% magrittr::set_names(ligand_color$ligand)
  receptor_color = circos_df %>% dplyr::distinct(receptor,color_receptor_type)
  grid_receptor_color = receptor_color$color_receptor_type %>% magrittr::set_names(receptor_color$receptor)
  grid_col =c(grid_ligand_color,grid_receptor_color)
  
  # give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
  transparency = circos_df %>% dplyr::mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% dplyr::mutate(transparency = 1-(weight*0.8)) %>% .$transparency
  
  # Define order of the ligands and receptors and the gaps
  ligand_order = circos_df$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    ligands = circos_df %>% dplyr::filter(sender == sender_oi) %>%  dplyr::arrange(ligand) %>% dplyr::distinct(ligand)
  }) %>% unlist()
  
  receptor_order = circos_df$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    receptors = circos_df %>% dplyr::filter(receiver == receiver_oi) %>%  dplyr::arrange(receptor) %>% dplyr::distinct(receptor)
  }) %>% unlist()
  
  order = c(ligand_order,receptor_order)
  
  width_same_cell_same_ligand_type = 0.275
  width_different_cell = 3
  width_ligand_receptor = 9
  width_same_cell_same_receptor_type = 0.275
  
  sender_gaps = circos_df$sender %>% unique() %>% sort() %>% lapply(function(sender_oi){
    sector = rep(width_same_cell_same_ligand_type, times = (circos_df %>% dplyr::filter(sender == sender_oi) %>% dplyr::distinct(ligand) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  sender_gaps = sender_gaps[-length(sender_gaps)]
  
  receiver_gaps = circos_df$receiver %>% unique() %>% sort() %>% lapply(function(receiver_oi){
    sector = rep(width_same_cell_same_receptor_type, times = (circos_df %>% dplyr::filter(receiver == receiver_oi) %>% dplyr::distinct(receptor) %>% nrow() -1))
    gap = width_different_cell
    return(c(sector,gap))
  }) %>% unlist()
  receiver_gaps = receiver_gaps[-length(receiver_gaps)]
  
  gaps = c(sender_gaps, width_ligand_receptor, receiver_gaps, width_ligand_receptor)
  
  if(length(gaps) != length(union(circos_df$ligand, circos_df$receptor) %>% unique())){
    warning("Specified gaps have different length than combined total of ligands and receptors - This is probably due to duplicates in ligand-receptor names")
  }
  
  links_circle$weight[links_circle$weight == 0] = 0.01
  circos.clear()
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle,
               directional = 1,
               order=order,
               link.sort = TRUE,
               link.decreasing = TRUE,
               grid.col = grid_col,
               transparency = transparency,
               diffHeight = 0.0075,
               direction.type = c("diffHeight", "arrows"),
               link.visible = links_circle$weight > 0.01,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.175),
               link.arr.length = 0.05, link.arr.type = "big.arrow",  link.lwd = 1.25, link.lty = 1,
               reduce = 0,
               scale = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  }, bg.border = NA) #
  
  title(title)
  p_circos = recordPlot()
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  # grid_col_all = c(colors_receiver, colors_sender)
  legend = ComplexHeatmap::Legend(at = circos_df$receiver %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_receiver[circos_df$receiver %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Receiver")
  ComplexHeatmap::draw(legend, just = c("left", "bottom"))
  
  legend = ComplexHeatmap::Legend(at = circos_df$sender %>% unique() %>% sort(),
                                  type = "grid",
                                  legend_gp = grid::gpar(fill = colors_sender[circos_df$sender %>% unique() %>% sort()]),
                                  title_position = "topleft",
                                  title = "Sender")
  ComplexHeatmap::draw(legend, just = c("left", "top"))
  
  p_legend = grDevices::recordPlot()
  
  circos_plot <- list(circos = p_circos,
                      legend = p_legend)
  
  return(circos_plot)
}


# Function to annotate amino acid with its basic properties from amino acid sequence
# Sequence must be in the format LYS 16, TYR 72 etc.

# Lookup table for amino acid properties
amino_acid_properties <- list(
  "ALA" = "Hydrophobic",
  "VAL" = "Hydrophobic",
  "ILE" = "Hydrophobic",
  "LEU" = "Hydrophobic",
  "MET" = "Hydrophobic",
  "PHE" = "Hydrophobic",
  "TYR" = "Hydrophobic",
  "TRP" = "Hydrophobic",
  "ASP" = "Electrically charged",
  "GLU" = "Electrically charged",
  "LYS" = "Electrically charged",
  "ARG" = "Electrically charged",
  "HIS" = "Electrically charged",
  "ASN" = "Polar uncharged",
  "GLN" = "Polar uncharged",
  "SER" = "Polar uncharged",
  "THR" = "Polar uncharged",
  "CYS" = "Special case",
  "GLY" = "Special case",
  "PRO" = "Special case"
)


annotateAminoAcid <- function(residue) {
  # Extract the amino acid code (first 3 characters)
  amino_acid <- substr(residue, 1, 3)
  
  # Lookup the property in the list
  return(amino_acid_properties[[amino_acid]])
}



# Function to process multimer clustering
extractMultimerClusters <-function(input_files){
  
  # Read in and process data
  multimer_data <- input_files %>% 
    map(read_delim, col_names = c("cluster_name", "cluster_member"))
  
  names(multimer_data) <- basename(input_files)
  
  return(multimer_data)
  
}


# Function to create networks for different species
generate_cluster_network <- function(distance_data, 
                                     TM_norm_threshold = 0.1,
                                     query_cluster_size_threshold = 3,
                                     target_cluster_size_threshold = 3,
                                     host_species = "none",
                                     file_name = "cluster_network.png",
                                     custom_colors = NULL) {
  
  # Filter Data based on thresholds
  distance_data <- distance_data %>%
    filter(
      TM_norm_query < TM_norm_threshold,
      query_cluster_size > query_cluster_size_threshold,
      target_cluster_size > target_cluster_size_threshold
    )
  
  # Filter by host species if specified
  if (host_species != "none") {
    distance_data <- distance_data %>%
      filter(
        query_host_species == host_species & 
        target_host_species == host_species
      )
  }

  # Create node and edge data
  nodes <- data.frame(
    name = unique(c(distance_data$query, distance_data$target)),
    stringsAsFactors = FALSE
  )

  # Annotate nodes with species and cluster size
  nodes <- nodes %>%
    rowwise() %>%
    mutate(
      species = {
        match_idx <- match(name, c(distance_data$query, distance_data$target))
        c(distance_data$query_lepto_species, distance_data$target_lepto_species)[match_idx]
      },
      size = {
        match_idx <- match(name, c(distance_data$query, distance_data$target))
        c(distance_data$query_cluster_size, distance_data$target_cluster_size)[match_idx]
      }
    )

  edges <- distance_data %>%
    select(from = query, to = target, weight = TM_norm_query)

  # Create tidygraph object
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

  # Define species colors
  species_colors <- if (is.null(custom_colors)) {
    project_palettes$nick_species
  } else {
    custom_colors
  }

  # Create ggraph plot
  p <- ggraph(graph, layout = "fr") +  # Fruchterman-Reingold layout
    geom_edge_link(aes(width = 1 - weight), alpha = 0.4, color = "grey40") +
    geom_node_point(aes(size = size, color = species)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_edge_width_continuous(range=c(0.25, 1.5), name = "Similarity (TM score)") +
    scale_size(range = c(3, 10), name = "Cluster Size") +
    scale_color_manual(values = species_colors, name = "Species") +
    theme_void(base_size=16) +
    theme(legend.margin = margin(t = 5, r = 10, b = 5, l = 5))

  # Save plot
  ggsave(file_name, plot = p, width = 12, height = 8, dpi = 300)

  # Return the graph object
  return(graph)
}
