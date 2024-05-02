################################################################################
##
## OP-TESTS
## Analyze data from smolt migration test
##
## Author: 
## Eneko Aspillaga (aspillaga@imedea.uib-csic.es)
## Instituto Mediterraneo de Estudios Avanzados (IMEDEA, CSIC)
##
## Last update: November 3, 2023
##
################################################################################

# Load libraries
library(data.table)
library(survival)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# 1. LOAD DATA =================================================================

# Load receiver locations
receivers <- fread("./data/smolt_migration/smolt_receiver_deployments.csv")
# Load fish metadata
fish_ref <- fread("./data/smolt_migration/smolt_fish_metadata.csv")
# Load detection data
detect <- fread("./data/smolt_migration/smolt_detection_data.csv.gz")

# Assign the array id to each detection
detect[, array := receivers$array[match(receiver_id, receivers$receiver_id)]]

# Assign species, release location and tag protocol to each detection
indx_f <- match(detect$acoustic_tag_id, fish_ref$acoustic_tag_id)
detect[, species := fish_ref$scientific_name[indx_f]]
detect[, release_location := fish_ref$release_location[indx_f]]
detect[, protocol := fish_ref$protocol[indx_f]]

# Convert acoustic_tag_id to factor to include all fish id as levels
detect[, acoustic_tag_id := factor(acoustic_tag_id, 
                                   levels = fish_ref$acoustic_tag_id)]
length(unique(detect$acoustic_tag_id))
length(levels(detect$acoustic_tag_id))
nrow(fish_ref)

# Count number of released smolts per location and species
n_smolts <- as.data.table(table(fish_ref$release_location, 
                                fish_ref$scientific_name,
                                fish_ref$protocol))
colnames(n_smolts) <- c("release_location", "species", "protocol", "n")
n_smolts                          

# Plot receiver locations
plot(latitude ~ longitude, data = receivers, col = factor(array),
     xlim = range(c(receivers$longitude, fish_ref$release_longitude)),
     ylim = range(c(receivers$latitude, fish_ref$release_latitude)),
     pch = 16)
# Add fish release locations
points(release_latitude ~ release_longitude, data = fish_ref, 
       col = factor(release_location), pch = 2, cex = 0.5)


# 2. CALCULATE NUMBER OF SMOLTS DETECTED IN EACH ARRAY =========================

# Filter detections. Remove single detections in each receivers
det_tab <- data.frame(table(detect$acoustic_tag_id, detect$receiver_id))
rm_det <- det_tab[det_tab$Freq <= 1, ]
rm_indx <- !(paste(detect$acoustic_tag_id, detect$receiver_id) %in%
               paste(rm_det[, 1], rm_det[, 2]))
detect <- detect[rm_indx, ]

# Total number of detected smolts by species
detect[, .(n = length(unique(acoustic_tag_id))), by = species]


# Calculate number of smolts detected in each array
smolt_det <- detect[, .(n_fish = length(unique(acoustic_tag_id)),
                        n_detections = sum(length(date_time))), 
                    by = list(species, release_location, array, protocol)]
smolt_det <- smolt_det[order(species, release_location, array)]

# For fish released in location R2, remove detections in upstream arrays (A1 
# and A2)
smolt_det <- smolt_det[!(release_location == "R2" & array %in% c("A1", "A2"))]

smolt_det


# 3. ANALYSIS OF MIGRATION SUCCESS =============================================

# Migration success will be analyzed by applying a survival analysis method.
# For the analysis, fish will be considered as "dead" (not having migrated) 
# in the first downstream array in which they were not longer detected.

arrays <- sort(unique(detect$array))
species <- c("Salmo trutta", "Salmo salar")

# Extract the arrays in which each fish was last detected
migration <- tapply(detect$array, detect$acoustic_tag_id, function(x) max(x))
migration <- data.table(tag_id = names(migration), 
                        last_array = as.character(migration))

# Add fish metadata
indx_m <- match(migration$tag_id, fish_ref$acoustic_tag_id)
migration[, species := fish_ref$scientific_name[indx_m]]
migration[, release_location := fish_ref$release_location[indx_m]]
migration[, protocol := fish_ref$protocol[indx_m]]

# Add fish status data. For each fish, we will assign a "status" value of 
# 1 ("dead") at the next array after their last detection ("status_array"). 
# Fish that were detected at the last array (A9) will be assigned a value of 0 
# (alive) at the same array.
# Add "status_array" (array after the last detection)
indx_a <- match(migration$last_array, arrays)
indx_a <- ifelse(indx_a != length(arrays), indx_a + 1, indx_a)
migration[, status_array := arrays[indx_a]]

# Fish that were never detected will be assigned with the first array after 
# their release location (R1 -> A1; R2 -> A3)
migration$status_array[is.na(migration$status_array) &
                               migration$release_location == "R1"] <- "A1"
migration$status_array[is.na(migration$status_array) &
                         migration$release_location == "R2"] <- "A3"
migration[, status_array := factor(status_array)]
migration[, array_indx := as.numeric(status_array)]

# Add fish status (1 = dead, 0 = alive at the last array (A9))
migration$status <- 1
migration$status[which(migration$last_array == "A9")] <- 0

# Fit survival models
migr_fit <- lapply(species, function(sp) {
  
  migr_sub <- migration[species == sp, ]
  
  migr_loc <- lapply(unique(migr_sub$release_location), function(loc) {
    
    # Subset data
    data <- migr_sub[release_location == loc]
    # Assign origin depending on the release location
    origin <- ifelse(loc == "R1", 0, 2)
    
    # Create a survival object
    surv <- Surv(data$array_indx, data$status, origin = origin)
    
    # Fit a Kaplan-Meier survival curve
    fit <- survfit(surv ~ data$protocol)
    
    # Log rank test to compare the survival curves of different protocols and
    # extract a p-value
    log_rank <- survdiff(surv ~ data$protocol)
    pval <- log_rank$pvalue
    
    return(list(fit = fit, data = data, log_rank = log_rank, pval = pval))
  })
  names(migr_loc) <- unique(migr_sub$release_location)
  return(migr_loc)
})
names(migr_fit) <- species


# 4. PLOT RESULTS ==============================================================

# Select colors for the different protocols
col <- c(OPi = "olivedrab1", R64K = "skyblue")
col2 <- c(OPi = "olivedrab4", R64K = "skyblue4")

# Create plot directory
dir.create("./plots/smolt_migration/", showWarnings = FALSE)

jpeg("./plots/smolt_migration/smolt_migration.jpeg", width = 2000, 
     height = 2700, res = 300, pointsize = 13)
par(mfrow = c(3, 2), mar = c(4.1, 4.3, 2.6, 1))

# Panel letter index
let_c <- 1

for (sp in species) {
  
  det_sub <- smolt_det[species == sp]
  migr_sub <- migr_fit[[sp]]
  comm <- ifelse(sp == "Salmo trutta", "Trout", "Salmon")
  
  for (loc in unique(det_sub$release_location)) {
    
    # Plot number of smolts detected in each array
    det_tmp <- det_sub[release_location == loc]
    det_tmp[, array := factor(det_tmp$array)]
    
    plot(1, xlim = range(as.numeric(det_tmp$array)), type = "n", 
         ylim = c(8, 28), axes = FALSE, ann = FALSE)
    axis(1, at = 1:length(levels(det_tmp$array)), 
         labels = levels(det_tmp$array), lwd = 0, lwd.ticks = 1)
    axis(2, lwd = 0, lwd.ticks = 1)
    box()
    for (prot in unique(det_sub$protocol)) {
      lines(n_fish ~ array, data = det_tmp[protocol == prot], col = col2[prot])
    }
    points(n_fish ~ array, data = det_tmp, pch = 21, col = col2[protocol], 
           bg = col[protocol])
    # Add number of released individuals to legend
    released_n <- n_smolts[species == sp & release_location == loc]
    released_n <- released_n$n[match(c("OPi", "R64K"), released_n$protocol)]
    legend("bottomleft",  bty = "n", inset = c(0.02, 0.01),
           legend = paste0(c("OPi", "R64K"), " (n = ", released_n, ")"), 
           pch = 21, pt.bg = col, col = col2, )
    title(ylab = "Number of detected individuals", line = 2.4)
    title(xlab = "Receiver array", line = 2.4)
    title(paste0(comm, " smolts (", loc, ")"))
    mtext(LETTERS[let_c], side = 2, las = 1, padj = -7.5, adj = 3, font = 2, 
          cex = 1.2)
    
    let_c <- let_c + 1
    
    # Plot migration success
    migr_loc <- migr_sub[[loc]][["fit"]]
    prot_tmp <- gsub(".*=(.*)", "\\1", names(migr_loc$strata))
    migr_loc$protocol <- unlist(sapply(1:length(prot_tmp), 
                                       function(p) rep(prot_tmp[p], 
                                                       migr_loc$strata[p])))
    
    plot(1, xlim = range(as.numeric(det_tmp$array)) + c(0, 0.5), type = "n", 
         ylim = c(0.1, 1), axes = FALSE, ann = FALSE)
    axis(1, at = 1:length(levels(det_tmp$array)), 
         labels = levels(det_tmp$array), lwd = 0, lwd.ticks = 1)
    axis(2, lwd = 0, lwd.ticks = 1)
    box()
    
    # Confidence intervals
    for (prot in unique(migr_loc$protocol)) {
      indx_p <- which(migr_loc$protocol == prot)
      n <- length(indx_p)
      
      x <- c(rep(migr_loc$time[indx_p], each = 2), par("usr")[2])[-1]
      y <- c(rep(migr_loc$upper[indx_p], each = 2), 
             rev(rep(migr_loc$lower[indx_p], each = 2)))
      polygon(x = c(x, rev(x)), y, col = adjustcolor(col[prot], 0.4), 
              border = "transparent")
    }
    
    # Mean survival
    for (prot in unique(migr_loc$protocol)) {
      indx_p <- which(migr_loc$protocol == prot)
      n <- length(indx_p)
      
      segments(x0 = migr_loc$time[indx_p], 
               x1 = c(migr_loc$time[indx_p][-1], par("usr")[2]),
               y0 = migr_loc$surv[indx_p], col = col2[prot], lwd = 2)
      segments(x0 = migr_loc$time[indx_p], y0 = c(1, migr_loc$surv[indx_p][-n]),
               y1 = migr_loc$surv[indx_p], col = col2[prot], lwd = 2)
    }
    
    legend("bottomleft", legend = c("OPi", "R64K"), bty = "n",
           lwd = 2, inset = c(0.02, 0.01), col = col2)
    title(ylab = "Migration success", line = 2.4)
    title(xlab = "Receiver array", line = 2.4)
    text(x = par("usr")[2], par("usr")[4], adj = c(1.2, 2.2), cex = 1.2,
         labels = paste0("p = ", round(migr_sub[[loc]][["pval"]], 2)))
    title(paste0(comm, " smolts (", loc, ")"))
    mtext(LETTERS[let_c], side = 2, las = 1, padj = -7.5, adj = 3, font = 2,
          cex = 1.2)
    let_c <- let_c + 1

  }
}

dev.off()

