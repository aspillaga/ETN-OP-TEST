################################################################################
##
## OP-TESTS
## Prepare data for range tests
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
library(lubridate)
library(sp)
library(geodist)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# 1. PREPARE DATA FOR THE ANALYSIS =============================================

# Load data files (downloaded from the ETN database)
path <- "./data/range_tests/etn_data/"
receivers <- fread(paste0(path, "etn_receiver_deployments.csv"))
rec_metadata <- fread(paste0(path, "etn_receiver_metadata.csv"))
tags <- fread(paste0(path, "etn_transmitter_deployments.csv"))
tag_metadata <- fread(paste0(path, "etn_transmitter_metadata.csv"))
detect <- fread(paste0(path, "rec_detection_data.csv.gz"))

# Function to rename manufacturer name
renameMfs <- function(x) {
  sub("THELMA BIOTEL|THELMA", "Thelma", 
      sub("LOTEK", "Lotek",
          sub("SONOTRONICS", "Sonotronics",
              sub("VEMCO|INNOVASEA", "Innovasea", x))))
}

# Function to rename acoustic signal protocols
renameProt <- function(x) {
  sub("OPI", "OPi", 
      sub("OPS", "OPs", 
          sub("A69-1303|SR2", "R64K", # Innovasea's A69-1303 and SR2 == R64K
              sub("A69-1105", "S256", x)))) # Innovasea's A69-1105 == S256
}


# 1.1. Prepare receiver deployment metadata ------------------------------------

# Add receiver model and manufacturer from receiver metadata
indx_r <- match(receivers$receiver_id, rec_metadata$receiver_id)
receivers$receiver_model <- rec_metadata$receiver_model[indx_r]
receivers$receiver_serial_number <- rec_metadata$receiver_serial_number[indx_r]
receivers$manufacturer <- renameMfs(rec_metadata$manufacturer[indx_r])

# Add code map information to Innovasea's receiver manufacturer label
inno <- which(receivers$manufacturer == "Innovasea")
map <- ifelse(grepl("MAP115", receivers$station_name[inno]), "MAP115", "MAP114")
receivers$manufacturer[inno] <- paste0(receivers$manufacturer[inno], "-", map)

# Homogenize station names
unique(receivers$station_name)
receivers[, station_name := gsub("(.*_.*)_.*", "\\1", station_name)]
receivers[, station_name := gsub("IJzer", "Yser", station_name)]
receivers[, station_name := gsub("Palma_R1", "Palma_0m", station_name)]
receivers[, station_name := gsub("Palma_R2", "Palma_450m", station_name)]
receivers[, station_name := gsub("Mooring_1", "Formosa_0m", station_name)]
receivers[, station_name := gsub("Mooring_4", "Formosa_300m", station_name)]


# Add new column to separate different tests (OPi and OPs tests)
receivers$test_id <- NA

# C-Power Wind Farm tests, OPi and OPs (Belgium):
receivers$test_id[grepl("CP", receivers$station_name) &
                    receivers$deploy_date_time < ymd("2021-02-01")] <- "CP-OPi"
receivers$test_id[grepl("CP", receivers$station_name) &
                    receivers$deploy_date_time > ymd("2021-02-01")] <- "CP-OPs"

# Yser river tests, OPi and OPs (Belgium):
receivers$test_id[grepl("Yser", receivers$station_name) &
                    receivers$deploy_date_time <
                    ymd("2021-05-06")] <- "Yser-OPi"
receivers$test_id[grepl("Yser", receivers$station_name) &
                    receivers$deploy_date_time >
                    ymd("2021-05-06")] <- "Yser-OPs"

# Palma Bay test, only OPi (Spain):
receivers$test_id[grepl("Palma", receivers$station_name)] <- "Palma-OPi"

# Ria Formosa tests, OPi and OPs (Portugal):
receivers$test_id[grepl("Formosa", receivers$station_name) &
                    receivers$deploy_date_time <
                    ymd("2021-07-28")] <- "Formosa-OPi"
receivers$test_id[grepl("Formosa", receivers$station_name) &
                    receivers$deploy_date_time >
                    ymd("2021-07-28")] <- "Formosa-OPs"

# Check dataframe
head(receivers)


# 1.2. Prepare transmitter deployment metadata ---------------------------------

# Add transmitter model, transmission interval and manufacturer from 
# transmitter metadata
indx_t <- match(tags$acoustic_tag_id, tag_metadata$acoustic_tag_id)
tags$model <- tag_metadata$model[indx_t]
tags$manufacturer <- renameMfs(tag_metadata$manufacturer[indx_t])
tags$power <- tag_metadata$step1_power[indx_t]
tags$min_delay <- tag_metadata$step1_min_delay[indx_t]
tags$max_delay <- tag_metadata$step1_max_delay[indx_t]
tags$avg_delay <- (tags$min_delay + tags$max_delay) / 2

# Rename tag protocols and add a separate column
tags$acoustic_tag_id <- renameProt(tags$acoustic_tag_id)
tags$protocol <- sub("(.*-?.*)-.*", "\\1", tags$acoustic_tag_id)

# Homogenize station names
unique(tags$release_location)
tags[, release_location := gsub("_base", "", release_location)]
tags[, release_location := gsub("IJzer", "Yser", release_location)]
tags[, release_location := gsub("Palma_T1", "Palma_0m", release_location)]
tags[, release_location := gsub("Palma_T2", "Palma_150m", release_location)]
tags[, release_location := gsub("Palma_T3", "Palma_225m", release_location)]
tags[, release_location := gsub("Mooring_1", "Formosa_0m", release_location)]
tags[, release_location := gsub("Mooring_2", "Formosa_100m", release_location)]
tags[, release_location := gsub("Mooring_3", "Formosa_150m", release_location)]

# Separate transmitter deployments by tests
tags[, test_id := gsub("(.*)_(.*)_.*", "\\2-\\1", release_location)]

# Check dataframe
head(tags)


# 1.3. Prepare detection data --------------------------------------------------

# Order detections by date and time
detect <- detect[order(detect$date_time), ]

# Rename protocol names
detect$acoustic_tag_id <- renameProt(detect$acoustic_tag_id)

# Extract tag protocol
detect$tag_protocol <- sub("(.*-?.*)-.*", "\\1", detect$acoustic_tag_id)

# Remove duplicated detections
rm <- duplicated(detect[, c("date_time", "acoustic_tag_id", "receiver_id")])
sum(rm)
detect <- detect[!rm, ]

# Assign the test ID and station to each detection depending on the deployment
# date and time
detect$test_id <- NA
detect$station <- NA

for (i in 1:nrow(receivers)) {
  t <- receivers$test_id[i]
  
  # Extract the initial and end date of tests from transmitter deployments
  date_ini <- max(tags$release_date_time[tags$test_id == t])
  date_ini <- ceiling_date(date_ini + 2*3600, "hours")
  date_end <- min(tags$recapture_date_time[tags$test_id == t])
  date_end <- floor_date(date_end - 2*3600, "hours")
  
  # Assign test ID to detections
  indx <- which(detect$receiver_id == receivers$receiver_id[i] &
                  detect$date_time >= date_ini & detect$date_time < date_end)
  detect$test_id[indx] <- t
  detect$station[indx] <- receivers$station_name[i]
}

# Check dataframe
head(detect)

detect <- detect[, c("test_id", "date_time", "station", "receiver_id",
                     "acoustic_tag_id", "tag_protocol")]


# 1.4. Export data -------------------------------------------------------------

# Remove data outside the test periods
detect <- detect[!is.na(detect$test_id)]

# Remove Lotek receivers from tests in Ria Formosa (they run out of battery
# after few days)
receivers <- receivers[!(grepl("Formosa", test_id) & manufacturer == "Lotek")]
detect <- detect[!(grepl("Formosa", test_id) & grepl("WHS6K", receiver_id))]

# Remove two malfunctioning transmitters (did not emit signals properly)
rm_tags <- c("OPi-10124", "A69-9007-16321")
tags <- tags[!acoustic_tag_id %in% rm_tags]

# Export data
fwrite(receivers, "./data/range_tests/receiver_deployments.csv")
fwrite(tags, "./data/range_tests/transmitter_deployments.csv")
fwrite(detect, "./data/range_tests/detection_data.csv.gz")


# 2. PREPARE DETECTION PROBABILITY DATA ========================================

# 2.1. Calculate distances between receivers and transmitters ------------------

test_ids <- unique(receivers$test_id)

distances <- lapply(test_ids, function(t) {
  cat(t, "\n")
  rec_sub <- receivers[receivers$test_id == t, ]
  tag_sub <- tags[test_id == t]
  setnames(rec_sub,  c("deploy_longitude", "deploy_latitude"), 
           c("longitude", "latitude"))
  setnames(tag_sub,   c("release_longitude", "release_latitude"), 
           c("longitude", "latitude"))
  dist_sub <- geodist(x = rec_sub[, c("longitude", "latitude")],
                      y = tag_sub[, c("longitude", "latitude")], 
                      measure = "geodesic")
  rownames(dist_sub) <- rec_sub$receiver_id
  colnames(dist_sub) <- tag_sub$acoustic_tag_id
  return(dist_sub)
})
names(distances) <- unique(receivers$test_id)
distances


# 2.2. Calculate number of detections at 2 h intervals -------------------------

detect_2h <- lapply(test_ids, function(t) {
  
  # Define time interval in hours
  time_int <- 2
  
  # Subset tag and receiver IDs for the specific test
  tag_id_sub <- tags$acoustic_tag_id[tags$test_id == t]
  rec_id_sub <- receivers$receiver_id[receivers$test_id == t]
  
  # Subset detections for the specific test
  detect_sub <- detect[test_id == t & acoustic_tag_id %in% tag_id_sub]

  # Create receiver and transmitter factors
  detect_sub[, receiver_id := factor(receiver_id, levels = rec_id_sub)]
  detect_sub[, acoustic_tag_id := factor(acoustic_tag_id, levels = tag_id_sub)]
  
  # Create time factor (divided in time bins)
  detect_sub[, hour := floor_date(date_time, paste0(time_int, "hours"))]
  detect_sub[, hour := factor(format(hour, "%Y-%m-%d %H:%M"))]
  
  # Calculate number of detections in each interval time interval
  detect_n <- table(detect_sub$receiver_id, detect_sub$acoustic_tag_id, 
                    detect_sub$hour, 
                    dnn = c("receiver_id", "acoustic_tag_id", "date_time"))
  detect_n <- data.table(as.data.frame(detect_n, responseName = "n"))
  
  # Format date time as POSIXct
  detect_n[ , date_time := ymd_hm(date_time)]
  
  # Add test id (location and main tested protocol, OPi or OPs)
  detect_n[, test_id := t]
  
  # Add test location (without tested protocol)
  detect_n[, test_location := sub("(.*)-.*", "\\1", t)]
  
  # Add receiver manufacturer
  indx_rec <- match(detect_n$receiver_id, receivers$receiver_id)
  detect_n[, receiver_manufacturer := receivers$manufacturer[indx_rec]]
  
  # Add tag manufacturer
  indx_tag <- match(paste(t, detect_n$acoustic_tag_id), 
                    paste(tags$test_id, tags$acoustic_tag_id))
  detect_n[, tag_manufacturer := paste0(tags$manufacturer, "-", 
                                        tags$protocol)[indx_tag]]
  
  # Assign expected number of signals (it has to be equal or larger than the 
  # detected ones)
  detect_n[ , expected_n := round(time_int*3600 / tags$avg_delay[indx_tag])]
  detect_n[ , expected_n := ifelse(n > expected_n, n, expected_n)]
  
  # Calculate detection probability
  detect_n[ , detect_prob := round(n / expected_n, 3)]
  
  # Assign distance to each detection count
  detect_n$distance <- round(sapply(1:nrow(detect_n), function(i) {
    distances[[t]][as.character(detect_n$receiver_id[i]), 
                   as.character(detect_n$acoustic_tag_id[i])]
  }))
  
  # Select columns of interest
  detect_n <- detect_n[, c("test_location", "test_id", "date_time", 
                           "receiver_manufacturer", "receiver_id", 
                           "tag_manufacturer", "acoustic_tag_id",
                           "distance", "n", "expected_n", "detect_prob")]
  
  return(detect_n)
  
})

detect_2h <- rbindlist(detect_2h)
detect_2h

# Export data
fwrite(detect_2h, file = "./data/range_tests/detections_2h_intervals.csv")


# 3. CHECK PLOTS ===============================================================

# Create directory for results
dir.create("./plots/range_tests/data_check/", showWarnings = FALSE, 
           recursive = TRUE)


# 3.1. Plot deployments --------------------------------------------------------

jpeg("./plots/range_tests/data_check/deployment_maps.jpeg", width = 1800, 
     height = 900, res = 200)
par(mfrow = c(2, 4), mar = c(3.1, 3.1, 3.1, 1.1))

for (t in unique(receivers$test_id)) {
  
  # Extract receiver positions
  r_pos <- receivers[test_id == t, c("deploy_longitude", "deploy_latitude")]
  r_pos <- r_pos[!duplicated(r_pos)][, type := "Receiver"]
  
  # Extract transmitter positions
  t_pos <- tags[test_id == t, c("release_longitude", "release_latitude")]
  t_pos <- t_pos[!duplicated(t_pos)][, type := "Tag"]
  colnames(t_pos) <- colnames(r_pos)
  
  # Merge receivers and transmitters and generate SpatialPoints
  pos <- rbind(r_pos, t_pos)
  pos$type <- factor(pos$type, levels = c("Receiver", "Tag"))
  pos <- SpatialPointsDataFrame(pos[, 1:2], data = pos[, 3],
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Plot positions
  plot(pos, pch = 16, col = c(1,3)[pos$type], cex = c(1.2, 0.8)[pos$type])
  axis(1, lwd = 0, lwd.ticks = 1)
  axis(2, lwd = 0, lwd.ticks = 1)
  
  # Add legend and title
  legend("bottomleft", legend = levels(pos$type), col = c(1, 3), 
         pch = 16, pt.cex = c(1.2, 0.8), inset = c(0.02, 0.04))
  title(main = t)
  
  box()
}

dev.off()


# 3.2. Plot time coverage for each test ----------------------------------------

jpeg(paste0("./plots/range_tests/data_check/temporal_coverage.jpeg"), 
     width = 2000, height = 1200, res = 150)
par(mfrow = c(3, 3), mar = c(4.1, 10, 3.1, 1.1))

for (t in unique(receivers$test_id)) {
  cat(t, "\n")
  
  # Extract deployments, transmitters and detections corresponding to each test
  r_tmp <- receivers[test_id == t]
  t_tmp <- tags[test_id == t]
  det_tmp <- detect[test_id == t]
  
  # Convert receiver_id into a factor
  det_tmp$receiver_id <- factor(det_tmp$receiver_id, 
                                levels = sort(r_tmp$receiver_id))
  
  # Extract temporal range of each test (deployment and retrieval dates)
  date_ini <- max(tags$release_date_time[tags$test_id == t])
  date_ini <- ceiling_date(date_ini + 2*3600, "hours")
  date_end <- min(tags$recapture_date_time[tags$test_id == t])
  date_end <- floor_date(date_end - 2*3600, "hours")
  time_range <- c(date_ini, date_end)
  
  # Plot ranges
  plot(det_tmp$date_time, det_tmp$receiver_id, pch = 16, cex = 0.6, 
       xlim = time_range, ann = FALSE, axes = FALSE)
  abline(v = time_range, col = "red")
  axis(2, at = seq_along(levels(det_tmp$receiver_id)), 
       labels = levels(det_tmp$receiver_id), las = 1, lwd = 0, lwd.ticks = 1)
  axis(1, at = pretty(time_range), lwd = 0, lwd.ticks = 1,
       labels = format(pretty(time_range), "%b-%d"))
  title(main = t)
  
}

dev.off()


# 3.3 Detection probability vs distance boxplots -------------------------------

# Names of receiver and transmitter manufacturers
rec_man <- c("Thelma", "Lotek", "Sonotronics", 
             "Innovasea-MAP114", "Innovasea-MAP115")
tag_man <- c("Thelma-OPi", "Thelma-OPs", "Thelma-R64K", 
             "Lotek-OPi", "Lotek-OPs", "Sonotronics-OPi", "Sonotronics-OPs",
             "Innovasea-R64K", "Innovasea-S256","Innovasea-A69-1602", 
             "Innovasea-A69-9007")

pdf(paste0("./plots/range_tests/data_check/prob_vs_distance.pdf"), 
    width = 12, height = 12, pointsize = 12)

for (t in unique(receivers$test_id)) {
  
  # Subset data to plot
  detect_sub <- detect_2h[test_id == t]
  
  # Extract distances
  distance_list <- sort(unique(detect_sub$dist))
  
  # Extract receiver and tag manufacturers
  rec_man_sub <- rec_man[rec_man %in% unique(detect_sub$receiver_manufacturer)]
  tag_man_sub <- tag_man[tag_man %in% unique(detect_sub$tag_manufacturer)]
  
  # Create plot
  par(mfrow = c(length(rec_man_sub), length(tag_man_sub)), 
      mar = c(4.1, 4.1, 4.1, 0.6), oma = c(0, 2, 3, 2))
  
  for (r in rec_man_sub) {
    for (g in tag_man_sub) {
      
      # Subset detections to plot
      det_tmp <- detect_sub[receiver_manufacturer == r & tag_manufacturer == g]
      indx <- which(det_tmp$n > 0)
      
      
      if(length(indx) > 0) {
        
        # If there are detections, plot boxplots
        time_range <- range(det_tmp$date_time[indx])
        det_tmp <- det_tmp[date_time >= time_range[1] & 
                             date_time < time_range[2]]
        dist_sub <- sort(unique(det_tmp$dist))
        bx <- boxplot(det_tmp$detect_prob ~ det_tmp$distance, ann = FALSE, 
                      axes = FALSE, xlim = c(0, length(dist_sub)) + 0.5, 
                      ylim = c(0, 1.1), at = which(dist_sub %in% distance_list), 
                      col = "#6aa3db")
        
      } else {
        
        # If there are no detections, create an empty plot
        dist_sub <- distance_list
        plot(1, axes = FALSE, ann = FALSE,type = "n",
             ylim = c(0, 1.1), xlim = c(0, length(distance_list)) + 0.5)
      }
      
      # Add axes and titles
      axis(1, at = seq_along(dist_sub), labels = dist_sub, lwd = 0, 
           lwd.ticks = 1)
      axis(2, las = 1, lwd = 0, lwd.ticks = 1)
      title(xlab = "Distance (m)", line = 2.4)
      title(ylab = "Detection prob.", line = 2.6)
      title(main = paste0("Receivers: ", r, "\n", "Tags: ", g))
      
      box()
    }
  }
  # Add main title with the test name
  title(main = t, outer = TRUE, cex.main = 2)
}

dev.off()


# 3.4 Detection probability vs time plots --------------------------------------

# List of receiver and transmitter manufacturers
rec_man <- c("Thelma", "Lotek", "Sonotronics", 
             "Innovasea-MAP114", "Innovasea-MAP115")
tag_man <- c("Thelma-OPi", "Thelma-OPs", "Thelma-R64K", 
             "Lotek-OPi", "Lotek-OPs", "Sonotronics-OPi", "Sonotronics-OPs",
             "Innovasea-R64K", "Innovasea-S256","Innovasea-A69-1602", 
             "Innovasea-A69-9007")

pdf(paste0("./plots/range_tests/data_check/detection_prob_vs_time.pdf"), 
    width = 16, height = 16, pointsize = 12)

for (t in unique(receivers$test_id)) {
  
  # Subset data to plot
  detect_sub <- detect_2h[test_id == t]
  
  # Extract distances
  distance_list <- sort(unique(detect_sub$dist))
  
  # Extract receiver and tag manufacturers
  rec_man_sub <- rec_man[rec_man %in% unique(detect_sub$receiver_manufacturer)]
  tag_man_sub <- tag_man[tag_man %in% unique(detect_sub$tag_manufacturer)]
  
  # Create plot
  par(mfrow = c(length(rec_man_sub), length(tag_man_sub)), 
      mar = c(4.1, 4.1, 4.1, 0.6), oma = c(0, 0, 3,4))
  
  for (r in rec_man_sub) {
    for (g in tag_man_sub) {

      # Subset specific detections for the plot
      det_tmp <- detect_sub[receiver_manufacturer == r & tag_manufacturer == g]
      
      # Extract distances to plot
      dist_sub <- sort(unique(det_tmp$dist))
      
      # Create empty plot
      plot(detect_prob ~ date_time, data = det_tmp, type = "n", 
           ylim = c(-0.25, 1.05), axes = FALSE, ann = FALSE)
      
      # Add detection probability vs time lines for each transmitter/
      # manufacturer combination
      for (k in seq_along(dist_sub)) {
        lines(detect_prob ~ date_time, data = det_tmp[distance == dist_sub[k]], 
              col = k)
      }
      
      # Add axes and titles
      axis(2, las = 1, lwd = 0, lwd.ticks = 1)
      dates <- sort(unique(floor_date(det_tmp$date_time, "days")))
      axis.POSIXct(1, at = dates, labels = "", lwd = 0, lwd.ticks = 1)
      axis.POSIXct(1, at = pretty(dates, 5), lwd = 0, lwd.ticks = 1,
                   tck = -0.04)
      title(xlab = "Date", line = 2.8)
      title(ylab = "Detection prob.", line = 2.6)
      title(main = paste0("Receivers: ", r, "\n", "Tags: ", g))
      
      # Add legend
      legend(x = "bottomright", lwd = 1, col = seq_along(dist_sub), 
             legend = dist_sub, horiz = TRUE, title = "Distance (m)",
             cex = 0.9, y.intersp = 0.8, x.intersp = 0.8)
      box()
    }
  }
  # Add main title
  title(main = t, outer = TRUE, cex.main = 2)
}

dev.off()

