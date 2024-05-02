################################################################################
##
## OP-TESTS
## Download acoustic-range data from ETN database
##
## Author: 
## Eneko Aspillaga (aspillaga@imedea.uib-csic.es)
## Instituto Mediterraneo de Estudios Avanzados (IMEDEA, CSIC)
##
## Last update: November 3, 2023
##
################################################################################

# NOTE: All this script should be run within the European Tracking Network
#       R studio server (https://rstudio.lifewatch.be/), because it is the only 
#       way to access the ETN database.

# Load libraries
library(data.table)
library(etn)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# 1. DOWNLOAD DATA FROM THE ETN DATABASE =======================================

# Database connection
usr <- "my_user" # Personal user name
pass <- "password" # Password
con <- connect_to_etn(usr, pass)

# Select acoustic telemetry project
proj <- "OP-Test"


# 1.1. Download receiver deployment metadata ---------------------------------

receivers <- get_acoustic_deployments(con, acoustic_project_code = proj)
receivers <- data.table(receivers)

# Select columns of interest
receivers <- receivers[, c("deployment_id", "receiver_id", "station_name",
                           "deploy_date_time", "deploy_latitude",
                           "deploy_longitude", "recover_date_time",
                           "download_file_name")]


# 1.2. Download receiver metadata --------------------------------------------

# Extract unique receiver IDs and donwload metadata
rec_ids <- unique(receivers$receiver_id)
rec_metadata <- get_acoustic_receivers(con, receiver_id = rec_ids)
rec_metadata <- data.table(rec_metadata)

# Select columns of interest
rec_metadata <- rec_metadata[, c("receiver_id", "manufacturer",
                                 "receiver_model", "receiver_serial_number")]


# 1.3. Download transmitter deployment metadata ------------------------------

# Test transmitters are stored as "animals" in the ETN database
tags <- get_animals(con, animal_project_code = proj)
tags <- data.table(tags)

# Select columns of interest
tags <- tags[, c("tag_serial_number", "acoustic_tag_id", "scientific_name",
                 "release_date_time", "release_location", "release_latitude",
                 "release_longitude", "recapture_date_time")]


# 1.4. Download transmitter metadata -----------------------------------------

# Extract unique tag IDs and download metadata
tag_ids <- unique(tags$acoustic_tag_id)
tag_metadata <- get_tags(con, acoustic_tag_id = tag_ids)
tag_metadata <- data.table(tag_metadata)

# Select columns of interest
tag_metadata <- tag_metadata[, c("tag_serial_number", "acoustic_tag_id",
                                 "manufacturer", "model", "frequency",
                                 "step1_min_delay", "step1_max_delay",
                                 "step1_power")]


# 1.5. Download detection data -----------------------------------------------

detect <- get_acoustic_detections(con, acoustic_project_code = proj)
detect <- data.table(detect)

# Select columns of interest
detect <- detect[, c("detection_id", "date_time", "tag_serial_number",
                     "acoustic_tag_id", "receiver_id", "source_file",
                     "deployment_id")]


# 1.6. Export files to load them offline -------------------------------------

fwrite(receivers, "etn_receiver_deployments.csv")
fwrite(rec_metadata, "etn_receiver_metadata.csv")
fwrite(tags, "etn_transmitter_deployments.csv")
fwrite(tag_metadata, "etn_transmitter_metadata.csv")
fwrite(detect, "rec_detection_data.csv.gz")
