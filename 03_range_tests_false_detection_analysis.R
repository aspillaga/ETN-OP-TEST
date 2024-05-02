################################################################################
##
## OP-TESTS
## Analysis of false detections
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

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
receivers <- fread("./data/range_tests/receiver_deployments.csv")
tags <- fread("./data/range_tests/transmitter_deployments.csv")
detections <- fread("./data/range_tests/detection_data.csv.gz")

# Lists of receiver manufacturers
rec_man <- c("Thelma", "Lotek", "Sonotronics",
             "Innovasea-MAP114", "Innovasea-MAP115")

# List of tag protocols
prot_lev <- c("OPi", "OPs", "R64K", "S256", "A69-1602", "A69-9007", "Other")


# 1. IDENTIFY AND COUNT FALSE DETECTIONS =======================================

# IDs of tags known to be in the study area (test tags + others)
known_tags <- c(unique(tags$acoustic_tag_id), "OPi-101", "OPi-102",
                "R64K-13745", "S256-73", "OPi-10124")

# Identify false detections
detections[, false_det := !acoustic_tag_id %in% known_tags]

# Protocol factor (ordered)
detections[, tag_protocol := factor(tag_protocol, levels = prot_lev)]
detections$tag_protocol[is.na(detections$tag_protocol)] <- "Other"

# Calculate number of false detections per test, day, receiver and protocol
false_det <- lapply(unique(detections$test_id), function(t) {
  
  # Subset detections
  det_sub <- detections[test_id == t]
  
  # Extract test days
  det_sub[, date := as.Date(date_time)]
  date_range <- seq(min(det_sub$date), max(det_sub$date), by = "days")
  det_sub[, date := factor(date, levels = date_range)]
  
  # Receiver factor
  det_sub[, receiver_id := factor(receiver_id)]
  
  # Count false detections for each protocol
  f_sub <- det_sub[which(false_det), ]
  false_tab <- table(f_sub$date, f_sub$tag_protocol, f_sub$receiver_id,
                     dnn = c("date", "protocol", "receiver_id"))
  false_tab <- data.table(false_tab)
  colnames(false_tab)[4] <- "n_false"
  
  # Count total number of detections for each protocol
  prot_tab <- table(det_sub$date, det_sub$tag_protocol, det_sub$receiver_id, 
                    dnn = c("date", "protocol", "receiver_id"))
  prot_tab <- data.table(prot_tab)
  false_tab$n_protocol <- prot_tab$N
  
  # Count total number of detections per receiver and day
  all_tab <- table(det_sub$date, det_sub$receiver_id, 
                   dnn = c("date", "receiver_id"))
  all_tab <- data.table(all_tab)
  indx <- match(paste(false_tab$date, false_tab$receiver_id),
                paste(all_tab$date, all_tab$receiver_id))
  false_tab$n_detect <- all_tab$N[indx]
  
  # Add test id and station id
  false_tab$test_id <- t
  false_tab[, station := f_sub$station[match(receiver_id, f_sub$receiver_id)]]
  
  return(false_tab)
  
})
false_det <- rbindlist(false_det)

# Delete days with no detections
false_det <- false_det[n_detect > 0]

# Format the date
false_det[, date := as.Date(date)]

# Add protocol and manufacturer information
false_det[, protocol := factor(protocol, levels = prot_lev)]
indx_rec <- match(false_det$receiver_id, receivers$receiver_id)
false_det[, manufacturer := factor(receivers$manufacturer[indx_rec], 
                                   levels = rec_man)]

# Add test location
false_det[, test_loc := factor(gsub("(.*)-.*", "\\1", test_id))]

# Calculate proportion of false detections
false_det[, prop_false := 100 * n_false / n_detect]

head(false_det)

# Summary per test, protocol and receiver manufacturer 
false_summ <- false_det[, .(total_n = sum(n_false),
                            total_detect = sum(n_detect),
                            mean_n = mean(n_false, na.rm = TRUE),
                            median_n = quantile(n_false, 0.5, na.rm = TRUE),
                            min_n = quantile(n_false, 0, na.rm = TRUE),
                            max_n = quantile(n_false, 1, na.rm = TRUE),
                            mean_prop = mean(prop_false, na.rm = TRUE),
                            median_prop = quantile(prop_false, 0.5, 
                                                   na.rm = TRUE),
                            min_prop = quantile(prop_false, 0.025, 
                                                na.rm = TRUE),
                            max_prop = quantile(prop_false, 0.975, 
                                                na.rm = TRUE)),
                        by = list(test_loc, manufacturer, protocol)]
false_summ <- false_summ[order(test_loc, manufacturer, protocol)]

# Remove non-compatible combinations
compatible <- tapply(false_det$n_protocol, 
                     list(false_det$manufacturer, false_det$protocol), sum)
compatible <- apply(compatible, 1, function(x) return(names(which(x == 0))))
indx_rm <- unlist(lapply(seq_along(compatible), function(i) {
  which(false_summ$manufacturer == names(compatible)[i] & 
          false_summ$protocol %in% compatible[[i]])
}))
false_summ <- false_summ[-indx_rm]

head(false_summ)


# 2. RESULTS FOR THE MANUSCRIPT ================================================

# Number and percentage of false detections in the entire dataset:
sum(detections$false_det)
round(100 * sum(detections$false_det) / nrow(detections), 2)

# Number ad percentage of false detections by protocol
n_protocol <- table(detections$tag_protocol[detections$false_det])
prop_protocol <- 100 * n_protocol / sum(detections$false_det)
prop_detect <- 100 * n_protocol / nrow(detections)

n_protocol
round(prop_detect, 2)
round(prop_protocol, 2)

sum(n_protocol[c("OPi", "OPs")]) # OPi and OPs
round(100 * sum(n_protocol[c("OPi", "OPs")])/nrow(detections), 2) # OPi and OPs
round(sum(prop_protocol[c("OPi", "OPs")]), 1) # OPi and OPs


# False detection occurrence rate per receiver and day (sum the protocols):
false_all <- false_det[, .(n_false = sum(n_false), n_detect = mean(n_detect)),
                       by = list(date, receiver_id, test_id, station, 
                                 manufacturer, test_loc)]
false_all[, false_prop := 100 * n_false / n_detect]

# False detection occurrence per receiver:
false_rec <- false_all[, .(n_false = sum(n_false), n_detect = sum(n_detect)),
                       by = list(receiver_id, test_id, station, manufacturer, 
                                 test_loc)]
false_rec[, false_prop := 100 * n_false / n_detect]

# Mean and 95% range of daily false detection proportions per deployment
round(mean(false_all$false_prop), 2)
round(quantile(false_all$false_prop, c(0.025, 0.975)), 2)

# Maximum daily false detection rates per environment
round(tapply(false_all$false_prop, 
             list(false_all$test_loc, false_all$manufacturer), max), 2)

# Check most recurrent false detection codes
n_code <- sort(table(detections$acoustic_tag_id[detections$false_det]), 
               decreasing = TRUE)
length(n_code)

# Proportion of IDs detected less than 10 times
sum(n_code < 10)
round(100 * sum(n_code < 10) / length(n_code), 1)

# Proportion of detections of codes detected less than 10 times
round(100 * sum(n_code[n_code < 10]) / sum(n_code), 1)

# Proportion of IDs detected 10 or more times
n_code[n_code >= 10]
round(100 * n_code[n_code >= 10] / sum(n_code), 2)
sum(n_code >= 10)
round(100 * sum(n_code >= 10) / length(n_code), 1)

# Proportion of detections of codes detected 10 or more times
round(100 * sum(n_code[n_code >= 10]) / sum(n_code), 1)

n_code[grepl("OP", names(n_code))][1:10]
round(100 * n_code[grepl("OP", names(n_code))] / sum(n_code), 2)[1:10]
length(n_code[grepl("OP", names(n_code))])

indx <- grepl("OP", names(n_code)) & n_code < 5
length(n_code[indx])
100 * sum(n_code[indx]) / sum(n_code)
n_code[][1:10]
round(100 * n_code[grepl("OP", names(n_code))] / sum(n_code), 2)[1:10]
length(n_code[grepl("OP", names(n_code))])


# 2. PLOT FALSE DETECION RATES =================================================

dir.create("./plots/range_tests/false_detections/", showWarnings = FALSE)

# Test labels
test_lab <- c(CP = "Open sea (BE)", Yser = "River (BE)", 
              Formosa = "Coastal lagoon (PT)", Palma = "Coastal habitat (SP)")

# Protocol colors
col <- c(OPi = "#4c76ba", OPs = "#71c2eb", R64K = "#a3c73e", S256 = "#e8d946", 
         "A69-1602" = "#ec6667", "A69-9007" = "#ef8c48", Others = "#816aac")
col_l <- c(OPi = "#153168", OPs = "#28809b", R64K = "#6c922d", S256 = "#9b9113", 
           "A69-1602" = "#af2322", "A69-9007" = "#b24818", Others = "#45277b")


# Plot number of false detections per day
jpeg("./plots/range_tests/false_detections/false_detection_rate.jpeg", 
     width = 2200, height = 2100, res = 300, pointsize = 12)
layout(matrix(c(1, 1, 5, 2, 2, 5, 3, 4, 4), ncol = 3, byrow = T), 
       widths = c(5, 3, 2))

for (i in seq_along(test_lab)) {
  
  data_sub <- false_summ[test_loc == names(test_lab[i])]
  
  sep <- 1.2
  indx_at <- c(TRUE, (data_sub$manufacturer[-1] == 
                        data_sub$manufacturer[-nrow(data_sub)]))
  data_sub[, at := cumsum(ifelse(indx_at, 0.7, sep))]
  xlim <- c(range(data_sub$at)) + c(-sep/2, sep/2)
  ylim <- c(0, max(data_sub$max_n))
  
  l = 0.2
  par(mar = c(4.1, 4.6, 2.6, 1.1))
  plot(max_n ~ at, data = data_sub, xlim = xlim, ylim = ylim, xaxs = "i", 
       type = "n", ann = FALSE, axes = FALSE)
  arrows(x0 = data_sub$at, y0 = data_sub$min_n,
         y1 = data_sub$max_n, col = col_l[data_sub$protocol],
         code = 3, angle = 90, length = 0.05)
  
  points(mean_n ~ at, data = data_sub, pch = 21, cex = 1.4,
         col = col_l[data_sub$protocol], bg = col[data_sub$protocol])
  
  abline(v = data_sub$at[!indx_at] - (sep/2), col = "gray50")
  axis(2, at = pretty(ylim, 4), las = 1, lwd = 0, lwd.ticks = 1)
  par(xpd = TRUE)
  rect(xleft = par("usr")[1], xright = par("usr")[2], 
       ytop = par("usr")[3],
       ybottom = par("usr")[3] - 1*diff(par("usr")[3:4])/5, col = "gray90")
  par(xpd = FALSE)
  
  axis(1, at = c(par("usr")[1], data_sub$at[!indx_at] - (sep/2), par("usr")[2]),
       labels = FALSE, lwd = 0, lwd.ticks = 1, tcl = -2.3)
  axis_at <- tapply(data_sub$at, data_sub$manufacturer, mean)
  axis(1, at = axis_at, labels = gsub("-", "\n", names(axis_at)), 
       lwd = 0, line = -0.7, padj = 0.5)
  
  title(ylab = "No. of false detections / day", line = 2.4)
  title(xlab = "Receiver manufacturer", line = 2.6, font.lab = 2)
  box()
  
  title(main = test_lab[i], adj = 0, line = 0.5)
  mtext(side = 2, text = LETTERS[i], font = 2, padj = -5.7, line = 2.2, las = 1, 
        cex = 1.2)
}

par(mar = c(1.8, 0, 0, 1.1))
plot(1, type = "n", axes = FALSE, ann = FALSE)
legend("left", legend = prot_lev, col = col_l, pt.bg = col, pch = 21, 
       pt.cex = 1.4, box.col = "gray50", inset = c(0.1, 0), cex = 1.4, 
       title = "Tag protocol", title.font = 2)

dev.off()


# Plot proportion of false detections per day

jpeg("./plots/range_tests/false_detections/false_detection_proportion.jpeg", 
     width = 2200, height = 2100, res = 300, pointsize = 12)
layout(matrix(c(1, 1, 5, 2, 2, 5, 3, 4, 4), ncol = 3, byrow = T), 
       widths = c(5, 3, 2))

for (i in seq_along(test_lab)) {
  
  data_sub <- false_summ[test_loc == names(test_lab[i])]
  
  sep <- 1.2
  indx_at <- c(TRUE, (data_sub$manufacturer[-1] == 
                        data_sub$manufacturer[-nrow(data_sub)]))
  data_sub[, at := cumsum(ifelse(indx_at, 0.7, sep))]
  xlim <- c(range(data_sub$at)) + c(-sep/2, sep/2)
  ylim <- c(0, max(data_sub$max_prop))
  
  l = 0.2
  par(mar = c(4.1, 4.6, 2.6, 1.1))
  plot(max_prop ~ at, data = data_sub, xlim = xlim, ylim = ylim, xaxs = "i", 
       type = "n", ann = FALSE, axes = FALSE)
  arrows(x0 = data_sub$at, y0 = data_sub$min_prop,
         y1 = data_sub$max_prop, col = col_l[data_sub$protocol],
         code = 3, angle = 90, length = 0.05)
  
  points(mean_prop ~ at, data = data_sub, pch = 21, cex = 1.4,
         col = col_l[data_sub$protocol], bg = col[data_sub$protocol])
  
  abline(v = data_sub$at[!indx_at] - (sep/2), col = "gray50")
  axis(2, at = pretty(ylim, 4), las = 1, lwd = 0, lwd.ticks = 1)
  par(xpd = TRUE)
  rect(xleft = par("usr")[1], xright = par("usr")[2], 
       ytop = par("usr")[3],
       ybottom = par("usr")[3] - 1*diff(par("usr")[3:4])/5, col = "gray90")
  par(xpd = FALSE)
  
  axis(1, at = c(par("usr")[1], data_sub$at[!indx_at] - (sep/2), par("usr")[2]),
       labels = FALSE, lwd = 0, lwd.ticks = 1, tcl = -2.3)
  axis_at <- tapply(data_sub$at, data_sub$manufacturer, mean)
  axis(1, at = axis_at, labels = gsub("-", "\n", names(axis_at)), 
       lwd = 0, line = -0.7, padj = 0.5)
  
  title(ylab = "% of detections / day", line = 2.6)
  title(xlab = "Receiver manufacturer", line = 2.6, font.lab = 2)
  box()
  
  title(main = test_lab[i], adj = 0, line = 0.5)
  mtext(side = 2, text = LETTERS[i], font = 2, padj = -5.7, line = 2.2, las = 1, 
        cex = 1.2)
}

par(mar = c(1.8, 0, 0, 1.1))
plot(1, type = "n", axes = FALSE, ann = FALSE)
legend("left", legend = prot_lev, col = col_l, pt.bg = col, pch = 21, 
       pt.cex = 1.4, box.col = "gray50", inset = c(0.1, 0), cex = 1.4, 
       title = "Tag protocol", title.font = 2)

dev.off()


