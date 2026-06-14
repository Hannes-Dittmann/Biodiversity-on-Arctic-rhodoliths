


#load packages
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)



#1.load the raw data



raw <- read.csv(
  "table_raw_data.CSV",
  sep       = ";",
  header    = TRUE,
  check.names = FALSE
)



depth_df <- read.csv(
  "depth_of_station.CSV",
  sep      = ";",
  header   = TRUE,
  stringsAsFactors = FALSE
)



#2.prepare data for analysis



#separate metadata from raw counts
station_row   <- raw[raw$group == "site",      ]
count_data    <- raw[!raw$group %in% c("site", "rhodolith"), ]



#create a vector with all sample-IDs
sample_ids <- colnames(raw)[-1]



#create dataframe "meta" from the metadata
meta <- data.frame(
  sample    = sample_ids,
  station   = as.character(station_row[1, -1]),
  stringsAsFactors = FALSE
)



#create vector with group names
group_names <- count_data$group



#create a final data_frame for analysis
counts_mat <- as.data.frame(lapply(count_data[, -1], as.numeric))
rownames(counts_mat) <- group_names
counts_mat <- t(counts_mat)
rownames(counts_mat) <- sample_ids



#add the mean depth of the depth range of each station to the depth-table
depth_df <- depth_df %>%
  mutate(
    depth_mid = sapply(strsplit(gsub("m", "", depth), "-"), function(x) {
      mean(as.numeric(trimws(x)))
    })
  )



#3.calculate shannon-wiener-index of each sample



shannon_sample <- diversity(counts_mat, index = "shannon")



#add the shannon-wiener index of each sample to the metadata-dataframe "meta"
meta$shannon_H <- shannon_sample



#4.calculate beta-diversity with whittaker´s turnover index for each station



beta_station <- meta %>%
  group_by(station) %>%
  group_modify(~ {
    
    #subset count matrix to the samples each station
    subs <- counts_mat[.x$sample, , drop = FALSE]
    
    #calculate the groups observed at least once in each station as gamma
    gamma <- sum(colSums(subs) > 0)
    
    #calculate alpha per sample as the number of groups with a count greater 0 to calculate mean alpha
    alpha_per_sample <- rowSums(subs > 0)
    mean_alpha <- mean(alpha_per_sample)
    
    #calculate whittaker´s turnover index
    beta_W <- (gamma / mean_alpha) - 1
    
    data.frame(
      n_samples  = nrow(subs),
      gamma_div  = gamma,
      mean_alpha = round(mean_alpha, 3),
      beta_W     = round(beta_W, 3)
    )
  }) %>%
  ungroup()



#add depth information to beta_station
beta_station <- beta_station %>%
  left_join(depth_df, by = "station") %>%
  arrange(depth_mid)



#5.calculate shannon-wiener diversity index per station



shannon_station_list <- list()

for (i in unique(meta$station)) {
  subs <- counts_mat[meta$station == i, , drop = FALSE]
  mean_counts <- colMeans(subs)
  H_station <- diversity(matrix(mean_counts, nrow = 1), index = "shannon")
  shannon_station_list[[i]] <- H_station[1]
}



#create dataframe "station_diversity" to add metadata and station alpha-diversity data together
station_diversity <- data.frame(
  station   = names(shannon_station_list),
  shannon_H = unlist(shannon_station_list),
  stringsAsFactors = FALSE
) %>%
  left_join(depth_df, by = "station") %>%
  arrange(depth_mid)



#calculate the median H´per station
station_medians <- aggregate(shannon_H ~ station, data = meta, FUN = median)
print(station_medians)



#6.plot the data



#boxplot of sample-level shannon-wiener indexes ordered by station and depth



station_order <- station_diversity$station

meta_plot <- meta %>%
  left_join(depth_df, by = "station") %>%
  mutate(station = factor(station, levels = station_order))

p1 <- ggplot(meta_plot, aes(x = station, y = shannon_H, fill = depth_mid)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, colour = "grey30") +
  scale_fill_gradient(low = "#cce5ff", high = "#003f7f",
                      name = "Mid-depth (m)") +
  labs(
    title    = "Sample-level Shannon-Wiener Index by Station",
    x        = "     Station  
    depth range (m)",
    y        = "Shannon-Wiener Index H'",
  ) +
  scale_x_discrete(
    labels = setNames(
      paste0(station_order, "\n(", station_diversity$depth, ")"),
      station_order
    )
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p1)



#plot of Station-level Shannon-Wiener Index vs water depth



p2 <- ggplot(station_diversity,
             aes(x = depth_mid, y = shannon_H, label = station)) +
  geom_line(colour = "steelblue", linewidth = 0.8) +
  geom_point(size = 4, colour = "steelblue") +
  ggrepel::geom_label_repel(size = 3, nudge_y = 0.02,
                            box.padding = 0.3) +
  labs(
    title = "Station-level Shannon-Wiener Index vs water depth",
    x     = "Mid depth (m)",
    y     = "Shannon-Wiener Index H' (from mean counts)"
  ) +
  theme_bw(base_size = 12)

print(p2)



#plot Beta-Diversity (Whittaker) per Station vs water depth



p3 <- ggplot(beta_station,
             aes(x = depth_mid, y = beta_W, label = station)) +
  geom_line(colour = "darkorange", linewidth = 0.8) +
  geom_point(size = 4, colour = "darkorange") +
  ggrepel::geom_label_repel(size = 3, nudge_y = 0.02,
                            box.padding = 0.3) +
  labs(
    title = "Beta-Diversity (Whittaker) per Station vs water depth",
    x     = "Mid-point depth (m)",
    y     = expression(paste("Whittaker's  ", beta[W]))
  ) +
  theme_bw(base_size = 12)

print(p3)



#7.statistical analysis



#shapiro test with results stored in table "normality_results"
normality_results <- meta %>%
  group_by(station) %>%
  summarise(
    n        = n(),
    mean_H   = round(mean(shannon_H), 3),
    sd_H     = round(sd(shannon_H),   3),
    SW_W     = ifelse(n >= 3,
                      round(shapiro.test(shannon_H)$statistic, 4),
                      NA_real_),
    SW_p     = ifelse(n >= 3,
                      round(shapiro.test(shannon_H)$p.value,   4),
                      NA_real_),
    normal   = case_when(
      is.na(SW_p)  ~ "too few samples",
      SW_p > 0.05  ~ "yes (p > 0.05)",
      TRUE         ~ "no (p ≤ 0.05)"
    )
  ) %>%
  left_join(depth_df %>% select(station, depth, depth_mid), by = "station") %>%
  arrange(depth_mid)



#kruskal-wallis test



kw_result <- kruskal.test(shannon_H ~ factor(station), data = meta)
print(kw_result)


