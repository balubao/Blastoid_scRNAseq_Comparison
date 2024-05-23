source("script/Initiation.R")


plotting.data <- readRDS("results/plotting_data_revision_2.RDS")
plotting.data$celltype <- plotting.data$annot_major

library(scales)
hex3 <- hue_pal()(3)
celltype_color <- c(hex3, "lightgrey")
names(celltype_color) <- c("EPI_ICM", "PE", "TE", "unknown")

ggplot(plotting.data, aes(x = umaprpca_1, y = umaprpca_2, color = celltype)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = celltype_color) +
    theme_classic()

# tmp_data <- plotting.data[plotting.data$rpca_clusters == 32, ]

## sample data for visualization (balanced)
prob_ct <- c(1, 1, 1, 1) / 4 # table(plotting.data$annot_major)/nrow(plotting.data)
names(prob_ct) <- c("EPI_ICM", "PE", "TE", "unknown")
prob_ct_vec <- prob_ct[match(plotting.data$annot_major, names(prob_ct))]
select.c <- sample(x = seq(nrow(plotting.data)), size = 30, prob = prob_ct_vec)
tmp_data <- plotting.data[select.c, ]

# tmp_data <- plotting.data[plotting.data$rpca_clusters == 32, ]
tmp_data$idu <- seq(nrow(tmp_data))

## annot sub ---
tmp.df <- tmp_data %>%
    select(idu, annot_scina, annot_sctype_clust, annot_snglr_clust_mrk) %>%
    pivot_longer(cols = c("annot_scina", "annot_sctype_clust", "annot_snglr_clust_mrk"))

p1 <- ggplot(tmp.df, aes(x = idu, y = name, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = celltype_color) +
    theme_minimal() +
    # theme_classic() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank()) +
    coord_fixed() +
    ylab("") +
    xlab("")

ggsave("03_annotation/schematic/annot_sub_tile.svg")

## consensus ---
tmp2.df <- tmp.df %>%
    group_by(idu, value) %>%
    summarise(total = n())
tmp2.df$total <- factor(tmp2.df$total, ordered = TRUE)

better_col_palette <- RColorBrewer::brewer.pal(4, name = "PuRd")

p2 <- ggplot(tmp2.df, aes(x = idu, y = value, fill = total)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = better_col_palette) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank()) +
    coord_fixed() +
    ylab("") +
    xlab("")

ggsave("03_annotation/schematic/consensus_tile.svg")

## annot major ---
tmp3.df <- tmp_data %>%
    select(idu, annot_major) %>%
    pivot_longer(cols = c("annot_major"))

p3 <- ggplot(tmp3.df, aes(x = idu, y = name, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = celltype_color) +
    theme_minimal() +
    # theme_classic() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank()) +
    coord_fixed() +
    ylab("") +
    xlab("")

ggsave("03_annotation/schematic/annot_major_tile.svg")

## Stack figures ----
(p1 + theme(plot.margin = unit(c(0, 0, 70, 0), "pt"))) /
    (p2 + theme(plot.margin = unit(c(0, 0, 70, 0), "pt"))) / p3
ggsave("03_annotation/schematic/stacked_schematic.svg")
