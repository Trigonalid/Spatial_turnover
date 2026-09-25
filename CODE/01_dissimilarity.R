
#      Spatial turnover amplifies with trophic level
#                in hyperdiverse food webs
#
#              Spatial turnover of species
#
# Computes pairwise spatial dissimilarity between sites for
#   caterpillar and parasitoid communities using three indices
#   (Bray-Curtis, Chao-Sorensen, Sorensen), with subsampling of
#   caterpillars to parasitoid sample size.
#
# Inputs:  DATA/MASTER.xlsx
# Outputs: dissimilarity summary table, Fig 1C, Fig S4
#
# Structure: STATISTICS first, then all FIGURES at the end
#   (below the FIGURES divider).
#
#----------------------------------------------------------#
# 0. Libraries -----
#----------------------------------------------------------#
library(dplyr)       # data wrangling
library(tidyr)       # pivot_wider / pivot_longer
library(purrr)       # map functions for downsampling loop
library(tibble)      # tibble, column_to_rownames
library(readxl)      # read_excel
library(here)        # here() for file paths
library(reshape2)    # melt() for matrices
library(vegan)       # vegdist()
library(betapart)    # beta.pair() for Sorensen
library(CommEcol)    # dis.chao() for Chao-Sorensen
library(ggplot2)     # plotting
library(ggpubr)      # ggarrange(), stat_compare_means()
library(flextable)   # flextable() for Word table export
library(officer)     # read_docx() for Word export

#----------------------------------------------------------#
# 1. Global parameters -----                    # -> 00_setup
#----------------------------------------------------------#
set.seed(1234)
n_rand <- 1000   # subsampling iterations (matches Methods)

#----------------------------------------------------------#
# 2. Load data -----                         # -> 01_load_data
#----------------------------------------------------------#
# MASTER: one row = one reared caterpillar.
#   CAT_sp = caterpillar species; PAR_sp = parasitoid species
#   (empty if not parasitised). guild: "CAT" = unparasitised,
#   "PAR" = parasitised host.
# Ohu2 excluded (temporal replicate, handled in separate script).
MASTER <- read_excel(here("DATA/MASTER.xlsx")) %>%
  as_tibble() %>%
  filter(locality != "Ohu2")

#----------------------------------------------------------#
# 3. Community matrices -----             # -> 02_prepare_data
#----------------------------------------------------------#
# --- Parasitoids ---
# table() drops empty PAR_sp automatically, so this counts
#   parasitoids only.
para_matrix    <- as.matrix(table(MASTER$locality, MASTER$PAR_sp))
para_matrix_PA <- (para_matrix > 0) * 1L   # presence/absence for Sorensen

# --- Caterpillars ---
# Built from the full MASTER, i.e. the whole caterpillar
#   community (both CAT and PAR host individuals).
cat_matrix    <- as.matrix(table(MASTER$locality, MASTER$CAT_sp))
cat_matrix_PA <- (cat_matrix > 0) * 1L     # presence/absence for Sorensen

#----------------------------------------------------------#
# 4. Helper: melt dissimilarity matrix to long -----  # -> functions
#    Converts a square dissimilarity matrix to a long tibble
#    and removes self-comparisons (diagonal zeros).
#----------------------------------------------------------#
matrix_to_long <- function(mat, guild_name, index_name) {
  reshape2::melt(mat, varnames = c("Locality_A", "Locality_B"),
                 value.name = "hodnota") %>%
    as_tibble() %>%
    mutate(
      Locality_A = as.character(Locality_A),
      Locality_B = as.character(Locality_B),
      guild      = guild_name,
      indexy     = index_name
    ) %>%
    filter(hodnota > 0)   # removes self-pairs (diagonal = 0)
}

#----------------------------------------------------------#
# 5. Helper: dissimilarity summary row -----          # -> functions
#----------------------------------------------------------#
summary_row <- function(index_name, guild_name, values_vec,
                        dataset = "All data") {
  data.frame(
    index     = index_name,
    Community = guild_name,
    mean      = round(mean(values_vec, na.rm = TRUE), 3),
    sd        = round(sd(values_vec,   na.rm = TRUE), 3),
    Dataset   = dataset,
    stringsAsFactors = FALSE
  )
}

#==========================================================#
#                                                          #
#                      STATISTICS                          #
#                                                          #
#==========================================================#

#----------------------------------------------------------#
# 7. Bray-Curtis dissimilarity -----
#    Abundance-based. {vegan} vegdist(method = "bray").
#----------------------------------------------------------#
bc_para <- vegan::vegdist(para_matrix, method = "bray") %>% as.matrix()
bc_cat  <- vegan::vegdist(cat_matrix,  method = "bray") %>% as.matrix()

BC_par <- matrix_to_long(bc_para, "Parasitoids",  "Bray-Curtis")
BC_cat <- matrix_to_long(bc_cat,  "Caterpillars", "Bray-Curtis")

result_bc_para <- summary_row("Bray-Curtis", "Parasitoids",  BC_par$hodnota)
result_bc_cat  <- summary_row("Bray-Curtis", "Caterpillars", BC_cat$hodnota)

wilcox_bc <- wilcox.test(BC_par$hodnota, BC_cat$hodnota)
print(wilcox_bc)

#----------------------------------------------------------#
# 8. Chao-Sorensen dissimilarity -----
#    Corrects for unsampled rare species. {CommEcol} dis.chao().
#    56% of parasitoid species were singletons or doubletons.
#----------------------------------------------------------#
cs_para <- CommEcol::dis.chao(para_matrix, index = "sorensen",
                              version = "rare") %>% as.matrix()
cs_cat  <- CommEcol::dis.chao(cat_matrix,  index = "sorensen",
                              version = "rare") %>% as.matrix()

CS_par <- matrix_to_long(cs_para, "Parasitoids",  "Chao-Sorensen")
CS_cat <- matrix_to_long(cs_cat,  "Caterpillars", "Chao-Sorensen")

result_cs_para <- summary_row("Chao-Sorensen", "Parasitoids",  CS_par$hodnota)
result_cs_cat  <- summary_row("Chao-Sorensen", "Caterpillars", CS_cat$hodnota)

wilcox_cs <- wilcox.test(CS_par$hodnota, CS_cat$hodnota)
print(wilcox_cs)

#----------------------------------------------------------#
# 9. Sorensen dissimilarity -----
#    Presence/absence based. {betapart} beta.pair().
#    Baseline comparison; less robust with many rare species.
#----------------------------------------------------------#
sor_para <- betapart::beta.pair(para_matrix_PA)$beta.sor %>% as.matrix()
sor_cat  <- betapart::beta.pair(cat_matrix_PA)$beta.sor  %>% as.matrix()

SOR_par <- matrix_to_long(sor_para, "Parasitoids",  "Sorensen")
SOR_cat <- matrix_to_long(sor_cat,  "Caterpillars", "Sorensen")

result_sor_para <- summary_row("Sorensen", "Parasitoids",  SOR_par$hodnota)
result_sor_cat  <- summary_row("Sorensen", "Caterpillars", SOR_cat$hodnota)

wilcox_sor <- wilcox.test(SOR_par$hodnota, SOR_cat$hodnota)
print(wilcox_sor)

#----------------------------------------------------------#
# 10. Downsampling -- caterpillars to parasitoid sample size --
#    At each site, the full caterpillar community (both CAT and
#    PAR host individuals) is randomly resampled WITH replacement
#    to match the number of parasitoid rearings at that site.
#    Repeated n_rand times; mean dissimilarity reported.
#----------------------------------------------------------#
# Number of parasitoid rearings per site (target sample size)
n_par_per_site <- MASTER %>%
  filter(guild == "PAR") %>%
  drop_na(PAR_sp) %>%
  group_by(locality) %>%
  summarise(N = n(), .groups = "drop")

site_levels <- unique(MASTER$locality)

# Pre-allocate result lists for all three indices
res_bc  <- vector("list", n_rand)
res_cs  <- vector("list", n_rand)
res_sor <- vector("list", n_rand)

for (i in seq_len(n_rand)) {
  
  # Resample the full caterpillar community at each site to match
  #   the parasitoid N (both CAT and PAR rows are eligible).
  resampled <- purrr::map_df(site_levels, function(loc) {
    
    N_target <- n_par_per_site %>%
      filter(locality == loc) %>%
      pull(N)
    
    if (length(N_target) == 0 || is.na(N_target) || N_target == 0) {
      return(tibble())
    }
    
    MASTER %>%
      filter(locality == loc) %>%
      sample_n(size = N_target, replace = TRUE)   # with replacement
  })
  
  # Build wide abundance matrix from resampled data
  df_wide <- resampled %>%
    select(locality, CAT_sp) %>%
    group_by(locality, CAT_sp) %>%
    summarise(N = n(), .groups = "drop") %>%
    pivot_wider(names_from = CAT_sp, values_from = N, values_fill = 0) %>%
    column_to_rownames("locality")
  
  df_pa <- (df_wide > 0) * 1L   # presence/absence for Sorensen
  
  res_bc[[i]]  <- vegan::vegdist(df_wide, method = "bray") %>% as.matrix()
  res_cs[[i]]  <- CommEcol::dis.chao(df_wide, index = "sorensen",
                                     version = "rare") %>% as.matrix()
  res_sor[[i]] <- betapart::beta.pair(df_pa)$beta.sor %>% as.matrix()
}

# Average across iterations (element-wise mean of 3D array)
mean_matrix <- function(res_list) {
  arr <- array(
    do.call(cbind, res_list),
    dim = c(dim(res_list[[1]]), length(res_list))
  )
  apply(arr, c(1, 2), mean, na.rm = TRUE)
}

bc_sub_mean  <- mean_matrix(res_bc)
cs_sub_mean  <- mean_matrix(res_cs)
sor_sub_mean <- mean_matrix(res_sor)

# Assign locality names to averaged matrices
loc_names <- rownames(para_matrix)
dimnames(bc_sub_mean)  <- list(loc_names, loc_names)
dimnames(cs_sub_mean)  <- list(loc_names, loc_names)
dimnames(sor_sub_mean) <- list(loc_names, loc_names)

# Convert to long format for plotting
BC_sub  <- matrix_to_long(bc_sub_mean,  "Caterpillars-subsampled", "Bray-Curtis")
CS_sub  <- matrix_to_long(cs_sub_mean,  "Caterpillars-subsampled", "Chao-Sorensen")
SOR_sub <- matrix_to_long(sor_sub_mean, "Caterpillars-subsampled", "Sorensen")

result_bc_sub  <- summary_row("Bray-Curtis",   "Caterpillars-subsampled",
                              BC_sub$hodnota,  "Subsampled")
result_cs_sub  <- summary_row("Chao-Sorensen", "Caterpillars-subsampled",
                              CS_sub$hodnota,  "Subsampled")
result_sor_sub <- summary_row("Sorensen",      "Caterpillars-subsampled",
                              SOR_sub$hodnota, "Subsampled")

#----------------------------------------------------------#
# 11. Summary table -- dissimilarity indices -----
#     Exported as Word document (Table Sxx -- final number TBD).
#----------------------------------------------------------#
all_results <- bind_rows(
  result_bc_para,  result_bc_cat,  result_bc_sub,
  result_cs_para,  result_cs_cat,  result_cs_sub,
  result_sor_para, result_sor_cat, result_sor_sub
) %>%
  select(index, Community, mean, sd, Dataset)
print(all_results)

ft <- flextable(all_results)
doc <- officer::read_docx() %>%
  officer::body_add_par("Table Sxx: Summary of dissimilarity indices",
                        style = "heading 1") %>%
  flextable::body_add_flextable(ft)
print(doc, target = "output/Table_Sxx_dissimilarity_summary.docx")

#==========================================================#
#                                                          #
#                        FIGURES                           #
#                                                          #
#==========================================================#

#----------------------------------------------------------#
# F0. Plotting parameters -- palette & theme -----  # -> 00_setup
#     (consistent across all scripts; used by figures only)
#----------------------------------------------------------#
# Factor ordering used across figures
desired_order <- c("Caterpillars", "Caterpillars-subsampled", "Parasitoids")

col_parasitoids       <- "#EDB439"
col_caterpillars      <- "#F25E0D"
col_caterpillars_sub  <- "#F018DA"  # subsampled caterpillars

fill_cols <- c(
  "Parasitoids"             = col_parasitoids,
  "Caterpillars"            = col_caterpillars,
  "Caterpillars-subsampled" = col_caterpillars_sub
)
shape_vals <- c(
  "Caterpillars"            = 21,   # filled circle with outline
  "Caterpillars-subsampled" = 22,   # filled square with outline
  "Parasitoids"             = 24    # filled triangle with outline
)
x_labels <- c(
  "Caterpillars",
  "Caterpillars\nsubsampled",
  "Parasitoids"
)

theme_pub_black <- theme_classic(base_size = 20) +
  theme(
    text                  = element_text(color = "black"),
    axis.line             = element_line(color = "black", linewidth = 0.8),
    axis.ticks            = element_line(color = "black", linewidth = 0.8),
    axis.text             = element_text(color = "black"),
    axis.title            = element_text(color = "black"),
    legend.text           = element_text(color = "black"),
    legend.title          = element_text(color = "black"),
    legend.box.background = element_rect(color = NA, fill = NA),
    plot.tag              = element_text(size = 28, face = "bold", color = "black")
  )

#----------------------------------------------------------#
# F1. Figure S4 -- violin plots of all three indices -----
#     Each panel shows dissimilarity distributions for
#     caterpillars, subsampled caterpillars, and parasitoids.
#----------------------------------------------------------#
figure_data <- bind_rows(
  BC_par,  BC_cat,  BC_sub,
  CS_par,  CS_cat,  CS_sub,
  SOR_par, SOR_cat, SOR_sub
) %>%
  mutate(guild = factor(guild, levels = desired_order))

my_comparisons <- list(
  c("Caterpillars", "Parasitoids"),
  c("Caterpillars", "Caterpillars-subsampled"),
  c("Parasitoids",  "Caterpillars-subsampled")
)

# Helper: one violin panel per index
plot_beta_panel <- function(index_name, y_label, tag_letter, ylim_upper = 1.01) {
  
  df <- figure_data %>% filter(indexy == index_name)
  
  ggplot(df, aes(x = guild, y = hodnota)) +
    geom_violin(aes(fill = guild), trim = TRUE, color = "black") +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(
      aes(shape = guild, fill = I(alpha("black", 0.4))),
      position    = position_jitter(0.2),
      size        = 2,
      color       = "black",
      stroke      = 0.6,
      show.legend = FALSE
    ) +
    stat_compare_means(
      comparisons = my_comparisons,
      method      = "wilcox.test",
      label       = "p.signif",
      color       = "black",
      tip.length  = 0.01,
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols   = c("***", "**", "*", "ns")
      )
    ) +
    scale_fill_manual(values = fill_cols) +
    scale_shape_manual(values = shape_vals) +
    scale_x_discrete(labels = x_labels) +
    coord_cartesian(ylim = c(0, ylim_upper)) +
    labs(x = "", y = y_label, tag = tag_letter) +
    guides(
      fill  = guide_legend(title = "Community"),
      shape = guide_legend(title = "Community")
    ) +
    theme_pub_black
}

BC_FIG  <- plot_beta_panel("Bray-Curtis",   "Bray\u2013Curtis dissimilarity",        "a", 1.01)
CS_FIG  <- plot_beta_panel("Chao-Sorensen", "Chao\u2013S\u00f8rensen dissimilarity", "b", 1.00)
SOR_FIG <- plot_beta_panel("Sorensen",      "S\u00f8rensen dissimilarity",           "c", 1.00)

Fig_S4 <- ggpubr::ggarrange(
  BC_FIG, CS_FIG, SOR_FIG,
  common.legend = TRUE,
  legend        = "bottom",
  nrow          = 1
); Fig_S4

ggsave("output/fig/Figure_S4_beta_diversity.pdf",
       Fig_S4, width = 16, height = 10, units = "in")
ggsave("output/fig/Figure_S4_beta_diversity.png",
       Fig_S4, width = 16, height = 10, units = "in", dpi = 300)
ggsave("output/fig/Figure_S4_beta_diversity.svg",
       Fig_S4, width = 16, height = 10, units = "in", dpi = 300)

#----------------------------------------------------------#
# F2. Figure 1C -- Bray-Curtis only -----
#     Caterpillars vs. parasitoids, two-group comparison
#     without subsampling.
#----------------------------------------------------------#
comparison_1c <- list(c("Caterpillars", "Parasitoids"))

fig1c_data <- bind_rows(BC_par, BC_cat) %>%
  mutate(guild = factor(guild, levels = c("Caterpillars", "Parasitoids")))

Fig_1C <- ggplot(fig1c_data, aes(x = guild, y = hodnota, fill = guild)) +
  geom_violin(trim = TRUE, color = "black") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(
    aes(shape = guild),
    position    = position_jitter(0.2),
    size        = 2,
    color       = "black",
    fill        = alpha("black", 0.4),
    stroke      = 0.6,
    show.legend = FALSE
  ) +
  stat_compare_means(
    comparisons = comparison_1c,
    method      = "wilcox.test",
    label       = "p.signif",
    color       = "black",
    tip.length  = 0.01,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols   = c("***", "**", "*", "ns")
    )
  ) +
  scale_fill_manual(values = c(
    "Caterpillars" = col_caterpillars,
    "Parasitoids"  = col_parasitoids
  )) +
  scale_shape_manual(values = c(
    "Caterpillars" = 21,
    "Parasitoids"  = 24
  )) +
  coord_cartesian(ylim = c(0, 1.01)) +
  labs(x = "", y = "Bray\u2013Curtis dissimilarity") +
  theme_pub_black +
  theme(legend.position = "none")

ggsave("output/fig/Figure_1C_BrayCurtis.pdf",
       Fig_1C, width = 6, height = 8, units = "in")
ggsave("output/fig/Figure_1C_BrayCurtis.png",
       Fig_1C, width = 6, height = 8, units = "in", dpi = 300)
ggsave("output/fig/Figure_1C_BrayCurtis.svg",
       Fig_1C, width = 6, height = 8, units = "in")