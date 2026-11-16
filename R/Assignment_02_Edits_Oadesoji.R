## Assignment 02 - Edits by Opeyemi Adesoji (oadesoji)
## Edited from: assignment_01_group10_scriptfile.R (Group 10)
## 2026-11-15

## Edit 1: summarise_te() function replaces the duplicated summary
##         code; data paths fixed from "../data/" to "data/"


## Edit 2: chi-square tests of TE composition between species, at order
##         level and subclass level (simulated p-value for small counts)



##***************************
## Assignment 01 - Group 10 
## 2026-11-02
##
##***************************
## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
library(ggrepel)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

rm(list = ls())


##Import TE annotation datasets -----


## Import T. patagoniensis dataset
setwd("../illustrating_worflow/")

df_patagoniensis <- read_csv("data/GCA_001447655.1.csv")
df_patagoniensis 


## Import T. spiralis dataset

df_spiralis = read_csv("data/trichinella_spiralis_te.csv")

df_spiralis 

## Clean dataset columns -----

df_patagoniensis =
  df_patagoniensis |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", "")),
    no_bp_in_complement = as.numeric(str_replace_all(no_bp_in_complement, "\\(|\\)", ""))
  )

df_spiralis =
  df_spiralis |> 
  mutate(
    qry_left = as.numeric(str_replace_all(qry_left, "\\(|\\)", "")),
    in_repeat_end = as.numeric(str_replace_all(in_repeat_end, "\\(|\\)", ""))
  )

## Define non-TE classes to exclude -----


ls_non_te = c("simple_repeat", "Simple_repeat",
              "Low_complexity", "Satellite",
              "rRNA", "snRNA", "ARTEFACT",
              "Unknown")
ls_non_te

## Function to summarise TE annotations for one species -----
# Takes a RepeatMasker table and a species label, removes non-TE classes,
# and returns the number of hits and total bp for each TE class.

summarise_te <- function(df, species_name, exclude = ls_non_te) {
  df |>
    filter(!matching_class %in% exclude) |>   # keep TE classes only
    group_by(matching_class) |>
    summarise(
      annotation_count = n(),
      total_length = sum(qry_width, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(species = species_name) |>
    arrange(desc(annotation_count))}

## Summarize TE annotations for each species -----
# Uses summarise_te() defined above, so the same steps are applied
# to both species and any future species with one line each.

df_spiralis_summary = summarise_te(df_spiralis, "T. spiralis")
df_spiralis_summary

df_patagoniensis_summary = summarise_te(df_patagoniensis, "T. patagoniensis")
df_patagoniensis_summary






## Combine species summaries -----

df_comparison =
  bind_rows(df_spiralis_summary,
            df_patagoniensis_summary)
df_comparison

## Statistical test to check if the two species differ in TE composition? -----
# Group subclasses (e.g. "LTR/Gypsy", "LTR/Pao") into main TE orders ("LTR"),
# because many subclasses have fewer than 5 hits, which makes a chi-square    # test unreliable.

df_order_counts =
  df_comparison |>
  mutate(te_order = str_extract(matching_class, "^[^/]+")) |>
  group_by(species, te_order) |>
  summarise(annotation_count = sum(annotation_count), .groups = "drop")

df_order_counts

# RC and SINE are still too rare (expected counts below 5), so orders with
# fewer than 20 hits in total (both species combined) are pooled into "Other".

df_order_counts =
  df_order_counts |>
  group_by(te_order) |>
  mutate(order_total = sum(annotation_count)) |>
  ungroup() |>
  mutate(te_order = if_else(order_total < 20, "Other", te_order)) |>
  group_by(species, te_order) |>
  summarise(annotation_count = sum(annotation_count), .groups = "drop")

df_order_counts

# Reshape into a contingency table: rows = TE order, columns = species
mat_order_counts =
  df_order_counts |>
  pivot_wider(names_from = species, values_from = annotation_count) |>
  column_to_rownames("te_order") |>
  as.matrix()

mat_order_counts

# Chi-square test of independence:
# H0 = both species have the same proportions of TE orders
chisq_order = chisq.test(mat_order_counts)
chisq_order

# Check the test's assumption: all expected counts should be >= 5
round(chisq_order$expected, 1)

# Subclass-level test: most expected counts are < 5, so the chi-square
# approximation is unreliable. A Monte Carlo (simulated) p-value is valid
# for small counts.
mat_class =
  df_comparison |>
  select(matching_class, species, annotation_count) |>
  pivot_wider(names_from = species, values_from = annotation_count, values_fill = 0) |>
  column_to_rownames("matching_class") |>
  as.matrix()

set.seed(123)   # makes the simulated p-value reproducible
chisq_class = chisq.test(mat_class, simulate.p.value = TRUE, B = 10000)
chisq_class

# Which subclasses drive the difference? Standardised residuals:
# values beyond about +/-2 show classes with more/fewer hits than expected
round(chisq_class$stdres, 1)

sqrt(chisq_order$statistic / sum(mat_order_counts))

#==> Order level: no significant difference in TE composition
#==>   (chi-square = 4.07, df = 4, p = 0.40; all expected counts >= 5).
#==> Subclass level: significant difference (chi-square = 43.90,
#==>   simulated p < 0.001, 10,000 replicates), driven mainly by LTR/Gypsy
#==>   (more in T. patagoniensis) and LTR/Pao and DNA/PIF-Harbinger
#==>   (more in T. spiralis). Grouping into orders hides these opposing shifts


## Combine raw TE annotation datasets -----

df_TE_combined =
  bind_rows(
    df_spiralis |> mutate(species = "T. spiralis"),
    df_patagoniensis |> mutate(species = "T. patagoniensis")
  )
df_TE_combined


## Rank TE classes by annotation count -----

df_comparison |>
  group_by(species) |>
  arrange(desc(annotation_count), .by_group = TRUE)


## Rank TE classes by total annotated length -----

df_comparison |>
  group_by(species) |>
  arrange(desc(total_length), .by_group = TRUE)


# Top 3 TE classes by total length for each species
df_comparison |> 
  group_by(species) |> 
  arrange(desc(total_length), .by_group = TRUE) |> 
  slice_head(n = 3)

## Plot distribution of TE annotation lengths -----

# Edit 3a: log-scaled x-axis so short and long TEs are both visible;
# density (not counts) so species with different numbers of hits are
# comparable; species overlaid in colour-blind-friendly viridis colours.

df_TE_combined |>
  filter(!matching_class %in% ls_non_te) |>
  ggplot(aes(x = qry_width, fill = species)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, alpha = 0.5, position = "identity") +
  scale_x_log10() +
  scale_fill_viridis_d(end = 0.8) +
  labs(
    title = "Distribution of TE Annotation Lengths",
    x = "Annotated TE length (bp, log scale)",
    y = "Density",
    fill = "Species")


## Plot annotation counts against total annotated lengths -----

# Edit 3b: removed geom_smooth() (too few points; its curve went below zero,
# which is impossible for lengths); log scales spread out the points;
# colour shows TE order; labels name the main classes (50+ hits).

df_comparison |>
  mutate(te_order = str_extract(matching_class, "^[^/]+")) |>
  ggplot(aes(x = annotation_count, y = total_length, colour = te_order)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(
    aes(label = if_else(annotation_count >= 50, matching_class, "")),
    size = 3, show.legend = FALSE,
    min.segment.length = 0, max.overlaps = Inf
  ) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.20))) +
  scale_colour_viridis_d(end = 0.9) +
  facet_grid(~ species) +
  labs(
    title = "Comparison of TE Abundance & Length",
    x = "Number of annotations (log scale)",
    y = "Total annotated length (bp, log scale)",
    colour = "TE order"
  ) +
  theme(legend.position = "bottom")



## Plot average TE length by class -----
# Edit 3c: classes sorted by average length; species side by side for direct
# comparison; classes with fewer than 5 hits in a species are left out
# because averages of 1-4 hits are unreliable.

df_comparison |>
  mutate(mean_length = total_length / annotation_count) |>
  filter(annotation_count >= 5) |>
  ggplot(aes(x = reorder(matching_class, mean_length),
             y = mean_length, fill = species)) +
  geom_col(position = position_dodge(preserve = "single")) +
  coord_flip() +
  scale_fill_viridis_d(end = 0.8) +
  labs(
    title = "Average TE Length by Class",
    x = "TE class",
    y = "Average annotated length (bp)",
    fill = "Species",
    caption = "Only classes with at least 5 annotations in a species are shown"
  ) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))

## Edit 3: improved plots - histogram on log scale with density and species overlaid; scatter plot without misleading smooth line, log scales,TE order colours and ggrepel labels; average length plot sorted, species side by side, classes with < 5 hits excluded. Note: requires the ggrepel package.



