# Generate the claw-fieldbook demo data (run from anywhere: Rscript demo/make_demo.R).
# Writes demo/fieldbook_database_export.csv: simulated observations for the
# 90-plot alpha design in demo/design_field_book.csv, in the Field Book app
# "database" export layout (one row per observation), with deliberate errors
# that the QC must catch (listed in demo/README.md). Seed 2026.

demo_dir <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[[1L]]))) else "demo"
})
set.seed(2026)
d <- utils::read.csv(file.path(demo_dir, "design_field_book.csv"), stringsAsFactors = FALSE)

ent <- sort(unique(d$entry))
g_y <- stats::setNames(stats::rnorm(length(ent), 0, 300), ent)
g_h <- stats::setNames(stats::rnorm(length(ent), 0, 8), ent)
g_f <- stats::setNames(round(stats::rnorm(length(ent), 0, 3)), ent)
trend <- 150 * (d$col - mean(d$col)) / sd(d$col) - 120 * (d$row - mean(d$row)) / sd(d$row)
yield <- round(3500 + g_y[d$entry] + trend + stats::rnorm(nrow(d), 0, 250))
height <- round(85 + g_h[d$entry] + stats::rnorm(nrow(d), 0, 5), 1)
flower <- as.Date("2026-07-10") + g_f[d$entry]

stamp <- function(day, i) sprintf("%s %02d:%02d:%02d.000-05:00", day, 9 + i %/% 60, i %% 60, (7 * i) %% 60)
mk <- function(trait, value, day) {
  data.frame(plot_id = d$plot_id, trait = trait, value = as.character(value),
             timestamp = stamp(day, seq_len(nrow(d))), person = "demo user",
             location = "", number = 1L, stringsAsFactors = FALSE)
}
obs <- rbind(mk("yield_kg_ha", yield, "2026-10-02"),
             mk("plant_height_cm", height, "2026-08-20"),
             mk("flowering_date", format(flower), as.character(flower)))
# flowering is recorded on the day it happens
fl <- obs$trait == "flowering_date"
obs$timestamp[fl] <- paste0(obs$value[fl], " 11:00:00.000-05:00")

# --- deliberate errors -------------------------------------------------------
obs <- obs[!obs$plot_id %in% c("DEMO-17", "DEMO-58"), ]                     # 2 plots never observed
obs$value[obs$plot_id == "DEMO-5" & obs$trait == "yield_kg_ha"] <- "NA"     # missing value
i <- which(obs$plot_id == "DEMO-12" & obs$trait == "yield_kg_ha")
obs$value[i] <- as.character(as.numeric(obs$value[i]) + 3000)               # digit typo (3xxx -> 6xxx)
obs$value[obs$plot_id == "DEMO-33" & obs$trait == "plant_height_cm"] <- "-8.5"   # impossible value
obs$value[obs$plot_id == "DEMO-40" & obs$trait == "plant_height_cm"] <- "35,5"   # decimal comma
obs$value[obs$plot_id == "DEMO-61" & obs$trait == "flowering_date"] <- "2026-04-30" # before planting
obs$timestamp[obs$plot_id == "DEMO-70" & obs$trait == "yield_kg_ha"] <- "2099-10-02 10:00:00.000-05:00" # future
dup <- obs[obs$plot_id == "DEMO-22" & obs$trait == "yield_kg_ha", ]
dup$value <- as.character(as.numeric(dup$value) + 400)                      # conflicting re-entry
obs$value[obs$plot_id == "DEMO-48" & obs$trait == "yield_kg_ha"] <- "35000"   # decimal shift (x10)
dup$timestamp <- "2026-10-03 08:00:00.000-05:00"
same <- obs[obs$plot_id == "DEMO-23" & obs$trait == "plant_height_cm", ]   # identical re-entry
unknown <- obs[obs$plot_id == "DEMO-30" & obs$trait == "yield_kg_ha", ]
unknown$plot_id <- "DEMO-999"                                              # plot not in the design
obs <- rbind(obs, dup, same, unknown)

utils::write.csv(obs, file.path(demo_dir, "fieldbook_database_export.csv"), row.names = FALSE, na = "")
utils::write.csv(data.frame(trait = c("yield_kg_ha", "plant_height_cm"), min = c(0, 20), max = c(8000, 150)),
                 file.path(demo_dir, "limits.csv"), row.names = FALSE)
cat("Wrote", nrow(obs), "observations.\n")
