# Tests for claw-fieldbook. Run from anywhere: Rscript tests/test_fieldbook.R
# A synthetic design and synthetic collected data with planted errors in the
# Field Book "database" layout, a wide "table" layout and two merged files;
# every planted error must be caught with the right issue code and action.
# Also: export format rules, robust z known answer, duplicate policies, CLI demo.

test_file <- local({
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[[1L]]), winslash = "/") else NA_character_
})
if (is.na(test_file)) stop("Run this file with Rscript.")
SKILL_TEST_DIR <- dirname(dirname(test_file))
source(file.path(SKILL_TEST_DIR, "scripts", "fieldbook.R"))  # loads functions; main() does not run

n_checks <- 0L
check <- function(ok, what) {
  if (!isTRUE(ok)) stop("FAIL: ", what, call. = FALSE)
  n_checks <<- n_checks + 1L
  invisible(TRUE)
}
expect_error <- function(expr, pattern, what) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  check(!is.na(msg) && grepl(pattern, msg), paste(what, "| got:", msg))
}
run_cli <- function(args) {
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- file.path(SKILL_TEST_DIR, "scripts", "fieldbook.R")
  out <- suppressWarnings(system2(rscript, c(shQuote(script), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, output = out)
}
O <- function(...) utils::modifyList(DEFAULTS, list(...))
has_flag <- function(fl, plot, issue, trait = NULL) {
  any(fl$plot_id == plot & fl$issue == issue & (is.null(trait) | fl$trait == if (is.null(trait)) "" else trait))
}

tmp <- tempfile("claw_fb_")
dir.create(tmp)

# Synthetic design: 5 rows x 8 columns, 2 FILLER plots
design <- data.frame(trial = "T", plot_id = paste0("T-", 1:40), plot = 1:40,
                     row = rep(1:5, each = 8), col = rep(1:8, 5), rep = rep(1:2, each = 20),
                     entry = c(paste0("E", 1:19), "FILLER", paste0("E", 1:19), "FILLER"),
                     is_check = FALSE, is_filler = rep(c(rep(FALSE, 19), TRUE), 2),
                     `seed lot [g]` = "L1", check.names = FALSE, stringsAsFactors = FALSE)
design_path <- file.path(tmp, "design.csv")
utils::write.csv(design, design_path, row.names = FALSE)
d <- read_design(design_path, "plot_id")
check(sum(d$is_filler) == 2, "design read with FILLER flags")

# 1. Export to the Field Book import format
ex <- export_fieldbook(d)
check(identical(names(ex)[1:3], c("plot_id", "row", "col")), "export: unique ID, primary and secondary order come first")
check("seed_lot_g" %in% names(ex) && !any(grepl("[][`\"' ]", names(ex))), "export: headers stripped of [ ] ` \" ' and spaces")
check(!anyDuplicated(ex$plot_id) && nrow(ex) == 40, "export: one row per plot, unique IDs")
check(all(diff(ex$row * 100 + ex$col) > 0), "export: rows ordered by primary then secondary order")
bad <- d; bad$plot_id[3] <- "T/3"
expect_error(export_fieldbook(bad), "rejects '/'", "export: '/' in IDs is refused")
dupd <- design; dupd$plot_id[2] <- "T-1"
utils::write.csv(dupd, file.path(tmp, "dup_design.csv"), row.names = FALSE)
expect_error(read_design(file.path(tmp, "dup_design.csv"), "plot_id"), "not unique.*T-1", "design with duplicated plot IDs is refused")
same_pos <- d; same_pos$col[2] <- 1
expect_error(export_fieldbook(same_pos), "share the same row and col", "export: two plots at one position is refused")

# 2. Robust z: known answer
z <- robust_z(c(1:9, 100))
check(abs(z[10] - 94.5 / (1.4826 * 2.5)) < 1e-6, "robust z: (x - median) / (1.4826 MAD) = 25.5 for the planted 100")
check(all(is.na(robust_z(rep(5, 12)))), "robust z: undefined (NA) when MAD = 0")

# 3. Database-format data with planted errors
set.seed(42)
real <- d[!d$is_filler, ]
base_y <- round(3000 + rnorm(nrow(real), 0, 200))
db <- rbind(
  data.frame(plot_id = real$plot_id, trait = "yield_kg_ha", value = as.character(base_y),
             timestamp = "2026-09-30 10:00:00.000-05:00", person = "a", location = "", number = 1),
  data.frame(plot_id = real$plot_id, trait = "height_cm", value = as.character(round(80 + rnorm(nrow(real), 0, 4), 1)),
             timestamp = "2026-08-01 10:00:00.000-05:00", person = "a", location = "", number = 1),
  data.frame(plot_id = real$plot_id, trait = "flower_date", value = "2026-07-10",
             timestamp = "2026-07-10 10:00:00.000-05:00", person = "a", location = "", number = 1))
db <- db[db$plot_id != "T-5", ]                                                    # plot never observed
db$value[db$plot_id == "T-6" & db$trait == "yield_kg_ha"] <- "NA"                  # missing value
db$value[db$plot_id == "T-7" & db$trait == "yield_kg_ha"] <- "900000"              # out of range
db$value[db$plot_id == "T-8" & db$trait == "yield_kg_ha"] <- "6500"                # robust outlier
db$value[db$plot_id == "T-9" & db$trait == "height_cm"] <- "45,2"                  # decimal comma
db$value[db$plot_id == "T-10" & db$trait == "flower_date"] <- "2026-13-45"         # impossible date
db$value[db$plot_id == "T-11" & db$trait == "flower_date"] <- "2026-04-01"         # before planting
db$timestamp[db$plot_id == "T-12" & db$trait == "height_cm"] <- "2031-08-01 10:00:00.000-05:00"  # future
db$timestamp[db$plot_id == "T-13" & db$trait == "height_cm"] <- "2026-10-10 10:00:00.000-05:00"  # after season end, before today
db$timestamp[db$plot_id == "T-16" & db$trait == "height_cm"] <- "2026-01-10 10:00:00.000-05:00"  # before planting
db$timestamp[db$plot_id == "T-14" & db$trait == "height_cm"] <- "yesterday"                      # bad timestamp
db$value[db$plot_id == "T-15" & db$trait == "flower_date"] <- "2031-07-01"         # date value in the future
v3 <- db[db$plot_id == "T-3" & db$trait == "yield_kg_ha", ]
db$value[db$plot_id == "T-3" & db$trait == "yield_kg_ha"] <- "3000"
v3$value <- "3500"; v3$timestamp <- "2026-10-01 09:00:00.000-05:00"               # conflicting duplicate (later)
v4 <- db[db$plot_id == "T-4" & db$trait == "height_cm", ]                         # identical duplicate
unk <- v4; unk$plot_id <- "T-999"                                                  # unknown plot
fil <- v4; fil$plot_id <- "T-20"; fil$value <- "70"                                # data on a FILLER plot
db <- rbind(db, v3, v4, unk, fil)
db_path <- file.path(tmp, "fb_database.csv")
utils::write.csv(db, db_path, row.names = FALSE, na = "")

opt <- O(design = design_path, data = db_path, limits = "yield_kg_ha=0:10000,height_cm=20:200",
         planting_date = "2026-05-01", season_start = "2026-04-15", season_end = "2026-10-03", today = "2026-10-15")
res <- import_fieldbook(d, db_path, opt)
fl <- res$flags
cl <- res$clean
check(identical(unname(res$formats), "database"), "database format detected from trait/value columns")
check(identical(unname(res$types[c("yield_kg_ha", "height_cm", "flower_date")]), c("numeric", "numeric", "date")), "trait types inferred: numeric, numeric, date")
check(has_flag(fl, "T-999", "unknown_plot_id"), "caught: unknown plot ID")
check(has_flag(fl, "T-20", "data_on_filler", "height_cm") && !"T-20" %in% cl$plot_id, "caught: data on FILLER plot; fillers dropped from clean table")
check(has_flag(fl, "T-5", "plot_not_observed") && all(is.na(cl[cl$plot_id == "T-5", c("yield_kg_ha", "height_cm")])), "caught: plot never observed")
check(has_flag(fl, "T-6", "missing_value", "yield_kg_ha") && is.na(cl$yield_kg_ha[cl$plot_id == "T-6"]), "caught: value recorded as NA")
check(has_flag(fl, "T-7", "out_of_range", "yield_kg_ha") && is.na(cl$yield_kg_ha[cl$plot_id == "T-7"]), "caught: out-of-range value, set to NA")
check(has_flag(fl, "T-8", "mad_outlier", "yield_kg_ha") && cl$yield_kg_ha[cl$plot_id == "T-8"] == 6500, "caught: robust outlier, flagged but kept")
check(has_flag(fl, "T-9", "non_numeric_value", "height_cm") && grepl("decimal comma", fl$detail[fl$plot_id == "T-9" & fl$issue == "non_numeric_value"]), "caught: decimal comma in a numeric trait")
check(has_flag(fl, "T-10", "unparseable_date", "flower_date"), "caught: impossible date value")
check(has_flag(fl, "T-11", "before_planting", "flower_date") && is.na(cl$flower_date[cl$plot_id == "T-11"]), "caught: date value before planting, set to NA")
check(has_flag(fl, "T-12", "future_date", "height_cm"), "caught: observation timestamp in the future")
check(has_flag(fl, "T-13", "outside_season", "height_cm"), "caught: observation after the season window")
check(has_flag(fl, "T-16", "before_planting", "height_cm") && !has_flag(fl, "T-16", "outside_season", "height_cm"),
      "caught: observation before planting (reported once, not also as outside the season)")
check(has_flag(fl, "T-14", "bad_timestamp", "height_cm"), "caught: unreadable timestamp")
check(has_flag(fl, "T-15", "future_date", "flower_date") && is.na(cl$flower_date[cl$plot_id == "T-15"]), "caught: date value in the future, set to NA")
check(has_flag(fl, "T-3", "duplicate_conflict", "yield_kg_ha") && is.na(cl$yield_kg_ha[cl$plot_id == "T-3"]), "caught: conflicting duplicate, set to NA under the default policy")
check(has_flag(fl, "T-4", "duplicate_identical", "height_cm") && !is.na(cl$height_cm[cl$plot_id == "T-4"]), "identical duplicate collapsed, value kept")
check(nrow(cl) == 38 && !anyDuplicated(cl$plot_id), "clean table: one row per non-filler design plot")
check(cl$flower_date_dap[cl$plot_id == "T-1"] == as.integer(as.Date("2026-07-10") - as.Date("2026-05-01")), "days after planting computed for date traits")
check(grepl("yield_kg_ha:out_of_range", cl$qc_issues[cl$plot_id == "T-7"]), "qc_issues column names the flag for each plot")
untouched <- setdiff(real$plot_id, paste0("T-", c(3:16)))
check(all(cl$yield_kg_ha[match(untouched, cl$plot_id)] == as.numeric(db$value[match(paste(untouched, "yield_kg_ha"), paste(db$plot_id, db$trait))])),
      "values without problems pass through unchanged")
s <- res$summary
check(s$out_of_range[s$trait == "yield_kg_ha"] == 1 && s$mad_outliers[s$trait == "yield_kg_ha"] == 1, "summary counts the yield flags")

# 4. Duplicate policies and keep-out-of-range
r_last <- import_fieldbook(d, db_path, utils::modifyList(opt, list(duplicates = "last")))
check(r_last$clean$yield_kg_ha[r_last$clean$plot_id == "T-3"] == 3500, "--duplicates last keeps the most recent record")
r_first <- import_fieldbook(d, db_path, utils::modifyList(opt, list(duplicates = "first")))
check(r_first$clean$yield_kg_ha[r_first$clean$plot_id == "T-3"] == 3000, "--duplicates first keeps the earliest record")
r_mean <- import_fieldbook(d, db_path, utils::modifyList(opt, list(duplicates = "mean")))
check(r_mean$clean$yield_kg_ha[r_mean$clean$plot_id == "T-3"] == 3250, "--duplicates mean averages numeric records")
r_keep <- import_fieldbook(d, db_path, utils::modifyList(opt, list(keep_out_of_range = TRUE)))
check(r_keep$clean$yield_kg_ha[r_keep$clean$plot_id == "T-7"] == 900000, "--keep-out-of-range keeps the value (still flagged)")

# 5. Wide data (Field Book table export or spreadsheet) with an unusual ID column
wide <- data.frame(`Plot ID` = real$plot_id, yield_kg_ha = as.character(base_y), check.names = FALSE, stringsAsFactors = FALSE)
wide$yield_kg_ha[wide$`Plot ID` == "T-2"] <- "lodged"
wide <- rbind(wide, data.frame(`Plot ID` = c("T-1", "X-77"), yield_kg_ha = c("9999", "3000"), check.names = FALSE))
wide_path <- file.path(tmp, "wide.csv")
utils::write.csv(wide, wide_path, row.names = FALSE)
rw <- import_fieldbook(d, wide_path, O(today = "2026-10-15"))
check(identical(unname(rw$id_cols), "Plot ID"), "wide: plot ID column found by matching design IDs")
check(has_flag(rw$flags, "T-2", "non_numeric_value", "yield_kg_ha"), "wide: text in a numeric column caught")
check(has_flag(rw$flags, "T-1", "duplicate_conflict", "yield_kg_ha"), "wide: duplicated plot row caught")
check(has_flag(rw$flags, "X-77", "unknown_plot_id"), "wide: unknown plot caught")

# 6. Two Field Book exports (two tablets) merged; a repeat across files is caught
p1 <- file.path(tmp, "tablet1.csv"); p2 <- file.path(tmp, "tablet2.csv")
y <- db[db$trait == "yield_kg_ha" & !db$plot_id %in% c("T-999", "T-3"), ]
utils::write.csv(y[1:20, ], p1, row.names = FALSE, na = "")
y2 <- y[20:nrow(y), ]; y2$value[1] <- "1234"
utils::write.csv(y2, p2, row.names = FALSE, na = "")
rm2 <- import_fieldbook(d, c(p1, p2), O(today = "2026-10-15"))
check(has_flag(rm2$flags, y$plot_id[20], "duplicate_conflict", "yield_kg_ha") && sum(!is.na(rm2$clean$yield_kg_ha)) >= 30,
      "two files merged; conflicting record across files caught")

# 7. Bad inputs fail loudly
expect_error(import_fieldbook(d, db_path, O(traits = "biomass_g")), "biomass_g", "unknown trait named in the error")
expect_error(parse_limits("yield=10"), "trait=min:max", "malformed --limits is refused")
expect_error(import_fieldbook(d, db_path, O(planting_date = "15/05/2026")), "YYYY-MM-DD", "non-ISO planting date is refused")

# 8. End to end on the demo data
out <- file.path(tmp, "demo_out")
res <- run_cli(c("--demo", "--out", shQuote(out)))
check(res$status == 0L, paste("demo run exits cleanly:", paste(res$output, collapse = "\n")))
for (f in c("report.md", "tables/fieldbook_import.csv", "tables/clean_data.csv", "tables/qc_flags.csv",
            "tables/qc_summary.csv", "figures/qc_field_maps.png", "reproducibility/commands.sh", "reproducibility/sessionInfo.txt")) {
  check(file.exists(file.path(out, f)), paste("demo writes", f))
}
dfl <- utils::read.csv(file.path(out, "tables", "qc_flags.csv"), stringsAsFactors = FALSE)
planted <- list(c("DEMO-17", "plot_not_observed"), c("DEMO-58", "plot_not_observed"), c("DEMO-5", "missing_value"),
                c("DEMO-12", "mad_outlier"), c("DEMO-48", "out_of_range"), c("DEMO-33", "out_of_range"),
                c("DEMO-40", "non_numeric_value"), c("DEMO-61", "before_planting"), c("DEMO-70", "future_date"),
                c("DEMO-22", "duplicate_conflict"), c("DEMO-23", "duplicate_identical"), c("DEMO-999", "unknown_plot_id"))
for (p in planted) check(any(dfl$plot_id == p[1] & dfl$issue == p[2]), paste("demo: planted error caught:", p[1], p[2]))
# Perfectly clean data: no flags, run completes, clean table equals the input
clean_path <- file.path(tmp, "clean_wide.csv")
utils::write.csv(data.frame(plot_id = real$plot_id, yield_kg_ha = base_y), clean_path, row.names = FALSE)
cr <- run_cli(c("--mode", "import", "--design", shQuote(design_path), "--data", shQuote(clean_path), "--out", shQuote(file.path(tmp, "clean_out"))))
check(cr$status == 0L, paste("import of clean data exits cleanly:", paste(cr$output, collapse = "
")))
cfl <- utils::read.csv(file.path(tmp, "clean_out", "tables", "qc_flags.csv"))
ccl <- utils::read.csv(file.path(tmp, "clean_out", "tables", "clean_data.csv"))
check(nrow(cfl) == 0 && isTRUE(all.equal(as.numeric(ccl$yield_kg_ha), as.numeric(base_y))),
      sprintf("clean data: no flags (%d), values unchanged", nrow(cfl)))
bad <- run_cli(c("--mode", "import", "--design", shQuote(design_path), "--out", shQuote(file.path(tmp, "no_data"))))
check(bad$status != 0L && any(grepl("--data", bad$output)) && !dir.exists(file.path(tmp, "no_data")),
      "import without --data exits with a helpful error before creating any output")

unlink(tmp, recursive = TRUE)
cat(sprintf("PASS: %d checks (%s)\n", n_checks, SKILL))
