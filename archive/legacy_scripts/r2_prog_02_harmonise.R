# Harmonise TCGA-LUAD and the six external cohorts onto one gene space.
#
# The archived signature failed external validation (pooled HR 1.18, P = 0.25,
# I2 = 0.69) and three re-derivations produced near-disjoint gene sets. Two design
# faults are visible in the archive and are fixed here:
#
#   1. Only 4 of 6 signature genes were measurable in the validation cohorts. A
#      gene that cannot be measured where the model must be tested should never
#      have been eligible, so selection is restricted up front to the intersection
#      of all seven datasets.
#   2. GSE31210 -- the one cohort that reached significance -- carries an explicit
#      "exclude for prognosis analysis" flag and normal-lung samples, and GSE30219
#      carries non-tumour lung. Tumour-only filtering and the exclusion flag are
#      applied here.
#
# Endpoint is overall survival throughout; cohorts reporting only relapse are not
# substituted with relapse-free survival.
suppressMessages({library(Biobase)})
P <- "/data/s01011/arid1a-nsclc/data/processed"
OUT <- "/data/s01011/arid1a_prognostic_redo"
dir.create(OUT, showWarnings = FALSE)

msg <- function(...) cat(..., "\n", sep = "")

# ---------------------------------------------------------------- TCGA
expr <- readRDS(file.path(P, "04c_tcga_luad_vst_expr_dedup.rds"))
cl <- readRDS(file.path(P, "tcga-luad_clinical.rds"))
cl$bcr <- substr(rownames(cl), 1, 12)
if (!"submitter_id" %in% colnames(cl)) cl$submitter_id <- cl$bcr
sid <- substr(colnames(expr), 1, 12)
# primary tumour only: TCGA sample code 01
keep <- substr(colnames(expr), 14, 15) == "01"
expr <- expr[, keep]; sid <- sid[keep]

vs <- as.character(cl$vital_status)
dd <- suppressWarnings(as.numeric(as.character(cl$days_to_death)))
df <- suppressWarnings(as.numeric(as.character(cl$days_to_last_follow_up)))
ev <- ifelse(grepl("dead", vs, ignore.case = TRUE), 1, 0)
tm <- ifelse(ev == 1, dd, df)
key <- if ("submitter_id" %in% colnames(cl)) as.character(cl$submitter_id) else cl$bcr
surv <- data.frame(id = key, time = tm, event = ev, stringsAsFactors = FALSE)
surv <- surv[!is.na(surv$time) & surv$time > 0 & !is.na(surv$event), ]
surv <- surv[!duplicated(surv$id), ]

m <- match(sid, surv$id)
ok <- !is.na(m)
tcga <- list(x = expr[, ok, drop = FALSE],
             time = surv$time[m[ok]], event = surv$event[m[ok]])
msg("TCGA: ", ncol(tcga$x), " tumours with OS, ", nrow(tcga$x), " genes, ",
    sum(tcga$event), " deaths")

# ---------------------------------------------------------------- externals
# per-cohort field map; each entry gives the pData columns and the time unit
SPEC <- list(
  gse31210  = list(status = "death:ch1", dead = "dead",
                   time = "days before death/censor:ch1", unit = 1/30.44,
                   tissue = "tissue:ch1", tumour = "primary lung tumor",
                   drop_flag = "exclude for prognosis analysis due to incomplete resection or adjuvant therapy:ch1",
                   drop_val = "exclude", hist = "LUAD"),
  gse37745  = list(status = "dead:ch1", dead = "yes",
                   time = "days to determined death status:ch1", unit = 1/30.44,
                   hist = "Mixed"),
  gse50081  = list(status = "status:ch1", dead = "dead",
                   time = "survival time:ch1", unit = 12, hist = "Mixed"),
  gse72094  = list(status = "vital_status:ch1", dead = "Dead",
                   time = "survival_time_in_days:ch1", unit = 1/30.44, hist = "LUAD"),
  gse30219  = list(status = "status:ch1", dead = "DEAD",
                   time = "follow-up time (months):ch1", unit = 1,
                   tissue = "histology:ch1", tumour_not = "NTL", hist = "Mixed"),
  gse157009 = list(status = "os_event:ch1", dead = "1",
                   time = "os_mo:ch1", unit = 1, hist = "SqCC")
)

collapse_to_symbol <- function(e, fd) {
  # column naming differs across the GEO platforms ("Gene symbol", "Gene Symbol",
  # "GeneSymbol"), so match on the normalised name rather than a fixed string
  norm <- tolower(gsub("[^a-z]", "", tolower(colnames(fd))))
  sc <- colnames(fd)[norm == "genesymbol"]
  if (!length(sc)) sc <- colnames(fd)[grepl("^gene.?symbol$", colnames(fd), ignore.case = TRUE)]
  stopifnot(length(sc) >= 1)
  sc <- sc[1]
  sym <- as.character(fd[[sc]])
  sym <- sub("\\s*///.*$", "", sym)          # multi-mapping probes: take the first
  good <- !is.na(sym) & sym != "" & sym != "---"
  e <- e[good, , drop = FALSE]; sym <- sym[good]
  o <- order(sym, -rowMeans(e, na.rm = TRUE))  # keep the highest-expressed probe
  e <- e[o, , drop = FALSE]; sym <- sym[o]
  e <- e[!duplicated(sym), , drop = FALSE]
  rownames(e) <- sym[!duplicated(sym)]
  e
}

ext <- list()
for (nm in names(SPEC)) {
  s <- SPEC[[nm]]
  g <- readRDS(file.path(P, paste0(nm, "_gsets.rds")))[[1]]
  pd <- pData(g); e <- exprs(g)
  # `$` partial-matches on lists, so s$tumour silently returned s$tumour_not and the
  # two tissue filters cancelled each other out, leaving zero samples. Match exactly.
  gv <- function(k) if (k %in% names(s)) s[[k]] else NULL
  keep <- rep(TRUE, ncol(e))
  if (!is.null(gv("tissue")) && !is.null(gv("tumour")))
    keep <- keep & grepl(gv("tumour"), as.character(pd[[gv("tissue")]]), ignore.case = TRUE)
  if (!is.null(gv("tissue")) && !is.null(gv("tumour_not")))
    keep <- keep & !grepl(gv("tumour_not"), as.character(pd[[gv("tissue")]]), ignore.case = TRUE)
  if (!is.null(gv("drop_flag")) && gv("drop_flag") %in% colnames(pd))
    keep <- keep & !grepl(gv("drop_val"), as.character(pd[[gv("drop_flag")]]), ignore.case = TRUE)
  st <- as.character(pd[[gv("status")]])
  ev <- as.integer(grepl(paste0("^", gv("dead"), "$"), st, ignore.case = TRUE))
  tm <- suppressWarnings(as.numeric(as.character(pd[[gv("time")]]))) * gv("unit")
  keep <- keep & !is.na(tm) & tm > 0 & !is.na(ev)
  e <- collapse_to_symbol(e[, keep, drop = FALSE], fData(g))
  ext[[nm]] <- list(x = e, time = tm[keep], event = ev[keep], hist = s[["hist"]])
  msg(sprintf("%-10s %4d tumours, %5d genes, %3d deaths (%s)",
              nm, ncol(e), nrow(e), sum(ev[keep]), s[["hist"]]))
}

common <- Reduce(intersect, c(list(rownames(tcga$x)), lapply(ext, function(z) rownames(z$x))))
msg("\ngenes measurable in TCGA and ALL six cohorts: ", length(common))
saveRDS(list(tcga = tcga, ext = ext, common = common),
        file.path(OUT, "harmonised.rds"))
msg("saved ", file.path(OUT, "harmonised.rds"))
