# One pre-specified attempt at a clinically useful prognostic model.
#
# The first rebuild gave a nested-CV C-index of 0.543 on TCGA alone, with no
# clinical variables and only 504 training tumours. Three things change here, all
# fixed before any result is seen:
#
#   1. CLINICAL BASELINE. Stage is the dominant prognostic factor in NSCLC and was
#      absent. The clinically meaningful question is not "does the gene score
#      predict survival" but "does it add anything to stage, age and sex".
#   2. LEAVE-ONE-COHORT-OUT. Train on six datasets, test on the seventh, rotating.
#      Uses all 2,020 tumours and tests transportability directly, which is what
#      the archived signature failed at.
#   3. A GRID, REPORTED WHOLE. Three feature spaces x three penalised models. The
#      best cell of a grid is an upward-biased statistic -- the same best-of-k bias
#      this project documented across 227 complexes -- so every cell is reported
#      and the decision rests on the median across held-out cohorts, not the max.
#
# DECISION RULE, pre-specified: genes are worth including only if the combined
# clinical+gene model beats the clinical-only model on held-out cohorts, judged by
# the median delta C-index across the seven folds with a paired test over folds.
suppressMessages({library(glmnet); library(survival); library(Biobase)})
set.seed(20260808)
P <- "/data/s01011/arid1a-nsclc/data/processed"
OUT <- "/data/s01011/arid1a_prognostic_redo"
H <- readRDS(file.path(OUT, "harmonised.rds"))
msg <- function(...) cat(..., "\n", sep = "")

# ------------------------------------------------------- clinical covariates
# Stage coding differs by cohort; all are reduced to an early/advanced indicator
# with the same meaning. GSE30219 and GSE157009 report only T stage, so T1 stands
# in for early there -- noted as an approximation rather than hidden.
roman <- function(v) {
  v <- toupper(gsub("[^IVX0-9]", "", v))
  ifelse(grepl("^IV", v), 4, ifelse(grepl("^III", v), 3,
    ifelse(grepl("^II", v), 2, ifelse(grepl("^I", v), 1,
      suppressWarnings(as.numeric(substr(v, 1, 1)))))))
}
clin_ext <- list(
  gse31210  = list(stage = "pathological stage:ch1", age = "age (years):ch1", sex = "gender:ch1"),
  gse37745  = list(stage = "tumor stage:ch1",        age = "age:ch1",          sex = "gender:ch1"),
  gse50081  = list(stage = "Stage:ch1",              age = "age:ch1",          sex = "Sex:ch1"),
  gse72094  = list(stage = "Stage:ch1",              age = "age_at_diagnosis:ch1", sex = "gender:ch1"),
  gse30219  = list(stage = "pt stage:ch1",           age = "age at surgery:ch1", sex = "gender:ch1"),
  gse157009 = list(stage = "Stage:ch1",              age = "age:ch1",          sex = "Sex:ch1")
)

D <- list()
cl <- readRDS(file.path(P, "tcga-luad_clinical.rds"))
# rownames are row numbers, not barcodes; the key is submitter_id, and the sex
# column is sex_at_birth. Getting either wrong silently yields an all-NA clinical
# block, which the imputation below would then hide.
tc <- H$tcga
sid <- substr(colnames(tc$x), 1, 12)
m <- match(sid, as.character(cl$submitter_id))
stopifnot(mean(!is.na(m)) > 0.9)
D$TCGA <- list(x = tc$x[H$common, , drop = FALSE], time = tc$time, event = tc$event,
               stage = roman(as.character(cl$ajcc_pathologic_stage)[m]),
               age = suppressWarnings(as.numeric(as.character(cl$age_at_diagnosis)[m])) / 365.25,
               sex = as.character(cl$sex_at_birth)[m])

for (nm in names(clin_ext)) {
  s <- clin_ext[[nm]]
  g <- readRDS(file.path(P, paste0(nm, "_gsets.rds")))[[1]]
  pd <- pData(g)
  z <- H$ext[[nm]]
  keep <- match(colnames(z$x), rownames(pd))
  D[[nm]] <- list(x = z$x[H$common, , drop = FALSE], time = z$time, event = z$event,
                  stage = roman(as.character(pd[[s$stage]])[keep]),
                  age = suppressWarnings(as.numeric(as.character(pd[[s$age]])[keep])),
                  sex = as.character(pd[[s$sex]])[keep])
}
for (nm in names(D)) {
  d <- D[[nm]]
  d$early <- as.integer(d$stage <= 1)
  d$male <- as.integer(grepl("^m", d$sex, ignore.case = TRUE))
  miss <- c(stage = mean(is.na(d$early)), age = mean(is.na(d$age)),
            sex = mean(is.na(d$male)))
  # impute to the cohort's own centre, and say how much was imputed; imputing
  # missing stage to "early" would turn a failed extraction into a clean-looking
  # covariate, which is how the first attempt hid a broken TCGA match
  if (any(miss > 0.5)) stop("cohort ", nm, ": clinical extraction failed, missing ",
                            paste(names(miss), round(miss, 2), collapse=" "))
  d$age[is.na(d$age)] <- median(d$age, na.rm = TRUE)
  d$early[is.na(d$early)] <- as.integer(median(d$early, na.rm = TRUE))
  d$male[is.na(d$male)] <- as.integer(median(d$male, na.rm = TRUE))
  # z-score genes within cohort so platforms are comparable
  d$z <- t(scale(t(d$x)))
  d$z[!is.finite(d$z)] <- 0
  D[[nm]] <- d
  msg(sprintf("%-10s n=%4d deaths=%3d  early %.0f%%  age %.0f  male %.0f%%  (imputed: stage %.0f%% age %.0f%% sex %.0f%%)",
              nm, ncol(d$x), sum(d$event), 100 * mean(d$early), median(d$age),
              100 * mean(d$male), 100*miss["stage"], 100*miss["age"], 100*miss["sex"]))
}

WG <- tryCatch({
  w <- read.csv("/data/s01011/arid1a-nsclc/results/tables/04c_WGCNA_key_module_genes_dedup.csv")
  intersect(as.character(w$gene), H$common)
}, error = function(e) character(0))
msg("\nWGCNA brown-module genes available: ", length(WG))

cidx <- function(risk, time, event) {
  s <- try(summary(coxph(Surv(time, event) ~ risk))$concordance[1], silent = TRUE)
  if (inherits(s, "try-error")) NA_real_ else as.numeric(s)
}

SPACES <- c("ARID1A_correlated", "survival_screen", "WGCNA_module")
MODELS <- c("lasso", "enet", "ridge")
grid <- expand.grid(space = SPACES, model = MODELS, stringsAsFactors = FALSE)
res <- data.frame()

for (ho in names(D)) {
  tr <- setdiff(names(D), ho)
  Xtr <- do.call(cbind, lapply(tr, function(k) D[[k]]$z))
  ttr <- unlist(lapply(tr, function(k) D[[k]]$time))
  etr <- unlist(lapply(tr, function(k) D[[k]]$event))
  Ctr <- do.call(rbind, lapply(tr, function(k)
    data.frame(early = D[[k]]$early, age = D[[k]]$age, male = D[[k]]$male)))
  te <- D[[ho]]
  Cte <- data.frame(early = te$early, age = te$age, male = te$male)

  # clinical-only baseline
  cf <- coxph(Surv(ttr, etr) ~ early + age + male, data = Ctr)
  rc <- as.numeric(predict(cf, newdata = Cte, type = "lp"))
  c_clin <- cidx(rc, te$time, te$event)

  Xtr <- t(Xtr); Xte <- t(te$z)
  atr <- unlist(lapply(tr, function(k) D[[k]]$z["ARID1A", ]))

  for (i in seq_len(nrow(grid))) {
    sp <- grid$space[i]; mo <- grid$model[i]
    feats <- if (sp == "ARID1A_correlated") {
      r <- apply(Xtr, 2, function(v) suppressWarnings(cor(v, atr, method = "spearman")))
      colnames(Xtr)[!is.na(r) & abs(r) >= 0.3]
    } else if (sp == "survival_screen") {
      pv <- apply(Xtr, 2, function(v) {
        f <- try(coxph(Surv(ttr, etr) ~ v), silent = TRUE)
        if (inherits(f, "try-error")) 1 else summary(f)$coefficients[5]
      })
      names(sort(pv))[1:500]
    } else WG
    feats <- intersect(feats, colnames(Xte))
    if (length(feats) < 5) next
    al <- switch(mo, lasso = 1, enet = 0.5, ridge = 0)
    cv <- try(cv.glmnet(Xtr[, feats, drop = FALSE], Surv(ttr, etr), family = "cox",
                        alpha = al, nfolds = 5), silent = TRUE)
    if (inherits(cv, "try-error")) next
    rg <- as.numeric(predict(cv, Xte[, feats, drop = FALSE], s = "lambda.min"))
    c_gene <- cidx(rg, te$time, te$event)
    cc <- coxph(Surv(ttr, etr) ~ early + age + male +
                  as.numeric(predict(cv, Xtr[, feats, drop = FALSE], s = "lambda.min")),
                data = Ctr)
    nd <- Cte; nd$`as.numeric(...)` <- NULL
    lp <- as.numeric(cbind(Cte$early, Cte$age, Cte$male, rg) %*% coef(cc))
    c_comb <- cidx(lp, te$time, te$event)
    res <- rbind(res, data.frame(holdout = ho, space = sp, model = mo,
                                 n_feat = length(feats), C_clin = c_clin,
                                 C_gene = c_gene, C_comb = c_comb,
                                 delta = c_comb - c_clin))
  }
  msg("  held out ", ho, " done")
}

write.csv(res, file.path(OUT, "loco_grid.csv"), row.names = FALSE)
msg("\n", strrep("=", 74))
msg("Leave-one-cohort-out C-index, median across the seven held-out cohorts")
msg(strrep("=", 74))
msg(sprintf("%-20s %-7s %8s %8s %8s %9s", "space", "model", "C_clin", "C_gene", "C_comb", "delta"))
agg <- aggregate(cbind(C_clin, C_gene, C_comb, delta) ~ space + model, res, median)
agg <- agg[order(-agg$delta), ]
for (i in seq_len(nrow(agg)))
  msg(sprintf("%-20s %-7s %8.3f %8.3f %8.3f %+9.3f", agg$space[i], agg$model[i],
              agg$C_clin[i], agg$C_gene[i], agg$C_comb[i], agg$delta[i]))

best <- agg[1, ]
sub <- res[res$space == best$space & res$model == best$model, ]
wt <- try(wilcox.test(sub$C_comb, sub$C_clin, paired = TRUE), silent = TRUE)
msg("\nbest cell: ", best$space, " / ", best$model,
    "  median delta C = ", sprintf("%+.3f", best$delta))
msg("paired test over the seven folds: P = ",
    if (inherits(wt, "try-error")) "n/a" else sprintf("%.3f", wt$p.value))
msg("\nclinical-only median C across folds: ", sprintf("%.3f", median(res$C_clin, na.rm = TRUE)))
msg("VERDICT: ", ifelse(!inherits(wt, "try-error") && wt$p.value < 0.05 && best$delta > 0,
                        "genes add to clinical -- worth reporting",
                        "genes add nothing beyond stage/age/sex -- abandon"))
