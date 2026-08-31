# Rebuild the ARID1A-associated prognostic signature, properly.
#
# PRE-SPECIFIED before any result is seen:
#   Candidate pool : genes measurable in TCGA and all six validation cohorts,
#                    restricted to those associated with ARID1A expression in TCGA
#                    (the paper's own definition of the ARID1A-associated space).
#   Model          : LASSO-Cox, lambda by 10-fold CV.
#   Honest accuracy: nested CV -- lambda chosen inside each training fold only, so
#                    the reported C-index is not inflated by having seen the test
#                    fold. The archived analysis reports a bootstrap-corrected
#                    R2 of 0.041, which suggests little to find; nested CV settles it.
#   Stability      : 1,000 bootstrap resamples of TCGA, recording how often each
#                    gene is selected. A signature whose members change between
#                    resamples is not a signature.
#   Validation     : model fitted once on all of TCGA, applied unchanged to each
#                    cohort, random-effects meta-analysis of the per-cohort HRs.
#
#   SUCCESS is declared only if ALL THREE hold:
#     (a) nested-CV C-index >= 0.60 on TCGA
#     (b) pooled external HR with 95% CI excluding 1
#     (c) between-cohort I2 < 50%
#   Anything less is reported as a negative result.
suppressMessages({library(glmnet); library(survival); library(meta)})
set.seed(20260808)
OUT <- "/data/s01011/arid1a_prognostic_redo"
H <- readRDS(file.path(OUT, "harmonised.rds"))
msg <- function(...) cat(..., "\n", sep = "")

tc <- H$tcga; common <- H$common
X <- t(tc$x[common, , drop = FALSE])
y <- Surv(tc$time, tc$event)
msg("TCGA ", nrow(X), " tumours x ", ncol(X), " candidate genes, ",
    sum(tc$event), " deaths")

# ---- candidate pool: association with ARID1A expression, pre-specified cut
a <- tc$x["ARID1A", ]
r <- apply(X, 2, function(v) suppressWarnings(cor(v, a, method = "spearman")))
pool <- colnames(X)[!is.na(r) & abs(r) >= 0.3]
msg("genes with |rho| >= 0.3 against ARID1A expression: ", length(pool))
Xp <- X[, pool, drop = FALSE]

# ---- nested CV for an honest C-index
folds <- sample(rep(1:5, length.out = nrow(Xp)))
pred <- rep(NA_real_, nrow(Xp))
for (k in 1:5) {
  tr <- folds != k; te <- !tr
  cv <- cv.glmnet(Xp[tr, , drop = FALSE], y[tr], family = "cox", alpha = 1, nfolds = 10)
  fit <- glmnet(Xp[tr, , drop = FALSE], y[tr], family = "cox", alpha = 1,
                lambda = cv$lambda.min)
  pred[te] <- as.numeric(predict(fit, Xp[te, , drop = FALSE]))
}
cidx <- summary(coxph(y ~ pred))$concordance
msg(sprintf("\nnested-CV C-index on TCGA: %.3f (SE %.3f)", cidx[1], cidx[2]))

# ---- stability selection
B <- 1000
sel <- integer(ncol(Xp)); names(sel) <- colnames(Xp)
nsel <- integer(B)
for (b in 1:B) {
  i <- sample(nrow(Xp), replace = TRUE)
  if (sum(tc$event[i]) < 10) next
  cv <- try(cv.glmnet(Xp[i, , drop = FALSE], y[i], family = "cox", alpha = 1, nfolds = 5),
            silent = TRUE)
  if (inherits(cv, "try-error")) next
  co <- as.matrix(coef(cv, s = "lambda.min"))
  g <- rownames(co)[co[, 1] != 0]
  sel[g] <- sel[g] + 1L
  nsel[b] <- length(g)
}
freq <- sort(sel / B, decreasing = TRUE)
# persist the expensive result BEFORE any formatting: the first run lost 1,000
# completed bootstrap resamples to a %d/%f mismatch on the line below
saveRDS(list(freq = freq, nsel = nsel, cindex = cidx),
        file.path(OUT, "bootstrap_stability.rds"))
msg(sprintf("bootstrap: median %.0f genes selected per resample", median(nsel[nsel > 0])))
msg("selection frequency, top 15:")
for (i in seq_len(min(15, length(freq))))
  msg(sprintf("   %-12s %.3f", names(freq)[i], freq[i]))
msg(sprintf("genes selected in >50%% of resamples: %d", sum(freq > 0.5)))

# ---- final model on all of TCGA
cvf <- cv.glmnet(Xp, y, family = "cox", alpha = 1, nfolds = 10)
cof <- as.matrix(coef(cvf, s = "lambda.min"))
sig <- cof[cof[, 1] != 0, , drop = FALSE]
msg(sprintf("\nfinal signature: %d genes", nrow(sig)))
print(round(sig, 4))
saveRDS(list(sig = sig, freq = freq, cindex = cidx, pool = pool),
        file.path(OUT, "model.rds"))

# ---- external validation, model applied unchanged
res <- data.frame()
for (nm in names(H$ext)) {
  z <- H$ext[[nm]]
  gsub_ <- intersect(rownames(sig), rownames(z$x))
  if (length(gsub_) < nrow(sig)) msg("  NOTE ", nm, ": ", length(gsub_), "/", nrow(sig), " genes present")
  if (!length(gsub_)) next
  # score on within-cohort z-scores so platform scale cannot drive the result
  em <- t(scale(t(z$x[gsub_, , drop = FALSE])))
  score <- as.numeric(t(em) %*% sig[gsub_, 1])
  cx <- coxph(Surv(z$time, z$event) ~ scale(score))
  s <- summary(cx)
  res <- rbind(res, data.frame(cohort = nm, hist = z$hist, n = ncol(z$x),
                               deaths = sum(z$event),
                               HR = s$conf.int[1], lo = s$conf.int[3],
                               hi = s$conf.int[4], P = s$coefficients[5],
                               C = s$concordance[1]))
}
msg("\nexternal validation (HR per SD of risk score):")
print(format(res, digits = 3), row.names = FALSE)

mm <- metagen(TE = log(res$HR), seTE = (log(res$hi) - log(res$lo)) / (2 * 1.96),
              studlab = res$cohort, sm = "HR", random = TRUE, common = FALSE)
msg(sprintf("\npooled HR %.3f (%.3f-%.3f), P = %.3f, I2 = %.0f%%",
            exp(mm$TE.random), exp(mm$lower.random), exp(mm$upper.random),
            mm$pval.random, mm$I2 * 100))

pass_a <- cidx[1] >= 0.60
pass_b <- exp(mm$lower.random) > 1 || exp(mm$upper.random) < 1
pass_c <- mm$I2 < 0.5
msg("\n", strrep("=", 62))
msg(sprintf("(a) nested-CV C-index >= 0.60 : %s (%.3f)", ifelse(pass_a, "PASS", "FAIL"), cidx[1]))
msg(sprintf("(b) pooled CI excludes 1      : %s", ifelse(pass_b, "PASS", "FAIL")))
msg(sprintf("(c) I2 < 50%%                  : %s (%.0f%%)", ifelse(pass_c, "PASS", "FAIL"), mm$I2 * 100))
msg("VERDICT: ", ifelse(pass_a && pass_b && pass_c,
                        "signature validated", "NOT validated -- report as negative"))
msg(strrep("=", 62))
write.csv(res, file.path(OUT, "external_validation_redo.csv"), row.names = FALSE)
