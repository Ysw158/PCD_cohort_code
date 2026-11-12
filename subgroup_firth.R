# Run using Rscript subgroup_firth.R private/subgroups.csv subgroup_pairs.csv out.csv
args <- commandArgs(trailingOnly = TRUE)
library(logistf)
patients <- read.csv(args[1], check.names = FALSE, na.strings = c("", "na", "none"))
pairs <- read.csv(args[2], check.names = FALSE, strip.white = FALSE)
pairs <- subset(pairs, (family == "S4" & phenotype != "Maternal miscarriage") | (family == "S6" & phenotype == "Infertility"))
rows <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  k <- pairs[i, ]
  dat <- data.frame(Outcome = patients[[k$phenotype]], Target_gene = as.integer(patients[[k$family]] == k$key),
                    Gender_new = factor(patients$sex), Enrollment.age_new = patients$age)
  dat <- droplevels(na.omit(dat))
  cells <- c(sum(dat$Target_gene == 1 & dat$Outcome == 1), sum(dat$Target_gene == 1 & dat$Outcome == 0),
             sum(dat$Target_gene == 0 & dat$Outcome == 1), sum(dat$Target_gene == 0 & dat$Outcome == 0))
  result <- c()
  for (adjusted in c(FALSE, TRUE)) {
    formula <- if (adjusted) Outcome ~ Target_gene + Gender_new + Enrollment.age_new else Outcome ~ Target_gene
    # Same estimator and contrasts as 5.Firth_regression.R. Retain the labelling of both arms.
    model <- logistf(formula, data = dat, plconf = 2, control = logistf.control(collapse = FALSE, maxit = 1000, maxhs = 0, maxstep = 5), plcontrol = logistpl.control(maxit = 300))
    if (any(model$pl.iter[2,1:2] >= 300)) {
      model <- logistf(formula, data=dat, plconf=2, control=logistf.control(collapse=FALSE, maxit=1000),
                       plcontrol=logistpl.control(maxit=1000,maxstep=0.25,maxhs=10))
      if (model$pl.iter[2,1] >= 1000) model$ci.lower[2] <- NA_real_
      if (model$pl.iter[2,2] >= 1000) model$ci.upper[2] <- NA_real_
    }
    est <- c(exp(coef(model)[2]), exp(model$ci.lower[2]), exp(model$ci.upper[2]), model$prob[2])
    if (adjusted && k$family == "S4" && k$phenotype == "Bronchiectasis" && k$key == "het") est[] <- NA_real_
    result <- c(result, unname(est))
  }
  rows[[i]] <- c(nrow(dat), cells, result)
}
res <- as.data.frame(do.call(rbind, rows))
names(res) <- c("N","yes_target","no_target","yes_other","no_other","or0","lo0","hi0","p0","or1","lo1","hi1","p1")
res <- cbind(pairs,res)
for (p in c("0","1")) {
  key <- paste0("p",p)
  res[[paste0("q",p)]] <- ave(res[[key]],res$family,res$phenotype,
    FUN=function(x){good <- !is.na(x); x[good] <- p.adjust(x[good], method="BH"); x})
}
write.csv(res,args[3],row.names=FALSE,na="")
