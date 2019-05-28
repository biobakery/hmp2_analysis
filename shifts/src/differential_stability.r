
source("./common/load_bugs.r")
source("./common/nonredundant_samples.r")
source("./common/disease_colors.r")

bugs.pcl <- bugs.pcl %>% pcl.filter.s(keep=!is.na(bugs.pcl$meta$subject))

bugs.pcl.species <- bugs.pcl %>% pcl.only(rank="s") %>%
    pcl.filter.f(mean(x>0.0001) > 0.2) %>%
    pcl.nicenames %>% nonredundant_samples

differences_nonibd <- matrix(0, nrow=0, ncol=bugs.pcl.species$nf)
differences_uc <- differences_nonibd
differences_cd <- differences_nonibd

for (subject in unique(bugs.pcl.species$meta$subject)) {
    ssp <- bugs.pcl.species %>%
        pcl.filter.s(keep=bugs.pcl.species$meta$subject == subject) %>%
        pcl.sort.s(collection)

    diffs <- diff(ssp$x)
    if (ssp$meta$diagnosis[1] == "nonIBD") {
        differences_nonibd <- rbind(differences_nonibd, diffs)
    } else if (ssp$meta$diagnosis[1] == "UC") {
        differences_uc <- rbind(differences_uc, diffs)
    } else {
        differences_cd <- rbind(differences_cd, diffs)
    }
}

for (i in seq_along(colnames(differences_nonibd))) {
    pdf(sprintf("./shifts/difference_ecdf-%s.pdf", colnames(differences_nonibd)[i]), 6, 5)
    ggp <- ggplot(data.frame(
        x=c(differences_nonibd[,i], differences_uc[,i], differences_cd[,i]),
        cls=c(rep("nonIBD", nrow(differences_nonibd)),
              rep("UC", nrow(differences_uc)),
              rep("CD", nrow(differences_cd)))), aes(x, colour=cls)) +
        stat_ecdf() + scale_color_manual(values=hmp2_disease_colors)
    print(ggp)
    dev.off()

    pdf(sprintf("./shifts/difference_density-%s.pdf", colnames(differences_nonibd)[i]), 6, 5)
    ggp <- ggplot(data.frame(
        x=c(differences_nonibd[,i], differences_uc[,i], differences_cd[,i]),
        cls=c(rep("nonIBD", nrow(differences_nonibd)),
              rep("UC", nrow(differences_uc)),
              rep("CD", nrow(differences_cd)))), aes(x, colour=cls)) +
        stat_density() + scale_color_manual(values=hmp2_disease_colors)
    print(ggp)
    dev.off()
}


df <- data.frame()

for (i in seq_along(colnames(differences_nonibd))) {
    dotest <- function(set1, set2) {
        dat1 <- set1[,i]
        dat2 <- set2[,i]
        #dat1 <- dat1[dat1>0]
        #dat2 <- dat2[dat2>0]
        dat1 <- abs(dat1)
        dat2 <- abs(dat2)
        return (ks.test(dat1, dat2, exact=T))
    }
    uct <- dotest(differences_nonibd, differences_uc)
    cdt <- dotest(differences_nonibd, differences_cd)
    mergedt <- dotest(differences_nonibd, rbind(differences_cd, differences_uc))

    df <- rbind(df, data.frame(
        feature = colnames(differences_nonibd)[i],
        stat_uc = uct$statistic,
        p.value_uc = uct$p.value,
        fdr_uc = NA,
        stat_cd = cdt$statistic,
        p.value_cd = cdt$p.value,
        fdr_cd = NA,
        stat_mg = mergedt$statistic,
        p.value_mg = mergedt$p.value,
        fdr_mg = NA
    ))
}

df$fdr_uc <- p.adjust(df$p.value_uc)
df$fdr_cd <- p.adjust(df$p.value_cd)
df$fdr_min <- pmin(df$fdr_uc, df$fdr_cd)

df$fdr_mg <- p.adjust(df$p.value_mg)

df <- df[order(df$fdr_min),]
df <- df[order(df$fdr_mg),]

write.table(df, "./shifts/differential_stability.tsv",
            quote=F, row.names=F, sep="\t")

