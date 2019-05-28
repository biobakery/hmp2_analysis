
source("./common/load_bugs.r")
source("./common/disease_colors.r")
source("./common/theme_nature.r")

test_stabi_differences <- function(pcl, what, min_thresh=0.1) {
    subj_max <- split(pcl$x[,what], pcl$meta$subject) %>% sapply(max)

    copri_people <- names(which(subj_max > min_thresh))

    pcl.ppl <- pcl.filter.s(pcl, subject %in% copri_people) %>%
        pcl.sort.s(collection)
    cat(sprintf("N=%d samples from %d subjects\n", pcl.ppl$ns, length(copri_people)))

    ts <- pcl.ppl$x[,what]

    mask <- pcl.ppl$meta$diagnosis == "nonIBD"
    delta_nonIBD <- split(ts[mask], pcl.ppl$meta$subject[mask]) %>% lapply(diff) %>% unlist
    mask <- pcl.ppl$meta$diagnosis == "UC"
    delta_UC <- split(ts[mask], pcl.ppl$meta$subject[mask]) %>% lapply(diff) %>% unlist
    mask <- pcl.ppl$meta$diagnosis == "CD"
    delta_CD <- split(ts[mask], pcl.ppl$meta$subject[mask]) %>% lapply(diff) %>% unlist

    cat("non-IBD -- UC")
    print(wilcox.test(abs(delta_nonIBD), abs(delta_UC)))
    cat("non-IBD -- CD")
    print(wilcox.test(abs(delta_nonIBD), abs(delta_CD)))
    cat("UC -- CD")
    print(wilcox.test(abs(delta_UC), abs(delta_CD)))
}

stabi_plots <- function(pcl, what, min_thresh=0.1, units="relative abundance", plotscale=1) {
    subj_max <- split(pcl$x[,what], pcl$meta$subject) %>% sapply(max)

    copri_people <- names(which(subj_max > min_thresh))

    pcl.ppl <- pcl.filter.s(pcl, subject %in% copri_people)

    library(ggplot2)
    ggp <-ggplot(data=data.frame(pcl.ppl$meta, y = pcl.ppl$x[,what]/plotscale),
           aes(x=week_num, y=y, color=diagnosis)) +
        geom_line(aes(group=subject), size=0.25) +
        geom_point(size=1, stroke=0) +
        facet_grid(diagnosis ~ .) +
        guides(color="none") +
        scale_color_manual(values=hmp2_disease_colors) +
        theme_nature() +
        ylab(sprintf("%s (%s)", what, units)) +
        xlab("Week")

    return (ggp)
}

# P copri timeseries plot for paper
species.pcl <- bugs.pcl %>% pcl.only(rank="s") %>% pcl.nicenames

pdf("./shifts/stability_pcopri.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Prevotella copri", 0.1))
dev.off()

test_stabi_differences(species.pcl, "Prevotella copri", 0.1)
# non-IBD -- UC
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_nonIBD) and abs(delta_UC)
# W = 5876, p-value = 4.202e-06
# alternative hypothesis: true location shift is not equal to 0
#
# non-IBD -- CD
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_nonIBD) and abs(delta_CD)
# W = 4589, p-value = 0.000108
# alternative hypothesis: true location shift is not equal to 0
#
# UC -- CD
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_UC) and abs(delta_CD)
# W = 2712, p-value = 0.2367
# alternative hypothesis: true location shift is not equal to 0


# other stuff

pdf("./shifts/stability_rhominis.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Roseburia hominis", 0.1))
dev.off()

pdf("./shifts/stability_fprau.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Faecalibacterium prausnitzii", 0.1))
dev.off()

pdf("./shifts/stability_ecoli.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Escherichia coli", 0.0))
dev.off()
test_stabi_differences(species.pcl, "Escherichia coli", -1)
# non-IBD -- UC
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_nonIBD) and abs(delta_UC)
# W = 70300, p-value = 0.0002224
# alternative hypothesis: true location shift is not equal to 0
#
# non-IBD -- CD
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_nonIBD) and abs(delta_CD)
# W = 124400, p-value = 0.02872
# alternative hypothesis: true location shift is not equal to 0
#
# UC -- CD
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_UC) and abs(delta_CD)
# W = 148490, p-value = 0.096
# alternative hypothesis: true location shift is not equal to 0

print(stabi_plots(species.pcl, "Escherichia coli", 0.1))
test_stabi_differences(species.pcl, "Escherichia coli", 0.1)

pdf("./shifts/stability_ecoli_larger.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Escherichia coli", 0.0))
dev.off()


pdf("./shifts/stability_buniformis.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Bacteroides uniformis", 0.1))
dev.off()

pdf("./shifts/stability_bfragilis.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Bacteroides fragilis", 0.1))
dev.off()

pdf("./shifts/stability_bvulgatus.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Bacteroides vulgatus", 0.1))
dev.off()

pdf("./shifts/stability_erectale.pdf", 2, 3.7/2)
print(stabi_plots(species.pcl, "Eubacterium rectale", 0.1))
dev.off()

source("./common/load_metabolites.r")
pdf("./shifts/stability_methylimidazole_acetic_acid.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.named.pcl, "methylimidazole acetic acid", 10))
dev.off()

pdf("./shifts/stability_phenylalanine.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.named.pcl, "phenylalanine", 10))
dev.off()

pdf("./shifts/stability_N-acetylhistamine.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.named.pcl, "N-acetylhistamine", 10))
dev.off()


pdf("./shifts/fig3_examples.pdf", 1.970*2, 1.774, onefile=F)
library(egg)
ggarrange(ncol=2,
          stabi_plots(species.pcl, "Prevotella copri", 0.1),
          stabi_plots(species.pcl, "Escherichia coli", -1))
dev.off()


# First set of non-targetted stuff shifting (caused by crazy outliers)

pdf("./shifts/stability_HILp_QI17524.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.pcl.nrm, "HILp_QI17524", 1e-7))
dev.off()

pdf("./shifts/stability_HILp_QI22820.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.pcl.nrm, "HILp_QI22820", 1e-7))
dev.off()

pdf("./shifts/stability_HILn_QI10991.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.pcl.nrm, "HILn_QI10991", 1e-7))
dev.off()

pdf("./shifts/stability_HILn_QI8869.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.pcl.nrm, "HILn_QI8869", 1e-7))
dev.off()


# Second set of shifters

pdf("./shifts/stability_methylimidazole_acetic_acid.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.named.pcl, "methylimidazole acetic acid", 10, units="abundance"))
dev.off()

pdf("./shifts/stability_HILp_QI20578.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.pcl.nrm, "HILp_QI20578", 1e-7))
print(stabi_plots(metabolites.pcl, "HILp_QI20578", 10, units="abundance"))
dev.off()

pdf("./shifts/stability_HILp_QI22918.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.pcl.nrm, "HILp_QI22918", 1e-7))
print(stabi_plots(metabolites.pcl, "HILp_QI22918", 10, units="abundance"))
dev.off()

pdf("./shifts/stability_urate.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.named.pcl, "urate", 10, units="abundance"))
dev.off()
test_stabi_differences(metabolites.named.pcl, "urate", 10)
# non-IBD -- UC
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_nonIBD) and abs(delta_UC)
# W = 7901, p-value = 0.001218
# alternative hypothesis: true location shift is not equal to 0
#
# non-IBD -- CD
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_nonIBD) and abs(delta_CD)
# W = 13325, p-value = 0.04368
# alternative hypothesis: true location shift is not equal to 0
#
# UC -- CD
# Wilcoxon rank sum test with continuity correction
#
# data:  abs(delta_UC) and abs(delta_CD)
# W = 11075, p-value = 0.0932
# alternative hypothesis: true location shift is not equal to 0


pdf("./shifts/stability_urobilin.pdf", 2, 3.7/2)
print(stabi_plots(metabolites.named.pcl, "urobilin", 10, units="abundance"))
dev.off()


foo <- c()
for (i in seq_along(foo)) {
    pdf("./shifts/stability_HILp_QI22918.pdf", 2, 3.7/2)
    print(stabi_plots(metabolites.pcl.nrm, "HILp_QI22918", 1e-7))
    print(stabi_plots(metabolites.pcl, "HILp_QI22918", 10, units="abundance"))
    dev.off()
}




