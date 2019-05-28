
source("./common/pcl_utils.r")
source("./env_config.r")
source("./common/disease_colors.r")
source("./common/match_datasets.r")
source("./common/merge_metadata.r")
source("./common/disease_activity.r")

library(dplyr)
library(ggplot2)

source("./common/load_bugs.r")
source("./common/load_biopsies.r")
source("./common/load_serology.r")

overview_battery <- function(pcl, alltaxa, ord, name) {
    if (!("diagnosis" %in% colnames(pcl$meta))) {
        pdf(sprintf("./overview/%s.pdf", name), 6, 5)
        print(pcl.ordplot(pcl, ord, size_abs=1.5))
        dev.off()
    } else {
        pdf(sprintf("./overview/%s_diagnosis_lines.pdf", name), 6, 5)
        print(pcl.ordplot(pcl, ord, colour="diagnosis", size_abs=1.5, connect="subject", sequence="week_num", colour_override=hmp2_disease_colors))
        dev.off()

        thm <- theme(plot.title = element_text(hjust = 0.5, face="bold"))
        sz <- 2
        make_plots <- function(pcos) list(
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Bacteroidetes",     size_abs=sz) + ggtitle("Bacteroidetes") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Firmicutes",        size_abs=sz) + ggtitle("Firmicutes") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Prevotella",        size_abs=sz) + ggtitle("Prevotella") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Enterobacteriaceae",size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Enterobacteriaceae") + thm,

            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Bacteroides",       size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Bacteroides") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Faecalibacterium",  size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Faecalibacterium") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Veillonellaceae",   size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Veillonellaceae") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Roseburia",         size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Roseburia") + thm,

            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Eubacterium",       size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Eubacterium") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Alistipes",         size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Alistipes") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Subdoligranulum",   size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Subdoligranulum") + thm,
            pcl.ordplot(pcos=pcos, alltaxa, ord, colour="Clostridium",       size_abs=sz, colour_log=T, colour_log_pseudocount=1e-3) + ggtitle("Clostridium") + thm
        )

        library(cowplot)
        pdf(sprintf("./overview/%s_taxonomy.pdf", name), 4*6.5, 3*6)
        print(plot_grid(plotlist=make_plots(2), nrow=3, ncol=4))
        dev.off()

        if (length(ord$ordnames) >= 4) {
            pdf(sprintf("./overview/%s_taxonomy_pc34.pdf", name), 4*6.5, 4*6)
            print(plot_grid(plotlist=make_plots(4), nrow=4, ncol=4))
            dev.off()
        }

        # Add some metadata
        pcl <- merge_metadata(pcl, "PDO.Number", "metagenomics")
        # Match and merge in serology
        pcls <- match_datasets(list(pcl, serology.eu.pcl), lenience=2)
        pcl$meta <- cbind(pcl$meta, pcls[[2]]$x[match(rownames(pcl$x), rownames(pcls[[1]]$x)),,drop=F])
        # Disease activity
        if (!("active" %in% colnames(pcl$meta))) {
            pcl <- merge_disease_activity(pcl)
        }

        make_plots <- function(pcos) list(
            pcl.ordplot(pcos=pcos, pcl, ord, colour="diagnosis", size_abs=sz, colour_override=hmp2_disease_colors) + ggtitle("Diagnosis") + thm,
            pcl.ordplot(pcos=pcos, pcl %>% pcl.filter.s(!is.na(hbi)), ord, colour="hbi", size_abs=sz) + ggtitle("HBI") + thm,
            pcl.ordplot(pcos=pcos, pcl %>% pcl.filter.s(!is.na(sccai)), ord, colour="sccai", size_abs=sz) + ggtitle("SCCAI") + thm,
            pcl.ordplot(pcos=pcos, pcl %>% pcl.filter.s(!is.na(fecalcal)), ord, colour="fecalcal", size_abs=sz) + ggtitle("FecalCal") + thm,

            pcl.ordplot(pcos=pcos, pcl, ord, colour="pilot",             size_abs=sz) + ggtitle("In Pilot") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="filtered_reads",    size_abs=sz, colour_log=T) + ggtitle("Filtered Read Depth") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="consent_age",       size_abs=sz) + ggtitle("Age at Consent") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="site_name",         size_abs=sz) + ggtitle("Recruitment Site") + thm,

            pcl.ordplot(pcos=pcos, pcl, ord, colour="IgA ASCA EU",       size_abs=sz) + ggtitle("Serology: IgA ASCA") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="IgG ASCA EU",       size_abs=sz) + ggtitle("Serology: IgG ASCA") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="OmpC. EU",          size_abs=sz) + ggtitle("Serology: OmpC") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="Cbir1 EU",          size_abs=sz) + ggtitle("Serology: Cbir1") + thm,

            pcl.ordplot(pcos=pcos, pcl, ord, colour="ANCA EU",           size_abs=sz) + ggtitle("Serology: ANCA") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="PDO.Number",        size_abs=sz) + ggtitle("PDO Number") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="active",            size_abs=sz) + ggtitle("Active Disease") + thm,
            pcl.ordplot(pcos=pcos, pcl, ord, colour="activity_index",    size_abs=sz) + ggtitle("Disease Activity Index") + thm
        )

        library(cowplot)
        pdf(sprintf("./overview/%s_metadata.pdf", name), 4*6.5, 4*6)
        print(plot_grid(plotlist=make_plots(2), nrow=4, ncol=4))
        dev.off()

        if (length(ord$ordnames) >= 4) {
            pdf(sprintf("./overview/%s_metadata_pc34.pdf", name), 4*6.5, 4*6)
            print(plot_grid(plotlist=make_plots(4), nrow=4, ncol=4))
            dev.off()
        }
    }
}

ordination_battery <- function(pcl, name) {
    # Heatmaps

    library(viridis)

    make_heatmap <- function(pcl, ...) {
        meta <- intersect(c("diagnosis"), colnames(pcl$meta))
        pcl.heatmap(pcl,
                    meta=meta, annotation_colors=list(diagnosis=hmp2_disease_colors),
                    show_rownames=T, color=viridis(100), logspace=T, minx=1e-4, zerospecial=NA, ...)
    }

    alltaxa <- pcl %>% pcl.nicenames

    do_plots <- function(pcl, subname) {
        pcl <- pcl.normalize(pcl)

        pdf(sprintf("./overview/%s_%s_heatmap_top30_ab.pdf", name, subname), 18, 7, onefile=F)
        make_heatmap(pcl %>% pcl.top.f(mean(x), n=30))
        dev.off()

        pdf(sprintf("./overview/%s_%s_heatmap_top30_abwhenpres.pdf", name, subname), 18, 7, onefile=F)
        make_heatmap(pcl %>% pcl.top.f((sum(x>0.001) > 50) * mean(x) / mean(x>0.001), n=30))
        dev.off()

        if (all("diagnosis", "site_sub_coll") %in% colnames(pcl$meta)) {
            pdf(sprintf("./overview/%s_%s_heatmap_top30_ab_nohclust.pdf", name, subname), 18, 9, onefile=F)
            make_heatmap(pcl %>% pcl.top.f(mean(x), n=30) %>%
                             pcl.sort.s(sprintf("%g%s", as.numeric(diagnosis), site_sub_coll)), cluster_cols=F)
            dev.off()
        }

        pdf(sprintf("./overview/%s_%s_heatmap_top30_pr.pdf", name, subname), 18, 9, onefile=F)
        make_heatmap(pcl %>% pcl.top.f(mean(x>0.001), n=30))
        dev.off()

        if (all("diagnosis", "site_sub_coll") %in% colnames(pcl$meta)) {
            pdf(sprintf("./overview/%s_%s_heatmap_top30_pr_nohclust.pdf", name, subname), 18, 9, onefile=F)
            make_heatmap(pcl %>% pcl.top.f(mean(x>0.001), n=30) %>%
                             pcl.sort.s(sprintf("%g%s", as.numeric(diagnosis), site_sub_coll)), cluster_cols=F)
            dev.off()
        }

        pdf(sprintf("./overview/%s_%s_heatmap_top30_abwhenpres.pdf", name, subname), 18, 9, onefile=F)
        make_heatmap(pcl %>% pcl.top.f((sum(x>0.001) > 50) * mean(x) / mean(x>0.001), n=30))
        dev.off()

        # Ordinations
        pcoa.ord <- pcl.pcoa(pcl, k=4)
        overview_battery(pcl, alltaxa, pcoa.ord, sprintf("%s_%s_pcoa_bc", name, subname))

        tsne.ord <- pcl.tsne(pcl)
        overview_battery(pcl, alltaxa, tsne.ord, sprintf("%s_%s_tsne_bc", name, subname))

        dmap.ord <- pcl.dmap(pcl, k=4)
        overview_battery(pcl, alltaxa, dmap.ord, sprintf("%s_%s_dmap_bc", name, subname))

    }

    # Ordinations
    if (pcl.only(pcl, rank="s")$nf > 0) {
        do_plots(pcl %>% pcl.only(rank="s") %>% pcl.nicenames, "species")
        do_plots(pcl %>% pcl.only(rank="g") %>% pcl.nicenames, "genus")
    } else {
        do_plots(pcl, "otu")
    }
}

ordination_battery(bugs.pcl %>% pcl.filter.f(!grepl("k__Viruses", Name)), "bugs_novirus")
ordination_battery(bugs.pcl, "bugs")
# for (site in levels(bugs.pcl$meta$site_name)) {
#     ordination_battery(pcl.filter.s(bugs.pcl, keep=bugs.pcl$meta$site_name==site), tolower(gsub("[^A-Za-z]", "_", site)))
# }
ordination_battery(biopsy_16s.pcl, "biopsy_16s")


# Color compression function to increase contrast.
compress_colors <- function(colors, compression = 0, bias=0) {
    x <- seq(0, 1, length=length(colors))
    y <- pbeta(x, exp(compression+bias), exp(compression-bias))

    newcols <- colorRamp(colors)(y) / 255
    return (rgb(newcols[,1], newcols[,2], newcols[,3]))
}


# Ordinations for main figure
library(viridis)
library(RColorBrewer)
source("./common/disease_colors.r")
species.pcl <- bugs.pcl %>% pcl.only(rank="s") %>% pcl.nicenames %>% pcl.normalize
species.ord <- pcl.pcoa(species.pcl)
species.pcl$meta$ginisimpson <- species.pcl %>% pcl.apply.s(1-sum(x^2))
ptsz <- 3
shrinkFactor <- 2.5
so.main <- pcl.ordplot(bugs.pcl, species.ord, colour="diagnosis", colour_override=hmp2_disease_colors, size_abs=ptsz) +
    guides(fill="none") + ggtitle("Taxonomy") +
    theme(text=element_text(size=15))
so.bact <- pcl.ordplot(bugs.pcl %>% pcl.nicenames, species.ord, colour="Bacteroidetes", size_abs=ptsz/shrinkFactor, pointoutline=F, sortby="Bacteroidetes") +
    guides(colour="none") + xlab(NULL) + ylab(NULL) + ggtitle("Bacteroidetes") +
    theme(text=element_text(size=15))
so.firm <- pcl.ordplot(bugs.pcl %>% pcl.nicenames, species.ord, colour="Firmicutes", size_abs=ptsz/shrinkFactor, pointoutline=F, sortby="Firmicutes") +
    guides(colour="none") + xlab(NULL) + ylab(NULL) + ggtitle("Firmicutes") +
    theme(text=element_text(size=15))
so.ent <-  pcl.ordplot(bugs.pcl %>% pcl.nicenames, species.ord, colour="Enterobacteriaceae", size_abs=ptsz/shrinkFactor, pointoutline=F, sortby="Enterobacteriaceae") +
    guides(colour="none") + xlab(NULL) + ylab(NULL) + ggtitle("Enterobacteriaceae") +
    theme(text=element_text(size=15))
gini_color_scale <- compress_colors(colorRampPalette(rev(brewer.pal(n=7, name="BuPu")))(101), 0.5, 0.2)
so.gini <-  pcl.ordplot(species.pcl, species.ord, colour="ginisimpson", size_abs=ptsz/shrinkFactor, pointoutline=F, sortby="ginisimpson", decreasing=T,
                        colour_override = gini_color_scale) +
    guides(colour="none") + xlab(NULL) + ylab(NULL) + ggtitle("Gini-Simpson") +
    theme(text=element_text(size=15))
library(cowplot)
w <- 0.23
sq <- 0.043
pdf("./overview/overview_taxonomy_ordinations.pdf", 7.5, 6)
print(ggdraw() +
          draw_plot(so.main,   0,   0, 1-w, 1) +
          draw_plot(so.gini,  1-w,  sq,   w, (1-sq)/3) +
          draw_plot(so.firm, 1-w, sq+(1-sq)/3,   w, (1-sq)/3) +
          draw_plot(so.bact, 1-w, sq+2*(1-sq)/3,   w, (1-sq)/3))
dev.off()
svg("./overview/overview_taxonomy_ordinations.svg", 7.5, 6)
print(ggdraw() +
    draw_plot(so.main,   0,   0, 1-w, 1) +
    draw_plot(so.gini,  1-w,  sq,   w, (1-sq)/3) +
    draw_plot(so.firm, 1-w, sq+(1-sq)/3,   w, (1-sq)/3) +
    draw_plot(so.bact, 1-w, sq+2*(1-sq)/3,   w, (1-sq)/3))
dev.off()
# for legends
pdf("./overview/overview_taxonomy_ordinations_legends.pdf", 7, 5)
pcl.ordplot(bugs.pcl, species.ord, colour="diagnosis", colour_override=hmp2_disease_colors)
pcl.ordplot(bugs.pcl %>% pcl.nicenames, species.ord, colour="Bacteroidetes")
pcl.ordplot(bugs.pcl %>% pcl.nicenames, species.ord, colour="Bacteroidetes", colour_override=gini_color_scale) +
    guides(fill=guide_colourbar())
dev.off()


# KS test of difference between PCo2's
pco2 <- species.ord$points[,2]
pco2_nibd <- pco2[species.pcl$meta$diagnosis == "nonIBD"]
pco2_uc <- pco2[species.pcl$meta$diagnosis == "UC"]
pco2_cd <- pco2[species.pcl$meta$diagnosis == "CD"]
ks.test(pco2_nibd, pco2_uc)
ks.test(pco2_nibd, pco2_cd)

# Other ordinations for Fig 1
source("./common/load_metabolites.r")
metabolites.pcoa.ord <- pcl.pcoa(metabolites.pcl.nrm)
source("./common/load_proteins.r")
proteins.pcoa.ord <- pcl.pcoa(pcl.normalize(proteins.pcl))
source("./common/load_kos.r")
ko.rna.unstrat.pcoa.ord <- pcl.pcoa(pcl.normalize(ko.rna.unstrat.pcl))
source("./common/load_biopsies.r")
biopsy_16s.pcoa.ord <- pcl.pcoa(pcl.normalize(biopsy_16s.pcl) %>% pcl.filter.s(!is.na(diagnosis)))

source("./common/theme_nature.r")
library(egg)
pdf("./overview/overview_other_ordinations.pdf", 2.417, 2.231)
ord_theme <- list(
    guides(fill="none"),
    theme_nature(),
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
)
ggarrange(
    pcl.ordplot(metabolites.pcl, metabolites.pcoa.ord, size_abs=1.3, outline_size=0.25,
                colour="diagnosis", colour_override=hmp2_disease_colors) +
        ord_theme + ggtitle("Metabolites"),
    pcl.ordplot(biopsy_16s.pcl, biopsy_16s.pcoa.ord, size_abs=1.3, outline_size=0.25,
                colour="diagnosis", colour_override=hmp2_disease_colors) +
        ord_theme + ggtitle("Taxonomy (Biopsies)"),
    pcl.ordplot(ko.rna.unstrat.pcl, ko.rna.unstrat.pcoa.ord, size_abs=1.3, outline_size=0.25,
                colour="diagnosis", colour_override=hmp2_disease_colors) +
        ord_theme + ggtitle("Transcripts"),
    pcl.ordplot(proteins.pcl, proteins.pcoa.ord, size_abs=1.3, outline_size=0.25,
                colour="diagnosis", colour_override=hmp2_disease_colors) +
        ord_theme + ggtitle("Proteins"),
    ncol=2)
dev.off()


# Sample dissimilarities over time
# These are an older version of the BC-over-time plots - newer versions are in the shifts folder
plot_dynamics <- function(pcl, name, index="bray/curtis",
                        binarize_thresh=NA, merge_biweeks=T) {

    library(labdsv)
    x <- if (!is.na(binarize_thresh)) {
        pcl$x >= binarize_thresh
    } else {
        pcl$x
    }
    if (index == "jaccard") {
        D <- as.matrix(dsvdis(x, index="bray/curtis"))
        D <- 1 - 2*D / (1 + D)
    } else {
        D <- as.matrix(dsvdis(x, index=index))
    }
    dw <- outer(pcl$meta$week_num, pcl$meta$week_num, FUN="-")
    dsub <- outer(as.numeric(pcl$meta$subject), as.numeric(pcl$meta$subject), FUN="-")
    allcd <- outer(pcl$meta$diagnosis=="CD", pcl$meta$diagnosis=="CD", FUN="&")
    alluc <- outer(pcl$meta$diagnosis=="UC", pcl$meta$diagnosis=="UC", FUN="&")
    allnonibd <- outer(pcl$meta$diagnosis=="nonIBD", pcl$meta$diagnosis=="nonIBD", FUN="&")

    dgn <- matrix("", nrow(allcd), ncol(allcd))
    dgn[allcd] <- "CD"
    dgn[alluc] <- "UC"
    dgn[allnonibd] <- "nonIBD"
    take <- (dsub==0) & (dw>0) & (allcd | alluc | allnonibd)
    dists <- D[take]
    dists_dgn <- dgn[take]
    dws <- dw[take]
    if (merge_biweeks)
        dws <- dws + (dws %% 2 == 1)
    ipdists <- D[(dsub!=0) & (dw>0) & (allcd | alluc | allnonibd)]
    ipdists_dgn <- dgn[(dsub!=0) & (dw>0) & (allcd | alluc | allnonibd)]

    index_name <- c("bray/curtis" = "Bray-Curtis Dissimilarity",
                    jaccard = "Jaccard Similarity")
    library(ggplot2)
    library(cowplot)
    df <- data.frame(posdws = factor(c(dws, rep(max(dws)+1, length(ipdists)))),
                     dists = c(dists, ipdists),
                     diagnosis = factor(c(dists_dgn, ipdists_dgn), levels=c("nonIBD", "UC", "CD")),
                     type = c(rep("Same person", length(dists)), rep("Between people", length(ipdists)))
    )
    ggp <- ggplot(data=df) + theme_cowplot() +
        geom_boxplot(aes(x=posdws, y=dists, fill=type), size=0.5, outlier.size=1) +
        facet_grid(. ~ diagnosis) +
        xlab(if (merge_biweeks) {"Time difference (week, odd weeks rounded up)"} else {"Time difference (week)"}) +
        ylab(if (index %in% names(index_name)) {index_name[index]} else {index}) +
        theme(axis.text.x=element_text(size=7))

    pdf(sprintf("./overview/community_dynamics_%s.pdf", name), 13, 5)
    print(ggp)
    dev.off()
}

plot_dynamics(pcl.only(bugs.pcl, rank="s"), "bc_species", index="bray/curtis")
plot_dynamics(pcl.only(bugs.pcl, rank="g"), "bc_genus", index="bray/curtis")
plot_dynamics(pcl.only(bugs.pcl, rank="s"), "bc_week_species", index="bray/curtis", merge_biweeks=F)
plot_dynamics(pcl.only(bugs.pcl, rank="g"), "bc_week_genus", index="bray/curtis", merge_biweeks=F)
plot_dynamics(pcl.only(bugs.pcl, rank="s"), "jaccard_.001_species", index="jaccard", binarize_thresh=0.001)
plot_dynamics(pcl.only(bugs.pcl, rank="g"), "jaccard_.001_genus", index="jaccard", binarize_thresh=0.001)
plot_dynamics(pcl.only(bugs.pcl, rank="s"), "jaccard_.0001_species", index="jaccard", binarize_thresh=0.0001)
plot_dynamics(pcl.only(bugs.pcl, rank="g"), "jaccard_.0001_genus", index="jaccard", binarize_thresh=0.0001)

# Viral dynamics
plot_dynamics(viruses.pcl, "jaccard_virus", index="jaccard", binarize_thresh=1)


# Alpha diversities

plot_alpha_divs <- function(pcl, name) {
    gs <- pcl.apply.s(pcl, 1-sum(x^2))
    is <- 1/(1-gs)
    sh <- pcl.apply.s(pcl, -sum(x[x>0]*log2(x[x>0])))

    df <- cbind(pcl$meta, data.frame(gs=gs, is=is, sh=sh))

    ggp_base <- ggplot(data=df, aes(x=diagnosis, fill=diagnosis)) +
        theme_cowplot() + xlab(NULL) + scale_fill_manual(values=hmp2_disease_colors)

    pdf(sprintf("./overview/alphadiversity_%s.pdf", name), 4.5, 3.5)
    print(ggp_base + geom_boxplot(aes(y=gs), size=0.5, outlier.size=1) +
              ylab("Gini-Simpson"))
    print(ggp_base + geom_boxplot(aes(y=is), size=0.5, outlier.size=1) +
              ylab("Inverse Simpson"))
    print(ggp_base + geom_boxplot(aes(y=sh), size=0.5, outlier.size=1) +
              ylab("Shannon Index (bits)"))
    dev.off()

}

plot_alpha_divs(pcl.only(bugs.pcl, rank="s"), "species")
plot_alpha_divs(pcl.only(bugs.pcl, rank="g"), "genus")



# CLR-based ordination
load(file.path(HMP2_data, "mgx", "species_CLR_overall.RData"))
species_CLR_overall.ord <- pcl.pcoa(D = dist(species_CLR_overall.pcl$x))
pcl.ordplot(species_CLR_overall.pcl, species_CLR_overall.ord, colour="diagnosis",
            colour_override=hmp2_disease_colors, size_abs=2)
overview_battery(bugs.pcl, pcl.nicenames(bugs.pcl), species_CLR_overall.ord, "species_CLR_overall")

load(file.path(HMP2_data, "mgx", "species_CLR_perfeature.RData"))
species_CLR_perfeature.ord <- pcl.pcoa(D = dist(species_CLR_perfeature.pcl$x))
pcl.ordplot(species_CLR_perfeature.pcl, species_CLR_perfeature.ord, colour="diagnosis",
            colour_override=hmp2_disease_colors, size_abs=2)
overview_battery(bugs.pcl, pcl.nicenames(bugs.pcl), species_CLR_perfeature.ord, "species_CLR_perfeature")

load(file.path(HMP2_data, "mgx", "species_CLR_persample.RData"))
species_CLR_persample.ord <- pcl.pcoa(D = dist(species_CLR_persample.pcl$x))
pcl.ordplot(species_CLR_persample.pcl, species_CLR_persample.ord, colour="diagnosis",
            colour_override=hmp2_disease_colors, size_abs=2)
overview_battery(bugs.pcl, pcl.nicenames(bugs.pcl), species_CLR_persample.ord, "species_CLR_persample")



# BC-based ordination using the bugs used in the CLR-based ordination
bugs.pcl.clrfilter <- bugs.pcl %>% pcl.only(rank="s") %>% pcl.nicenames
bugs.pcl.clrfilter$x <- bugs.pcl.clrfilter$x[,colnames(pcl.nicenames(species_CLR.pcl)$x)]
bugs.pcl.clrfilter <- pcl.normalize(bugs.pcl.clrfilter)

species_BC_CLR.ord <- pcl.pcoa(bugs.pcl.clrfilter)
pcl.ordplot(bugs.pcl.clrfilter, species_BC_CLR.ord, colour="diagnosis",
            colour_override=hmp2_disease_colors, size_abs=2)

overview_battery(bugs.pcl, pcl.nicenames(bugs.pcl), species_BC_CLR.ord, "species_CLRfiltered_BC")




