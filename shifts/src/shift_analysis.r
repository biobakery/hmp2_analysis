
library(plyr)
library(dplyr)

source("./common/nonredundant_samples.r")
source("./common/disease_colors.r")
source("./common/theme_nature.r")
source("./common/disease_activity.r")

library(cowplot)

shift_thresholds <- function(D, di, dw, samesubject, nonibd) {
    shift_definition <- function(D_within_nonibd, D_interpersonal_nonibd) {
        if (length(D_within_nonibd) < 20) {
            return (NA)
        }

        dens_adj_nonibd <- density(D_within_nonibd, n=401, from=0, to=1)
        dens_diff_nonibd <- density(D_interpersonal_nonibd, n=401, from=0, to=1)

        return (dens_diff_nonibd$x[min(which(dens_diff_nonibd$y-0.0001 > dens_adj_nonibd$y))])
    }

    D_interpersonal_nonibd <- D[(!samesubject) & (di>0) & (nonibd)]
    thresholds <- c()
    threshold_x <- c()
    minthresh <- 0
    for (wd in seq(2, 40, by=2)) {
        D_within_nonibd <- D[samesubject & (di>0) & (dw>=wd-1) & (dw <= wd) & nonibd]
        thr <- shift_definition(D_within_nonibd, D_interpersonal_nonibd)
        minthresh <- max(minthresh, thr)
        thr <- max(minthresh, thr)
        thresholds <- c(thresholds, thr)
        threshold_x <- c(threshold_x, wd)
    }
    thresholds <- cummax(thresholds)

    return (list(thresholds=thresholds, threshold_x=threshold_x))
}

find_shifts <- function(pcl, nicename, name, onethresh=F, zerothresh=F, italic=F) {
    face <- ifelse(italic, "italic", "plain")

    # Remove replicates
    pcl <- nonredundant_samples(pcl)

    if (!("active" %in% colnames(pcl$meta))) {
        pcl <- merge_disease_activity(pcl, 4)
    }

    # Get distances
    library(vegan)
    D <- as.matrix(vegdist(pcl$x, method="bray"))
    dw <- outer(pcl$meta$week_num, pcl$meta$week_num, FUN="-")
    di <- outer(seq_along(pcl$meta$week_num), seq_along(pcl$meta$week_num), FUN="-")
    samesubject <- outer(pcl$meta$subject, pcl$meta$subject, FUN="==")
    nonibd <- outer(pcl$meta$diagnosis == "nonIBD", pcl$meta$diagnosis == "nonIBD", FUN="&")

    if (zerothresh) {
        threshold_x <- seq(2, 20, by=2)
        thresholds <- 0 * threshold_x
    } else if (!onethresh) {
        threshs <- shift_thresholds(D, di, dw, samesubject, nonibd)
        threshold_x <- threshs$threshold_x
        thresholds <- threshs$thresholds
    } else {
        threshs <- shift_thresholds(D, di, 1+0*dw, samesubject, nonibd)
        threshold_x <- threshs$threshold_x
        thresholds <- rep(threshs$thresholds[1], length(threshs$threshold_x))
    }

    # shift definition plots
    if (onethresh) {
        dens_t12_nonibd <- density(D[(di>=1) & samesubject & nonibd], n=401, from=0, to=1)
    } else {
        dens_t12_nonibd <- density(D[(dw>=1 & dw<=2) & samesubject & nonibd], n=401, from=0, to=1)
    }
    dens_diff_nonibd <- density(D[!samesubject & nonibd], n=401, from=0, to=1)

    pdf(sprintf("./shifts/%s_bc_thresholds.pdf", name), 6, 4.5)
    print(ggplot(data=data.frame(x=threshold_x, y=thresholds)) +
              geom_line(aes(group=1, x=x, y=y), size=1.2) +
              theme_cowplot() + ylab("Bray-Curtis Dissimilarity") + xlab("Time difference (weeks)"))
    dev.off()

    pdf(sprintf("./shifts/%s_shift_definition_n.pdf", name), 3, 2)
    ss_name <- sprintf("Same subject (n=%d)", dens_t12_nonibd$n)
    ds_name <- sprintf("Different subject (n=%d)", dens_diff_nonibd$n)
    colors <- c("Same subject"="slateblue", "Different subject"="firebrick")
    names(colors) <- c(ss_name, ds_name)
    ggp <- ggplot(data=data.frame(x=c(dens_t12_nonibd$x, dens_diff_nonibd$x),
                                  y=c(dens_t12_nonibd$y, dens_diff_nonibd$y),
                                  type=c(
                                      rep(ss_name, length(dens_t12_nonibd$y)),
                                      rep(ds_name, length(dens_diff_nonibd$y))))) +
        geom_path(aes(x=x, y=y, color=type), size=0.4) +
        geom_vline(xintercept=thresholds[threshold_x==2], color="black", size=0.8) +
        scale_color_manual(values=colors) +
        xlab("Bray-Curtis dissimilarity") + ylab("Density") +
        theme_nature() + theme(legend.position="bottom")
    print(ggp + guides(color="none"))
    print(ggp)
    dev.off()

    # Split distances into different subjects, and within-subject adjacent samples
    diffsubjD_nonibd <- D[!samesubject & lower.tri(D) & nonibd]
    samesubjpairs <- which(samesubject & (dw>=0) & !diag(T, nrow=nrow(D), ncol=ncol(D)), arr.ind=T)
    samesubjD <- data.frame(i1 = samesubjpairs[,2], i2 = samesubjpairs[,1])
    # Order by week sup
    samesubjD <- samesubjD[order(pcl$meta$collection[samesubjD$i1] + pcl$meta$collection[samesubjD$i2]),]
    # Remove non-adjacent pairs
    samesubjD <- samesubjD[!duplicated(samesubjD$i1),]
    # Merge in other info
    samesubjD$D <- D[cbind(samesubjD$i1, samesubjD$i2)]
    samesubjD$dw <- pcl$meta$week_num[samesubjD$i2] - pcl$meta$week_num[samesubjD$i1]
    samesubjD$subject <- pcl$meta$subject[samesubjD$i1]
    samesubjD$diagnosis <- pcl$meta$diagnosis[samesubjD$i1]
    samesubjD$coll1 <- pcl$meta$collection[samesubjD$i1]
    # Remove one side of same-week pairs
    samesubjD <- samesubjD[!((samesubjD$dw==0) & (samesubjD$i1>samesubjD$i2)),]
    # Sort by subject+collection
    samesubjD <- samesubjD[order(as.numeric(samesubjD$subject) + 0.01 * pcl$meta$collection[samesubjD$i1]),]

    # Dump how many potential shifts coult have been observed
    write.table(file=sprintf("./shifts/%s_possible_shifts.txt", name),
                quote=F, row.names=F, col.names=T, sep="\t",
                data.frame(nonIBD=sum(samesubjD$diagnosis=="nonIBD"),
                           UC=sum(samesubjD$diagnosis=="UC"),
                           CD=sum(samesubjD$diagnosis=="CD"),
                           total=nrow(samesubjD)))

    shiftdat <- list()
    shiftdat$thresholds <- thresholds
    shiftdat$threshold_x <- threshold_x

    # pdf(sprintf("./shifts/%s_bc_distributions.pdf", name), 6, 4.5)
    # print(ggplot(data=data.frame(x=c(shiftdat$shifted_x, shiftdat$shifted_x),
    #                              type=c(rep("Different Person", length(shiftdat$shifted_x)), rep("Same Person (adjacent sampled timepoints)", length(shiftdat$shifted_x))),
    #                              y=c(shiftdat$dens_diff_nonibd$y, shiftdat$dens_adj_nonibd$y))) +
    #           geom_line(aes(group=type, x=x, y=y, color=type), size=1.2) +
    #           theme_cowplot() + xlab("Bray-Curtis Dissimilarity") + ylab("KDE") +
    #           theme(legend.position="bottom") + scale_color_discrete(name=NULL) +
    #           ggtitle(sprintf("%s Bray-Curtis dissimilarities\namong non-IBD subjects (Intersect at %.3f)", nicename, shiftdat$shifted_adj_thresh)))
    # dev.off()

    # Identify and characterize shifts
    shifts <- samesubjD[samesubjD$D > thresholds[if(onethresh){1}else{max(1,ceiling(samesubjD$dw/2))}],]
    shifts$main_contributor <- NA
    shifts$delta_abundance <- NA
    shifts$return_bc <- NA
    shifts$weeks_left <- NA
    shifts$delta_dysbiosis <- pcl$meta$active[shifts$i2] - pcl$meta$active[shifts$i1]
    shifts$week_num <- pcl$meta$week_num[shifts$i2]
    shifts$collection <- pcl$meta$collection[shifts$i2]
    shift_profiles <- data.frame()
    shift_deltas <- pcl$x[c(),]
    shift_deltas_meta <- pcl$meta[c(),]
    for (i in seq_along(shifts$main_contributor)) {
        # Who is the main contributor?
        delta_ab <- pcl$x[shifts$i2[i],] - pcl$x[shifts$i1[i],]
        maxi <- which.max(abs(delta_ab))
        shifts$main_contributor[i] <- colnames(pcl$x)[maxi]
        shifts$delta_abundance[i] <- delta_ab[maxi]

        # Did they return to their pre-shift space?
        other_timepoints <- which(pcl$meta$subject==shifts$subject[i])
        other_timepoints <- other_timepoints[order(pcl$meta$collection[other_timepoints])]
        later_timepoints <- other_timepoints[pcl$meta$collection[other_timepoints] >= pcl$meta$collection[shifts$i2[i]]]
        shifts$return_bc[i] <- min(D[shifts$i1[i], later_timepoints])
        shifts$weeks_left[i] <- max(pcl$meta$week_num[later_timepoints]) - pcl$meta$week_num[shifts$i2[i]]

        # Asymmetry?
        main_ab <- pcl$x[other_timepoints, maxi]

        # Save the dissimilarity profile for plotting later
        shift_profiles <- rbind(shift_profiles, data.frame(
            D = D[shifts$i1[i], other_timepoints],
            shift_i = rep(i, length(other_timepoints)),
            diagnosis = rep(shifts$diagnosis[i], length(other_timepoints)),
            dw = pcl$meta$week_num[other_timepoints] - pcl$meta$week_num[shifts$i2[i]],
            main_ab = main_ab,
            subject = rep(shifts$subject[i], length(other_timepoints)),
            main_contributor = rep(shifts$main_contributor[i], length(other_timepoints))
        ))

        # Save the actual shift
        shift_deltas <- rbind(shift_deltas, delta_ab)
        shift_deltas_meta <- rbind(shift_deltas_meta, pcl$meta[shifts$i2[i],])
    }
    rownames(shift_deltas) <- rownames(shift_deltas_meta)
    shifts$class <- ifelse(shifts$delta_abundance>0, "bloom", "crash")

    # Dump the list of shifts
    write.table(shifts[,c("D", "subject", "diagnosis", "collection", "week_num", "main_contributor", "delta_abundance", "weeks_left", "delta_dysbiosis")],
                quote=F, row.names=F, col.names=T, sep="\t",
                file=sprintf("./shifts/%s_shift_list.txt", name))

    # Dump entries into/leaving dysbiosis
    shifts_nonIBD <- shifts[shifts$diagnosis=="nonIBD",]
    shifts_UC <- shifts[shifts$diagnosis=="UC",]
    shifts_CD <- shifts[shifts$diagnosis=="CD",]
    write.table(data.frame(row.names=c("nonIBD", "UC", "CD", "total"),
                           entry=c(sum(shifts_nonIBD$delta_dysbiosis==1), sum(shifts_UC$delta_dysbiosis==1), sum(shifts_CD$delta_dysbiosis==1), sum(shifts$delta_dysbiosis==1)),
                           exit=c(sum(shifts_nonIBD$delta_dysbiosis==-1), sum(shifts_UC$delta_dysbiosis==-1), sum(shifts_CD$delta_dysbiosis==-1), sum(shifts$delta_dysbiosis==-1))),
        quote=F, row.names=T, col.names=T, sep="\t",
        file=sprintf("./shifts/%s_shift_dysbiosis.txt", name)
    )

    # Count the total observation time by diagnosis so we can transform shift counts into shift rates
    observation_weeks <- rep(0, length(levels(pcl$meta$diagnosis)))
    names(observation_weeks) <- levels(pcl$meta$diagnosis)
    for (subject in levels(pcl$meta$subject)) {
        subji <- which(pcl$meta$subject == subject)
        if (length(subji) > 1) {
            dg <- pcl$meta$diagnosis[subji[1]]
            observation_weeks[dg] <- observation_weeks[dg] +
                max(pcl$meta$week_num[subji]) - min(pcl$meta$week_num[subji])
        }
    }
    observation_years <- observation_weeks / 52

    # Produce summary tabes
    dump_shift_rates <- function(tbl, file) {
        df <- as.data.frame(tbl)
        df$rate <- df$n / observation_years[df$diagnosis]
        write.table(df, sep="\t", quote=F, row.names=F, file=file)
    }
    dump_shift_rates(dplyr::count(shifts, diagnosis),
                     sprintf("./shifts/%s_shifts_diagnosis.txt", name))
    dump_shift_rates(dplyr::count(shifts, diagnosis, main_contributor),
                     sprintf("./shifts/%s_shifts_diagnosis_contributor.txt", name))

    # Make a cleaner version of shifts to remove shifts with main contributors
    # that only occur once
    shiftMCN <- dplyr::count(shifts, main_contributor)
    ord_MC <- shiftMCN$main_contributor[order(shiftMCN$n, decreasing=T)]
    shifts_clean <- shifts[shifts$main_contributor %in% ord_MC[1:10],]
    shifts_clean$main_contributor <- factor(shifts_clean$main_contributor, levels=rev(ord_MC[1:10]))

    # Barplot of shifts by main contributor
    library(tidyr)
    shiftN <- dplyr::count(shifts_clean, diagnosis, main_contributor) %>%
        complete(diagnosis, main_contributor, fill=list(n=0))
    shiftN$rate <- shiftN$n / observation_years[shiftN$diagnosis]
    #shiftN$main_contributor <- factor(shiftN$main_contributor, levels=shiftN$main_contributor[shiftN$diagnosis=="CD"][order(shiftN$rate[shiftN$diagnosis=="CD"])])

    pdf(sprintf("./shifts/%s_shift_contributors.pdf", name), 8, 7.5)
    print(ggplot(data=shiftN) +
        geom_bar(aes(x=main_contributor, y=rate, fill=diagnosis), stat="identity", position="dodge") +
        scale_fill_manual(values=hmp2_disease_colors, name=NULL) +
        theme_cowplot() + coord_flip() + xlab(NULL) + ylab("Shift rate (per observed year)") +
        theme(legend.position="bottom") +
        ggtitle(sprintf("%s shifts by primary contributor", nicename)))
    dev.off()

    pdf(sprintf("./shifts/%s_shift_contributors_figurequality.pdf", name), 1.823, 2.233)
    shiftdat$contributors_paperqual <- ggplot(data=shiftN) +
              geom_bar(aes(x=main_contributor, y=rate, fill=diagnosis), stat="identity", position="dodge") +
              scale_fill_manual(values=hmp2_disease_colors, name=NULL) +
              theme_cowplot() + coord_flip() + xlab(NULL) + ylab("Shift frequency (per year)") +
              ggtitle("Primary contributors to shifts") +
              theme_nature() + theme(axis.text.y=element_text(face=face)) +
              theme(legend.position="bottom")
    print(shiftdat$contributors_paperqual)
    dev.off()

    # Barplot of shifts by main contributor, split between bloom/crash
    shiftNbc <- dplyr::count(shifts_clean, diagnosis, main_contributor, class) %>%
        complete(diagnosis, main_contributor, class, fill=list(n=0))
    shiftNbc$rate <- shiftNbc$n / observation_years[shiftNbc$diagnosis]
    shiftNbc$main_contributor <- factor(shiftNbc$main_contributor, levels=levels(shiftN$main_contributor))

    pdf(sprintf("./shifts/%s_shift_contributors_bloomcrash.pdf", name), 8, 7.5)
    print(ggplot(data=shiftNbc) +
              geom_bar(aes(x=main_contributor, y=rate, fill=diagnosis), stat="identity", position="dodge") +
              scale_fill_manual(values=hmp2_disease_colors, name=NULL) +
              facet_grid(. ~ class) +
              theme_cowplot() + coord_flip() + xlab(NULL) + ylab("Shift rate (per observed year)") +
              theme(legend.position="bottom") +
              ggtitle(sprintf("%s shifts by primary contributor", nicename)))
    dev.off()

    pdf(sprintf("./shifts/%s_shift_returned_weekstoend.pdf", name), 7, 5.7)
    print(ggplot(data=shifts, aes(x=weeks_left, y=return_bc, fill=diagnosis)) +
        geom_point(aes(shape=class), size=3) +
        scale_shape_manual(values=c(bloom=24, crash=25)) +
        scale_fill_manual(values=hmp2_disease_colors, name=NULL) +
        scale_color_manual(values=hmp2_disease_colors) +
        guides(color="none") +
        theme_cowplot() + theme(legend.position="bottom") +
        xlab("Weeks until end of timeseries") +
        ylab("Minimum Bray-Curtis dissimilarity between the\ntimepoint before the shift and all future timepoints") +
        ggtitle(sprintf("%s shift recovery", nicename)))
    dev.off()

    pdf(sprintf("./shifts/%s_shift_returned.pdf", name), 7, 6)
    print(ggplot(data=shifts[shifts$weeks_left > 10,], aes(x=diagnosis, y=return_bc, fill=diagnosis)) +
        geom_boxplot() +
        scale_fill_manual(values=hmp2_disease_colors, name=NULL) +
        guides(color="none") +
        theme_cowplot() + theme(legend.position="none") +
        xlab(NULL) + ylab("Minimum Bray-Curtis dissimilarity between the\ntimepoint before the shift and all future timepoints\n(minimum 10 weeks until end of timeseries)") +
        ggtitle(sprintf("%s shift recovery after 10 weeks", nicename)))
    dev.off()

    pdf(sprintf("./shifts/%s_shift_returned_bloomcrash.pdf", name), 7, 6)
    print(ggplot(data=shifts[shifts$weeks_left > 10,], aes(x=diagnosis, y=return_bc, fill=diagnosis)) +
              geom_boxplot() +
              scale_fill_manual(values=hmp2_disease_colors, name=NULL) +
              guides(color="none") +
              facet_grid(. ~ class) +
              theme_cowplot() + theme(legend.position="none") +
              xlab(NULL) + ylab("Minimum Bray-Curtis dissimilarity between the\ntimepoint before the shift and all future timepoints\n(minimum 10 weeks until end of timeseries)") +
              ggtitle(sprintf("%s shift recovery after 10 weeks", nicename)))
    dev.off()

    pdf(sprintf("./shifts/%s_shift_profiles_dist_birdsnest.pdf", name), 10, 7.5)
    print(ggplot(data=shift_profiles) +
        geom_line(aes(x=dw, y=D, group=shift_i, color=diagnosis)) +
        scale_color_manual(values=hmp2_disease_colors, name=NULL) +
        theme_cowplot() + theme(legend.position="bottom") +
        ylab("Bray-Curtis distance from timepoint before shift") +
        xlab("Weeks since shift") +
        ggtitle(sprintf("%s shift profiles", nicename)))
    dev.off()

    shift_profiles_clean <- shift_profiles[shift_profiles$main_contributor %in% levels(shiftN$main_contributor),]
    shift_profiles_clean$main_contributor <- factor(shift_profiles_clean$main_contributor, levels=levels(shiftN$main_contributor))
    pdf(sprintf("./shifts/%s_shift_profiles_main_split.pdf", name), length(unique(shift_profiles_clean$main_contributor))*3.2, 7)
    print(ggplot(data=shift_profiles_clean, aes(x=dw, y=main_ab, color=subject)) +
        geom_point(size=1) +
        geom_line(aes(group=shift_i)) +
        theme_cowplot() +
        facet_grid(diagnosis ~ main_contributor) +
        ylab("Relative abundance") +
        xlab("Weeks since shift") +
        ggtitle(sprintf("%s main contributor shift profiles", nicename)))
    dev.off()

    # Shift ordinations
    deltas.pcl <- pcl.make(shift_deltas, meta=shift_deltas_meta)
    pdf(sprintf("./shifts/%s_shift_pca_symmetric.pdf", name), 6, 5)
    deltas.pca <- pcl.pcoa(deltas.pcl, D=dist(deltas.pcl$x))
    print(pcl.ordplot(deltas.pcl, deltas.pca,
                      colour="diagnosis", colour_override=hmp2_disease_colors))
    dev.off()
    adeltas.pcl <- deltas.pcl
    adeltas.pcl$x <- abs(adeltas.pcl$x)
    centroid <- NA
    if ("Escherichia coli" %in% colnames(adeltas.pcl$x)) {
        centroid <- c("Escherichia coli",
                                     "Prevotella copri",
                                     "Bacteroides uniformis",
                                     "Roseburia intestinalis",
                                     "Bacteroides stercoris",
                                     "Bacteroides fragilis",
                                     "Bacteroides ovatus",
                                     "Bacteroides vulgatus",
                                     "Dialister invisus",
                      "Faecalibacterium prausnitzii")
        centroid <- centroid[centroid %in% colnames(adeltas.pcl$x)]
    }
    printshiftord <- function(pcl, ord, pco) {
        print(pcl.ordplot(pcl, ord, pco=pco,
                          colour="diagnosis", colour_override=hmp2_disease_colors,
                          centroid=centroid))
    }
    pdf(sprintf("./shifts/%s_shift_pcoa_bc_abs.pdf", name), 6, 5)
    adeltas.pco <- pcl.pcoa(adeltas.pcl, k=4)
    adeltas.pcl.norm <- adeltas.pcl %>% pcl.normalize
    adeltas.pco.norm <- pcl.pcoa(adeltas.pcl.norm, k=4)
    printshiftord(adeltas.pcl, adeltas.pco, 2)
    printshiftord(adeltas.pcl, adeltas.pco.norm, 2)
    printshiftord(adeltas.pcl, adeltas.pco, 4)
    printshiftord(adeltas.pcl, adeltas.pco.norm, 4)
    dev.off()

    # PERMANOVA for differences by disease group
    library(vegan)
    Dshift <- vegdist(adeltas.pcl.norm$x, method="bray")
    sink(sprintf("./shifts/%s_shift_absdelta_adonis.txt", name))
    print(adonis(Dshift ~ ., data=adeltas.pcl$meta[,"diagnosis",drop=F], permutations=9999))
    sink()

    # Drivers of the ordination
    xy <- t(adeltas.pco.norm$points[,1:2]) %*% adeltas.pcl.norm$x
    drivers <- data.frame(feature=colnames(adeltas.pcl.norm$x), r=colSums(xy^2))
    drivers <- drivers[order(drivers$r, decreasing = T),,drop=F]
    write.table(drivers, sep="\t", quote=F, row.names=F,
                file=sprintf("./shifts/%s_shift_pcoa_bc_abs_drivers.txt", name))
    pdf(sprintf("./shifts/%s_shift_pcoa_bc_abs_drivers_paper.pdf", name), 0.634, 0.634)
    for (i in 1:20) {
        ggp <- pcl.ordplot(adeltas.pcl.norm, adeltas.pco.norm, colour=drivers$feature[i],
                           pointoutline=F, sortby=drivers$feature[i], size_abs=0.6) +
            ggtitle(drivers$feature[i]) +
            xlab(NULL) + ylab(NULL) + theme_nature() +
            theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
                  axis.ticks = element_blank(), plot.title = element_text(size=3))
        print(ggp)
        print(ggp + ggtitle(NULL) + guides(fill="none", color="none"))
    }
    dev.off()

    # Shift clustering
    ann_colors <- list(diagnosis=hmp2_disease_colors)
    pdf(sprintf("./shifts/%s_shift_heatmap.pdf", name), 10, 4.5)
    pcl.heatmap(deltas.pcl %>% pcl.top.f(max(abs(x)), 20), meta="diagnosis", minx=-Inf, annotation_colors = ann_colors)
    pcl.heatmap(adeltas.pcl %>% pcl.top.f(max(abs(x)), 20), meta="diagnosis", annotation_colors = ann_colors)
    pcl.heatmap(adeltas.pcl %>% pcl.normalize %>% pcl.top.f(max(abs(x)), 20), meta="diagnosis", annotation_colors = ann_colors)
    dev.off()

    pdf(sprintf("./shifts/%s_shift_heatmap_top10.pdf", name), 9, 3)
    pcl.heatmap(deltas.pcl %>% pcl.filter.f(Name %in% ord_MC[1:10]), meta="diagnosis", minx=-Inf, annotation_colors = ann_colors,
                clustering_method="average", color=colorRampPalette(c("steelblue4","black","goldenrod2"))(101))
    pcl.heatmap(adeltas.pcl %>% pcl.filter.f(Name %in% ord_MC[1:10]), meta="diagnosis", annotation_colors = ann_colors)
    pcl.heatmap(adeltas.pcl %>% pcl.normalize %>% pcl.filter.f(Name %in% ord_MC[1:10]), meta="diagnosis", annotation_colors = ann_colors)
    dev.off()

    shiftdat$shifts <- shifts
    return (shiftdat)
}



source("./common/load_bugs.r")
source("./common/load_kos.r")

species_shifts <- find_shifts(bugs.pcl %>% pcl.only(rank="s") %>% pcl.nicenames, "Species", "species", italic=T)
#genera_shifts <- find_shifts(bugs.pcl %>% pcl.only(rank="g") %>% pcl.nicenames, "Genus", "genera", italic=T)

speciesRNA_shifts <- find_shifts(bugs.fromko.rna.pcl, "Species (RNA)", "speciesRNA")

# source("./common/load_pathways.r")
#
# invisible(find_shifts(pwy.dna.unstrat.pcl, "Pathway-DNA", "pathwayDNA"))
# invisible(find_shifts(pwy.rna.unstrat.pcl, "Pathway-RNA", "pathwayRNA"))
#
# source("./common/load_ecs.r")
#
# invisible(find_shifts(ec.dna.unstrat.pcl, "EC-DNA", "ecDNA"))
# invisible(find_shifts(ec.rna.unstrat.pcl, "EC-RNA", "ecRNA"))
# invisible(find_shifts(bug.rna.pcl, "Bug-EC-RNA", "bug-ecRNA"))
#
source("./common/load_metabolites.r")
colnames(metabolites.named.pcl$x)[colnames(metabolites.named.pcl$x)=="sphingosine-isomer3"] <- "sphingosine isomer"

#metabolites_shifts <- find_shifts(metabolites.named.pcl, "Metabolite", "metabolites", onethresh=T)
metabolites.pcl.nrm.names <- metabolites.pcl.nrm
nm <- colnames(metabolites.pcl.nrm$x)
nm[metabolites.pcl.nrm.names$metaf$Metabolite!=""] <- metabolites.pcl.nrm.names$metaf$Metabolite[metabolites.pcl.nrm.names$metaf$Metabolite!=""]
colnames(metabolites.pcl.nrm.names$x) <- nm
metabolites_shifts <- find_shifts(metabolites.pcl.nrm.names, "Metabolite", "metabolites", onethresh=T)


source("./common/match_datasets.r")
b_mt <- match_datasets(list(bugs.pcl, metabolites.pcl.nrm.names), lenience=2)[[1]]
species_mbxsamp_shifts <- find_shifts(b_mt %>% pcl.only(rank="s") %>% pcl.nicenames, "Species (Metabolite Samples)", "species_mbxsamp", italic=T, onethresh=T)

#
# source("./common/load_proteins.r")
#
# invisible(find_shifts(proteins.pcl %>% pcl.normalize, "Protein", "proteins"))


# Shift contributors paper figure
library(egg)
pdf("./shifts/fig3_contributors.pdf", 3.220, 2*2.154, onefile=F)
ggarrange(ncol=1,
          species_shifts$contributors_paperqual +
              guides(fill="none") + ylab(NULL) +
              scale_y_continuous(limits=c(0, 0.6)),
          metabolites_shifts$contributors_paperqual +
              guides(fill="none") + ggtitle(NULL) +
              scale_y_continuous(limits=c(0, 0.6)))
dev.off()



