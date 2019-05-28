
source("./common/load_bugs.r")
source("./common/load_metabolites.r")
source("./common/load_proteins.r")
source("./common/load_kos.r")
source("./common/theme_nature.r")
source("./common/disease_activity.r")


shift_thresholds <- function(D, di, dw, samesubject, nonibd) {
    shift_definition <- function(D_within_nonibd, D_interpersonal_nonibd) {
        if (length(D_within_nonibd) < 20) {
            return (NA)
        }

        dens_adj_nonibd <- density(D_within_nonibd, n=401, from=0, to=1)
        dens_diff_nonibd <- density(D_interpersonal_nonibd, n=401, from=0, to=1)

        return (dens_diff_nonibd$x[min(which(dens_diff_nonibd$y > dens_adj_nonibd$y))])
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

ftest <- function(sse_c, df_c, sse_r, df_r, n) {
    df1 <- df_c - df_r
    df2 <- n - df_c - 1
    F_stat <- ((sse_r - sse_c) / df1) / (sse_c / df2)
    pval <- pf(F_stat, df1, df2, lower.tail=F)

    return (list(p.value=pval, statistic=F_stat))
}

powerlaw_fit <- function(x, y, xt) {
    # initial guesses
    C <- mean(y[x>=quantile(x,0.8)])
    a <- 4*(mean(y[x<=quantile(x,0.1)])-C)
    g <- -1.5

    # optimization
    library(nloptr)
    optparams <- neldermead(x0=c(C, a, g), fn=function(theta) {
        sum((y - (theta[1] + theta[2] * x^theta[3]))^2) },
        control=list(maxeval=100000))
    C <- optparams$par[1]
    a <- optparams$par[2]
    g <- optparams$par[3]

    yt <- a * xt^g + C

    # how did we do?
    # ggplot(data=data.frame(x=x, y=y), aes(x=x, y=y)) +
    #     geom_point() +
    #     geom_path(data=data.frame(x=xt, y=yt), size=2, color="blue") +
    #     scale_x_log10() + scale_y_log10()

    # F-test against flat
    sse_c <- sum((y - (C + a*x^g))^2)
    sse_r <- sum((y - mean(y))^2)

    ftest_flat <- ftest(sse_c, 3, sse_r, 1, length(x))

    return (list(xt=xt, yt=yt, C=C, a=a, g=g, x=x, y=y,
                 sse=sse_c, df_c=3, isflat=ftest_flat))
}

powerlaw_common_cmp <- function(fit1, fit2) {
    x <- c(fit1$x, fit2$x)
    y <- c(fit1$y, fit2$y)

    common_fit <- powerlaw_fit(x, y, fit1$xt)

    # F-test split model versus common model
    sse_r <- common_fit$sse
    df_r <- common_fit$df

    ftest_common <- ftest(fit1$sse + fit2$sse, fit1$df + fit2$df, sse_r, df_r, length(x))

    return (list(common_fit=common_fit, test=ftest_common, fit1=fit1, fit2=fit2))
}

# Sample dissimilarities over time
plot_dynamics <- function(pcl, index="bray/curtis",
                          binarize_thresh=NA, merge_to_weeks=2,
                          output_cmps=NA, plot_threshold=T) {

    pcl <- pcl.filter.s(pcl, !is.na(subject))
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
    dw[is.na(dw)] <- 1
    dcol <- outer(pcl$meta$collection, pcl$meta$collection, FUN="-")
    dcol[is.na(dcol)] <- 1
    dsub <- outer(as.numeric(pcl$meta$subject), as.numeric(pcl$meta$subject), FUN="-")
    dsub[is.na(dsub)] <- 1
    di <- outer(seq_along(pcl$meta$subject), seq_along(pcl$meta$subject), FUN="-")
    allcd <- outer(pcl$meta$diagnosis=="CD", pcl$meta$diagnosis=="CD", FUN="&")
    allcd[is.na(allcd)] <- F
    alluc <- outer(pcl$meta$diagnosis=="UC", pcl$meta$diagnosis=="UC", FUN="&")
    alluc[is.na(alluc)] <- F
    allnonibd <- outer(pcl$meta$diagnosis=="nonIBD", pcl$meta$diagnosis=="nonIBD", FUN="&")
    allnonibd[is.na(allnonibd)] <- F

    thr <- shift_thresholds(D, di, dw, dsub==0, allnonibd)

    dgn <- matrix("", nrow(allcd), ncol(allcd))
    dgn[allcd] <- "CD"
    dgn[alluc] <- "UC"
    dgn[allnonibd] <- "nonIBD"

    # tech rep distances
    take <- (di != 0) & (dsub==0) & (dcol==0) & (allcd | alluc | allnonibd)
    trdists <- D[take]
    trdists_dgn <- dgn[take]

    # intra-person over time distances
    take <- (dsub==0) & (dcol>0) & (dw <= 40) & (allcd | alluc | allnonibd)
    dists <- D[take]
    dists_dgn <- dgn[take]
    dws <- dw[take]
    dws[dws==0] <- 1

    # power-law fits
    plx <- seq(merge_to_weeks/2, 40 + merge_to_weeks/2, length=201)
    mask <- dists_dgn == "nonIBD"
    pl_nonIBD <- powerlaw_fit(dws[mask], dists[mask], plx)
    mask <- dists_dgn == "UC"
    pl_UC <- powerlaw_fit(dws[mask], dists[mask], plx)
    mask <- dists_dgn == "CD"
    pl_CD <- powerlaw_fit(dws[mask], dists[mask], plx)

    # difference of parameters for IBD vs non-IBD
    if (!is.na(output_cmps)) {
        plcmp_nonIBD_UC <- powerlaw_common_cmp(pl_nonIBD, pl_UC)
        plcmp_nonIBD_CD <- powerlaw_common_cmp(pl_nonIBD, pl_CD)

        sink(sprintf("./shifts/powerlawfits_%s.txt", output_cmps))
        printfit <- function(name, fit, test) {
            cat("========================================\n")
            cat(sprintf("    %s\n\n", name))
            cat(sprintf("y = %.3g + %.3g * x ^ %.3g\n", fit$C, fit$a, fit$g))
            cat(sprintf("SSE: %.3g, n: %d, F: %.3g\n", fit$sse, length(fit$x), test$statistic))
            cat(sprintf("P: %.3g\n", test$p.value))
            cat("\n")
            cat("\n")
        }
        printfit("non-IBD fit", pl_nonIBD, pl_nonIBD$isflat)
        printfit("UC fit", pl_UC, pl_UC$isflat)
        printfit("CD fit", pl_CD, pl_CD$isflat)
        printfit("non-IBD != UC", plcmp_nonIBD_UC$common_fit, plcmp_nonIBD_UC$test)
        printfit("non-IBD != CD", plcmp_nonIBD_CD$common_fit, plcmp_nonIBD_CD$test)
        sink()
    }

    # round times up
    dws <- ceiling(dws/merge_to_weeks) * merge_to_weeks

    # inter-personal distances
    ipdists <- D[(dsub!=0) & (dw>0) & (allcd | alluc | allnonibd)]
    ipdists_dgn <- dgn[(dsub!=0) & (dw>0) & (allcd | alluc | allnonibd)]

    # prepare dataframes
    index_name <- c("bray/curtis" = "Bray-Curtis Dissimilarity",
                    jaccard = "Jaccard Similarity")
    library(ggplot2)
    library(cowplot)
    df <- data.frame(posdws = c(rep(0, length(trdists)), dws, rep(max(dws)+merge_to_weeks, length(ipdists))),
                     dists = c(trdists, dists, ipdists),
                     diagnosis = factor(c(trdists_dgn, dists_dgn, ipdists_dgn), levels=c("nonIBD", "UC", "CD")),
                     type = c(rep("Tech rep", length(trdists)), rep("Same person", length(dists)), rep("Different people", length(ipdists)))
    )
    if (plot_threshold) {
        dfthr <- data.frame(x=c(thr$threshold_x, thr$threshold_x[length(thr$threshold_x)]+merge_to_weeks), y=c(thr$thresholds, thr$thresholds[length(thr$thresholds)]))
        dfthr$diagnosis <- factor("nonIBD", levels=c("nonIBD", "UC", "CD"))
        plot_threshold <- geom_step(data=dfthr, aes(x=x-merge_to_weeks/2, y=y), color="orange", size=0.4)
    } else {
        plot_threshold <- list()
    }
    dfpl <- data.frame(x=rep(plx, 3),
                       y=c(pl_nonIBD$yt, pl_UC$yt, pl_CD$yt),
                       diagnosis=factor(c(rep("nonIBD", length(plx)), rep("UC", length(plx)), rep("CD", length(plx))), levels=c("nonIBD", "UC", "CD")))

    sink(sprintf("./shifts/samplesizes_%s.txt", output_cmps))
    print(t(t(table(df$posdws, df$diagnosis))))
    sink()

    ggp <- ggplot(data=df) + theme_cowplot() +
        geom_boxplot(aes(group=factor(posdws), x=posdws, y=dists, fill=type, color=type), size=0.25, width=merge_to_weeks*1.75/2, outlier.size=0.6, outlier.stroke=0) +
        plot_threshold +
        geom_line(data=dfpl, aes(x=x, y=y), color="blue", size=0.4) +
        facet_grid(. ~ diagnosis) +
        xlab("Time difference (week)") +
        ylab(if (index %in% names(index_name)) {index_name[index]} else {index}) +
        scale_fill_manual(values=c("Tech rep" = "forestgreen", "Same person" = "steelblue2", "Different people" = "firebrick2")) +
        scale_color_manual(values=c("Tech rep" = rgb(colorRamp(c("forestgreen","black"))(0.5), maxColorValue=255),
                                    "Same person" = rgb(colorRamp(c("steelblue2","black"))(0.5), maxColorValue=255),
                                    "Different people" = rgb(colorRamp(c("firebrick2","black"))(0.5), maxColorValue=255))) +
        theme_nature() + guides(fill="none", color="none")

    return (ggp)
}

bug_dyn <- plot_dynamics(pcl.only(bugs.pcl, rank="s"), output_cmps="bugs", merge_to_weeks=2)
#bug_dyn
#mbx_named_dyn <- plot_dynamics(pcl.normalize(metabolites.named.pcl), merge_to_weeks=4)
mbx_dyn <- plot_dynamics(metabolites.pcl.nrm, merge_to_weeks=4, output_cmps="mbx")
mtx_dyn <- plot_dynamics(pcl.normalize(ko.rna.unstrat.pcl), merge_to_weeks=3, output_cmps="mtx")
mpx_dyn <- plot_dynamics(pcl.normalize(proteins.kos.pcl), merge_to_weeks=4, output_cmps="mpx")
bugmtx_dyn <- plot_dynamics(pcl.normalize(bugs.fromko.rna.pcl), merge_to_weeks=3, output_cmps="mtxbugs")


library(egg)

pdf("./shifts/community_dynamics_mainfig.pdf", 2.8, 4.2)
ggarrange(ncol=1, bug_dyn, mtx_dyn, mbx_dyn)
dev.off()

pdf("./shifts/community_dynamics_edfig.pdf", 3.973, 2.956)
ggarrange(ncol=1, mpx_dyn, bugmtx_dyn)
dev.off()

pdf("./shifts/community_dynamics_bugs_mainfig.pdf", 2.7, 1.3)
print(bug_dyn)
dev.off()

pdf("./shifts/community_dynamics_mbx_mainfig.pdf", 2.7, 1.3)
print(mbx_dyn)
dev.off()

pdf("./shifts/community_dynamics_mtx_mainfig.pdf", 2.7, 1.3)
print(mtx_dyn)
dev.off()

pdf("./shifts/fig3_community_dynamics.pdf", 3.8, 2.35, onefile=F)
library(egg)
ggarrange(ncol=1,
          bug_dyn + xlab(NULL) + ylab(NULL) + theme(strip.background = element_blank(),
                                       strip.text.x = element_blank()),
          mbx_dyn + xlab(NULL) + theme(strip.background = element_blank(),
                                 strip.text.x = element_blank()),
          mtx_dyn + ylab(NULL) + theme(strip.background = element_blank(),
                          strip.text.x = element_blank()))
dev.off()

# pdf("./shifts/community_dynamics_mpx_mainfig.pdf", 2.7, 0.95)
# print(mpx_dyn)
# dev.off()


# Dynamics for never-active subjects

subj_anyactive <- sapply(split(bugs.pcl$meta$active, factor(bugs.pcl$meta$subject)), any)
subj_N <- sapply(split(bugs.pcl$meta$active, factor(bugs.pcl$meta$subject)), length)
subj_anyactive_diag <- bugs.pcl$meta$diagnosis[match(names(subj_anyactive), bugs.pcl$meta$subject)]
table(subj_anyactive, subj_anyactive_diag)
#               subj_anyactive_diag
# subj_anyactive nonIBD UC CD
#          FALSE     16 23 34
#          TRUE      11 15 31
table(subj_anyactive[subj_N>=10], subj_anyactive_diag[subj_N>=10])
table(subj_N, subj_anyactive_diag)
#       subj_anyactive_diag
# subj_N nonIBD UC CD
#     1       0  7 12
#     2       0  0  1
#     3       0  1  2
#     5       1  0  0
#     7       0  2  0
#     8       0  1  2
#     9       0  1  1
#     10      3  3  5
#     11      2  5  5
#     12      2  4  9
#     13      3  2  8
#     14      6  1  4
#     15      1  0  2
#     16      0  1  2
#     18      0  1  1
#     20      0  1  3
#     21      1  1  1
#     22      2  2  1
#     23      3  2  5
#     24      2  2  0
#     25      1  1  0
#     26      0  0  1

bugs_noactive.pcl <- bugs.pcl %>%
    pcl.filter.s(keep=!subj_anyactive[as.character(bugs.pcl$meta$subject)])
bug_neveractive_dyn <- plot_dynamics(pcl.only(bugs_noactive.pcl, rank="s"),
    output_cmps="bugs_neveractive", merge_to_weeks=2, plot_threshold=F)

pdf("./shifts/community_dynamics_bugs_neveractive.pdf", 3.973, 2.956/2)
print(bug_neveractive_dyn)
dev.off()

