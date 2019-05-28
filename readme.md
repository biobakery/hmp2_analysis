# HMP2 Analysis
The HMP2 analysis repository contains code to generate figures for the IBDMDB paper.
The directory structure is as follows:

* `common/` scripts for common operations including:
	* Loading of all measurement types
	* Loading and merging metadata in PCL tables
	* Lenient matching of datasets
	* Common plot format settings (colors, themes, etc..)
	* Dysbiosis score (called "activity" in the code)
* `differential_abundance/` differential abundance testing
* `disease_activity/` dysbiosis timing and plots
* `mechanistic_correlations/` Cross-measurement-type correlations
* `mgx-mtx-mpx_correlations/` Gene family correlations across 3 measurement types
* `overview/` Overview plots and analyses
* `shifts/` Microbial datatype shift analysis

# Figure Guide
The source for most figures can be found in the following scripts. These scripts are generally written to output figures and results into their parent folder, e.g. `overview/src/taxonomy_overview.r` outputs to `overview/`.

* 1C: `overview/src/overview_venn.r`
* 1D: `overview/src/taxonomy_overview.r`
* 1E: `overview/src/mantel_tests.r`
* 1F: `overview/src/omnibus_tests.r`
* 2A: `disease_activity/src/disease_differences_boxplots.r`
* 2B: `overview/src/severity_correlations.r`
* 2C: `disease_activity/src/define_disease_activity.r`
* 2D: `disease_activity/src/define_disease_activity.r`
* 2E: `disease_activity/src/define_disease_activity.r`
* 2F: `disease_activity/src/disease_differences_boxplots.r` and `overview/src/serology_biopsy16s_overviews.r`
* 3A: `shifts/src/community_dynamics.r`
* 3B: `shifts/src/shift_analysis.r`
* 3C: `shifts/src/pcopri_timeseries.r`
* 3D: `shifts/src/pcopri_timeseries.r`
* 3E: `shifts/src/pcopri_timeseries.r`
* ED 1A: `overview/src/sample_overview.r`
* ED 1B: `overview/src/overview_venn.r`
* ED 2A: `overview/src/ed_tsne_figures.r`
* ED 2B: `overview/src/mgx-mtx-mpx_correlation_plot.r`
* ED 3A: `overview/src/severity_correlations.r`
* ED 3B: `disease_activity/src/define_disease_activity.r`
* ED 3C: `disease_activity/src/define_disease_activity.r`
* ED 3E: `disease_activity/src/hmp1ii_concordance.r`
* ED 4A: `disease_activity/src/alpha_diversity_boxplots.R`
* ED 4B: `disease_activity/src/disease_differences_boxplots.r`
* ED 4C: `disease_activity/src/disease_differences_boxplots.r`
* ED 4D: `disease_activity/src/disease_differences_boxplots.r`
* ED 5A: `shifts/src/community_dynamics.r`
* ED 5B: `shifts/src/shift_analysis.r`
* ED 5C: `shifts/src/shift_analysis.r`
* ED 5D: `shifts/src/pcopri_timeseries.r`
* ED 5E: `shifts/src/shift_analysis.r`
* ED 5F: `shifts/src/shift_analysis.r`
* S1: `mechanistic_correlations/src/mtx-rxn-mbx_associations.r`
* S2: `overview/src/subject_timeseries.r`

Fig. 4C and ED Figs. 7, 8, and 9 were generated using Cytoscape. Instructions can be found in `mechanistic_correlations/network construction.txt`.

## Setup
A script called `env_config.r` is expected in the root folder, containing paths to data folders. You can copy and fill in the template `env_config.r.template`. The path pointed to by `HMP2_Dropbox` is expected to have the following structure:

* `data/`
	* `16s/`
		* `biopsy_16s_otus_Rfriendly.pcl.tsv` biopsy 16S data
	* `htx/`
		* `host_tx_counts.pcl.tsv` host transcriptomes
	* `mbx/`
		* `iHMP_metabolomics.pcl.csv.gz` metabolomics profiles
	* `metadata/` full sample metadata tables
	* `mgx/`
		* `ecs_relab.pcl.tsv.gz` MGX EC abundances with stratification
		* `ecs_relab.slim.pcl.tsv.gz` MGX EC abundances without stratification
		* `kos_relab.pcl.tsv.gz` MGX KO abundances with stratification
		* `kos_relab.slim.pcl.tsv.gz` MGX KO abundances with stratification
		* `kos_relab_bugsummary.pcl.tsv` Species-level MGX KO rollups
		* `pathabundance_relab.pcl.tsv.gz` MGX pathway abundances with stratification
		* `taxonomic_profiles.pcl.tsv` MetaPhlAn2 taxonomic profiles
	* `mpx/`
		* `hmp2_MPX_not_normalized_ecs.names.pcl.tsv` EC rollup of protein counts
		* `hmp2_MPX_not_normalized_kos.names.pcl.tsv` KO rollup of protein counts
		* `iHMP_all_proteomics_?pep-pro_1p_FDR_?ppm.pcl.tsv.gz` Raw protein count tables
	* `mtx/`
		* `ecs_relab.pcl.tsv.gz` MTX EC abundances with stratification
		* `ecs_relab.slim.pcl.tsv.gz` MTX EC abundances without stratification
		* `kos_relab.pcl.tsv.gz` MTX KO abundances with stratification
		* `kos_relab.slim.pcl.tsv.gz` MTX KO abundances without stratification
		* `kos_relab_bugsummary.pcl.tsv` Species-level MTX KO rollups
		* `pathabundance_relab.pcl.tsv.gz` MTX pathway abundances with stratification
	* `mvx/`
		* `HMP2.Virome.VirMAP.pcl.tsv` VirMAP viral profiles
		* `metaphlan2_taxonomic_profiles.pcl.tsv` MetaPhlAn2 viral profiles
	* `serology/`
		* `hmp2_serology_Compiled_ELISA_Data.pcl.tsv` Serological profiles

All data tables are in PCL format with metadata attached at the top.
Additionally, `common/install_dependencies.r` can be run to install dependencies for most scripts in this repository.
