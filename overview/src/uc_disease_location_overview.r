
source("./common/uc_disease_locations.r")

uc_participants <- unique(hmp2_sample_metadata$Participant.ID[hmp2_sample_metadata$diagnosis == "UC"])

locs_matched <- uc_disease_locations_master[match(uc_participants, names(uc_disease_locations_master))]

library(ggplot2)
pdf("./overview/uc_location_groups_barplot.pdf", 5, 3)
print(ggplot(data=data.frame(c=locs_matched)) +
    geom_bar(aes(c)) +
    xlab(NULL) +
    ggtitle("UC locations (All UC patients)"))


source("./common/load_bugs.r")
locs_matched <- uc_disease_locations_master[match(intersect(uc_participants, unique(bugs.pcl$meta$subject)), names(uc_disease_locations_master))]

library(ggplot2)
print(ggplot(data=data.frame(c=locs_matched)) +
          geom_bar(aes(c)) +
          xlab(NULL) +
          ggtitle("UC locations (UC patients with MGX)"))
dev.off()

