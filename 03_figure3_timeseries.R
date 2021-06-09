

# ============= "Reconciliation" MS Figure 3: Temporal trends in data =================

# dependencies and basedir
pacman::p_load("dplyr", "magrittr", "ggplot2")

# data
clover = read.csv("./CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv", stringsAsFactors = FALSE)

# set year
clover$Year = clover$PublicationYear
clover$Year[ is.na(clover$Year) ] = clover$ReleaseYear[ is.na(clover$Year) ]

# plot 1: overall temporal trends in associations
cl1 = clover %>%
  dplyr::group_by(Year, Database, Host, Virus) %>%
  dplyr::summarise(NReports = length(Year)) %>%
  dplyr::group_by(Year, Database) %>%
  dplyr::summarise(NReports=length(Year))
p1 = ggplot(cl1) + 
  geom_bar(aes(Year, NReports, fill=Database), stat="identity") + 
  theme_minimal() + 
  scale_fill_viridis_d(begin=0.05, end=0.95) + 
  xlab("Year") + 
  ylab("Number of associations") + 
  ggtitle("Total associations by source database") +
  scale_x_continuous(limits = c(1920, 2020), breaks=seq(1920, 2020, by=20)) + 
  theme(plot.title =  element_text(size=13.5, hjust=0.5, vjust=-4),
        legend.position=c(0.08, 0.5), legend.title=element_blank(),
        axis.text = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11))

# plot 2: temporal trends in novel associations
cl2 = clover %>%
  dplyr::mutate(HostVirus = paste(Host, Virus, sep="_")) %>%
  dplyr::group_by(Database, HostVirus) %>%
  dplyr::summarise(Year = min(Year, na.rm=TRUE)) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(NReports = n_distinct(HostVirus))
cl2$Year[ cl2$Year == Inf ] = NA
p2 = ggplot(cl2) + 
  geom_bar(aes(Year, NReports), stat="identity", fill="grey50", col=NA) + 
  theme_minimal() + 
  xlab("Year") + 
  ylab("Number of associations") + 
  ggtitle("Novel host-virus associations (previously unreported)") +
  scale_x_continuous(limits = c(1920, 2020), breaks=seq(1920, 2020, by=20)) + 
  theme(plot.title =  element_text(size=13.5, hjust=0.5, vjust=-4),
        axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12))

# combine
p_comb = gridExtra::grid.arrange(p1, p2, ncol=1)
#ggsave(p_comb, file="./Figure3_TemporalTrends.png", device="png", units="in", scale=0.95, width=9, height=6, dpi=300)



# effort by year
eff_plot = read.csv("./pubmed_viruscounts/PubMed_HostsEffort_PerYear_VirusRelated_19302020.csv") %>%
  dplyr::filter(NumPubs > 0 & Year <= 2020) %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(NumPubs = sum(NumPubs),
                   `Host species` = n_distinct(Host)) %>%
  ggplot() + 
  geom_smooth(aes(Year, NumPubs), method="gam", se=FALSE, alpha=0.9) +
  geom_point(aes(Year, NumPubs, size=`Host species`), pch=21, col="grey20", fill=viridis::viridis(200)[100], alpha=0.5) +
  theme_minimal() + xlab("Year") + ylab("Total virus-related publications") +
  ggtitle("Virus-related publications") +
  scale_x_continuous(limits=c(1930, 2020), breaks=seq(1940, 2020, by=20)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=11), legend.title=element_text(size=11),
        legend.position = c(0.18, 0.6),
        plot.title =  element_text(size=13.5, hjust=0.5, vjust=-4))



p_comb = gridExtra::grid.arrange(grobs = list(p1, p2, eff_plot), 
                                 ncol=2, nrow=4, layout_matrix=rbind(c(1, NA), c(1, 3),
                                                                     c(2, 3), c(2, NA)), 
                                 widths=c(0.8, 0.4), heights=c(0.25, 0.5, 0.5, 0.25))
p_comb = ggpubr::as_ggplot(p_comb)  +
  cowplot::draw_plot_label(label = c("a", "b", "c"),
                           fontface = "bold", size = 24,
                           x = c(0.05, 0.05, 0.68), y = c(0.95, 0.48, 0.82))
ggsave(p_comb, file="./Figure3_TemporalTrends.jpeg", units="in", scale=0.9, width=12, height=6.8, dpi=600)

