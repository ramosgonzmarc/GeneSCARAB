########################################################################
#                              GENESCARAB                              #
########################################################################
library(GeneSCARAB)

# Load example data
data("circa_table_genescarab")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))


###### QUICK START #####

# The basic pipeline implemented in GeneSCARAB consists of taking an annotation package,
# a table with the estimated phases for each gene, and a set of GO terms of interest,
# and the software will determine the parameters of the circular distributions for each process,
# as well as significance measures to identify rhythmic processes with respect to one or more
# diel time points.

# Select all GOs that appear in the annotation package
library(GO.db)
library(org.Otauri.eg.db)
library(circular)
library(CircMLE)

# Generate a subset of GOs using the annotation from the package
functional_data <- select(org.Otauri.eg.db, 
                          keys = keys(org.Otauri.eg.db,keytype="GID"), 
                          columns =  c("GID", "GO"))
complete_gos <- unique(functional_data$GO)
complete_gos <- complete_gos[!is.na(complete_gos)]

set.seed(2345)
subset_gos <- complete_gos[sample(1:length(complete_gos), 50, replace = F)]

# Generate GO terms list with genes per GO term (including ancestors of selected GO terms)
go.list.test <- create_gene_list_go(go_vector = subset_gos, org.package = "org.Otauri.eg.db",
                                    go_column = "GO", id_column = "GID")

# Generate GO terms list with phases per GO term for SD
phases.list.sd <- gene.list.to.phases(go.list.test, total_phases_table_sd)

# Remove GOs with no associated rhythmic genes
phases.list.sd.clean <- phases.list.sd[which(sapply(phases.list.sd, nrow) != 0)]
length(phases.list.sd.clean)

# Test preferential concentration of GO terms during the diel cycle
library(circular)

# If Rayleigh or Kuiper tests reject uniformity, in most cases it does not
# make sense to apply HR and Rao, which are tests with lower power in unimodal cases where
# Rayleigh has maximum power; this is controlled by the parameters force.hr,
# force.hr.th, force.rao, and force.rao.th. On the other hand, in multimodal scenarios,
# HR and Rao can be more sensitive. However, a very high number of phases associated 
# with a set almost ensures a departure from uniformity and greatly increases computation 
# time; therefore, the use of these tests for large sets can be controlled 
# using the parameters hr.on.large.sets, hr.on.large.sets.th, rao.on.large.sets, 
# and rao.on.large.sets.th
go.circa.table.sd <- complete_circular_table(phase.list = phases.list.sd.clean, 
                                             hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                                             force.hr = F, force.hr.th = 0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                             force.rao = F, force.rao.th = 0.05,
                                             iter.rao = 999, iter.hr = 999)

# Explore results
# It is  recommended to filter by n > 2 (or more), as HR may give significant results 
# for gene sets with a small number of associated rhythmic genes
go.circa.table.sd <- subset(go.circa.table.sd, n > 2)

#
sum(go.circa.table.sd$rayleigh_p_value_adj < 0.05)
sum(go.circa.table.sd$kuiper_p_value_adj < 0.05)
sum(go.circa.table.sd$hr_p_value_adj < 0.05, na.rm = T)
sum(go.circa.table.sd$rao_p_value_adj < 0.05, na.rm = T)

# Since the tests performed when using HR and Rao are typically only for negative results
# in Rayleigh and Kuiper, the expected distribution of p-values will favor high p-values,
# which can lead to an increase in false negatives when adjusting the
# p-value; therefore, it may be more informative to filter by p-value in these cases.
# In any case, this decision is left to the user, and you can always force the
# tests to be performed
sum(go.circa.table.sd$hr_p_value < 0.05, na.rm = T)
sum(go.circa.table.sd$rao_p_value < 0.05, na.rm = T)

# We can plot these results with 3 different visualizations, let's plot a sample of
# 8 biological processes
# First a circular boxplot. THis function exports a png containing the plot. For 
# colors, MetBrewer's palettes are used
ex_plots <- c("GO:0016833", "GO:0015698", "GO:0043603", "GO:0016790",
              "GO:0055085","GO:0019205", "GO:0004553", "GO:0008169")
plot_table_8 <- go.circa.table.sd[ex_plots,]
circular_boxplot(plot_table_8, filename = "boxplot_example.png", color.palette = "Tam")

# Second a circular dotplot
plot_phase_list_8 <- phases.list.sd.clean[ex_plots]
circular_dotplot(plot_phase_list_8, filename = "dotplot_example.png", color.palette = "Tam")

# And a circular histogram. nbins controls the number of breaks in the histogram
circular_histogram(plot_phase_list_8, filename = "histogram_example.png", color.palette = "Tam", 
                   nbins = 48)

# As we can see, there exists one GO that is significant following HR adjusted p-value but not
# that of Rayleigh or Kuiper, it is probably due to multimodality
multi_table <- go.circa.table.sd[which(go.circa.table.sd$hr_p_value < 0.05),]
multi_go <- rownames(multi_table)
plot(circular(phases.list.sd.clean[[multi_go]]$phase*pi/12))


# These kind of GOs may be further analyzed using GeneSCARAB. multimodal_analysis takes a phase 
# table from an individual process or gene set,which ideally has already been determined to be 
# non-uniform, and allows the user to identify its modality, identifying processes with periods
# of 24, 12, and 8 hours. It also groups the genes into each of the corresponding modes
phase_table_multi_go <- phases.list.sd.clean[[multi_go]]
multimodal_list <- multimodal_analysis(phase_table_multi_go)

# This GO presents significant bimodality, with circular medians in ZT5 and ZT17 and high rhos
# We can also use GeneSCARABA to plot these clusters. When creating the list, remember that it must 
# be a named list
circular_dotplot(phases.list = list(cluster1 = multimodal_list$cluster_1, 
                                    cluster2 = multimodal_list$cluster_2), 
                 filename = "multimodal_go.png", color.palette = "Tam", tr.height = 0.15)


##############################     GROUPED DATA    ########################### 
# Some software returns discrete phases rather than continuous ones; in this case, the pipeline
# is identical to the previous one, but uses the complete_circular_table_grouped function


############# OTHER WAYS OF CREATING GENE SETS LISTS ################

##### KEGG

# Just as in the previous case for GO terms, GeneSCARAB allows you to automatically generate 
# lists of genes associated with specific KEGG pathways based on an annotation package and a set
# of KO terms to include (or calculate it for all those appearing in the annotation package) 
# using the create_gene_list_kegg function instead of create_gene_list_go. In addition, it 
# allows you to filter pathways based on whether they belong to bacteria, archaea, plants,
# animals, fungi, or protists


##### CUSTOM GENE SETS
# GeneSCARAB also allows you to use custom gene lists, for example,
# based on associations found in other databases or previous analyses. All you need to do is define
# a list structure in which each element corresponds to a vector of genes (or elements)
setA <- unique(functional_data$GID)[1:50]
setB <- unique(functional_data$GID)[100:120]
setC <- unique(functional_data$GID)[400:450]

custom.list.test <- list(setA = setA, setB = setB, setC = setC)

# And, as in the previous cases, the gene.list.to.phases function creates the phase table 
# for the rhythmic genes contained in each of the custom gene sets
phases.list.custom <- gene.list.to.phases(custom.list.test, total_phases_table_ld)


###### COMPARISON OF THE DISTRIBUTIONS OF GENE SETS UNDER TWO DIFFERENT CONDITIONS #######

# Different biological processes or gene sets may exhibit significantly different circular distributions
# under different experimental conditions. GeneSCARAB allows for a direct statistical comparison using 
# the test_two_dist function. This is useful, for example, for determining advancement or delay in 
# the diel expression of different processes for different conditions or genotypes
phases.list.ld <- gene.list.to.phases(go.list.test, total_phases_table_ld)
phases.list.ld.clean <- phases.list.ld[which(sapply(phases.list.ld, nrow) != 0)]

go.circa.table.ld <- complete_circular_table(phase.list = phases.list.ld.clean, 
                                             hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                                             force.hr = F, force.hr.th = 0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                             force.rao = F, force.rao.th = 0.05,
                                             iter.rao = 999, iter.hr = 999)

# We select terms with 3 or more elements that are rhythmic according to any of the tests
rhythmic_gos_ld <- rownames(go.circa.table.ld)[go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05 | 
    go.circa.table.ld$hr_p_value < 0.05 | go.circa.table.ld$rao_p_value < 0.05]
rhythmic_gos_ld <- rhythmic_gos_ld[!is.na(rhythmic_gos_ld)]
rhythmic_gos_sd <- rownames(go.circa.table.sd)[go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05 | 
                                                 go.circa.table.sd$hr_p_value < 0.05 | go.circa.table.sd$rao_p_value < 0.05]
rhythmic_gos_sd <- rhythmic_gos_sd[!is.na(rhythmic_gos_sd)]

rhythmic_gos <- intersect(rhythmic_gos_ld, rhythmic_gos_sd)
length(rhythmic_gos)

diff_distributed_gos <- test_two_dist(phases.list.ld.clean[rhythmic_gos], phases.list.sd.clean[rhythmic_gos])
sum(diff_distributed_gos$p_value_adj < 0.05)

# 81 out of 125 rhythmic sets in both conditions significantly change their distribution between LD and SD

# We can plot one of these
go_circa_plot_diff <- rbind(go.circa.table.ld["GO:0008026",], 
                            go.circa.table.sd["GO:0008026",])
rownames(go_circa_plot_diff) <- c("LD", "SD")

circular_boxplot(go_circa_plot_diff, filename = "boxplot_diff_gos.png", color.palette = "Austria")

# Alternatively, we could plot several of them adding more rows to the table



###### IDENTIFICATION OF GENE SETS NOT FOLLOWING THE EXPERIMENTAL PHASE DISTRIBUTION ######

not_experimentally_distributed_gos <- test_against_gen_dist(phases.list.ld.clean[rhythmic_gos], total_phases_table_ld)
sum(not_experimentally_distributed_gos$gen_dist_p_value_adj < 0.05)



####### CONTRIBUTION OF EACH INDIVIDUAL ELEMENT TO THE COMPLETE DISTRIBUTION OF THE SET #######
processes_for_gene_contribution <- rhythmic_gos[100:105]
gene_contributions <- gene_contribution_to_set(circa_result = go.circa.table.ld[processes_for_gene_contribution,],
                                                phase_list = phases.list.ld.clean)

# Identify the genes that have the greatest impact on the mean of each process
ranked_mean <- sapply(gene_contributions, function(x) rownames(x)[1])
ranked_mean

# Identify the genes that have the greatest impact on the variance of each process
ranked_rho <- sapply(gene_contributions, function(x) names(which.min(x[,"Rho_rank"])))
ranked_rho

# And we can visualize these genes with respect to the complete distribution
gene_of_interest_mean <- subset(phases.list.ld.clean$`GO:1902494`, names == "ostta10g03010")
circular_dotplot(phases.list = list(ranked_gene = gene_of_interest_mean, 
                                    complete_go = phases.list.ld.clean$`GO:1902494`), 
                 filename = "gene_contribution.png", color.palette = "Austria", tr.height = 0.15)


