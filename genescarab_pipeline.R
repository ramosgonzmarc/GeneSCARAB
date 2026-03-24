
# Ejemplos de uso de GeneSCARAB
source("genescarab_package.R")

# Cargar datos
circa_table_genescarab <- read.table("../../../Marchantia/rna_marchantia/circa_marchantia_paper/tables/circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table_ld <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["ld.peak.time.hours"]]))
total_phases_table_sd <- data.frame(names=rownames(circa_table_genescarab), phase=as.numeric(circa_table_genescarab[["sd.peak.time.hours"]]))



###### INICIO RÁPIDO PARA IMPACIENTES #####

# La pipeline más básico implementado en GeneSCARAB consiste en tomar un paquete de 
# anotación, una tabla con las fases estimadas para cada gen y un grupo de términos de GO de interés,
# y el software determinará los parámetros de las distribuciones circulares de cada proceso,
# además de medidas de significancia para identificar procesos rítmicos respecto a uno o varios
# puntos temporales diurnos.

# Select all GOs that appear in the annotation package
library(GO.db)
library(org.Mpolymorpha71.eg.db)
library(circular)
library(CircMLE)

# Generate a subset of GOs using the annotation from the package
functional_data <- select(org.Mpolymorpha71.eg.db, 
                          keys = keys(org.Mpolymorpha71.eg.db,keytype="GID"), 
                          columns =  c("GID", "GO")) # load annotation data
complete_gos <- unique(functional_data$GO)
complete_gos <- complete_gos[!is.na(complete_gos)]
subset_gos <- complete_gos[sample(1:length(complete_gos), 30, replace = F)]

# Generate GO terms list with genes per GO term (including ancestors of selected GO terms)
go.list.test <- create_gene_list_go(go_vector = subset_gos, org.package = "org.Mpolymorpha71.eg.db",
                                    go_column = "GO", id_column = "GID")

# Generate GO terms list with phases per GO term for LD and SD
phases.list.ld <- gene.list.to.phases(go.list.test, total_phases_table_ld)

# Remove GOs with no associated rhythmic genes
phases.list.ld.clean <- phases.list.ld[which(sapply(phases.list.ld, nrow) != 0)]

# Test preferential concentration of GO terms during the diel cycle
library(circular)

go.circa.table.ld <- complete_circular_table(phase.list = phases.list.ld.clean, 
                                             hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                                             force.hr = F, force.hr.th = 0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                             force.rao = F, force.rao.th = 0.05,
                                             iter.rao = 999, iter.hr = 999)

# Explore results
sum(go.circa.table.ld$rayleigh_p_value_adj < 0.05)
sum(go.circa.table.ld$kuiper_p_value_adj < 0.05)
sum(go.circa.table.ld$hr_p_value_adj < 0.05, na.rm = T)
sum(go.circa.table.ld$rao_p_value_adj < 0.05, na.rm = T)


# HR positive results may be due to low number of associated genes
sum(go.circa.table.ld[go.circa.table.ld$n > 3,]$hr_p_value_adj < 0.05, na.rm = T)
plot(circular(phases.list.ld$`GO:0140053`$phase))



############### DATOS CIRCULARES AGRUPADOS ################ 




############# OTRAS MANERAS DE CREAR LISTAS DE GENES ################

##### KEGG

# Igual que en el caso anterior para los términos de GO, GeneSCARAB permite crear
# automáticamente las listas de genes asociados a rutas KEGG concretas a partir de 
# un paquete de anotación y un set de términos KO a incluir (o calcularlo para todos
# los que aparezcan en el paquete de anotación). Además, permite filtrar las rutas según pertenezcan a
# bacteria, archaea, plants, animals, fungi or protists.
library(KEGGREST)
kegg.list.test <- create_gene_list_kegg(ko_vector = "all", org.package = "org.Mpolymorpha71.eg.db", 
                                  ko_column = "KEGG", id_column = "GID", 
                                  ko_prefix = "map", species = "plants")

phases.list.kegg <- gene.list.to.phases(kegg.list.test, total_phases_table_ld)


##### CUSTOM GENE SETS

# GeneSCARAB también permite emplear listas de genes customizadas, por ejemplo,
# por asociaciones en otras bases de datos o análisis previos. Basta con definir
# una estructura de lista en la que cada elemento corresponda a un vector de genes
setA <- unique(functional_data$GID)[1:50]
setB <- unique(functional_data$GID)[100:120]
setC <- unique(functional_data$GID)[400:450]

custom.list.test <- list(setA = setA, setB = setB, setC = setC)

# Y, como en los casos anteriores, la función gene.list.to.phases crea la tabla de 
# fases para los genes rítmicos contenidos en cada uno de los custom gene sets
phases.list.custom <- gene.list.to.phases(custom.list.test, total_phases_table_ld)



###### COMPARACIÓN DE LA RITMICIDAD DE LOS GRUPOS DE GENES EN DOS CONDICIONES DISTINTAS #######

# Los distintos procesos biológicos o grupos de genes pueden mostrar distribuciones circulares
# significativamente diferentes en distintas condiciones experimentales. GeneSCARAB permite 
# realizar un cálculo estadístico directo mediante la función test_two_dist. Esto sirve, por ejemplo,
# para determinar adelantamientos o retrasos en la expresión diurna de los distintos procesos para
# diferentes condiciones o genotipos
phases.list.sd <- gene.list.to.phases(go.list.test, total_phases_table_sd)
phases.list.sd.clean <- phases.list.sd[which(sapply(phases.list.sd, nrow) != 0)]

go.circa.table.sd <- complete_circular_table(phase.list = phases.list.sd.clean, 
                                             hr.on.large.sets = F , hr.on.large.sets.th = 200, 
                                             force.hr = F, force.hr.th = 0.05, 
                                             rao.on.large.sets = F, rao.on.large.sets.th = 200, 
                                             force.rao = F, force.rao.th = 0.05,
                                             iter.rao = 999, iter.hr = 999)

# Seleccionamos los términos con 3 o más genes que sean rítmicos según alguno de los tests
rhythmic_gos_ld <- rownames(go.circa.table.ld)[which(go.circa.table.ld$n > 3 & (go.circa.table.ld$rayleigh_p_value_adj < 0.05 | go.circa.table.ld$kuiper_p_value_adj < 0.05 | 
    go.circa.table.ld$hr_p_value_adj < 0.05 | go.circa.table.ld$rao_p_value_adj < 0.05))]
rhythmic_gos_sd <- rownames(go.circa.table.sd)[which(go.circa.table.sd$n > 3 & (go.circa.table.sd$rayleigh_p_value_adj < 0.05 | go.circa.table.sd$kuiper_p_value_adj < 0.05 | 
                                                                                  go.circa.table.sd$hr_p_value_adj < 0.05 | go.circa.table.sd$rao_p_value_adj < 0.05))]

rhythmic_gos <- intersect(rhythmic_gos_ld, rhythmic_gos_sd)

diff_distributed_gos <- test_two_dist(phases.list.ld.clean[rhythmic_gos], phases.list.sd.clean[rhythmic_gos])
sum(diff_distributed_gos$p_value_adj < 0.05)




###### IDENTIFICACIÓN DE GRUPOS DE GENES QUE NO SIGUEN LA DINÁMICA GENERAL DE FASES ##########
not_experimentally_distributed_gos <- test_against_gen_dist(phases.list.ld.clean[rhythmic_gos], total_phases_table_ld)
sum(not_experimentally_distributed_gos$gen_dist_p_value_adj < 0.05)



####### CONTRIBUCIÓN DE CADA GEN INDIVIDUAL A LA DISTRIBUCIÓN CIRCULAR DEL GRUPO DE GENES #######
processes_for_gene_contribution <- rhythmic_gos[1:5]
gene_contributions <- gene_contribution_to_set(circa_result = go.circa.table.ld[processes_for_gene_contribution,],
                                                phase_list = phases.list.ld.clean)

# Comprobamos los genes que más impacto tienen en la media de cada proceso
sapply(gene_contributions, function(x) rownames(x)[1])

# Comprobamos los genes que más impacto tienen en la varianza de cada proceso
sapply(gene_contributions, function(x) names(which.min(x[,"Rho_rank"])))



######## DETECTAR Y CARACTERIZAR ESCENARIOS MULTIMODALES #########


######## REPRESENTAR LOS RESULTADOS ##########

###### BOXPLOTS

###### DOTPLOTS

###### HISTOGRAMAS

