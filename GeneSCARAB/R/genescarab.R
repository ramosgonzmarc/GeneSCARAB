
#' Critical value calculation for HR test
#'
#' A helper function that calculates the critical value for the HR test,
#' as implemented in the article "Grouped circular data
#' in biology: advice for effectively implementing statistical procedures"
#' https://doi.org/10.1007/s00265-020-02881-6.
#'
#' @param sample Data points in radians
#'
#' @returns The critical value for the test
#'
#' @examples
HermansRasson2T <- function(sample){
  n <- length(sample)
  total <- 0
  for (i in 1:n){
    for (j in 1:n){ total <- total + abs(abs(sample[i]-sample[j])-pi)-
      (pi/2)
    total <- total - (2.895*(abs(sin(sample[i]-sample[j]))-(2/pi)))}}
  T <- total/n
  return(T)}

#' HR test's p-value calculation for grouped data
#'
#' Helper function that calculates the p-value of an adaptation of the Hermans-Rasson test
#' for grouped data (discrete phases) as implemented in the article "Grouped circular data
#' in biology: advice for effectively implementing statistical procedures"
#' https://doi.org/10.1007/s00265-020-02881-6.
#'
#'
#' @param sample Data points in radians
#' @param m Number of bins
#' @param k Kappa for the error distribution
#' @param iter Number of iterations
#'
#' @returns The p-value for the test
#'
#' @examples
HermansRasson2PGroupedRad <- function(sample, m, k=1000, iter=9999){
  sample<-circular::circular(sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  n <- length(sample)
  univals <- iter
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- circular::rcircularuniform(n, control.circular=list(units="radians"))
    data1 <- trunc(data1*m/(2*pi))
    data1 <- data1*2*pi/m
    errorsamp <- circular::rvonmises(n, 0, k,control.circular=list(units="radians"))
    data1 <- data1+errorsamp
    data1 <- ifelse((data1>(2*pi)),(data1-(2*pi)), data1)
    testset[f] <- HermansRasson2T(data1)}
  errorsamp2 <- circular::rvonmises(n, 0, k,control.circular=list(units="radians"))
  sample<-sample+errorsamp2
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  Tsample <- HermansRasson2T(sample)
  counter <- 0
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}


#' Critical value calculation for Rao's spacing test
#'
#' Helper function that calculates the critical value for Rao's spacing test as
#' implemented in the article "Grouped circular data in biology: advice for effectively
#' implementing statistical procedures" https://doi.org/10.1007/s00265-020-02881-6.
#'
#' @param sample Data points in radians
#'
#' @returns The critical value for the test
#'
#' @examples
RaoTestValue <- function(sample){n <- length(sample)
f <- sort(sample)
fplus <- c(f[2:n],f[1])
T <- fplus - f
T[n] <- (2*pi) - f[n] + f[1]
abs_diff <- abs(T-(2*pi/n))
U <- 0.5*sum(abs_diff)
return(U)}

#' Rao's spacing test's p-value calculation for ungrouped data
#'
#'
#' Helper function that calculates the p-value of Rao's spacing test
#' for ungrouped data (continuous phases) as implemented in the article "Circular
#' statistics meets practical limitations: a simulation-based Rao’s spacing
#' test for non-continuous data" https://doi.org/10.1186/s40462-019-0160-x.
#'
#' @param sample Data points in radians
#' @param iter Number of iterations
#'
#' @returns The p-value for the test
#'
#' @examples
RaoTestUngroupedRad <- function(sample, iter=9999){
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  n <- length(sample)
  univals <- iter
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- circular::rcircularuniform(n, control.circular=list(units="radians"))
    testset[f] <- RaoTestValue(data1)}
  Tsample <- RaoTestValue(sample)
  counter <- 0
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}


#' Rao's spacing test's p-value calculation for grouped data
#'
#' Helper function that calculates the p-value of an adaptation of the Rao's spacing test
#' for grouped data (discrete phases) as implemented in the article "Grouped circular data in biology: advice for effectively
#' implementing statistical procedures" https://doi.org/10.1007/s00265-020-02881-6.
#'
#' @param sample Data points in radians
#' @param m Number of bins
#' @param k Kappa for the error distribution
#' @param iter Number of iterations
#'
#' @returns The p-value for the test
#'
#' @examples
RaoPGroupedRad <- function(sample, m, k=1000, iter=9999){
  sample<-circular::circular(sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  n <- length(sample)
  univals <- iter
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- circular::rcircularuniform(n, control.circular=list(units="radians"))
    data1 <- trunc(data1*m/(2*pi))
    data1 <- data1*2*pi/m
    errorsamp <- circular::rvonmises(n, 0, k,control.circular=list(units="radians"))
    data1 <- data1+errorsamp
    data1 <- ifelse((data1>(2*pi)),(data1-(2*pi)), data1)
    testset[f] <- RaoTestValue(data1)}
  errorsamp2 <- circular::rvonmises(n, 0, k,control.circular=list(units="radians"))
  sample<-sample+errorsamp2
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  Tsample <- RaoTestValue(sample)
  counter <- 0
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}



#' Creation of input lists for GO terms from an annotation package
#'
#' Function to create a list of genes associated with each GO term based on a
#' species annotation package. This includes all ancestor GO terms of the selected ones.
#'
#' @param go_vector GO terms vector to be included in the analysis, along with their ancestors. "all" takes the complete set from the annotation package
#' @param org.package Species annotation package
#' @param go_column GO term column name in the annotation package
#' @param id_column Gene ID column name in the annotation package
#'
#' @returns A list linking each GO term with the genes that share that annotation
#' @export
#'
#' @examples
create_gene_list_go <- function(go_vector = c("all"), org.package, go_column, id_column)
{

  functional_data <- select(get(org.package),
                            keys = keys(get(org.package),keytype=id_column),
                            columns =  c(id_column, go_column))
  if (length(go_vector) == 1 & go_vector[1] == "all")
    {
    go_vector <- unique(functional_data[[go_column]])
    go_vector <- go_vector[!is.na(go_vector)]
    }


  # Get ancestors of selected GOs for BP, MF and CC
  bp_ancestors <- BiocGenerics::mget(go_vector, GO.db::GOBPANCESTOR, ifnotfound=NA)
  mf_ancestors <- BiocGenerics::mget(go_vector, GO.db::GOMFANCESTOR, ifnotfound=NA)
  cc_ancestors <- BiocGenerics::mget(go_vector, GO.db::GOCCANCESTOR, ifnotfound=NA)

  l <- list(bp_ancestors, mf_ancestors, cc_ancestors)
  l_keys <- unique(unlist(lapply(l, names)))

  go_ancestors <- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, l_keys))), l_keys)
  go_ancestors <- lapply(go_ancestors, function(x) x[!is.na(x)])
  go_ancestors <- lapply(go_ancestors, function(x) x[x != "all"])
  go_vec_ancestors <- unlist(go_ancestors)
  names(go_vec_ancestors) <- NULL

  # Add go terms in study to the ancestors list to get successors
  query_vec_ancestors <- c(go_vec_ancestors, go_vector)
  query_vec_ancestors <- unique(query_vec_ancestors)

  # Get succesors for the complete list of current GO terms and their ancestors
  bp_offspring <- BiocGenerics::mget(query_vec_ancestors, GO.db::GOBPOFFSPRING, ifnotfound=NA)
  mf_offspring <- BiocGenerics::mget(query_vec_ancestors, GO.db::GOMFOFFSPRING, ifnotfound=NA)
  cc_offspring <- BiocGenerics::mget(query_vec_ancestors, GO.db::GOCCOFFSPRING, ifnotfound=NA)

  # Iterate and get genes for each parental based in its offspring GO
  l_offspring <- list(bp_offspring, mf_offspring, cc_offspring)
  l_offspring_keys <- unique(unlist(lapply(l_offspring, names)))
  go_offspring <- setNames(do.call(mapply, c(FUN=c, lapply(l_offspring, `[`, l_offspring_keys))), l_offspring_keys)

  go_offspring_clean <- lapply(go_offspring, function(x) x[!is.na(x)])

  # Add self to each list
  go_offspring_comp <- mapply(function(x, y) c(x,y), go_offspring_clean, names(go_offspring_clean))

  # Direct subset
  go.list <- lapply(go_offspring_comp, function(x) subset(functional_data, functional_data[[go_column]] %in% x)[[id_column]])
  go.list <- lapply(go.list, unique)

  return(go.list)

}


#' Creation of input lists for KEGG pathways from an annotation package
#'
#' Function to create a list of genes associated with each KEGG pathway based on a
#' species annotation package. It also allows filtering by species group.
#'
#' @param ko_vector KO vector to be included in the analysis, along with their pathways ancestors. "all" takes the complete set from the annotation package
#' @param org.package Species annotation package
#' @param ko_column KO column name in the annotation package
#' @param id_column Gene ID column name in the annotation package
#' @param ko_prefix Prefix used to denote KEGG pathways, either "map" or "ko", depending on package version
#' @param species Species group to filter final pathways ("all" for no filtering), either "bacteria", "archaea", "animals", "plants", "fungi" or "protists"
#'
#' @returns A list linking each KEGG pathway with the genes that share that annotation
#' @export
#'
#' @examples
create_gene_list_kegg <- function(ko_vector = c("all"), org.package, ko_column, id_column,
                                  ko_prefix = "map", species = "all")
{
  functional_data <- select(get(org.package),
                            keys = keys(get(org.package),keytype=id_column),
                            columns =  c(id_column, ko_column))

  if (length(ko_vector) == 1 & ko_vector[1] == "all")
  {

    complete_kos <- unique(functional_data[[ko_column]])
    ko_vector <- complete_kos[!is.na(complete_kos)]
  }

  # Get pathways for selected KOs
  {
  if (length(ko_vector) < 500)
  {
    present_pathways <- KEGGREST::keggLink("pathway", paste0(ko_vector, collapse = "+"))
  }
  else
  {
    present_pathways <- KEGGREST::keggLink("pathway", paste0(ko_vector[1:500], collapse = "+"))

    for (i in 2:(length(ko_vector)%/%500 + 1))
      {
      present_pathways <- append(present_pathways, KEGGREST::keggLink("pathway", paste0(ko_vector[(1+(500*(i-1))):(i*500)], collapse = "+")))
      }
  }
  }

  unique_pathways <- unique(unlist(present_pathways))

  map_pathways <- unique_pathways[grep(ko_prefix, unique_pathways)]
  list_pathways <- sapply(strsplit(map_pathways, ":"), function(x) x[[2]])


  {
    if (length(list_pathways) < 500)
    {
      kos_pathways <- KEGGREST::keggLink("ko", paste0(list_pathways, collapse = "+"))
    }
    else
    {
      kos_pathways <- KEGGREST::keggLink("ko", paste0(list_pathways[1:500], collapse = "+"))

      for (i in 2:(length(list_pathways)%/%500 + 1))
      {
        kos_pathways <- append(kos_pathways, KEGGREST::keggLink("ko", paste0(list_pathways[(1+(500*(i-1))):(i*500)], collapse = "+")))
      }
    }
  }

  kos_pathways <- tapply(kos_pathways, INDEX = names(kos_pathways), FUN = function(x) unlist(x), simplify = F)
  names(kos_pathways) <- sapply(strsplit(names(kos_pathways), split = "path:"), function(x) x[[2]])

  kos_pathways_clean <- lapply(kos_pathways, function(x) unique(x))
  kos_pathways_clean <- lapply(kos_pathways_clean, function(x) sapply(strsplit(x, ":"), function(y) y[[2]]))
  ko_list <- lapply(kos_pathways_clean, function(x) subset(functional_data, functional_data[[ko_column]] %in% x)[[id_column]])
  ko_list <- lapply(ko_list, unique)

  {
    if (species == "all")
    {
      return(ko_list)
    }
    else
    {
      kegg_to_keep <- kegg_pathways_per_species[[species]]
      kegg_to_keep <- paste0(ko_prefix, kegg_to_keep)
      ko_list <- ko_list[names(ko_list) %in% kegg_to_keep]
      return(ko_list)
    }
  }

}


#' Creation of phase table lists
#'
#' Function that takes a list of genes associated with specific biological processes (GO, KEGG) or gene sets
#' and a general phase table, and returns a list consisting of the phase tables specific to each set.
#'
#' @param gene.list Process or gene sets list (as the output of create_gene_list_go and create_gene_list_kegg)
#' @param phase.table Table showing acrophases for the complete gene set under study, stated in hours
#'
#' @returns A list containing phase tables for each process or gene set
#' @export
#'
#' @examples
gene.list.to.phases <- function(gene.list, phase.table)
{
  circa_reduced <- lapply(gene.list, function(x) subset(phase.table, names %in% x))
  circa_reduced_clean <- circa_reduced[which(sapply(circa_reduced, nrow) != 0)]
  return(circa_reduced_clean)

}


#' Creation of the circular table and deviation from the uniform distribution for a single set
#'
#' Function that takes a vector of gene phases (continuous measurements or parametric estimates)
#' and returns statistical measures of their circular distribution and their deviations from a
#' circular uniform distribution based on the Rayleigh, Kuiper, Hermans-Rasson and Rao tests.
#'
#' @param circa_vector Phases vector
#' @param hr.on.large.sets Whether to perform HR test on sets of a greater size than the number stated in hr.on.large.sets.th, logical
#' @param hr.on.large.sets.th Threshold indicating the maximum size of a set to perform HR test on
#' @param iter.hr Number of iterations for HR test
#' @param force.hr Whether to perform HR test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.hr.th Significance level to use for force.hr
#' @param rao.on.large.sets Whether to perform Rao test on sets of a greater size than the number stated in rao.on.large.sets.th, logical
#' @param rao.on.large.sets.th Threshold indicating the maximum size of a set to perform Rao test on
#' @param iter.rao Number of iterations for Rao test
#' @param force.rao Whether to perform Rao test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.rao.th Significance level to use for force.rao
#'
#' @returns A table indicating number of phases used for computation (n); mean, rho, variance, first, second and third quantiles of the circular distribution in hours; and p-values for each of the four tests
#' @export
#'
#' @examples
#'
create_circular_table <- function(circa_vector,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                  iter.hr = 999, force.hr=F, force.hr.th=0.05,
                                  rao.on.large.sets = F, rao.on.large.sets.th=400,
                                  iter.rao = 999, force.rao=F, force.rao.th=0.05)
{
  circa_radians <- circular::circular(circa_vector*pi/12)

  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"

  # Calculate circular variance and transform to hours again
  circa_var <- circular::var.circular(circa_radians)

  # Apply kupier's test
  circa_kupier <- Directional::kuiper(circa_radians, rads = T, R=1)$p.value

  # Apply Rayleigh's test
  circa_ray <- circular::rayleigh.test(circa_radians, mu=circular::circular(circa_summary["Mean"]))$p.value

  # Apply HR test
  if (iter.hr == 0)
  {
    circa_hr <- NA
  }

  else if (!(force.hr) & (circa_kupier <= force.hr.th | circa_ray <= force.hr.th))
  {
    circa_hr <- NA
  }

  else if (!(hr.on.large.sets) & (length(circa_radians) >= hr.on.large.sets.th))
  {
    circa_hr <- NA
  }

  else
  {
    circa_hr <- CircMLE::HR_test(circa_radians, original = F, iter = iter.hr)
  }

  # Apply Rao test
  if (iter.rao == 0)
  {
    circa_rao <- NA
  }

  else if (!(force.rao) & (circa_kupier <= force.rao.th | circa_ray <= force.rao.th))
  {
    circa_rao <- NA
  }

  else if (!(rao.on.large.sets) & (length(circa_radians) >= rao.on.large.sets.th))
  {
    circa_rao <- NA
  }

  else
  {
    circa_rao <- RaoTestUngroupedRad(circa_radians, iter = iter.rao)
  }

  # Create table
  circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi,
                            median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
                            third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"],
                            var=circa_var, kuiper_p_value=circa_kupier, rayleigh_p_value=circa_ray,
                            hr_p_value = circa_hr[2], rao_p_value = circa_rao )


  return(circa_table)

}


#' Creation of the circular table and deviation from the uniform distribution for a single set of grouped data
#'
#' Function that takes a vector of gene phases (discrete measurements or non-parametric estimates)
#' and returns statistical measures of their circular distribution and their deviations from a
#' circular uniform distribution based on adapted Rayleigh, Kuiper, Hermans-Rasson and Rao tests.
#'
#'
#' @param circa_vector Discrete phases vector
#' @param n_bins Number of bins
#' @param hr.on.large.sets Whether to perform HR test on sets of a greater size than the number stated in hr.on.large.sets.th, logical
#' @param hr.on.large.sets.th Threshold indicating the maximum size of a set to perform HR test on
#' @param iter.hr Number of iterations for HR test
#' @param force.hr Whether to perform HR test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.hr.th Significance level to use for force.hr
#' @param rao.on.large.sets Whether to perform Rao test on sets of a greater size than the number stated in rao.on.large.sets.th, logical
#' @param rao.on.large.sets.th Threshold indicating the maximum size of a set to perform Rao test on
#' @param iter.rao Number of iterations for Rao test
#' @param force.rao Whether to perform Rao test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.rao.th Significance level to use for force.rao
#' @param error_kappa Kappa for the error distribution for HR and Rao grouped tests
#'
#' @returns A data.frame indicating number of phases used for computation (n); mean, rho, variance, first, second and third quantiles of the circular distribution in hours; and p-values for each of the four tests
#' @export
#'
#' @examples
create_circular_table_grouped <- function(circa_vector, n_bins,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                  iter.hr = 999, force.hr=F, force.hr.th=0.05,
                                  rao.on.large.sets = F, rao.on.large.sets.th=400,
                                  iter.rao = 999, force.rao=F, force.rao.th=0.05, error_kappa=1000)
{

  circa_radians <- circular::circular(circa_vector*pi/12)

  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"

  # Calculate circular variance and transform to hours again
  circa_var <- circular::var.circular(circa_radians)

  # Apply kupier's test
  circa_kupier <- Directional::kuiper(circa_radians, rads = T, R=1)$p.value

  # Apply Rayleigh's test
  circa_ray <- circular::rayleigh.test(circa_radians, mu=circular::circular(circa_summary["Mean"]))$p.value

  # Apply HR test
  if (is.na(circa_ray))
  {
    circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi,
                              median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
                              third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"],
                              var=circa_var, kuiper_p_value=circa_kupier, rayleigh_p_value=1,
                              hr_p_value = 1, rao_p_value = 1 )


    return(circa_table)
  }

  if (iter.hr == 0)
  {
    circa_hr <- NA
  }

  else if (!(force.hr) & (circa_kupier <= force.hr.th | circa_ray <= force.hr.th))
  {
    circa_hr <- NA
  }

  else if (!(hr.on.large.sets) & (length(circa_radians) >= hr.on.large.sets.th))
  {
    circa_hr <- NA
  }

  else
  {
    circa_hr <- HermansRasson2PGroupedRad(circa_radians, m=n_bins, k=error_kappa, iter=iter.hr)
  }

  # Apply Rao test
  if (iter.rao == 0)
  {
    circa_rao <- NA
  }

  else if (!(force.rao) & (circa_kupier <= force.rao.th | circa_ray <= force.rao.th))
  {
    circa_rao <- NA
  }

  else if (!(rao.on.large.sets) & (length(circa_radians) >= rao.on.large.sets.th))
  {
    circa_rao <- NA
  }

  else
  {
    circa_rao <- RaoPGroupedRad(circa_radians, m=n_bins, k=error_kappa, iter=iter.rao)
  }

  # Create table
  circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi,
                            median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
                            third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"],
                            var=circa_var, kuiper_p_value=circa_kupier, rayleigh_p_value=circa_ray,
                            hr_p_value = circa_hr[2], rao_p_value = circa_rao )


  return(circa_table)

}


#' Creation of the circular table and deviation from the uniform distribution for a phases list
#'
#' Function that takes a list of continuous phases associated with each process or gene set (output of the
#' gene.list.to.phases function) and calculates the table of statistical measures and deviations from circular uniformity
#' for each one using the Rayleigh, Kuiper, Hermans-Rasson, and Rao tests.
#'
#' @param phase.list Phases list in the same format as the output of gene.list.to.phases
#' @param hr.on.large.sets Whether to perform HR test on sets of a greater size than the number stated in hr.on.large.sets.th, logical
#' @param hr.on.large.sets.th Threshold indicating the maximum size of a set to perform HR test on
#' @param iter.hr Number of iterations for HR test
#' @param force.hr Whether to perform HR test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.hr.th Significance level to use for force.hr
#' @param rao.on.large.sets Whether to perform Rao test on sets of a greater size than the number stated in rao.on.large.sets.th, logical
#' @param rao.on.large.sets.th Threshold indicating the maximum size of a set to perform Rao test on
#' @param iter.rao Number of iterations for Rao test
#' @param force.rao Whether to perform Rao test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.rao.th Significance level to use for force.rao
#'
#' @returns A data.frame indicating, for each process or gene set: number of phases used for computation (n); mean, rho, variance, first, second and third quantiles of the circular distribution in hours; p-values and adjusted p-values for each of the four tests
#' @export
#'
#' @examples
complete_circular_table <- function(phase.list, hr.on.large.sets = F, hr.on.large.sets.th=400,
                                    iter.hr = 999, force.hr=F, force.hr.th=0.05,
                                    rao.on.large.sets = F, rao.on.large.sets.th=400,
                                    iter.rao = 999, force.rao=F, force.rao.th=0.05)
{

  # Generate circular statistics table per GO
  go_circa_res <- lapply(phase.list, function(x) create_circular_table(x$phase, hr.on.large.sets = hr.on.large.sets,
                                                                       hr.on.large.sets.th=hr.on.large.sets.th, iter.hr=iter.hr,
                                                                       force.hr=force.hr, force.hr.th=force.hr.th,
                                                                       rao.on.large.sets = rao.on.large.sets,
                                                                       rao.on.large.sets.th=rao.on.large.sets.th, iter.rao=iter.rao,
                                                                       force.rao=force.rao, force.rao.th=force.rao.th))

  # Adjust p-values due to multiple testing
  go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
  go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
  kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
  go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
  ray_bh <- p.adjust(go_circa_ray, method = "BH")
  go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
  rh_bh <- p.adjust(go_circa_hr, method = "BH")
  go_circa_rao <- sapply(go_circa_res, function(x) unlist(x[11]))
  rao_bh <- p.adjust(go_circa_rao, method = "BH")

  # Return updated table
  go_circa_table <- data.frame(go_circa_num, kuiper_p_value=go_circa_kupier, kuiper_p_value_adj=kupier_bh,
                               rayleigh_p_value=go_circa_ray, rayleigh_p_value_adj=ray_bh,
                               hr_p_value=go_circa_hr, hr_p_value_adj=rh_bh, rao_p_value=go_circa_rao, rao_p_value_adj=rao_bh)
  rownames(go_circa_table) <- names(phase.list)


  # Transform negative means and quantiles to positive
  pos_res <- apply(go_circa_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
  go_circa_table[,2:5] <- pos_res

  return(go_circa_table)
}


#' Creation of the circular table and deviation from the uniform distribution for a discrete (grouped) phases list
#'
#' Function that takes a list of discrete phases associated with each process or gene set (output of the
#' gene.list.to.phases function) and calculates the table of statistical measures and deviations from circular uniformity
#' for each one using the adapted Rayleigh, Kuiper, Hermans-Rasson, and Rao tests.
#'
#' @param phase.list Phases list in the same format as the output of gene.list.to.phases
#' @param n_bins Number of bins
#' @param hr.on.large.sets Whether to perform HR test on sets of a greater size than the number stated in hr.on.large.sets.th, logical
#' @param hr.on.large.sets.th Threshold indicating the maximum size of a set to perform HR test on
#' @param iter.hr Number of iterations for HR test
#' @param force.hr Whether to perform HR test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.hr.th Significance level to use for force.hr
#' @param rao.on.large.sets Whether to perform Rao test on sets of a greater size than the number stated in rao.on.large.sets.th, logical
#' @param rao.on.large.sets.th Threshold indicating the maximum size of a set to perform Rao test on
#' @param iter.rao Number of iterations for Rao test
#' @param force.rao Whether to perform Rao test on sets that have rejected uniformity for Rayleigh or Kuiper tests at a significance level of force.hr.th, logical
#' @param force.rao.th Significance level to use for force.rao
#' @param error_kappa Kappa for the error distribution for HR and Rao grouped tests
#'
#' @returns A data.frame indicating, for each process or gene set: number of phases used for computation (n); mean, rho, variance, first, second and third quantiles of the circular distribution in hours; p-values and adjusted p-values for each of the four tests
#' @export
#'
#' @examples
complete_circular_table_grouped <- function(phase.list, n_bins,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                    iter.hr = 999, force.hr=F, force.hr.th=0.05,
                                    rao.on.large.sets = F, rao.on.large.sets.th=400,
                                    iter.rao = 999, force.rao=F, force.rao.th=0.05, error_kappa=1000)
{

  # Generate circular statistics table per GO
  go_circa_res <- lapply(phase.list, function(x) create_circular_table_grouped(x$phase, n_bins=n_bins, hr.on.large.sets = hr.on.large.sets,
                                                                       hr.on.large.sets.th=hr.on.large.sets.th, iter.hr=iter.hr,
                                                                       force.hr=force.hr, force.hr.th=force.hr.th,
                                                                       rao.on.large.sets = rao.on.large.sets,
                                                                       rao.on.large.sets.th=rao.on.large.sets.th, iter.rao=iter.rao,
                                                                       force.rao=force.rao, force.rao.th=force.rao.th, error_kappa=error_kappa))

  # Adjust p-values due to multiple testing
  go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
  go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
  kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
  go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
  ray_bh <- p.adjust(go_circa_ray, method = "BH")
  go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
  rh_bh <- p.adjust(go_circa_hr, method = "BH")
  go_circa_rao <- sapply(go_circa_res, function(x) unlist(x[11]))
  rao_bh <- p.adjust(go_circa_rao, method = "BH")

  # Return updated table
  go_circa_table <- data.frame(go_circa_num, kuiper_p_value=go_circa_kupier, kuiper_p_value_adj=kupier_bh,
                               rayleigh_p_value=go_circa_ray, rayleigh_p_value_adj=ray_bh,
                               hr_p_value=go_circa_hr, hr_p_value_adj=rh_bh, rao_p_value=go_circa_rao, rao_p_value_adj=rao_bh)
  rownames(go_circa_table) <- names(phase.list)


  # Transform negative means and quantiles to positive
  pos_res <- apply(go_circa_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
  go_circa_table[,2:5] <- pos_res

  return(go_circa_table)
}

#' Deviation from the total experimental distribution for a phases list
#'
#' Function that compares whether the distributions of specific processes or gene sets differ
#' significantly from the overall phase distribution of the entire study. Significance is calculated
#' based on a MANOVA analysis as described in https://doi.org/10.1186/s40462-022-00323-8.
#'
#' @param phase.list Phases list in the same format as the output of gene.list.to.phases
#' @param total.phase.table Table showing acrophases for the complete gene set under study, stated in hours
#'
#' @returns A data.frame with p-value and adjusted p-value (null hyphotesis: same distribution as the general one) for each process or gene set
#' @export
#'
#' @examples
test_against_gen_dist <- function(phase.list, total.phase.table)
{
  gen_dist <- total.phase.table$phase*pi/12
  len_gen_dist <- length(gen_dist)

  # p_values_gen_test <- sapply(phase.list, function(x) summary(manova(cbind(sin(c(gen_dist, x$phase*pi/12)),cos(c(gen_dist, x$phase*pi/12)))
  #                                                                    ~ c(rep("general", len_gen_dist), rep("specific", nrow(x)))))$stats[1,6])

  p_values_gen_test <- sapply(phase.list, function(x) {
    angles <- c(gen_dist, x$phase * pi / 12)

    if (length(unique(x$phase)) == 1) {
      angles <- angles + rnorm(length(angles), mean = 0, sd = 0.01)
    }

    p_values_gen_test <- summary(manova(cbind(sin(angles),cos(angles))
                                        ~ c(rep("general", len_gen_dist), rep("specific", nrow(x)))))$stats[1,6]
  })

  gen_dist_bh <- p.adjust(p_values_gen_test, method = "BH")

  gen_dist_table <- data.frame(gen_dist_p_value=p_values_gen_test,
                               gen_dist_p_value_adj=gen_dist_bh)

  return(gen_dist_table)
}


#' Comparison of circular distributions between two conditions
#'
#' Function that compares whether the distributions of specific processes or gene sets differ
#' significantly between two different conditions or genotypes. The two lists must have the same
#' number of elements, although the number of genes associated with each process or gene set may differ.
#' Significance is calculated based on a MANOVA analysis as described in https://doi.org/10.1186/s40462-022-00323-8.
#'
#'
#' @param phase.list.1 Phases list for condition 1, in the same format as the output of gene.list.to.phases
#' @param phase.list.2 Phases list for condition 2, in the same format as the output of gene.list.to.phases
#'
#' @returns A data.frame with p-value and adjusted p-value (null hyphotesis: same distribution in both conditions) for each process or gene set
#' @export
#'
#' @examples
test_two_dist <- function(phase.list.1, phase.list.2)
{
  if (length(phase.list.1) != length(phase.list.2))
  {
    stop("Error: length of lists differ.")
  }

  if (sum(names(phase.list.1) == names(phase.list.2)) != length(phase.list.1))
  {
    warning("The names of the sets do not match, check if they should.")
  }


  # p_values_diff_test <- mapply(function(x, y) summary(manova(cbind(sin(c(x$phase*pi/12, y$phase*pi/12)),cos(c(x$phase*pi/12, y$phase*pi/12)))
  #                                                            ~ c(rep("first", nrow(x)), rep("second", nrow(y)))))$stats[1,6],
  #                              phase.list.1, phase.list.2)

  p_values_diff_test <- mapply(function(x,y) {
    angles <- c(x$phase * pi / 12, y$phase * pi / 12)

    if (length(unique(x$phase)) == 1 | length(unique(y$phase)) == 1) {
      angles <- angles + rnorm(length(angles), mean = 0, sd = 0.01)
    }

    p_values_gen_test <- summary(manova(cbind(sin(angles),cos(angles))
                                        ~ c(rep("first", nrow(x)), rep("second", nrow(y)))))$stats[1,6]
  }, phase.list.1, phase.list.2)

  p_values_diff_bh <- p.adjust(unlist(p_values_diff_test), method = "BH")

  p_values_diff_table <- data.frame(p_value=unlist(p_values_diff_test),
                                    p_value_adj=p_values_diff_bh)

  return(p_values_diff_table)
}


#' Contribution of each gene in a set to the temporal cohesion of the whole set
#'
#' Function that measures the contribution to the mean and rho of the circular distribution of a
#' process or gene set for each of its constituent genes by removing them one by one. For sets
#' consisting on more than 100 genes, it takes a random sample of 100 gene phases to normalize
#' and ensure noticiable changes.
#'
#' @param circa_result Table with circular distribution results, as the output from complete_circular_table or complete_circular_table_grouped
#' @param phase_list Phases list in the same format as the output of gene.list.to.phases, containing the row names of circa_result as names of each element
#'
#' @returns A list whose elements are tables indicating the difference in circular mean and rho due to removing each gene for each process or gene set, and the rank of each difference
#' @export
#'
#' @examples
gene_contribution_to_set <- function(circa_result, phase_list)
{
  phase_list <- phase_list[rownames(circa_result)]
  res_lot <- list()

  for (x in names(phase_list))
  {
    new_list <- lapply(1:nrow(phase_list[[x]]), function(y) phase_list[[x]][-y,])
    new_phases <- lapply(new_list, function(y) circular::circular(y$phase*pi/12))
    new_phases <- lapply(new_phases, function(y) if(length(y) > 100) y[sample(1:length(y), size = 100, replace = F)] else y)
    new_summaries <- t(sapply(new_phases, function(y) summary(y)))[,c("Mean", "Rho")]
    rownames(new_summaries) <- phase_list[[x]]$names
    new_summaries[,"Mean"] <- new_summaries[,"Mean"]*12/pi
    pos_res <- ifelse(as.numeric(new_summaries[,"Mean"]) < 0, 24 + as.numeric(new_summaries[,"Mean"]) ,as.numeric(new_summaries[,"Mean"]))
    new_summaries[,"Mean"] <- pos_res - circa_result[x,"mean"]
    new_summaries[,"Rho"] <- new_summaries[,"Rho"] - circa_result[x,"rho"]
    # Rank impact
    new_summaries <- cbind(new_summaries, Mean_rank = rank(-abs(new_summaries[,"Mean"])))
    new_summaries <- cbind(new_summaries, Rho_rank = rank(-abs(new_summaries[,"Rho"])))

    res_lot[[x]] <- new_summaries
  }

  res_lot <- lapply(res_lot, function(x) x[order(x[,"Mean_rank"]),])
  res_lot <- lapply(res_lot, function(x) { colnames(x) <- c("Mean_diff", "Rho_diff", "Mean_rank", "Rho_rank"); x })

  return(res_lot)

}


#' Multimodality analysis
#'
#' Function to determine the most likely number of modes in the presence of f-fold symmetry (1, 2 or 3),
#' cluster genes in f clusters and determine the statistic summary for each cluster
#'
#' @param circa_table A phase table, as the individual elements of the list derived from gene.list.to.phases
#'
#' @returns A list containing the most likely mode number, the p-value for that number, the p-values corresponding to each mode number, the circular summary table for each cluster and the genes contained in each cluster
#' @export
#'
#' @examples
multimodal_analysis <- function(circa_table)
{
  circa_radians <- circular::circular(circa_table$phase*pi/12)
  modes_summary <- summary(circa_radians)
  modes_ray <- circular::rayleigh.test(circa_radians, mu=circular::circular(modes_summary["Mean"]))$p.value

  circa_radians_2 <- circular::circular(((circa_table$phase*2) %% 24)*pi/12)
  modes_summary_2 <- summary(circa_radians_2)
  modes_ray_2 <- circular::rayleigh.test(circa_radians_2, mu=circular::circular(modes_summary_2["Mean"]))$p.value

  circa_radians_3 <- circular::circular(((circa_table$phase*3) %% 24)*pi/12)
  modes_summary_3 <- summary(circa_radians_3)
  modes_ray_3 <- circular::rayleigh.test(circa_radians_3, mu=circular::circular(modes_summary_3["Mean"]))$p.value

  expected_modes <- which.min(c(modes_ray, modes_ray_2, modes_ray_3))
  expected_p_value <- min(c(modes_ray, modes_ray_2, modes_ray_3))

  result_FOCC <- OptCirClust::CirClust(circa_table$phase, expected_modes, 24, method = "FOCC")
  clusters_focc <- lapply(1:max(result_FOCC$cluster), function(x) circa_table[result_FOCC$cluster == x,])
  names(clusters_focc) <- paste0("cluster_", 1:max(result_FOCC$cluster))

  cluster_radians <- lapply(clusters_focc, function(x) circular::circular(x$phase*pi/12))
  cluster_summary <- sapply(cluster_radians, function(x) summary(x))

  cluster_table <- data.frame(n=cluster_summary["n",], first=cluster_summary["1st Qu.",]*12/pi,
                              median=cluster_summary["Median",]*12/pi, mean=cluster_summary["Mean",]*12/pi,
                              third=cluster_summary["3rd Qu.",]*12/pi, rho=cluster_summary["Rho",],
                              var=1-cluster_summary["Rho",])
  pos_cluster <- apply(cluster_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
  cluster_table[,2:5] <- pos_cluster

  cluster_list <- list(expected_modes = expected_modes, expected_p_value = expected_p_value,
                       total_p_values = c(modes_ray, modes_ray_2, modes_ray_3), summary = cluster_table)
  cluster_list_res <- append(cluster_list, clusters_focc)

  return(cluster_list_res)
}


#' Circular boxplot
#'
#' Function to plot a circular boxplot showing the circular distribution of
#' different processes or gene sets
#'
#' @param res.table A circular results table, as the output from complete_circular_table or complete_circular_table_grouped
#' @param filename PNG file name to export the plot
#' @param color.palette Color palette to use. Available palettes are those from MetBrewer
#' @param tr.height Height of each track. For more than 8 tracks, reduce this parameter
#' @param width PNG file width in px
#' @param height PNG file height in px
#' @param res PNG file resolution
#'
#' @returns A PNG file containing the plot
#' @export
#'
#' @examples
circular_boxplot <- function(res.table, filename, color.palette = "Tam",
                             tr.height = 0.1, width = 1600, height = 1600, res = 180)
{

  grDevices::png(filename, width = width, height = height, res = res)

  sem <- sqrt(res.table$var)/sqrt(res.table$n)
  quan.table <- data.frame(go=rownames(res.table), res.table, min=res.table$first-sem, max=res.table$third+sem)

  # Create axis and initialize plot
  circlize::circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0), gap.degree=0, track.margin =c(0.005, 0.005))
  circlize::circos.initialize("a", xlim = c(0, 24)) # a means that there is a sector, which will go from 0 to 24

  circlize::circos.track(ylim = c(0, 1.5), track.height = 0.001,
               bg.col = NA, bg.border=NA, panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circlize::circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks),
                             labels.cex = 1, lwd = 2)
               })


  my_color=MetBrewer::met.brewer(color.palette, n = nrow(quan.table))

  # Depending on the data, it may happen that the boxplot cuts the sunrise in its
  # lower segment, in the box itself, in its upper segment or not at all,
  # affecting the graph to be plotted. The for loop runs through all the GO
  # terms analyzed (if there are more than 8, the dimensions and the palette
  # used must be changed).


  for (i in 1:nrow(quan.table))
  {
    {
      # If the box does not cut the dawn time
      if ((quan.table$first[i] < quan.table$median[i]) && (quan.table$median[i] < quan.table$third[i]))
      {
        # If neither of the segments cut the sunrise, there is nothing to unfold
        if (quan.table$min[i] < quan.table$max[i])
        {
          circlize::circos.track(ylim = c(0, 1.5), track.height = tr.height,
                       bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                         xlim = circlize::CELL_META$xlim
                         # Estos primeros segmentos son para marcar los ZT
                         circlize::circos.segments(seq(0,24,4), 0,
                                         seq(0,24,4), 1.5,
                                         col = "gray", lwd = 3, lty=3)

                         circlize::circos.segments(quan.table[i,]$min, 0.75,
                                         quan.table[i,]$max, 0.75,
                                         col = my_color[i], lwd =1.5)

                         circlize::circos.rect(quan.table[i,]$first, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                     col = my_color[i], border = my_color[i], lwd=1.5)

                         circlize::circos.segments(c(quan.table[i,]$min, quan.table[i,]$max),
                                         c(0.75-0.2, 0.75-0.2),
                                         c(quan.table[i,]$min, quan.table[i,]$max),
                                         c(0.75+0.2, 0.75+0.2),
                                         col = my_color[i], lwd = 2)
                         circlize::circos.segments(quan.table[i,]$median, 0.75-0.4,
                                         quan.table[i,]$median, 0.75+0.4,
                                         col = "black", lwd = 2)

                       })
        }
        # If the lower or upper segment cut the dawn time, unfold the corresponding one
        else if(quan.table$min[i] > quan.table$first[i] | quan.table$max[i] < quan.table$third[i])
        {

          circlize::circos.track(ylim = c(0, 1.5), track.height = tr.height,
                       bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                         xlim = circlize::CELL_META$xlim

                         circlize::circos.segments(seq(0,24,4), 0,
                                         seq(0,24,4), 1.5,
                                         col = "gray", lwd = 3, lty=3)

                         circlize::circos.segments(quan.table[i,]$min, 0.75,
                                         24, 0.75,
                                         col = my_color[i], lwd =1.5)
                         circlize::circos.segments(0, 0.75,
                                         quan.table[i,]$max, 0.75,
                                         col = my_color[i], lwd =1.5)

                         circlize::circos.rect(quan.table[i,]$first, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                     col = my_color[i], border = my_color[i], lwd=1.5)

                         circlize::circos.segments(c(quan.table[i,]$min, quan.table[i,]$max),
                                         c(0.75-0.2, 0.75-0.2),
                                         c(quan.table[i,]$min, quan.table[i,]$max),
                                         c(0.75+0.2, 0.75+0.2),
                                         col = my_color[i], lwd = 2)
                         circlize::circos.segments(quan.table[i,]$median, 0.75-0.4,
                                         quan.table[i,]$median, 0.75+0.4,
                                         col = "black", lwd = 2)

                       })
        }

      }

      # If the box cuts the dawn time
      else
      {
        # In this case, the segment must be cut from min to 24 and from 0 to max.
        # Also, the box must be unfolded from first to 24 and from 0 to third
        circlize::circos.track(ylim = c(0, 1.5), track.height = tr.height,
                     bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                       xlim = circlize::CELL_META$xlim

                       circlize::circos.segments(seq(0,24,4), 0,
                                       seq(0,24,4), 1.5,
                                       col = "gray", lwd = 3, lty=3)

                       circlize::circos.segments(quan.table[i,]$min, 0.75,
                                       24, 0.75,
                                       col = my_color[i], lwd =1.5)
                       circlize::circos.segments(0, 0.75,
                                       quan.table[i,]$max, 0.75,
                                       col = my_color[i], lwd =1.5)

                       circlize::circos.rect(quan.table[i,]$first, 0.75 - 0.4, 24 , 0.75 + 0.4,
                                   col = my_color[i], border = my_color[i], lwd=1.5)
                       circlize::circos.rect(0, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                   col = my_color[i], border = my_color[i], lwd=1.5)

                       circlize::circos.segments(c(quan.table[i,]$min, quan.table[i,]$max),
                                       c(0.75-0.2, 0.75-0.2),
                                       c(quan.table[i,]$min, quan.table[i,]$max),
                                       c(0.75+0.2, 0.75+0.2),
                                       col = my_color[i], lwd = 2)
                       circlize::circos.segments(quan.table[i,]$median, 0.75-0.4,
                                       quan.table[i,]$median, 0.75+0.4,
                                       col = "black", lwd = 2)


                     })


      }
    }
  }

  # Legend, establish values, lines, colors and the adjust position in the plot
  my_legend <- ComplexHeatmap::Legend(at = quan.table$go,
                      legend_gp = grid::gpar(fill = my_color), title_position = "topleft",
                      title = "")
  ComplexHeatmap::draw(my_legend, x = grid::unit(1, "npc") - grid::unit(2, "mm"), y = grid::unit(10, "mm"),
       just = c("right", "bottom"))

  circlize::circos.clear()
  grDevices::dev.off()
}


#' Circular dotplot
#'
#' Function to plot a circular dotplot showing the circular distribution of
#' different processes or gene sets
#'
#' @param phases.list A named phase tables list in the same format as the output of gene.list.to.phases
#' @param filename PNG file name to export the plot
#' @param color.palette Color palette to use. Available palettes are those from MetBrewer
#' @param tr.height Height of each track. For more than 8 tracks, reduce this parameter
#' @param width PNG file width in px
#' @param height PNG file height in px
#' @param res PNG file resolution
#'
#' @returns A PNG file containing the plot
#' @export
#'
#' @examples
circular_dotplot <- function(phases.list, filename, color.palette = "Tam",
                             tr.height = 0.1, width = 1600, height = 1600, res = 180)
{

  grDevices::png(filename, width = width, height = height, res = res)

  # Create axis and initialize plot
  circlize::circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0), gap.degree=0, track.margin =c(0.005, 0.005))
  circlize::circos.initialize("a", xlim = c(0, 24)) # a means that there is a sector, which will go from 0 to 24

  circlize::circos.track(ylim = c(0, 1.5), track.height = 0.001,
               bg.col = NA, bg.border=NA, panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circlize::circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks),
                             labels.cex = 1, lwd = 2)
               })


  my_color=MetBrewer::met.brewer(color.palette, n = length(phases.list))
  my_transparent = sapply(my_color, function(x) circlize::add_transparency(x, 0.5))


  # Depending on the data, it may happen that the boxplot cuts the sunrise in its
  # lower segment, in the box itself, in its upper segment or not at all,
  # affecting the graph to be plotted. The for loop runs through all the GO
  # terms analyzed (if there are more than 8, the dimensions and the palette
  # used must be changed).


  for (i in 1:length(phases.list))
  {

    circlize::circos.track(ylim = c(0, 1.5), track.height = tr.height,
                 bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                   xlim = circlize::CELL_META$xlim
                   # Estos primeros segmentos son para marcar los ZT
                   circlize::circos.segments(seq(0,24,4), 0,
                                   seq(0,24,4), 1.5,
                                   col = "gray", lwd = 3, lty=3)

                   circlize::circos.points(phases.list[[i]]$phase, 0.75,
                                 col = my_color[i], pch=21, bg=my_transparent[i])

                 })

  }

  # Legend, establish values, lines, colors and the adjust position in the plot
  my_legend <- ComplexHeatmap::Legend(at = names(phases.list),
                      legend_gp = grid::gpar(fill = my_color), title_position = "topleft",
                      title = "")
  ComplexHeatmap::draw(my_legend, x = grid::unit(1, "npc") - grid::unit(2, "mm"), y = grid::unit(10, "mm"),
       just = c("right", "bottom"))

  circlize::circos.clear()
  grDevices::dev.off()
}



#' Circular histogram
#'
#' Function to plot a circular histogram showing the circular distribution of
#' different processes or gene sets
#'
#' @param phases.list A phase tables list in the same format as the output of gene.list.to.phases
#' @param filename PNG file name to export the plot
#' @param color.palette Color palette to use. Available palettes are those from MetBrewer
#' @param tr.height Height of each track. For more than 8 tracks, reduce this parameter
#' @param nbins Number of bins (breaks) to split the histogram
#' @param width PNG file width in px
#' @param height PNG file height in px
#' @param res PNG file resolution
#'
#' @returns A PNG file containing the plot
#' @export
#'
#' @examples
circular_histogram <- function(phases.list, filename,
                               color.palette = "Tam", tr.height = 0.08,
                               nbins = 24, width = 1600, height = 1600, res = 180)
{

  # Preparar fichero
  grDevices::png(filename, width = width, height = height, res = res)

  circlize::circos.par(start.degree = 90,
             cell.padding = c(0, 0, 0, 0),
             gap.degree = 0,
             track.margin = c(0.005, 0.005))
  circlize::circos.initialize(factors = "a", xlim = c(0, 24))

  circlize::circos.track(ylim = c(0, 1.5), track.height = 0.001, bg.col = NA, bg.border = NA,
               panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circlize::circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks),
                             labels.cex = 1, lwd = 2)
               })

  my_color <- MetBrewer::met.brewer(color.palette, n = length(phases.list))
  my_transparent <- sapply(my_color, function(x) circlize::add_transparency(x, 0.55))

  breaks_common <- seq(0, 24, length.out = nbins + 1)

  counts_list <- lapply(phases.list, function(x) {
    h <- hist(x$phase, breaks = breaks_common, plot = FALSE, right = FALSE)
    h$counts
  })

  for (i in seq_along(counts_list)) {
    counts <- counts_list[[i]]
    circlize::circos.track(ylim = c(0, max(counts_list[[i]])), track.height = tr.height,
                 bg.col = "#eaeded", bg.border = NA,
                 panel.fun = function(x, y) {
                   for (b in seq_along(counts)) {
                     xleft <- breaks_common[b]
                     xright <- breaks_common[b + 1]
                     ybottom <- 0
                     ytop <- counts[b]
                     circlize::circos.rect(xleft, ybottom, xright, ytop,
                                 sector.index = circlize::CELL_META$sector.index,
                                 col = my_color[i],
                                 border = my_color[i],
                                 lwd = 0.3)
                   }
                   circlize::circos.segments(seq(0,24,4), 0, seq(0,24,4), max(counts_list[[i]]),
                                   col = "gray", lwd = 1, lty = 3)
                 })
  }

  my_legend <- ComplexHeatmap::Legend(at = names(phases.list),
                      legend_gp = grid::gpar(fill = my_color),
                      title = "")
  ComplexHeatmap::draw(my_legend, x = grid::unit(1, "npc") - grid::unit(2, "mm"),
       y = grid::unit(10, "mm"), just = c("right", "bottom"))

  circlize::circos.clear()
  grDevices::dev.off()

}

#' plants
#'
#' KEGG pathways associated to plants
#'
#' @format Data frame with phase information in two different conditions for Marchantia polymorpha
"circa_table_genescarab"
