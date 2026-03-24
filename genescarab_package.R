# Functions for tests

#install.packages("Directional")
#install.packages("CircMLE")


# HR for grouped data, with TB as implemented in 
# "Grouped circular data in biology: advice for effectively implementing 
# statistical procedures"
# Critical value calculation
HermansRasson2T <- function(sample){
  n <- length(sample)
  total <- 0
  for (i in 1:n){
    for (j in 1:n){ total <- total + abs(abs(sample[i]-sample[j])-pi)-
      (pi/2)
    total <- total - (2.895*(abs(sin(sample[i]-sample[j]))-(2/pi)))}}
  T <- total/n
  return(T)}

# Function
HermansRasson2PGroupedRad <- function(sample, m, k=1000, iter=9999){
  sample<-circular(sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  n <- length(sample)
  univals <- iter
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- rcircularuniform(n, control.circular=list(units="radians"))
    data1 <- trunc(data1*m/(2*pi))
    data1 <- data1*2*pi/m
    errorsamp <- rvonmises(n, 0, k,control.circular=list(units="radians"))
    data1 <- data1+errorsamp
    data1 <- ifelse((data1>(2*pi)),(data1-(2*pi)), data1)
    testset[f] <- HermansRasson2T(data1)}
  errorsamp2 <- rvonmises(n, 0, k,control.circular=list(units="radians"))
  sample<-sample+errorsamp2
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  Tsample <- HermansRasson2T(sample)
  counter <- 0
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}

# Rao spacing test for continuous data
RaoTestUngroupedRad <- function(sample, iter=9999){
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  n <- length(sample)
  univals <- iter
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- rcircularuniform(n, control.circular=list(units="radians"))
    testset[f] <- RaoTestValue(data1)}
  Tsample <- RaoTestValue(sample)
  counter <- 0
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}

# Rao spacing test with TB for grouped data
# Critical function calculation
RaoTestValue <- function(sample){n <- length(sample)
f <- sort(sample)
fplus <- c(f[2:n],f[1])
T <- fplus - f
T[n] <- (2*pi) - f[n] + f[1]
abs_diff <- abs(T-(2*pi/n))
U <- 0.5*sum(abs_diff)
return(U)}

# Function
RaoPGroupedRad <- function(sample, m, k=1000, iter=9999){
  sample<-circular(sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  sample<- ifelse((sample<(0)),(sample+(2*pi)), sample)
  n <- length(sample)
  univals <- iter
  testset<- rep(0,univals)
  for (f in 1:univals){
    data1 <- rcircularuniform(n, control.circular=list(units="radians"))
    data1 <- trunc(data1*m/(2*pi))
    data1 <- data1*2*pi/m
    errorsamp <- rvonmises(n, 0, k,control.circular=list(units="radians"))
    data1 <- data1+errorsamp
    data1 <- ifelse((data1>(2*pi)),(data1-(2*pi)), data1)
    testset[f] <- RaoTestValue(data1)}
  errorsamp2 <- rvonmises(n, 0, k,control.circular=list(units="radians"))
  sample<-sample+errorsamp2
  sample<- ifelse((sample>(2*pi)),(sample-(2*pi)), sample)
  Tsample <- RaoTestValue(sample)
  counter <- 0
  for(j in 1:univals){if(testset[j]>=Tsample){counter <- counter+1}}
  p <- counter/(univals+1)
  return(p)}


# Function to create input lists for GO terms
create_gene_list_go <- function(go_vector = "all", org.package, go_column, id_column)
{
  if (go_vector == "all")
    {
    functional_data <- select(get(org.package), 
                                keys = keys(get(org.package),keytype=id_column), 
                                columns =  c(id_column, go_column))
    
    go_vector <- unique(functional_data[[go_column]])
    go_vector <- go_vector[!is.na(go_vector)] 
    }
  
  
  # Get ancestors of selected GOs for BP, MF and CC 
  bp_ancestors <- BiocGenerics::mget(go_vector, GOBPANCESTOR, ifnotfound=NA)
  mf_ancestors <- BiocGenerics::mget(go_vector, GOMFANCESTOR, ifnotfound=NA)
  cc_ancestors <- BiocGenerics::mget(go_vector, GOCCANCESTOR, ifnotfound=NA)
  
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
  bp_offspring <- BiocGenerics::mget(query_vec_ancestors, GOBPOFFSPRING, ifnotfound=NA)
  mf_offspring <- BiocGenerics::mget(query_vec_ancestors, GOMFOFFSPRING, ifnotfound=NA)
  cc_offspring <- BiocGenerics::mget(query_vec_ancestors, GOCCOFFSPRING, ifnotfound=NA)
  
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


# Function to create input lists for KEGG pathways
create_gene_list_kegg <- function(ko_vector = "all", org.package, ko_column, id_column, 
                                  ko_prefix = "map", species = "all")
{
  if (ko_vector == "all")
  {
    functional_data <- select(get(org.package), 
                              keys = keys(get(org.package),keytype=id_column), 
                              columns =  c(id_column, ko_column))
    
    complete_kos <- unique(functional_data[[ko_column]])
    ko_vector <- complete_kos[!is.na(complete_kos)] 
  }
  
  # Get pathways for selected KOs 
  {
  if (length(ko_vector) < 500)
  {
    present_pathways <- keggLink("pathway", paste0(ko_vector, collapse = "+"))
  }
  else
  {
    present_pathways <- keggLink("pathway", paste0(ko_vector[1:500], collapse = "+"))
    
    for (i in 2:(length(ko_vector)%/%500 + 1)) 
      {
      present_pathways <- append(present_pathways, keggLink("pathway", paste0(ko_vector[(1+(500*(i-1))):(i*500)], collapse = "+")))
      }
  }
  }
  
  unique_pathways <- unique(unlist(present_pathways))
  
  map_pathways <- unique_pathways[grep(ko_prefix, unique_pathways)] 
  list_pathways <- sapply(strsplit(map_pathways, ":"), function(x) x[[2]])
  
  
  {
    if (length(list_pathways) < 500)
    {
      kos_pathways <- keggLink("ko", paste0(list_pathways, collapse = "+"))
    }
    else
    {
      kos_pathways <- keggLink("ko", paste0(list_pathways[1:500], collapse = "+"))
      
      for (i in 2:(length(list_pathways)%/%500 + 1)) 
      {
        print(i)
        kos_pathways <- append(kos_pathways, keggLink("ko", paste0(list_pathways[(1+(500*(i-1))):(i*500)], collapse = "+")))
      }
    }
  }
  
  kos_pathways <- tapply(kos_pathways, INDEX = names(kos_pathways), FUN = function(x) unlist(x), simplify = F)
  names(kos_pathways) <- sapply(strsplit(names(kos_pathways), split = "path:"), function(x) x[[2]])
  
  kos_pathways_clean <- lapply(kos_pathways, function(x) unique(x))
  kos_pathways_clean <- lapply(kos_pathways_clean, function(x) sapply(strsplit(x, ":"), function(y) y[[2]]))
  ko_list <- lapply(kos_pathways_clean, function(x) subset(functional_data, functional_data[[ko_column]] %in% x)[[id_column]])
  ko_list <- lapply(ko_list, unique)
  
  # TODO filter paths by species
  {
    if (species == "all")
    {
      return(ko_list)
    }
    else
    {
      kegg_to_keep <- read.csv(paste0("pathways_", species, ".txt"), sep = "\t", header = T)$x
      kegg_to_keep <- paste0(ko_prefix, sapply(strsplit(kegg_to_keep, "ko"), function(x) x[[2]]))
      ko_list <- ko_list[names(ko_list) %in% kegg_to_keep]
      return(ko_list)
    }
  }
  
}


# Statistics and general workflow functions
# Function to create a phase table list from gene lists
gene.list.to.phases <- function(gene.list, phase.table)
{
  circa_reduced <- lapply(gene.list, function(x) subset(phase.table, names %in% x))
  circa_reduced_clean <- circa_reduced[which(sapply(circa_reduced, nrow) != 0)]
  return(circa_reduced_clean)
  
}

# Function to create the circular table and deviation from the uniform distribution for a single set
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
  circa_ray <- circular::rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
  
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

# Function to create the circular table and deviation from the uniform distribution for a single set of grouped data
create_circular_table_grouped <- function(circa_vector, n_bins,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                  iter.hr = 999, force.hr=F, force.hr.th=0.05, 
                                  rao.on.large.sets = F, rao.on.large.sets.th=400,
                                  iter.rao = 999, force.rao=F, force.rao.th=0.05, error_kappa=1000)
{
  
  circa_radians <- circular::circular(circa_vector*pi/12)
  
  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
  
  # Calculate circular variance and transform to hours again
  circa_var <- var.circular(circa_radians)
  
  # Apply kupier's test
  circa_kupier <- Directional::kuiper(circa_radians, rads = T, R=1)$p.value
  
  # Apply Rayleigh's test
  circa_ray <- circular::rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
  
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

# Function to create the circular table and deviation from the uniform distribution for a phases list
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

# Function to create the circular table and deviation from the uniform distribution for a grouped phases list
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

# Function to compute the deviation from the total experimental distribution for a phases list
test_against_gen_dist <- function(phase.list, total.phase.table)
{
  gen_dist <- total.phase.table$phase*pi/12
  len_gen_dist <- length(gen_dist)
  
  p_values_gen_test <- sapply(phase.list, function(x) summary(manova(cbind(sin(c(gen_dist, x$phase*pi/12)),cos(c(gen_dist, x$phase*pi/12))) 
                                                                     ~ c(rep("general", len_gen_dist), rep("specific", nrow(x)))))$stats[1,6])
  gen_dist_bh <- p.adjust(p_values_gen_test, method = "BH")
  
  gen_dist_table <- data.frame(gen_dist_p_value=p_values_gen_test, 
                               gen_dist_p_value_adj=gen_dist_bh)
  
  return(gen_dist_table)
}

# Function for comparing the distributions of to different phase lists
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
  
  p_values_diff_test <- mapply(function(x, y) summary(manova(cbind(sin(c(x$phase*pi/12, y$phase*pi/12)),cos(c(x$phase*pi/12, y$phase*pi/12))) 
                                                             ~ c(rep("first", nrow(x)), rep("second", nrow(y)))))$stats[1,6],
                               phase.list.1, phase.list.2)
  
  p_values_diff_bh <- p.adjust(unlist(p_values_diff_test), method = "BH")
  
  p_values_diff_table <- data.frame(p_value=unlist(p_values_diff_test), 
                                    p_value_adj=p_values_diff_bh)
  
  return(p_values_diff_table)
}

# Function to compute the contribution of each gene in a set to the temporal cohesion of the whole set
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

# Function to determine the most likely number of modes in the pressence of f-fold symmetry (2 or 3),
# cluster genes in f clusters and determine the circular tables for each cluster
multimodal_analysis <- function(circa_table)
{
  circa_radians <- circular::circular(circa_table$phase*pi/12)
  modes_summary <- summary(circa_radians)
  modes_ray <- circular::rayleigh.test(circa_radians, mu=circular(modes_summary["Mean"]))$p.value
  
  circa_radians_2 <- circular::circular(((circa_table$phase*2) %% 24)*pi/12)
  modes_summary_2 <- summary(circa_radians_2)
  modes_ray_2 <- circular::rayleigh.test(circa_radians_2, mu=circular(modes_summary_2["Mean"]))$p.value
  
  circa_radians_3 <- circular::circular(((circa_table$phase*3) %% 24)*pi/12)
  modes_summary_3 <- summary(circa_radians_3)
  modes_ray_3 <- circular::rayleigh.test(circa_radians_3, mu=circular(modes_summary_3["Mean"]))$p.value
  
  expected_modes <- which.min(c(modes_ray, modes_ray_2, modes_ray_3))
  expected_p_value <- min(c(modes_ray, modes_ray_2, modes_ray_3))
  
  result_FOCC <- CirClust::CirClust(circa_table$phase, expected_modes, 24, method = "FOCC")
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

# Functions for plotting results
# Circular boxplot function
circular_boxplot <- function(res.table, filename, color.palette = "Tam",
                             tr.height = 0.1, width = 1600, height = 1600, res = 180)
{
  
  png(filename, width = width, height = height, res = res)
  
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
                         xlim = CELL_META$xlim
                         # Estos primeros segmentos son para marcar los ZT
                         circlize::circos.segments(seq(0,24,4), 0,
                                         seq(0,24,4), 1.5,
                                         col = "gray", lwd = 3, lty=3)
                         
                         circos.segments(quan.table[i,]$min, 0.75,
                                         quan.table[i,]$max, 0.75,
                                         col = my_color[i], lwd =1.5)
                         
                         circos.rect(quan.table[i,]$first, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                     col = my_color[i], border = my_color[i], lwd=1.5)
                         
                         circos.segments(c(quan.table[i,]$min, quan.table[i,]$max), 
                                         c(0.75-0.2, 0.75-0.2),
                                         c(quan.table[i,]$min, quan.table[i,]$max), 
                                         c(0.75+0.2, 0.75+0.2),
                                         col = my_color[i], lwd = 2)
                         circos.segments(quan.table[i,]$median, 0.75-0.4,
                                         quan.table[i,]$median, 0.75+0.4,
                                         col = "black", lwd = 2)
                         
                       })
        }
        # If the lower or upper segment cut the dawn time, unfold the corresponding one
        else if(quan.table$min[i] > quan.table$first[i] | quan.table$max[i] < quan.table$third[i]) 
        {
          
          circos.track(ylim = c(0, 1.5), track.height = tr.height, 
                       bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                         xlim = CELL_META$xlim
                         
                         circos.segments(seq(0,24,4), 0,
                                         seq(0,24,4), 1.5,
                                         col = "gray", lwd = 3, lty=3)
                         
                         circos.segments(quan.table[i,]$min, 0.75,
                                         24, 0.75,
                                         col = my_color[i], lwd =1.5)
                         circos.segments(0, 0.75,
                                         quan.table[i,]$max, 0.75,
                                         col = my_color[i], lwd =1.5)
                         
                         circos.rect(quan.table[i,]$first, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                     col = my_color[i], border = my_color[i], lwd=1.5)
                         
                         circos.segments(c(quan.table[i,]$min, quan.table[i,]$max), 
                                         c(0.75-0.2, 0.75-0.2),
                                         c(quan.table[i,]$min, quan.table[i,]$max), 
                                         c(0.75+0.2, 0.75+0.2),
                                         col = my_color[i], lwd = 2)
                         circos.segments(quan.table[i,]$median, 0.75-0.4,
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
        circos.track(ylim = c(0, 1.5), track.height = tr.height, 
                     bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                       xlim = CELL_META$xlim
                       
                       circos.segments(seq(0,24,4), 0,
                                       seq(0,24,4), 1.5,
                                       col = "gray", lwd = 3, lty=3)
                       
                       circos.segments(quan.table[i,]$min, 0.75,
                                       24, 0.75,
                                       col = my_color[i], lwd =1.5)
                       circos.segments(0, 0.75,
                                       quan.table[i,]$max, 0.75,
                                       col = my_color[i], lwd =1.5)
                       
                       circos.rect(quan.table[i,]$first, 0.75 - 0.4, 24 , 0.75 + 0.4,
                                   col = my_color[i], border = my_color[i], lwd=1.5)
                       circos.rect(0, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                   col = my_color[i], border = my_color[i], lwd=1.5)
                       
                       circos.segments(c(quan.table[i,]$min, quan.table[i,]$max), 
                                       c(0.75-0.2, 0.75-0.2),
                                       c(quan.table[i,]$min, quan.table[i,]$max), 
                                       c(0.75+0.2, 0.75+0.2),
                                       col = my_color[i], lwd = 2)
                       circos.segments(quan.table[i,]$median, 0.75-0.4,
                                       quan.table[i,]$median, 0.75+0.4,
                                       col = "black", lwd = 2)
                       
                       
                     })
        
        
      }
    }
  }
  
  # Legend, establish values, lines, colors and the adjust position in the plot
  my_legend <- Legend(at = quan.table$go,
                      legend_gp = gpar(fill = my_color), title_position = "topleft", 
                      title = "")
  draw(my_legend, x = unit(1, "npc") - unit(2, "mm"), y = unit(10, "mm"), 
       just = c("right", "bottom"))
  
  circos.clear()
  dev.off()
}

# Circular dotplot function
circular_dotplot <- function(phases.list, filename, color.palette = "Tam",
                             tr.height = 0.1, width = 1600, height = 1600, res = 180)
{
  
  png(filename, width = width, height = height, res = res)
  
  # Create axis and initialize plot
  circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0), gap.degree=0, track.margin =c(0.005, 0.005))
  circos.initialize("a", xlim = c(0, 24)) # a means that there is a sector, which will go from 0 to 24
  
  circos.track(ylim = c(0, 1.5), track.height = 0.001, 
               bg.col = NA, bg.border=NA, panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks), 
                             labels.cex = 1, lwd = 2)
               })
  
  
  my_color=met.brewer(color.palette, n = length(phases.list))
  #my_colors_rgb=lapply(my_color, function(x) col2rgb(x))
  #my_transparent = sapply(my_colors_rgb, function(x) rgb(red = x[1,1], green = x[2,1], blue = x[3,1], alpha = 0.9, maxColorValue = 255))
  my_transparent = sapply(my_color, function(x) add_transparency(x, 0.5))
  
  
  # Depending on the data, it may happen that the boxplot cuts the sunrise in its
  # lower segment, in the box itself, in its upper segment or not at all, 
  # affecting the graph to be plotted. The for loop runs through all the GO 
  # terms analyzed (if there are more than 8, the dimensions and the palette 
  # used must be changed).
  
  
  for (i in 1:length(phases.list))
  {
    
    circos.track(ylim = c(0, 1.5), track.height = tr.height, 
                 bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                   xlim = CELL_META$xlim
                   # Estos primeros segmentos son para marcar los ZT
                   circos.segments(seq(0,24,4), 0,
                                   seq(0,24,4), 1.5,
                                   col = "gray", lwd = 3, lty=3)
                   
                   circos.points(phases.list[[i]]$phase, 0.75,
                                 col = my_color[i], pch=21, bg=my_transparent[i])
                   
                 })
    
  }
  
  # Legend, establish values, lines, colors and the adjust position in the plot
  my_legend <- Legend(at = names(phases.list),
                      legend_gp = gpar(fill = my_color), title_position = "topleft", 
                      title = "")
  draw(my_legend, x = unit(1, "npc") - unit(2, "mm"), y = unit(10, "mm"), 
       just = c("right", "bottom"))
  
  circos.clear()
  dev.off()
}

# Circular histogram function
circular_histogram <- function(phases.list, filename,
                               color.palette = "Tam", tr.height = 0.08,
                               nbins = 24, width = 1600, height = 1600, res = 180)
{
  # paquetes necesarios
  library(circlize)
  library(MetBrewer)
  library(grid)
  
  # helper para transparencia (si no la tienes, usa esta)
  add_transparency <- function(col, alpha = 0.5) {
    rgb_val <- grDevices::col2rgb(col)
    grDevices::rgb(rgb_val[1,]/255, rgb_val[2,]/255, rgb_val[3,]/255, alpha = alpha)
  }
  
  # Preparar fichero
  png(filename, width = width, height = height, res = res)
  
  # inicializar círculo
  circos.clear()
  circos.par(start.degree = 90,
             cell.padding = c(0, 0, 0, 0),
             gap.degree = 0,
             track.margin = c(0.005, 0.005))
  circos.initialize(factors = "a", xlim = c(0, 24))
  
  # Eje con ZT
  circos.track(ylim = c(0, 1.5), track.height = 0.001, bg.col = NA, bg.border = NA,
               panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks),
                             labels.cex = 1, lwd = 2)
               })
  
  # Paleta y transparencias
  my_color <- met.brewer(color.palette, n = length(phases.list))
  my_transparent <- sapply(my_color, function(x) add_transparency(x, 0.55))
  
  # Construir los breaks comunes (desde 0 a 24)
  breaks_common <- seq(0, 24, length.out = nbins + 1)
  
  # Calcular los histograms sin plotear y hallar el máximo para escala vertical
  counts_list <- lapply(phases.list, function(x) {
    # asumimos que phases.list[[i]]$phase es un vector numérico en [0,24)
    h <- hist(x$phase, breaks = breaks_common, plot = FALSE, right = FALSE)
    h$counts
  })
  #max_count <- max(sapply(counts_list, max), 1) # evitar 0
  
  # Dibujar un anillo por cada conjunto (con misma escala vertical)
  for (i in seq_along(counts_list)) {
    counts <- counts_list[[i]]
    circos.track(ylim = c(0, max(counts_list[[i]])), track.height = tr.height,
                 bg.col = "#eaeded", bg.border = NA,
                 panel.fun = function(x, y) {
                   # dentro del panel dibujamos rectángulos para cada bin
                   for (b in seq_along(counts)) {
                     xleft <- breaks_common[b]
                     xright <- breaks_common[b + 1]
                     ybottom <- 0
                     ytop <- counts[b]
                     # dibujar barra (rect) en coordenadas del sector "a"
                     circos.rect(xleft, ybottom, xright, ytop,
                                 sector.index = CELL_META$sector.index,
                                 #col = my_transparent[i],
                                 col = my_color[i],
                                 border = my_color[i],
                                 lwd = 0.3)
                   }
                   # opcional: dibujar líneas divisoras ZT
                   circos.segments(seq(0,24,4), 0, seq(0,24,4), max(counts_list[[i]]),
                                   col = "gray", lwd = 1, lty = 3)
                 })
  }
  
  # Leyenda (mantener nombres y colores)
  my_legend <- Legend(at = names(phases.list),
                      legend_gp = gpar(fill = my_color),
                      title = "")
  draw(my_legend, x = unit(1, "npc") - unit(2, "mm"),
       y = unit(10, "mm"), just = c("right", "bottom"))
  
  circos.clear()
  dev.off()
}
