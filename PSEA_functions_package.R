# Functions for tests

# Kuiper test with bootstrap, as implemented in the package Directional
kuiper <- function(u, rads = FALSE, R = 1) {
  ## u is a vector with circular data
  ## if data are in rads set it to TRUE
  ## R is for Monte Carlo estimate of the p-value
  if ( !rads )   u <- u / 180 * pi
  u <- Rfast::Sort(u) / (2 * pi)
  n <- length(u)
  i <- 1:n
  f <- sqrt(n)
  Vn <- f * ( max(u - (i - 1)/n ) + max(i/n - u) )
  
  if ( R == 1 ) {  ## asymptotic p-value is returned
    m <- (1:50)^2
    a1 <- 4 * m * Vn^2
    a2 <- exp( -2 * m * Vn^2 )
    b1 <- 2 * ( a1 - 1 ) * a2
    b2 <- 8 * Vn / ( 3 * f ) * m * (a1 - 3) * a2
    p.value <- sum(b1 - b2)
    
  } else {
    x <- matrix( Rfast2::Runif(n * R, 0, 2 * pi), ncol = R)
    x <- Rfast::colSort(x) / (2 * pi)
    bvn <- f * ( Rfast::colMaxs(x - (i - 1)/n, value = TRUE) + Rfast::colMaxs(i/n - x, value = TRUE) )
    p.value <- ( sum(bvn > Vn) + 1 ) / (R + 1)
  }
  
  p.value <- ifelse(p.value < 0, 0, p.value)
  parameter <- "NA"     ;   names(parameter) <- "df"
  statistic <- Vn  ;   names(statistic) <- "Test statistic"
  alternative <- "The distribution is not circular uniform"
  method <- "Kuiper test of uniformity with circular data"
  data.name <- c("data")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}

# Hermann-Rasson test with bootstrap, as implemented in the package CircMLE
HR_test <- function(data, original = F, iter = 9999){
  
  #Check input parameters
  if (is.null(data) | missing(data)) stop("Please provide input data vector")
  if (missing(original)) original = F else original = original
  if (!is.logical(original)) stop ("Please set \"original\" to TRUE (T) or FALSE (F)")
  if (missing(iter)) iter = 9999 else iter = iter
  if (length(iter) != 1 | !is.numeric(iter)) stop("Please set the number of iterations correctly for p-value estimation")
  if (iter <= 99) warning("WARNING: You set an awfully small number of iterations for p-value estimation, it is recommended to run at least 999 for accurate estimation.")
  
  
  # Get all circular attributes from data
  params = circularp(data)
  
  # Begin Hermans-Rasson test
  n <- length(data)
  # sample <- circular(data, units = "radians")  # not needed because of check_data() above
  testset<- rep(0, iter)
  for (f in 1:iter){
    data1 <- matrix(rcircularuniform(n, control.circular = params))
    if(original) testset[f] <- HermansRassonT(data1) else testset[f] <- HermansRasson2T(data1)
  }
  
  # Run appropriate method for T
  if(original) Tsample <- HermansRassonT(data) else Tsample <- HermansRasson2T(data)
  counter <- 1
  for(j in 1:iter){
    if(testset[j] > Tsample){
      counter <- counter + 1
    }
  }
  p <- counter/(iter + 1)
  out <- c(as.numeric(Tsample), p)
  names(out) <- c("Test statistic (T)", "p-value")
  return(out)
}

# Original Hermans-Rasson test formula, referred to as HR-infinity by Landler et al. 2019
HermansRassonT <- function(data){
  n <- length(data)
  total <- 0
  for (i in 1:n){
    for (j in 1:n){
      total <- total + abs(sin(data[i] - data[j]))
    }
  }
  T <- abs((n / pi) - (total / (2 * n)))
  return(T)
}

# Updated Hermans-Rasson test formula, referred to as HR by Landler et al. 2019
HermansRasson2T <- function(data){
  n <- length(data)
  total <- 0
  for (i in 1:n){
    for (j in 1:n){
      total <- total + abs(abs(data[i] - data[j]) - pi) - (pi / 2)
      total <- total - (2.895 * (abs(sin(data[i] - data[j])) - (2 / pi)))
    }
  }
  T <- total / n
  return(T)
}

# HR for grouped data, with TB
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
  # sample is the list of data points in radians measure, this can be a simple list 
  # or a list of class circular
  #sample <- rad(circular(sample)) for input in degrees
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



# Functions to create input lists for GO terms and KEGG pathways
create_gene_list_go <- function(go_vector, org.package, go_column, id_column)
{
  functional_data <- select(get(org.package), 
                            keys = keys(get(org.package),keytype=id_column), 
                            columns =  c(id_column, go_column))
  
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

# REVISAR CUÁL ES LA BUENA
create_gene_list_kegg <- function(ko_vector, org.package, ko_column, id_column)
{
  functional_data <- select(get(org.package), 
                            keys = keys(get(org.package),keytype=id_column), 
                            columns =  c(id_column, ko_column))
  
  # Get pathways for selected KOs 
  present_pathways <- sapply(ko_vector, function(x) keggLink("pathway", x))
  unique_pathways <- unique(unlist(present_pathways))
  map_pathways <- unique_pathways[grep("map", unique_pathways)] # maybe change this for ko in other versions
  list_pathways <- sapply(strsplit(map_pathways, ":"), function(x) x[[2]])
  
  
  kos_pathways <- lapply(list_pathways, function(x) keggLink("ko",x))
  names(kos_pathways) <- list_pathways
  kos_pathways_clean <- lapply(kos_pathways, function(x) unique(x))
  kos_pathways_clean <- lapply(kos_pathways_clean, function(x) sapply(strsplit(x, ":"), function(y) y[[2]]))
  ko_list <- lapply(kos_pathways_clean, function(x) subset(functional_data, functional_data[["KO"]] %in% x)[["GID"]])
  ko_list <- lapply(ko_list, unique)
  
  return(ko_list)
  
}
# TODO for more than 500 KO terms

create_gene_list_kegg <- function(ko_vector, org.package, ko_column, id_column)
{
  functional_data <- select(get(org.package), 
                            keys = keys(get(org.package),keytype=id_column), 
                            columns =  c(id_column, ko_column))
  
  # Get pathways for selected KOs 
  present_pathways <- sapply(ko_vector, function(x) keggLink("pathway", x))
  unique_pathways <- unique(unlist(present_pathways))
  map_pathways <- unique_pathways[grep("map", unique_pathways)] # maybe change this for ko in other versions
  list_pathways <- sapply(strsplit(map_pathways, ":"), function(x) x[[2]])
  
  
  kos_pathways <- lapply(list_pathways, function(x) keggLink("ko",x))
  names(kos_pathways) <- list_pathways
  kos_pathways_clean <- lapply(kos_pathways, function(x) unique(x))
  kos_pathways_clean <- lapply(kos_pathways_clean, function(x) sapply(strsplit(x, ":"), function(y) y[[2]]))
  ko_list <- lapply(kos_pathways_clean, function(x) subset(functional_data, grepl(x, functional_data[["KO"]]))[["GID"]])
  ko_list <- lapply(ko_list, unique)
  
  return(ko_list)
  
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
create_circular_table_bk <- function(circa_vector, hr.on.large.sets = F, hr.on.large.sets.th=400,
                                  iter = 999, force.hr=F, force.hr.th=0.05)
{
  
  circa_radians <- circular(circa_vector*pi/12)
  
  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
  
  # Calculate circular variance and transform to hours again
  circa_var <- var.circular(circa_radians)
  
  # Apply kupier's test
  circa_kupier <- kuiper(circa_radians, rads = T, R=1)$p.value
  
  # Apply Rayleigh's test
  circa_ray <- rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
  
  # Apply HR test
  if (iter == 0)
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
    circa_hr <- HR_test(circa_radians, original = F, iter = iter)
  }
  
  
  # Create table
  circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi, 
                            median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
                            third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"], 
                            var=circa_var, kuiper_p_value=circa_kupier, rayleigh_p_value=circa_ray, 
                            hr_p_value = circa_hr[2] )
  
  
  return(circa_table)
  
}

create_circular_table <- function(circa_vector,hr.on.large.sets = F, hr.on.large.sets.th=400,
                                  iter.hr = 999, force.hr=F, force.hr.th=0.05, 
                                  rao.on.large.sets = F, rao.on.large.sets.th=400,
                                  iter.rao = 999, force.rao=F, force.rao.th=0.05)
{
  
  circa_radians <- circular(circa_vector*pi/12)
  
  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
  
  # Calculate circular variance and transform to hours again
  circa_var <- var.circular(circa_radians)
  
  # Apply kupier's test
  circa_kupier <- kuiper(circa_radians, rads = T, R=1)$p.value
  
  # Apply Rayleigh's test
  circa_ray <- rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
  
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
    circa_hr <- HR_test(circa_radians, original = F, iter = iter.hr)
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
  
  circa_radians <- circular(circa_vector*pi/12)
  
  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
  
  # Calculate circular variance and transform to hours again
  circa_var <- var.circular(circa_radians)
  
  # Apply kupier's test
  circa_kupier <- kuiper(circa_radians, rads = T, R=1)$p.value
  
  # Apply Rayleigh's test
  circa_ray <- rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
  
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


# Fuction for comparing the distributions of to different phase lists
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
  
  p_values_diff_table <- data.frame(gen_dist_p_value=unlist(p_values_diff_test), 
                                    gen_dist_p_value_adj=p_values_diff_bh)
  
  return(p_values_diff_table)
}

# Function to compute the contribution of each gene in a set to the temporal cohesion of the whole set
gene_contribution_to_set <- function(circa_result, phase_list)
{
  new_phase_list <- phase_list[rownames(circa_result)]
  res_lot <- list()
  
  for (x in names(phase_list))
  {
    new_list <- lapply(1:nrow(phase_list[[x]]), function(y) phase_list[[x]][-y,])
    new_phases <- lapply(new_list, function(y) circular(y$phase*pi/12))
    new_summaries <- t(sapply(new_phases, function(y) summary(y)))[,c("Mean", "Rho")]
    rownames(new_summaries) <- phase_list[[x]]$names
    new_summaries[,"Mean"] <- new_summaries[,"Mean"]*12/pi
    pos_res <- ifelse(as.numeric(new_summaries[,"Mean"]) < 0, 24 + as.numeric(new_summaries[,"Mean"]) ,as.numeric(new_summaries[,"Mean"]))
    new_summaries[,"Mean"] <- pos_res - circa_result[x,"mean"]
    new_summaries[,"Rho"] <- new_summaries[,"Rho"] - circa_result[x,"rho"]
    res_lot[[x]] <- new_summaries
  }
  
  return(res_lot)
  
}

# Function to determine the most likely number of modes in the pressence of f-fold symmetry (2 or 3),
# cluster genes in f clusters and determine the circular tables for each cluster
multimodal_analysis <- function(circa_table)
{
  circa_radians <- circular(circa_table$phase*pi/12)
  modes_summary <- summary(circa_radians)
  modes_ray <- rayleigh.test(circa_radians, mu=circular(modes_summary["Mean"]))$p.value
  
  circa_radians_2 <- circular(((circa_table$phase*2) %% 24)*pi/12)
  modes_summary_2 <- summary(circa_radians_2)
  modes_ray_2 <- rayleigh.test(circa_radians_2, mu=circular(modes_summary_2["Mean"]))$p.value
  
  circa_radians_3 <- circular(((circa_table$phase*3) %% 24)*pi/12)
  modes_summary_3 <- summary(circa_radians_3)
  modes_ray_3 <- rayleigh.test(circa_radians_3, mu=circular(modes_summary_3["Mean"]))$p.value
  
  expected_modes <- which.min(c(modes_ray, modes_ray_2, modes_ray_3))
  expected_p_value <- min(c(modes_ray, modes_ray_2, modes_ray_3))
  
  result_FOCC <- CirClust(circa_table$phase, expected_modes, 24, method = "FOCC")
  clusters_focc <- lapply(1:max(result_FOCC$cluster), function(x) circa_table[result_FOCC$cluster == x,])
  names(clusters_focc) <- paste0("cluster_", 1:max(result_FOCC$cluster))
  
  cluster_radians <- lapply(clusters_focc, function(x) circular(x$phase*pi/12))
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
                             tr.height = 0.1)
{
  
  png(filename, width = 800, height = 800)
  
  sem <- sqrt(res.table$var)/sqrt(res.table$n)
  quan.table <- data.frame(go=rownames(res.table), res.table, min=res.table$first-sem, max=res.table$third+sem)
  
  # Create axis and initialize plot
  circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0), gap.degree=0, track.margin =c(0.005, 0.005))
  circos.initialize("a", xlim = c(0, 24)) # a means that there is a sector, which will go from 0 to 24
  
  circos.track(ylim = c(0, 1.5), track.height = 0.001, 
               bg.col = NA, bg.border=NA, panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks), 
                             labels.cex = 1, lwd = 2)
               })
  
  
  my_color=met.brewer(color.palette, n = nrow(quan.table))
  
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
          circos.track(ylim = c(0, 1.5), track.height = tr.height, 
                       bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                         xlim = CELL_META$xlim
                         # Estos primeros segmentos son para marcar los ZT
                         circos.segments(seq(0,24,4), 0,
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
                             tr.height = 0.1)
{
  
  png(filename, width = 1600, height = 1600, res=300)
  
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