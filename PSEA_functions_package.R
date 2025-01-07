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

# Main function: Individual circular table creation
# create_circular_table <- function(circa_vector, iter = 9999)
# {
#   
#   circa_radians <- circular(circa_vector*pi/12)
#   
#   # Calculate circular measures
#   circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
#   
#   # Calculate circular variance and transform to hours again
#   circa_var <- var.circular(circa_radians)
#   
#   # Apply kupier's test
#   circa_kupier <- kuiper(circa_radians, rads = T, R=1)$p.value
#   
#   # Apply Rayleigh's test
#   circa_ray <- rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
#   
#   # Apply HR test
#   circa_hr <- HR_test(circa_radians, original = F, iter = iter)
#   
#   # Create table
#   circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi, 
#                             median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
#                             third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"], 
#                             var=circa_var, kuiper_p_value=circa_kupier, rayleigh_p_value=circa_ray, 
#                             hr_p_value = circa_hr[2] )
#   
#   
#   return(circa_table)
#   
# }

# Circular table complete
create_circular_table <- function(circa_vector, hr.on.large.sets = F, hr.on.large.sets.th=400,
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


# Functions for GO analysis
# Function to generate GO terms list with ancestors
go_list_ancestors <- function(go_vector, functional_data, go_column="GO", id_column="GID")
{
  functional_data <- select(org.Mpolymorpha.eg.db, 
                            keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), 
                            columns =  c("GID", "GO"))
  functional_data <- select(org.Mpolymorpha.eg.db, keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), columns =  c("GID", "GO"))
  
  go.list.first <- lapply(go_vector, 
                          function(x) go_phase_table(go_column="GO", go_term=x, total_phases_table, id_column="GID", functional_data))
  names(go.list.first) <- go_vector
}

# Function to create tables filtered by go
go_phase_table <- function(go_column, go_term, total_phases_table, id_column, functional_data)
{
  
  functional_reduced <- subset(functional_data, functional_data[[go_column]] == go_term)
  circa_reduced <- subset(total_phases_table, total_phases_table$names %in% functional_reduced[[id_column]] )
  return(circa_reduced)
}

#Function to enter a vector of GOs and compute the resulting circular statistiscs table
gos_circular_table <- function(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db",
                               iter=9999)
{
  
  # Load annotation
  functional_data <- select(get(annot_package), keys = keys(get(annot_package),keytype=id_column), 
                            columns =  c(id_column, go_column))
  
  # Generate phases lists per GO
  go_phases_list <- lapply(go_vector, function(x) go_phase_table(go_column, go_term = x, total_phases_table, id_column, functional_data))
  
  # Generate circular statistics table per GO
  go_circa_res <- lapply(go_phases_list, function(x) create_circular_table(x$phase, iter=iter))
  
  # Adjust p-values due to multiple testing
  go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
  go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
  kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
  go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
  ray_bh <- p.adjust(go_circa_ray, method = "BH")
  go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
  rh_bh <- p.adjust(go_circa_hr, method = "BH")
  
  # Return updated table
  go_circa_table <- data.frame(go_circa_num, kuiper_p_value=kupier_bh, rayleigh_p_value=ray_bh, hr_p_value=rh_bh)
  rownames(go_circa_table) <- go_vector
  
  return(go_circa_table)
}

# Function to adjust resulting table to non-negative phases
non_negative_table <- function(res_table)
{
  pos_res <- apply(res_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
  res_table[,2:5] <- pos_res
  
  return(res_table)
}

# Same functions using data.table for speed
# Function to create tables filtered by go
go_phase_table_rapid <- function(go_column, go_term, total_phases_table, id_column, functional_data)
{
  functional_reduced <- functional_data[get(go_column) == go_term,]
  circa_reduced <- total_phases_table[names %in% unlist(functional_reduced[,.(get(id_column))])]
  return(circa_reduced)
}

#Function to enter a vector of GOs and compute the resulting circular statistiscs table
gos_circular_table_rapid <- function(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db",
                                     iter=9999)
{
  
  functional_data <- select(get(annot_package), keys = keys(get(annot_package),keytype=id_column), 
                            columns =  c(id_column, go_column))
  functional_data <- data.table(functional_data)
  total_phases_table <- data.table(total_phases_table)
  
  go_phases_list <- lapply(go_vector, function(x) go_phase_table_rapid(go_column, go_term = x, total_phases_table, id_column, functional_data))
  
  go_circa_res <- lapply(go_phases_list, function(x) create_circular_table(x$phase, iter=iter))
  
  go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
  go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
  kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
  go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
  ray_bh <- p.adjust(go_circa_ray, method = "BH")
  go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
  rh_bh <- p.adjust(go_circa_hr, method = "BH")
  
  go_circa_table <- data.frame(go_circa_num, kuiper_p_value=kupier_bh, rayleigh_p_value=ray_bh, hr_p_value=rh_bh)
  rownames(go_circa_table) <- go_vector
  
  return(go_circa_table)
}

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
