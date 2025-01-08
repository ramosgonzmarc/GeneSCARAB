.libPaths("/home/marcos/R/x86_64-pc-linux-gnu-library/4.1")

#install.packages("Directional")
# See also https://personality-project.org/r/html/cosinor.html
#library(CircMLE)
#library(Directional)
library(circular)
library(data.table)
library(org.Mpolymorpha.eg.db)

# Generate data from the uniform distribution on the circle.
data <- circular(runif(100, 0, 2*pi))
as.character(kuiper.test(data))
plot(data)
# Generate data from the von Mises distribution.
data <- rvonmises(n=3, mu=circular(0), kappa=3)
hola <- kuiper.test(data)
hola$statistic



# Mean resultant lenght
data <- circular(runif(100, 0, 2*pi))
rho.circular(data)

# Variance
var.circular(data)

# Median and mean
mean.circular(data)
median.circular(data)
hola <- summary(data)

# Try building the complete table
circa_table <- read.table("circacompare_intersection.tsv", header = T, sep = "\t")
head(circa_table)

# Function which takes one GO and the complete table and returns its circular statistics
library(org.Mpolymorpha.eg.db)

functional_data <- select(org.Mpolymorpha.eg.db, keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), columns =  c("GID", "GO"))
functional_reduced <- subset(functional_data, GO == "GO:0015979")

circa_reduced_ld <- circa_table[functional_reduced$GID, "ld.peak.time"]
circa_reduced_sd <- circa_table[functional_reduced$GID, "sd.peak.time"]

circa_reduced_ld <- circa_reduced_ld[!is.na(circa_reduced_ld)]
circa_reduced_sd <- circa_reduced_sd[!is.na(circa_reduced_sd)]

circa_vector <- circa_reduced_ld

# Individual function
create_circular_table_old <- function(circa_vector)
{
  
  circa_radians <- circular(circa_vector*pi/12)
  
  # Apply kupier's test
  kuiper_crit <- kuiper.test(circa_radians)$statistic
  kuiper_crit <- kuiper(circa_radians, rads = T)
  kuiper(circa_vector*pi/12, rads = T, R = 2)
  # P-value ranges using critical values
  
  {
  if (kuiper_crit > 2.001)
  {
    p_value <- "< 0.01"
  }
  else if (kuiper_crit > 1.862)
  {
    p_value <- "0.1 - 0.025"
  }
  else if (kuiper_crit > 1.747)
  {
    p_value <- "0.025 - 0.05"
  }
  else if (kuiper_crit > 1.62)
  {
    p_value <- "0.05 - 0.1"
  }
  else if (kuiper_crit > 1.537)
  {
    p_value <- "0.1 - 0.15"
  }
  else
  {
    p_value <- "> 0.15"
  }
  }
  
  
  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
  
  # Calculate circular variance and transform to hours again
  circa_var <- var.circular(circa_radians)
  circa_ray <- rayleigh.test(circa_radians, mu=circular(circa_summary["Median"]))$p.value
  
  circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi, 
                            median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
                            third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"], 
                            var=circa_var, kuiper_p_value=p_value, rayleigh_p_value=circa_ray )
  
  
  
  return(circa_table)
  
}

create_circular_table<- function(circa_vector)
{
  
  circa_radians <- circular(circa_vector*pi/12)
  
  # Apply kupier's test
  circa_kupier <- kuiper(circa_radians, rads = T, R=1)$p.value
  
  # Apply Rayleigh's test
  circa_ray <- rayleigh.test(circa_radians, mu=circular(circa_summary["Mean"]))$p.value
  
  # Apply HR test
  circa_hr <- HR_test(circa_radians, original = F, iter = 9999)
  
  # Calculate circular measures
  circa_summary <- summary(circa_radians) # "n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Rho"
  
  # Calculate circular variance and transform to hours again
  circa_var <- var.circular(circa_radians)
  
  circa_table <- data.frame(n=circa_summary["n"], first=circa_summary["1st Qu."]*12/pi, 
                            median=circa_summary["Median"]*12/pi, mean=circa_summary["Mean"]*12/pi,
                            third=circa_summary["3rd Qu."]*12/pi, rho=circa_summary["Rho"], 
                            var=circa_var, kuiper_p_value=circa_kupier, rayleigh_p_value=circa_ray, 
                            hr_p_value = circa_hr[2] )
  
  
  return(circa_table)
  
}
# Function to create tables filtered by go
go_phase_table <- function(go_column, go_term, total_phases_table, id_column, functional_data)
{
  
  functional_reduced <- subset(functional_data, functional_data[[go_column]] == go_term)
  circa_reduced <- subset(total_phases_table, total_phases_table$names %in% functional_reduced[[id_column]] )
  return(circa_reduced)
}

# Input must be in data.table format
go_phase_table_rapid <- function(go_column, go_term, total_phases_table, id_column, functional_data)
{
  functional_reduced <- functional_data[get(go_column) == go_term,]
  circa_reduced <- total_phases_table[names %in% unlist(functional_reduced[,.(get(id_column))])]
  return(circa_reduced)
}

# Function to enter a list of GOs and compute the resulting circular statistiscs table

gos_circular_table_rapid <- function(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db")
{
  
  functional_data <- select(get(annot_package), keys = keys(get(annot_package),keytype=id_column), 
                            columns =  c(id_column, go_column)) 
  functional_data <- data.table(functional_data)
  total_phases_table <- data.table(total_phases_table)
  
  go_phases_list <- lapply(go_vector, function(x) go_phase_table_rapid(go_column, go_term = x, total_phases_table, id_column, functional_data))
  
  go_circa_res <- sapply(go_phases_list, function(x) create_circular_table(x$phase))
  
  
  colnames(go_circa_res) <- go_vector
  
  return(go_circa_res)
}

gos_circular_table <- function(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db")
{
  
  functional_data <- select(get(annot_package), keys = keys(get(annot_package),keytype=id_column), 
                            columns =  c(id_column, go_column))
  
  go_phases_list <- lapply(go_vector, function(x) go_phase_table(go_column, go_term = x, total_phases_table, id_column, functional_data))
  
  go_circa_res <- lapply(go_phases_list, function(x) create_circular_table(x$phase))
  
  go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
  go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
  kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
  go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
  ray_bh <- p.adjust(go_circa_ray, method = "BH")
  go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
  rh_bh <- p.adjust(go_circa_hr, method = "BH")
  
  go_circa_table <- data.frame(go_circa_num, kuiper_p_value=kupier_bh, rayleigh_p_value=ray_bh, hr_p_value=hr_bh)
  rownames(go_circa_table) <- go_vector
  
  return(go_circa_table)
}

# Function to test if the aggregation between a point is significant
go_vector <- c("GO:0042254", "GO:0015979")
circa_table <- read.table("circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table <- data.frame(names=rownames(circa_table), phase=as.numeric(circa_table[["ld.peak.time"]]))
functional_data <- select(org.Mpolymorpha.eg.db, keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), columns =  c("GID", "GO"))

go_phase_table(go_column="GO", go_term="GO:0015979", total_phases_table, id_column="GID", functional_data)

res_table <- gos_circular_table(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db")


# Function to adjust resulting table to non-negative phases
non_negative_table <- function(res_table)
{
  pos_res <- apply(res_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x < 0), 24 + x , x))
  res_table[,2:5] <- pos_res
  
  return(res_table)
}


# Function to detect if a group is clustered on a certain point
x <- rvonmises(n=25, mu=circular(pi), kappa=2)
# General alternative
rayleigh.test(x)
# Specified alternative
hola <- rayleigh.test(x)
hola$p.value



watson(x, rads = TRUE)
x <- rvonmises(40, m = 2, k = 0)
kuiper(x, rads = TRUE)
watson(x, rads = TRUE)


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

x <- rvonmises(n = 40, m = 2, k = 10)
hola <- kuiper(u=x, rads = TRUE)

data <- circular(runif(100, 0, 2*pi))
kuiper(data, rads = T)
kuiper.test(data)
# Generate data from the von Mises distribution.
data <- rvonmises(n=100, mu=circular(0), kappa=3)
kuiper.test(data)


# Gráfico
circular_boxplot <- function(res.table, filename, color.palette = "Tam",
                             tr.height = 0.1)
{
  
  png(filename, width = 800, height = 800)
  
  sem <- sqrt(res.table$var)/sqrt(res.table$n)
  quan.table <- data.frame(go=rownames(res_table), res_table, min=res.table$first-sem, max=res.table$third+sem)
  
  # Primero creamos el eje e inicializamos el gráfico
  circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0), gap.degree=0, track.margin =c(0.005, 0.005))
  circos.initialize("a", xlim = c(0, 24)) # a significa que hay un sector, que irá de 0 a 24
  
  circos.track(ylim = c(0, 1.5), track.height = 0.001, 
               bg.col = NA, bg.border=NA, panel.fun = function(x, y) {
                 breaks = seq(0, 24, by = 4)
                 circos.axis(h = "top", major.at = breaks, labels = paste0("ZT", breaks), 
                             labels.cex = 1, lwd = 2)
               })
  
  
  my_color=met.brewer(color.palette, n = nrow(quan.table))
  
  # En función de los datos, puede ocurrir que el boxplot corte al amanecer en
  # su segmento inferior, en la propia caja, en su segmento superior o que no lo
  # corte, afectando al gráfico que se representará
  
  # El bucle for recorre todos los términos de GO analizados (si son más de 8,
  # hay que cambiar las dimensiones y la paleta empleada)
  
  
  
  for (i in 1:nrow(quan.table))
  {
    {
      # Si la caja no corta el amanecer, hay tres posibilidades
      if ((quan.table$first[i] < quan.table$median[i]) && (quan.table$median[i] < quan.table$third[i]))
      {
        # Si ninguno de los segmentos corta el amanecer, no hay que desdoblar nada
        if (quan.table$min[i] > 0 && quan.table$max[i] < 24)
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
        # Si el segmento inferior corta el amanecer (min será negativo)
        else if(quan.table$min[i] < 0)
        {
          # Habría que desdoblar el segmento desde min a 24 y desde 0 al max
          # sumándole 24 a min que será negativo
          circos.track(ylim = c(0, 1.5), track.height = tr.height, 
                       bg.col = "#eaeded", bg.border=NA, panel.fun = function(x, y) {
                         xlim = CELL_META$xlim
                         
                         circos.segments(seq(0,24,4), 0,
                                         seq(0,24,4), 1.5,
                                         col = "gray", lwd = 3, lty=3)
                         
                         circos.segments(24+quan.table[i,]$min, 0.75,
                                         24, 0.75,
                                         col = my_color[i], lwd =1.5)
                         circos.segments(0, 0.75,
                                         quan.table[i,]$max, 0.75,
                                         col = my_color[i], lwd =1.5)
                         
                         circos.rect(quan.table[i,]$first, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                     col = my_color[i], border = my_color[i], lwd=1.5)
                         
                         circos.segments(c(24+quan.table[i,]$min, quan.table[i,]$max), 
                                         c(0.75-0.2, 0.75-0.2),
                                         c(24+quan.table[i,]$min, quan.table[i,]$max), 
                                         c(0.75+0.2, 0.75+0.2),
                                         col = my_color[i], lwd = 2)
                         circos.segments(quan.table[i,]$median, 0.75-0.4,
                                         quan.table[i,]$median, 0.75+0.4,
                                         col = "black", lwd = 2)
                         
                       })
        }
        # El segmento superior corta al amanecer (max será mayor que 24)
        else if(quan.table$max[i] > 24)
        {
          # Habría que desdoblar el segmento desde min a 24 y desde 0 al max
          # restándole 24 a max que será mayor que 24
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
                                         quan.table[i,]$max-24, 0.75,
                                         col = my_color[i], lwd =1.5)
                         
                         circos.rect(quan.table[i,]$first, 0.75 - 0.4, quan.table[i,]$third , 0.75 + 0.4,
                                     col = my_color[i], border = my_color[i], lwd=1.5)
                         
                         circos.segments(c(quan.table[i,]$min, quan.table[i,]$max-24), 
                                         c(0.75-0.2, 0.75-0.2),
                                         c(quan.table[i,]$min, quan.table[i,]$max-24), 
                                         c(0.75+0.2, 0.75+0.2),
                                         col = my_color[i], lwd = 2)
                         circos.segments(quan.table[i,]$median, 0.75-0.4,
                                         quan.table[i,]$median, 0.75+0.4,
                                         col = "black", lwd = 2)
                         
                       })
        }
      }
      
      # Si la caja corta el amanecer
      else
      {
        # En este caso, habrá que cortar el segmento desde min a 24 y desde 0 a max
        # sin aplicar correcciones sobre los propios números, ya que min será positivo
        # y max será menor que 24. Además, habrá que desdoblar el rectángulo desde 
        # first a 24 y desde 0 a third
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
  
  
  # Ponemos una leyenda: establecemos los valores, que represente líneas,
  # definimos los colores y la posición de ajuste en el gráfico
  
  my_legend <- Legend(at = quan.table$go,
                      legend_gp = gpar(fill = my_color), title_position = "topleft", 
                      title = "")
  draw(my_legend, x = unit(1, "npc") - unit(2, "mm"), y = unit(10, "mm"), 
       just = c("right", "bottom"))
  
  circos.clear()
  dev.off()
}




# Pruebas
library(circular)
library(data.table)
library(org.Mpolymorpha.eg.db)
source("PSEA_functions.R")

go_vector <- c("GO:0042254", "GO:0015979")
# En un principio vamos a hacerlo para fotosíntesis, catabolismo de lípidos, síntesis de lípidos y síntesis de pigmentos
#go_vector <- c("GO:0015979", "GO:0016042", "GO:0008610","GO:0046148", "GO:0019684")

circa_table <- read.table("circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table <- data.frame(names=rownames(circa_table), phase=as.numeric(circa_table[["ld.peak.time"]]))
functional_data <- select(org.Mpolymorpha.eg.db, keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), columns =  c("GID", "GO"))

go_phase_table(go_column="GO", go_term="GO:0016117", total_phases_table, id_column="GID", functional_data)

table(functional_data$GO)

res_table <- gos_circular_table(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db")

nn_res_table <- non_negative_table(res_table)


library(circlize)
library(MetBrewer)
library(ComplexHeatmap)

circular_boxplot(nn_res_table, "prueba3.png", tr.height = 0.1, color.palette = "Austria")

gos_circular_table_rapid <- function(total_phases_table, go_vector, id_column="GID", go_column="GO", annot_package="org.Mpolymorpha.eg.db")
{
  
  functional_data <- select(get(annot_package), keys = keys(get(annot_package),keytype=id_column), 
                            columns =  c(id_column, go_column)) 
  functional_data <- data.table(functional_data)
  total_phases_table <- data.table(total_phases_table)
  
  go_phases_list <- lapply(go_vector, function(x) go_phase_table_rapid(go_column, go_term = x, total_phases_table, id_column, functional_data))
  
  go_circa_res <- sapply(go_phases_list, function(x) create_circular_table(x$phase))
  
  
  colnames(go_circa_res) <- go_vector
  
  return(go_circa_res)
}

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


# Comparación de dos grupos (MANOVA)


# Recrear su MANOVA
Distribution1_list = list(pigeons_c,Ants_1,bats_C)
Distribution2_list = list(pigeons_on,Ants_2,bats_M)

Distribution1<-Distribution1_list[[1]]
Distribution2<-Distribution2_list[[1]]

# Fases como tal

nv1<-length(Distribution1)
nv2<-length(Distribution2)

#Making treatment column
Group1 <-(as.matrix( rep((as.matrix((strrep(c("A"), 1)))),nv1),dim=c(nv1, 1)))
Group2 <- (as.matrix( rep((as.matrix((strrep(c("B"), 1)))),nv2),dim=c(nv2, 1)))
GroupsAB<-as.data.frame(c(Group1[,1], Group2[,1]))
# Grupo, una sola columna, pero en forma de data.frame, no vector a secas

DistGroups<-c(Distribution1, Distribution2) # concatenadas las fases en un vector
ManCirc<-summary(manova(cbind(cos(DistGroups),sin(DistGroups)) ~ GroupsAB[,1]))
Results_Examples$Man[i] <- ManCirc$stats[1,6]

# Sacar los que tengan un GO que esperamos que cambie
go_interest <- "GO:0015979"

# Try building the complete table
circa_table <- read.table("circacompare_intersection.tsv", header = T, sep = "\t")

library(org.Mpolymorpha.eg.db)

functional_data <- select(org.Mpolymorpha.eg.db, keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), columns =  c("GID", "GO"))
functional_reduced <- subset(functional_data, GO == go_interest)

circa_reduced_ld <- circa_table[functional_reduced$GID, "ld.peak.time"]
names(circa_reduced_ld) <- functional_reduced$GID
circa_reduced_sd <- circa_table[functional_reduced$GID, "sd.peak.time"]
names(circa_reduced_sd) <- functional_reduced$GID

circa_reduced_ld <- circa_reduced_ld[!is.na(circa_reduced_ld)]
circa_reduced_sd <- circa_reduced_sd[!is.na(circa_reduced_sd)]

circa_for_manova <- data.frame(names=names(circa_reduced_sd), ld_phase=circa_reduced_ld, sd_phase=circa_reduced_sd)

phase_vector <- c(circa_for_manova$ld_phase*pi/12, circa_for_manova$sd_phase*pi/12)
cond_vector <- c(rep("LD", nrow(circa_for_manova)), rep("SD", nrow(circa_for_manova)))
summary(manova(cbind(sin(phase_vector),cos(phase_vector)) ~ cond_vector))
plot.circular(circa_for_manova$ld_phase*pi/12)
plot.circular(circa_for_manova$sd_phase*pi/12)

mean.circular(circa_for_manova$ld_phase*pi/12)
var.circular(circa_for_manova$ld_phase*pi/12)
mean.circular(circa_for_manova$sd_phase*pi/12)
var.circular(circa_for_manova$sd_phase*pi/12)

# Linear discriminant analysis
library(MASS)
lda_res <- lda(cond_vector ~ cbind(sin(phase_vector),cos(phase_vector)), CV=F)
help("lda")

predict(lda_res, )


# Get ancestors
source("PSEA_functions.R")

go_vector <- c("GO:0042254", "GO:0015979", "GO:0016117", "GO:0006414")
circa_table <- read.table("circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table <- data.frame(names=rownames(circa_table), phase=as.numeric(circa_table[["ld.peak.time"]]))
functional_data <- select(org.Mpolymorpha.eg.db, 
                          keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), 
                          columns =  c("GID", "GO"))

go.list.first <- lapply(go_vector, 
       function(x) go_phase_table(go_column="GO", go_term=x, total_phases_table, id_column="GID", functional_data))
names(go.list.first) <- go_vector

# Mejor cambiar la de hacer las listas, para que sólo haga una de GOs con ids, luego ya hacer las fases
# De momento, lo solucionamos con:
go.list.sec <- lapply(go.list.first, function(x) as.vector(x[["names"]]))

library(GO.db)
bp_ancestors <- mget(go_vector, GOBPANCESTOR, ifnotfound=NA)
mf_ancestors <- mget(go_vector, GOMFANCESTOR, ifnotfound=NA)
cc_ancestors <- mget(go_vector, GOCCANCESTOR, ifnotfound=NA)

l <- list(bp_ancestors, mf_ancestors, cc_ancestors)
l_keys <- unique(unlist(lapply(l, names)))

go_ancestors <- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, l_keys))), l_keys)
go_ancestors <- lapply(go_ancestors, function(x) x[!is.na(x)])
go_ancestors <- lapply(go_ancestors, function(x) x[x != "all"])

# Usamos ahora go_ancestors y go.list.sec
# Añade los genes a cada uno nuevo si existe, si no, créalo


go.list <- go_ancestors
x <- l_keys[1]
y <- go_ancestors[[l_keys[1]]][1]

for (x in l_keys)
{
  for (y in go_ancestors[[x]])
  {
    if (y %in% names(go.list))
    {
      go.list[[y]] <- c(go.list[[y]], go.list.sec[[x]])
    }
    else
    {
      go.list[[y]] <- go.list.sec[[x]]
      
    }
    
  }
  
}

lapply(go.list, length)
go.list2 <- lapply(go.list, unique)
lapply(go.list2, length)


# With lapply

for (x in l_keys)
  {
    lapply(go_ancestors[[x]], function(y) 
      {if (y %in% names(go.list)) {go.list[[y]] <- c(go.list[[y]], go.list.sec[[x]])} else {go.list[[y]] <- go.list.sec[[x]]}})
  }
go.list    




# Get ALL successors for each parental
library(circular)
library(data.table)
library(org.Mpolymorpha.eg.db)
source("PSEA_functions_package.R")

go_vector <- c("GO:0042254", "GO:0015979", "GO:0016117", "GO:0006414")
circa_table <- read.table("circacompare_intersection.tsv", header = T, sep = "\t")
total_phases_table <- data.frame(names=rownames(circa_table), phase=as.numeric(circa_table[["ld.peak.time"]]))
functional_data <- select(org.Mpolymorpha.eg.db, 
                          keys = keys(org.Mpolymorpha.eg.db,keytype="GID"), 
                          columns =  c("GID", "GO"))

library(GO.db)
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
  go_column <- "GO"
  go.list <- lapply(go_offspring_comp, function(x) subset(functional_data, functional_data[[go_column]] %in% x)[["GID"]])
  go.list <- lapply(go.list, unique)
  
  return(go.list)
  
}

go.list.test <- create_gene_list_go(go_vector = go_vector, org.package = "org.Mpolymorpha.eg.db",
                                    go_column = "GO", id_column = "GID")

# Function to convert gene lists to phases table list
gene.list.to.phases <- function(gene.list, phase.table)
{
  circa_reduced <- lapply(gene.list, function(x) subset(phase.table, names %in% x))
  return(circa_reduced)

}
# TODO remove GOs with no associated genes

phases.list.test <- gene.list.to.phases(go.list.test, total_phases_table)

# Function to compute circular table for complete gene set
# complete_circular_table <- function(phase.list, iter)
# {
#   
#   # Generate circular statistics table per GO
#   go_circa_res <- lapply(phase.list, function(x) create_circular_table(x$phase, iter=iter))
#   
#   # Adjust p-values due to multiple testing
#   go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
#   go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
#   kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
#   go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
#   ray_bh <- p.adjust(go_circa_ray, method = "BH")
#   go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
#   rh_bh <- p.adjust(go_circa_hr, method = "BH")
#   
#   # Return updated table
#   go_circa_table <- data.frame(go_circa_num, kuiper_p_value=kupier_bh, rayleigh_p_value=ray_bh, hr_p_value=rh_bh)
#   rownames(go_circa_table) <- names(phase.list)
#   
#   # TODO incorporate raw p-values
#   
#   # Transform negative means and quantiles to positive
#   pos_res <- apply(go_circa_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
#   go_circa_table[,2:5] <- pos_res
#   
#   return(go_circa_table)
# }


# Highly sampled gene sets usually follows the general distribution for RNA accumulation
# The computation of HR metric in this cases takes too long and provide almost no information,
# so we recommend to perform it only if kupier and rayleigh don't report preferential direction
# and there is a strong reason to look for multimodality in these datasets
# For that, we split the sets by its n and implement another test against the general distribution
filter_high_sampled_sets <- function(phase.list)
{
  n.genes <- setNames(sapply(phase.list, nrow), NULL)
  phase.res <- phase.list[which(n.genes < 400)] # Adapt this based on experimental data
  return(phase.res)
}

phases.list.test.red <- filter_high_sampled_sets(phases.list.test)
length(phases.list.test)
length(phases.list.test.red)

circa_radians_test <- circular(phases.list.test$`GO:0008150`$phase*pi/12)
create_circular_table(circa_radians_test, 1)


# Paper Evaluating the power of a recent method for comparing two circular distributions: an alternative to the Watson U2 test
# WatsonU2, MANOVA o Kuiper
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

gen.dist.test <- test_against_gen_dist(phases.list.test, total_phases_table)

# Two or more distributions: MANOVA


# Adjust functions to compute circular table for complete gene set, not performing
# HR on all

# Individual circular table creation
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

create_circular_table(circa_vector = phases.list.test$`GO:0008150`$phase, hr.on.large.sets = F, hr.on.large.sets.th = 2000, iter = 999,
                      force.hr=F, force.hr.th=0)
create_circular_table(phases.list.test$`GO:0009987`$phase, hr.on.large.sets = F, threshold.hr = 400, iter = 999)
create_circular_table(phases.list.test$`GO:0006414`$phase, hr.on.large.sets = F, threshold.hr = 400, iter = 999)


complete_circular_table <- function(phase.list, hr.on.large.sets = F, hr.on.large.sets.th=400,
                                    iter = 999, force.hr=F, force.hr.th=0.05)
{
  
  # Generate circular statistics table per GO
  go_circa_res <- lapply(phase.list, function(x) create_circular_table(x$phase, hr.on.large.sets = hr.on.large.sets, 
                                                                       hr.on.large.sets.th=hr.on.large.sets.th, iter=iter,
                                                                       force.hr=force.hr, force.hr.th=force.hr.th))
  
  # Adjust p-values due to multiple testing
  go_circa_num <- t(sapply(go_circa_res, function(x) unlist(x[1:7])))
  go_circa_kupier <- sapply(go_circa_res, function(x) unlist(x[8]))
  kupier_bh <- p.adjust(go_circa_kupier, method = "BH")
  go_circa_ray <- sapply(go_circa_res, function(x) unlist(x[9]))
  ray_bh <- p.adjust(go_circa_ray, method = "BH")
  go_circa_hr <- sapply(go_circa_res, function(x) unlist(x[10]))
  rh_bh <- p.adjust(go_circa_hr, method = "BH")
  
  # Return updated table
  go_circa_table <- data.frame(go_circa_num, kuiper_p_value=go_circa_kupier, kuiper_p_value_adj=kupier_bh, 
                               rayleigh_p_value=go_circa_ray, rayleigh_p_value_adj=ray_bh,
                               hr_p_value=go_circa_hr, hr_p_value_adj=rh_bh)
  rownames(go_circa_table) <- names(phase.list)
  
  
  # Transform negative means and quantiles to positive
  pos_res <- apply(go_circa_table[,2:5], MARGIN=2, function(x) ifelse(as.numeric(x) < 0, 24 + as.numeric(x) , as.numeric(x)))
  go_circa_table[,2:5] <- pos_res
  
  return(go_circa_table)
}

res.table <- complete_circular_table(phases.list.test, iter=10, hr.on.large.sets = F)

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

library(circlize)
library(ComplexHeatmap)
library(MetBrewer)
circular_dotplot(phases.list = phases.list.test[1:8], "prueba_dot.png", tr.height = 0.1, color.palette = "Austria")


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

contributions_per_gene <- gene_contribution_to_set(res.table, phases.list.test)
contributions_per_gene$`GO:0006414`
circular_dotplot(phases.list = phases.list.test["GO:0006414"], "prueba_GO:0006414.png", tr.height = 0.1, color.palette = "Austria")

contributions_per_gene$`GO:0016117`
circular_dotplot(phases.list = phases.list.test["GO:0016117"], "prueba_GO:0016117.png", tr.height = 0.1, color.palette = "Austria")



# Determinar el número de modos
# GO:0006414 se puede considerar bimodal
# Rayleigh normal
res.table["GO:0006414",]

prueba_modos <- circular(phases.list.test$`GO:0006414`$phase*pi/12)
modos_summary <- summary(prueba_modos)
modos_ray <- rayleigh.test(prueba_modos, mu=circular(modos_summary["Mean"]))$p.value
modos_hr <- HR_test(prueba_modos, original = F, iter = 99999) #hr identifies multimodality

# Transformar en simetría de 2
prueba_2_modos <- (phases.list.test$`GO:0006414`$phase*2) %% 24
prueba_2_modos_radians <- circular(prueba_2_modos*pi/12)
modos_summary_2 <- summary(prueba_2_modos_radians)
modos_ray_2 <- rayleigh.test(prueba_2_modos_radians, mu=circular(modos_summary_2["Mean"]))$p.value # significant

# Transformar en simetría de 3
prueba_3_modos <- (phases.list.test$`GO:0006414`$phase*3) %% 24
prueba_3_modos_radians <- circular(prueba_3_modos*pi/12)
modos_summary_3 <- summary(prueba_3_modos_radians)
modos_ray_3 <- rayleigh.test(prueba_3_modos_radians, mu=circular(modos_summary_3["Mean"]))$p.value

# Transformar en simetría de 4 (esta ya no es frecuente)
prueba_4_modos <- (phases.list.test$`GO:0006414`$phase*4) %% 24
prueba_4_modos_radians <- circular(prueba_4_modos*pi/12)
modos_summary_4 <- summary(prueba_4_modos_radians)
modos_ray_4 <- rayleigh.test(prueba_4_modos_radians, mu=circular(modos_summary_4["Mean"]))$p.value

# Determine number of modes in the pressence of f-fold symmetry
number_of_modes <- function(circa_table)
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
  
  print(paste0("The most probable number of modes is ", expected_modes, 
               ", p-value: ", expected_p_value))
  
  return(expected_modes)
  
}

# This should be used in cases where HR but not Kuiper or Rayleigh reported significancy
aaaa <- number_of_modes(phases.list.test$`GO:0006414`)

# Compare two different distributions
test_two_dist <- function(phase.list.1, phase.list.2)
{
  if (length(phase.list.1) != length(phase.list.2))
  {
    stop("Error: length of lists differ.")
  }
  
  if (names(phase.list.1) != names(phase.list.2))
  {
    warning("The names of the sets do not match, check if they should.")
  }
  
  p_values_diff_test <- mapply(function(x, y) summary(manova(cbind(sin(c(x$phase*pi/12, y$phase*pi/12)),cos(c(x$phase*pi/12, y$phase*pi/12))) 
                                                        ~ c(rep("first", nrow(x)), rep("second", nrow(y)))))$stats[1,6],
                               phase.list.1, phase.list.2)
  
  p_values_diff_bh <- p.adjust(p_values_diff_test, method = "BH")
  
  p_values_diff_table <- data.frame(p_value=p_values_diff_test, 
                               p_value_adj=p_values_diff_bh)
  
  return(p_values_diff_table)
}


# Mirar dryR para sustituir los Venn y circacompare
# Leer también Circadian transcriptome processing and analysis: a workflow for muscle stem cells
# Este paper para dejar claro que lo que hablan de multimodalidad es absurdo si usan Kuiper en PSEA:
# Circular data in biology: advice for effectively implementing statistical procedures 
# A comparative study of algorithms detecting differential rhythmicity in transcriptomic data


# GeneSCARAB: Gene Set Cohesion Analysis for Rhythms Ascertainment in Bioprocesses
