### The BEGINNING ~~~~~
##
# NLSparrow--LD | First written by Filipe G. Vieira & Emma Fox with later slight modifications by George Pacheco.


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(optparse, tools, tidyverse, plyr, reshape2, RColorBrewer, extrafont, lemon)


# Imports extra fonts ~
loadfonts(device = "win", quiet = TRUE)


# Specifies options ~
options(width = 200) 


# Specifies arguments ~
option_list <- list(
  make_option(c('--ld_files'), action='store', type='character', default=NULL, help = 'File with list of LD files to fit and plot (if ommited, can be read from STDIN)'),
  make_option(c('--header'), action='store_true', type='logical', default=FALSE, help='Input file has header'),
  make_option(c('--col'), action='store', type='numeric', default=3, help='Which column is distance between sites? [%default]'),
  make_option(c('--ld'), action='store', type='character', default="r2", help='Which LD stats to plot (r2pear, D, Dp, r2) [%default]'),
  make_option(c('--n_ind'), action='store', type='numeric', default=0, help='How many individuals per-sample (for r^2 fitting correction)?'),
  make_option(c('-r', '--use_recomb_rate'), action='store_true', type='logical', default=FALSE, help='Assume constant recombination rate. [%default]'),
  make_option(c('--recomb_rate'), action='store', type='numeric', default=1, help='Recombination rate (or probability of recombination between adjacent sites in cM/Mb) to calculate genetic distances from physical distances. It is assumed to be constant throughout the whole dataset and, for human datasets, a common rule-of-thumb value is 1cM/Mb (1e-6). [%default]'),
  make_option(c('--max_kb_dist'), action='store', type='numeric', default=Inf, help='Maximum distance between SNPs (in kb) to include in the fitting analysis. [%default]'),
  make_option(c('--fit_boot'), action='store', type='numeric', default=0, help='Number of bootstrap replicates for fitting CI. [%default]'),
  make_option(c('--fit_bin_size'), action='store', type='numeric', default=250, help='Bin data into fixed-sized windows for fitting. [default %default bps]'),
  make_option(c('--fit_level'), action='store', type='numeric', default=1, help='Fitting level 0) no fitting, best of 1) Nelder-Mead, 2) and BFGS, 3) and L-BFGS-B). [%default]'),
  make_option(c('--plot_group'), action='store', type='character', default='File', help='Group variable'),
  make_option(c('--plot_data'), action='store_true', type='logical', default=FALSE, help='Also plot data points?'),
  make_option(c('--plot_bin_size'), action='store', type='numeric', default=0, help='Bin data into fixed-sized windows for plotting. [default %default bps]'),
  make_option(c('--plot_x_lim'), action='store', type='numeric', default=NULL, help='X-axis plot limit (in kb). [%default]'),
  make_option(c('--plot_y_lim'), action='store', type='numeric', default=NULL, help='Y-axis plot limit. [%default]'),
  make_option(c('--plot_axis_scales'), action='store', type='character', default='fixed', help='Plot axis scales: fixed (default), free, free_x or free_y'),
  make_option(c('--plot_size'), action='store', type='character', default='1,2', help='Plot size (height,width). [%default]'),
  make_option(c('--plot_scale'), action='store', type='numeric', default=1.5, help='Plot scale. [%default]'),
  make_option(c('--plot_wrap'), action='store', type='numeric', default=0, help='Plot in WRAP with X columns (default in GRID)'),
  make_option(c('--plot_no_legend'), action='store_true', type='logical', default=FALSE, help='Remove legend from plot'),
  make_option(c('--plot_shapes'), action='store_true', type='logical', default=FALSE, help='Use also shapes (apart from colors)'),
  make_option(c('--plot_line_smooth'), action='store', type='numeric', default=1000, help='LD decay curve smoothness'),
  make_option(c('--bin_quant'), action='store', type='numeric', default=0, help='Quantile to represent the bins (e.g. 0 = mean, 50 = median). [%default]'),
  make_option(c('-f','--plot_wrap_formula'), action='store', type='character', default=NULL, help='Plot formula for WRAP. [%default]'),
  make_option(c('-o','--out'), action='store', type='character', default=NULL, help='Output file'),
  make_option(c('--seed'), action='store', type='numeric', default=NULL, help='Seed for random number generator'),
  make_option(c('--debug'), action='store_true', type='logical', default=FALSE, help='Debug mode. Extra output...'))


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


opt$ld_files = "NLSparrow--LD.list"
opt$plot_x_lim = 100
opt$fit_level = 100
opt$plot_size = "2,5"
opt$max_kb_dist = 100
opt$fit_bin_size = 100
opt$col = 3
opt$bin_quant = 50
opt$fit_boot = 0
opt$plot_wrap = 1
opt$plot_wrap_formula = "File ~."
opt$plot_data = TRUE
opt$out_file = "Layka.pdf"


#opt$plot_axis_scales = "free_y"


# Set random seed ~
if(is.null(opt$seed))
  opt$seed <- as.integer(runif(1, 1, 100000))
cat("Random seed:", opt$seed, fill = TRUE)
set.seed(opt$seed)


# Parse input LD files ~
if(is.null(opt$ld_files)) 
  opt$ld_files <- file('stdin')
ld_files <- read.table(opt$ld_files, header=opt$header, stringsAsFactors=FALSE)
colnames(ld_files)[1] <- "File"
# Keep file list ordered
for (id in colnames(ld_files))
  ld_files[,id] <- factor(ld_files[,id], unique(ld_files[,id]), ordered=TRUE)
n_files <- nrow(ld_files)
if(opt$debug)
  print(ld_files)


# Parse LD stats to plot ~
if(!is.null(opt$ld))
opt$ld <- unlist(strsplit(opt$ld, ","))
if(!all(opt$ld %in% c("r2pear", "D", "Dp", "r2")))
stop("Invalid LD measure to plot", call. = opt$debug)
n_ld = length(opt$ld)


# Check if number of individuals was specified ~
if(opt$n_ind < 0)
stop("Number of individuals must be greater than zero", call. = opt$debug)
if(opt$n_ind > 0 && length(unique(ld_files$File)) > 1)
warning("Sample size for correction will be assumed equal for all samples!", call. = opt$debug)
if(!any(opt$ld %in% c("r2","r2pear")) && opt$n_ind)
stop("Number of individuals is only used for r^2 fitting", call. = opt$debug)
for(i in opt$ld)
cat("==> Fitting", i, "LD decay assuming a", ifelse(i == "Dp" || opt$n_ind == 0, "three (rate of decay, max LD and min LD)",
    "one (rate of decay)"), "parameter decay model", fill = TRUE)


# Check max_kb_dist parameter ~ 
if(opt$max_kb_dist < 50)
warning("Fitting of LD decay is highly unreliable at short distances (<50kb).", call. = opt$debug)


# Check binning sizes ~
if(opt$plot_bin_size > 0 && opt$plot_bin_size < opt$fit_bin_size)
stop("Ploting bin size must be greater than fiting bin size!", call.=opt$debug)


# Parse plot_size ~
opt$plot_size <- as.numeric(unlist(strsplit(opt$plot_size, ",")))
if(length(opt$plot_size) < 2)
  opt$plot_size <- c(opt$plot_size, opt$plot_size)


# Plot formula ~
if(!is.null(opt$plot_wrap_formula))
opt$plot_wrap_formula <- as.formula(opt$plot_wrap_formula)


# Set output file name (if not defined) ~
if(is.null(opt$out))
{if(is.null(opt$ld_files))
stop('Output file name required, when reading LD files from STDIN')
opt$out <- paste(basename(file_path_sans_ext(opt$ld_files)),".pdf", sep = "")}


# Load LD data ~
header <- c("Dist", "r2pear", "D", "Dp", "r2")
ld_data <- data.frame()
for (i in 1:n_files) {
     ld_file <- as.character(ld_files$File[i])

  
     # Reads point data ~
     if(opt$debug) cat("Reading file:", ld_file, fill = TRUE)
     tmp_data <- read.table(gzfile(ld_file), sep = "\t", quote = "\"", dec = ".")[-(1:(opt$col-1))]
     
     # Checks if file is valid ~
     if(ncol(tmp_data) < 5)
       stop('Invalid LD file format.\n', call. = opt$debug)
     
     
     # Add column labels ~
     colnames(tmp_data) <- header
     
     # Extract relevant columns ~
     tmp_data <- tmp_data[, which(names(tmp_data) %in% c("Dist", opt$ld))]
     
     # Filter by minimum distance ~
     tmp_data <- tmp_data[which(tmp_data$Dist < opt$max_kb_dist * 1000), ]
     
     # Convert all 'Inf' to NA ~
     tmp_data[mapply(is.infinite, tmp_data)] <- NA
     
     # Calculate genetic distances, according to Haldane's formula (assumes constant rate across all dataset) ~
     if(opt$use_recomb_rate && !is.null(opt$recomb_rate))
       tmp_data$Dist <- (1 - (1 - opt$recomb_rate*0.01/1e6)^(tmp_data$Dist)) / 2
     
     # Bin data ~
     if(opt$fit_bin_size > 1) {
       breaks <- seq(0, max(tmp_data$Dist)+opt$fit_bin_size, opt$fit_bin_size)
       tmp_data$Dist <- cut(tmp_data$Dist, breaks, head(breaks, -1))
       if(opt$bin_quant > 0) {
         tmp_data <- aggregate(. ~ Dist, data = tmp_data, quantile, probs=opt$bin_quant / 100, na.rm = TRUE)
       } else {
         tmp_data <- aggregate(. ~ Dist, data = tmp_data, mean, na.rm = TRUE)}}
     tmp_data$File <- ld_file
     ld_data <- rbind(ld_data, melt(tmp_data, c("File", "Dist"), variable.name = "LD", na.rm = TRUE))}


# Clean-up ~
rm(tmp_data)


# Removes factor in Dist ~
ld_data$Dist <- as.numeric(levels(ld_data$Dist))[ld_data$Dist]


# Sets maximum X-axis ~
if(is.null(opt$plot_x_lim)) {
  opt$plot_x_lim = max(ld_data$Dist)} else {
  opt$plot_x_lim = opt$plot_x_lim * 1000}


# Sets maximum Y-axis ~
if(!is.null(opt$plot_y_lim))
  opt$plot_y_lim <- c(0, opt$plot_y_lim)


# Adds extra info ~
ld_data <- merge(ld_files, ld_data, by="File", sort=FALSE)
# 
n_groups <- length(unique(ld_data[,opt$plot_group]))
n_plots <- n_files * n_ld / n_groups
if(opt$debug)
  print(head(ld_data))

# Fits decay ~

# Model function
ld_exp <- function(par, dist, ld_stat) {
  par <- as.numeric(par)
  if(ld_stat == "r2" || ld_stat == "r2pear") {
    C = par[1] * dist
    r2h = par[2]
    r2l = par[3]
    if(opt$n_ind){
      # LD decay curve adjusted for finite samples sizes
      ((10+C) / ((2+C)*(11+C))) * (1+((3+C)*(12+12*C+C^2))/(opt$n_ind*(2+C)*(11+C)))
    } else{
      # Theoretical expectation under to drift
      #1 / (1 + C)
      # Theoretical expectation with r2_high and r2_low
      (r2h - r2l) / (1 + C) + r2l}} else if(ld_stat == "Dp") {
    D0 = 1
    t = par[1]
    Dh = par[2]
    Dl = par[3]
    Dl + (Dh-Dl) * D0 * (1 - dist * opt$recomb_rate/1e6)^t}}


# Evaluates function ~
fit_eval <- function(par, obs_data) {
  if(length(unique(obs_data$LD)) != 1)
    stop("Invalid data.frame (several LD measures)", call.=opt$debug)
  model <- ld_exp(par, obs_data$Dist, obs_data$LD[1])
  eval <- sum((model - obs_data$value)^2)
  if(opt$debug) {
    str(par)
    print(eval)}
  eval}


# Fits function ~
fit_func <- function(x, fit_level) {
  ld_stat <- x$LD[1]
  # There is no fitting model for D
  if(ld_stat == "D") return(NULL)

  optim_tmp <- list()
  n_iter <- ifelse(fit_level>=10,fit_level,1)
  for(iter in 1:n_iter){
    # Fit LD model
    init_vals <- runif(3)
    if(ld_stat == 'Dp') {
      init_vals[1] = runif(1,10,20)
      par <- data.frame(init=init_vals, low_lim=c(0,0,0), up_lim=c(Inf,1,1))
    } else { # r2 and r2pear
      init_vals[1] = runif(1,0,0.1)
      par <- data.frame(init=init_vals, low_lim=c(0,0,0), up_lim=c(1,1,1))
    }
    if(opt$debug) str(par)
    
    optim_tmp <- append(optim_tmp, list("BFGS" = optim(par$init, fit_eval, obs_data=x, method="BFGS")) )
    if(fit_level > 1) optim_tmp <- append(optim_tmp, list("Nelder-Mead" = optim(par$init, fit_eval, obs_data=x, method="Nelder-Mead")) )
    if(fit_level > 2) optim_tmp <- append(optim_tmp, list("L-BFGS-B" = optim(par$init, fit_eval, obs_data=x, method="L-BFGS-B", lower=par$low_lim, upper=par$up_lim)) )
  }
  if(opt$debug) str(optim_tmp)

  # If not using the theoretical r2 decay curve (with r2h and r2l)
  if(opt$n_ind > 0)
    if(ld_stat != 'Dp') optim_tmp <- lapply(optim_tmp, function(x){x$par[2]=x$par[3]=0;x})
  
  # Filter out runs that not-converged and/or with out-of-bound parameters
  # Columns stand for: par1 (rate), par2 (high), par3 (low), score, counts, convergence, message
  optim_tmp <- Filter(function(x){x$convergence == 0 &
                                  x$par[1] >= par$low_lim[1] & x$par[1] <= par$up_lim[1] & 
                                  x$par[2] >= par$low_lim[2] & x$par[2] <= par$up_lim[2] &
                                  x$par[3] >= par$low_lim[3] & x$par[3] <= par$up_lim[3] &
                                  x$par[2] >= x$par[3]}, optim_tmp)
  
  # Pick best run
  optim_tmp <- optim_tmp[order(sapply(optim_tmp,'[[',2))[1]]
  if(length(optim_tmp[[1]]) == 0) stop("convergence analyses failed. Please try increasing the fit level (--fit_level)", call.=opt$debug)
  if(opt$debug) cat("Best fit for ", as.character(x$File[1]), " (",as.vector(ld_stat),"): ", names(optim_tmp), sep="", fill=TRUE)
  optim_tmp[[1]]$par
}

# Fit LD decay distribution
if(opt$fit_level > 0) {
  # Define line "resolution"
  smooth <- seq(1, opt$plot_x_lim, length=opt$plot_line_smooth)
  # Full data
  optim_fit <- ddply(ld_data, .(File,LD), fit_func, fit_level=opt$fit_level)

  # Bootstrap
  if(opt$fit_boot > 0) {
    boot_rep_fit <- c()
    for (b in 1:opt$fit_boot) {
      ld_bootdata <- ddply(ld_data, .(File, LD), function(x) {x[sample(nrow(x), size = nrow(x), replace = TRUE),]} )
      boot_rep_fit <- rbind(boot_rep_fit, data.frame(ddply(ld_bootdata, .(File,LD), fit_func, fit_level = opt$fit_level), Rep = b))}
    boot_fit <- as.data.frame(as.matrix(aggregate(cbind(V1, V2, V3) ~ File+LD, boot_rep_fit, quantile, probs = c(.025, .975), names = FALSE)), stringsAsFactors = FALSE)
    optim_fit <- merge(optim_fit, boot_fit, sort=FALSE)
    optim_data <- ddply(optim_fit, .(File,LD), function(x) data.frame(LD = x[, "LD"], Dist = smooth, value = ld_exp(x[, c("V1", "V2", "V3")],
                                                                      smooth, x[, "LD"]), ci_l = ld_exp(x[, c("V1.2", "V2.1", "V3.1")], smooth, x[, "LD"]), ci_u = ld_exp(x[, c("V1.1", "V2.2", "V3.2")], smooth, x[, "LD"])) )
  } else {
    optim_data <- ddply(optim_fit, .(File,LD), function(x) data.frame(LD = x[, "LD"], Dist = smooth, value = ld_exp(x[, c("V1", "V2", "V3")], smooth, x[, "LD"])))}
  
  
  
  
  
# Print best FIT parameters ~
optim_fit <- plyr::rename(optim_fit, c("V1" = "DecayRate", "V2" = "LDmax", "V3" = "LDmin", "V1.1" = "DecayRate_CI.u",
                                       "V1.2" = "DecayRate_CI.l", "V2.1" = "LDmax_CI.l", "V2.2" = "LDmax_CI.u",
                                        "V3.1" = "LDmin_CI.l", "V3.2" = "LDmin_CI.u"), warn_missing = FALSE)
  print(optim_fit)

  # Merge data together with extra info from input
  fit_data <- merge(ld_files, optim_data, sort = FALSE)
  if(opt$debug)
    print(head(fit_data, n = 10))}


# Corrects LD State (fit_data) ~
levels(fit_data$File <- sub("AllSamples_haplotypecaller.raw.vcf.FilteredUp.LD.gz", "No Pruning", fit_data$File))
levels(fit_data$File <- sub("AllSamples_haplotypecaller.raw.vcf.FilteredUp.Weight2.Pruned.LD.gz", "Pruning", fit_data$File))


# Corrects LD State (ld_data) ~
levels(ld_data$File <- sub("AllSamples_haplotypecaller.raw.vcf.FilteredUp.LD.gz", "No Pruning", ld_data$File))
levels(ld_data$File <- sub("AllSamples_haplotypecaller.raw.vcf.FilteredUp.Weight2.Pruned.LD.gz", "Pruning", ld_data$File))


# Reorders LD State ~
fit_data$File <- factor(fit_data$File, ordered = T,
                        levels = c("No Pruning", "Pruning"))


ld_data$File <- factor(ld_data$File, ordered = T,
                       levels = c("No Pruning", "Pruning"))


# Creates base plot ~
cat("==> Plotting data...", fill=TRUE)
plot <- ggplot() + 
  theme(panel.spacing = unit(1,"lines")) +
  coord_cartesian(xlim = c(0, opt$plot_x_lim), ylim = opt$plot_y_lim) +
  scale_colour_hue() +
  ylab("Linkage Disequilibrium") +
  xlab("Distance")


if(!is.null(opt$plot_wrap_formula)) {
  if(opt$plot_wrap) {
    plot <- plot + facet_rep_grid(opt$plot_wrap_formula, scales = opt$plot_axis_scales)
  } else {
    plot <- plot + facet_rep_grid(opt$plot_wrap_formula, scales = opt$plot_axis_scales)}}


# Adds LD decay fit CI ~
if(opt$fit_boot > 0) {
  grp <- NULL
  if(n_files == n_groups)
    grp <- opt$plot_group
  plot <- plot + 
          geom_ribbon(data = fit_data, aes_string(x = "Dist", ymin = "ci_l", ymax = "ci_u", group = opt$plot_group),
                      fill = "#000000", alpha = .15, show.legend = FALSE)}


# Adds data points ~
if(opt$plot_data){


# Checks format ~
  if(ncol(ld_data) < 4)
    stop("Invalid `ld_data` format.\n", call. = opt$debug)
  
  
# Bins data ~
  if(opt$plot_bin_size > 1) {
    breaks <- seq(0, max(ld_data$Dist) + opt$plot_bin_size, opt$plot_bin_size)
    ld_data$Dist <- cut(ld_data$Dist, breaks, head(breaks, -1))
    ld_data$Dist <- as.numeric(levels(ld_data$Dist))[ld_data$Dist]
    if(opt$bin_quant > 0) {
      ld_data <- aggregate(value ~ ., data = ld_data, quantile, probs = opt$bin_quant / 100)
    } else {
      ld_data <- aggregate(value ~ ., data = ld_data, mean)}}
  if(opt$debug)
    print(head(ld_data, n = 10))
  
  
# Add points ~
plot <- plot + 
        geom_point(data = ld_data, aes_string(x = "Dist", y = "value"),
                   size = .05, colour = "red", alpha = .25, show.legend = FALSE)}


# Adds LD decay best fit ~
if(length(opt$ld) > 0) {

  
# Selects fields that are variable ~
header <- names(which(lapply(lapply(fit_data, unique), length) > 1))


# Excludes non-relevant fields ~
grp <- header[!header %in% unique(c(as.character(opt$plot_wrap_formula), opt$plot_group, "Dist", "value", "File", "ci_l", "ci_u"))]
if(opt$debug) {
   print(grp)
   print(opt$plot_group)}
  

# Defines line type ~
  if(length(grp) == 0) grp <- NULL
  if(length(grp) > 1) stop("invalid number of linetype groups!")
  plot <- plot + 
          geom_line(data = fit_data, aes_string(x = "Dist", y = "value", linetype = grp), colour = "#000000", alpha = .5)
  

# If plotting data, adds a thin black line to help see the line ~
if(opt$plot_data)
plot <- plot +
        geom_line(data = fit_data, aes_string(x = "Dist", y = "value", colour = opt$plot_group), size = .1, alpha = .2,
                  colour = "#000000") +
        scale_x_continuous("Genomic Distance (Kb)",
                           breaks = c(5e+04, 1e+05, 1.5e+05, 2e+05),
                           labels = c("50", "100", "150", "200"),
                           limits = c(0, 100000),
                           expand = c(.005, .005)) +
        scale_y_continuous("Linkage Disequilibrium (r2)",
                           #breaks = c(.005, .01),
                           #labels = c("0.005", "0.01"),
                           #limits = c(0, .2),
                           expand = c(0, 0)) +
        theme(panel.background = element_rect(fill = "#ffffff"),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.title.x = element_text(size = 11, face = "bold", color = "#000000", margin = margin(t = 20, r = 0, b = 0, l = 0)),
              axis.title.y = element_text(size = 11, face = "bold", color = "#000000", margin = margin(t = 0, r = 20, b = 0, l = 0)),
              axis.text = element_text(color = "#000000", size = 8, face = "bold"),
              axis.line = element_line(colour = "#000000", size = .3),
              axis.ticks = element_line(color = "#000000", size = .3),
              strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = .3),
              strip.text = element_text(colour = "#000000", size = 9, face = "bold"),
              legend.position = "none",
              legend.background = element_blank(),
              legend.key = element_blank(),
              legend.spacing.y = unit(.4, "cm"),
              legend.key.height = unit(.45, "cm"),
              legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
              legend.box.margin = margin(t = 5, b = -20, r = 0, l = 30)) +
        guides(colour = guide_legend(title = "LD State:", title.theme = element_text(size = 12, face = "bold"),
                                     label.theme = element_text(size = 10),
                                     override.aes = list(size = 3, alpha = .9, starstroke = .15)))
  
  
# If plotting, apart from linetypes, also shapes (for B/W or color-blind printing) ~
if(opt$plot_shape) {
   smooth <- seq(1, opt$plot_x_lim, length=opt$plot_line_smooth)[seq(2,opt$plot_line_smooth,length=5)]
   sample_fit_data <- subset(fit_data, Dist %in% smooth)
   plot <- 
   plot + geom_point(data = sample_fit_data, aes_string(x = "Dist", y = "value", colour = opt$plot_group, shape = opt$plot_group),
                     show.legend = FALSE)}}


### Sets plot size ~
n_plots <- length(unique(ggplot_build(plot)$data[[1]]$PANEL))
if(!is.null(opt$plot_wrap_formula)) {
  par <- dcast(fit_data, opt$plot_wrap_formula, length, fill = 0)
  rownames(par) <- par[, 1]
  par <- par[,-1, drop = FALSE]
} else {
  par <- matrix(ncol = 1)}
if(opt$debug){
  cat(n_files, n_ld, n_groups, n_plots, fill = TRUE)
  cat(nrow(par), ncol(par), fill = TRUE)}
plot_height <- opt$plot_size[1] * nrow(par)
plot_width <- opt$plot_size[2] * ncol(par)


# Removes legend if plotting just a single variable ~
if(n_groups < 2 || opt$plot_no_legend) {
  plot <- plot + theme(legend.position = "none")
} else {
  plot_width = plot_width + 1}


# Saves plot ~
ggsave(plot = plot, file = "NLSparrow--LD_Median.pdf",
       device = cairo_pdf, scale = 1, width = 8, height = 8, dpi = 600, limitsize = FALSE)
x <- warnings()
