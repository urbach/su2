library("hadron")
library(optparse)

## Example use of script:
## cd /hiskp4/gross/masterthesis/analyse/flowresults/
## Rscript ../code/U1_analyse_potential/gflowanisotropy.R --beta 1.6521 --skip 1900 -r 16 -t 20 --betaone 1.6521 --xidiff -x 0.8 --datapath confL16.T20.b1.6521.x0.80 --plotpath results/

## TODO: replace measurement of xi, at the moment, we are calculating the expactation value of the ratio E_ts/E_ss,
## but we should be calculating the ratio of the expectation values
## to do this: read in E_ss as well, determine E_ts=E-E_ss,
## do this for all configs and then take expectation value of E_ts and E_ss

## read in parameters from command line
if (TRUE) {
    # set option list
option_list <- list(
    make_option(c("-b", "--beta"), type = "double", default = 1.7,
    help = "beta-Parameter of simulation [default %default]"),
    make_option(c("--skip"), type = "integer", default = 1000,
    help = "how many configurations are skipped when reading in [default %default]"),

    make_option(c("-r", "--Ns"), type = "integer", default = 16,
    help = "L of lattice (spatial extent) [default %default]"),

    make_option(c("-t", "--Nt"), type = "integer", default = 16,
    help = "T of lattice (temporal extent) [default %default]"),
    make_option(c("--betaone"), type = "double", default = 0,
    help = "input beta at corresponding xi = 1, appears in summarytable of results [default %default]"),
    make_option(c("-x", "--xi"), type = "double", default = 0,
    help = "xi used in lattice, only used if xidiff = TRUE, else
            xi is assumed to be L/T [default %default]"),

    make_option(c("--xidiff"), action = "store_true", default = FALSE,
    help = "Is xi different to L/T? [default %default]"),
    make_option(c("--datapath"), type = "character", default = "./",
    help = "path to where the datafiles are stored [default %default]"),
    make_option(c("--plotpath"), type = "character", default = "./",
    help = "path to where the plots are stored [default %default]"),
    make_option(c("--basename"), type = "character", default = "gradient_flow",
    help = "basename of resultfiles [default %default]"),

    make_option(c("--myfunctions"), type = "character",
        default = "/hiskp4/gross/masterthesis/analyse/code/U1_analyse_potential/",
    help = "path to where this script is stored,
            relative to folder where script is executed,
            only needed for printing git hash [default %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
}

xiin <- opt$Ns / opt$Nt
if (opt$xidiff) {
    xiin <- opt$xi
}

# determine git commit hash
printgitcommit <- function(pathtogit) {
    #get git commit hash of myfunctions.R, should be the same as of any script
    cwd <- setwd(pathtogit)
    githash <- try(system("git rev-parse --short HEAD", intern = TRUE))
    setwd(cwd)
    print(paste("## on git commit", githash))
    return(githash)
}
try(githash <- printgitcommit(opt$myfunctions))

## as determied in https://arxiv.org/pdf/2212.09627.pdf
## take c0 to reproduce the result from the potential measurements
c0 <- 1.628e-3
dc0 <- 9.1e-5


## get all filenames
filelist <- getorderedfilelist(path = opt$datapath, basename = opt$basename, last.digits = 6, ending = "")
## prepare emtpy containers
xilistrange <- c()
timesrange <- c()
xilistmin <- c()
timesmin <- c()

## read in the data to determine the times,
## only read in selected columns to save time
data <- read.table(file = filelist[1], header = F, skip = 1,
    # colClasses = c("numeric", "numeric", "NULL", "NULL", "NULL", "numeric", rep("NULL", 3)),
    colClasses = c("numeric", "numeric", "NULL", "NULL", "numeric", rep("NULL", 4)),
    # col.names = c("t", "xi", NA, NA, NA, "E", rep(NA, 3)))
    col.names = c("t", "xi", NA, NA, "E", rep(NA, 4)))

timesteps <- length(data$t)
timelist <- data$t

## empty arrays in which results for all times for xi, tsqE are stored
resultxi <- array(data = rep(NA, timesteps * (length(filelist) - opt$skip)),
        dim = c(timesteps, length(filelist) - opt$skip))

resulttsqE <- resultxi

for (index in seq(opt$skip + 1, length(filelist))) {
    ## for each file, read in necessary columns: t, xi, E, skip header
    data <- read.table(file = filelist[index], header = F, skip = 1,
    # colClasses = c("numeric", "numeric", "NULL", "NULL", "NULL", "numeric", rep("NULL", 3)),
    colClasses = c("numeric", "numeric", "NULL", "NULL", "numeric", rep("NULL", 4)),
    # col.names = c("t", "xi", NA, NA, NA, "E", rep(NA, 3)))
    col.names = c("t", "xi", NA, NA, "E", rep(NA, 4)))
    ## determine t^2E, take 2*E due to Lüschers formula
    data$tsqE <- 2 * data$t * data$t * data$E
    ## determine indices for which t^2E in c0 +/- dc0
    match <- which(abs(data$tsqE - c0) < dc0)
    ## select anisotropies corresponding to these indices
    xilistrange <- append(xilistrange, data$xi[match])
    timesrange <- append(timesrange, data$t[match])
    ## determine indices for which abs(t^2E-c0) is minimal
    match <- which(abs(data$tsqE - c0) == min(abs(data$tsqE - c0)))
    # print(match)
    ## select anisotropies corresponding to these indices
    xilistmin <- append(xilistmin, data$xi[match])
    timesmin <- append(timesmin, data$t[match])
    ## TODO: determine anisotropy, Energy by interpolating between the t^2Es that are closest to c0
    ## save all values for plot
    resultxi[, index - opt$skip] <- data$xi
    resulttsqE[, index - opt$skip] <- 2 * data$t^2 * data$E
}

## get mean and error
## if all times are the same, uwerr does not work
uwerrxirange <- uwerrprimary(xilistrange)
try(uwerrtrange <- uwerrprimary(timesrange))
uwerrximin <- uwerrprimary(xilistmin)
try(uwerrtmin <- uwerrprimary(timesmin))


## save result
## nom: number of occurences which fulfill abs(t^2E-c0) < dc0
## noc: number of different starting configurations
err <- try(result <- data.frame(beta = opt$beta, L = opt$Ns, T = opt$Nt, xiin = xiin,
xirange = uwerrxirange$value, dxirange = uwerrxirange$dvalue, timerange = uwerrtrange$value, dtimerange = uwerrtrange$dvalue,
ximin = uwerrximin$value, dximin = uwerrximin$dvalue, timemin = uwerrtmin$value, dtimemin = uwerrtmin$dvalue,
c0 = c0, dc0 = dc0, githash = githash,
noc = length(filelist) - opt$skip, nom = length(xilistrange)))

## if taking the means did not work, assume the time is the same everywhere
if(inherits(err, "try-error")){
 result <- data.frame(beta = opt$beta, L = opt$Ns, T = opt$Nt, xiin = xiin,
xirange = uwerrxirange$value, dxirange = uwerrxirange$dvalue, timerange = times[1], dtimerange = 0,
ximin = uwerrximin$value, dximin = uwerrximin$dvalue, timemin = times[1], dtimemin = 0,
c0 = c0, dc0 = dc0, githash = githash,
noc = length(filelist) - opt$skip, nom = length(xilistrange))
}
print(result)

## save results
filename <- sprintf("%ssummaryanisotropygflowb%.3f.csv",
                    opt$plotpath, opt$betaone)
columnnames <- FALSE
if (!file.exists(filename)) {
    columnnames <- TRUE
}
write.table(result, filename,
        append = TRUE, row.names = FALSE, col.names = columnnames)

## prepare list of xi(t), E(t), plot
## for each time, determine mean value and sd of Energy and anisotropy
## Maybe replace sd by different function to take autocorrelation into account?
## plot and save the results

xitime <- apply(resultxi, MARGIN = 1, FUN = mean)
dxitime <- apply(resultxi, MARGIN = 1, FUN = sd)
Etime <- apply(resulttsqE, MARGIN = 1, FUN = mean)
dEtime <- apply(resulttsqE, MARGIN = 1, FUN = sd)

timeresult <- data.frame(t = timelist, xi = xitime, dxi = dxitime, tsqE = Etime, dtsqE = dEtime)

filename <- sprintf("%s/gflowb%fx%fL%dT%d", opt$plotpath, opt$beta, xiin, opt$Ns, opt$Nt)
write.table(x  = timeresult, file = paste(filename, ".csv", sep = ""), row.names = F, col.names = T)

pdf(paste(filename, ".pdf", sep = ""), title = "")

plotwitherror(x = timeresult$t, y = timeresult$tsqE, dy = timeresult$dtsqE,
main = "Energy", xlab = "t", ylab = "t^2E")

plotwitherror(x = timeresult$t, y = timeresult$xi, dy = timeresult$dxi,
main = "Anisotropy", xlab = "t", ylab = "xi")

