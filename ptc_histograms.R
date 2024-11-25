library(dplyr)
library(ggplot2)

############################################################################################
plot_histogram_read_counts_per_target <- function (x=vector(), count=numeric(), id=character()) {
    #function to plot histogram of per target coverage per sample
    #validation of input vector
    if (length(x)==0) {
        stop("Input vector 'x' is empty. Please provide data.")
    }

    #print histogram for vector
    data_v <- data.frame( values = x)
    output_file <- paste0(id, "_histogram_of_read_counts_per_target_", count, ".pdf")
    pdf(output_file)
    ggplot(data = data_v, aes(x=values)) + geom_histogram(binwidth = 100, fill = "skyblue", color = "black", alpha = 0.7) + labs(title = paste("Read counts per target, req. sample id:", id), x = "Coverage per target", y = "Frequency") + theme_minimal()
    dev.off()

    # Check if the file is created and is non-zero size
    if (file.info(output_file)$size == 0) {
        stop("Error: PDF file is empty, something went wrong with the plot creation.")
    }
}
##########################################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Please provide dir path to per target coverage metrics files")
}
#fetch path to files from cmd arguments
files_path <- args[1]
cat("Path to metrics files provided:", files_path)

#List files and read files in a list of data frames
data.list <- list.files(files_path, pattern="*.metrics", full.names=T )
print(paste0("My files to read:", data.list))

#read files into list of dataframes
data.frames <- lapply(data.list, read.csv, sep="\t", header=T, row.names=NULL)
#generate list of req sample ids
#get list of lists after splitting with "/"
names.data.frames <- lapply(data.list, strsplit, "/")

#get last element after "/"
names.data.files <- lapply(names.data.frames, function(l){tail(tail(l[[1]],1),1)})

#using file names to extract first entry after splitting with "." -> this is the request sample id
req.sample.ids <- sapply(names.data.files, function (x){head(unlist(strsplit(x, "\\."))[[1]],1)})

#Get read counts in a list of vectors
read.counts <- lapply(data.frames, function(df) {df$read_count})
#add req sample ids as names
read.counts <- setNames(read.counts, req.sample.ids)

#plot histograms of read counts per target
c <- 1
for (id in names(read.counts)) {
    print(paste("Plotting histogram for sample ID:", id))
    plot_histogram_read_counts_per_target(read.counts[[id]], c, id)
    c <- c+1
}

