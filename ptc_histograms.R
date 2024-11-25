library(ggplot2)

############################################################################################
plot_histogram_read_counts_per_target <- function (x=vector(), count=numeric(), id=character(), max_value=numeric()) {
    #function to plot histogram of per target coverage per sample
    #validation of input vector
    if (length(x)==0) {
        stop("Input vector 'x' is empty. Please provide data.")
    }

    #print histogram for vector
    suppressWarnings(data_v <- data.frame( values = as.numeric(x) ))
    #remove NAs
    data_v <- data_v[!is.na(data_v$values),, drop=FALSE]
    #if values larger than 2500 replace them with 2500
    data_v$values[data_v$values>as.numeric(max_value)] <- as.numeric(max_value) #to avoid extremely large values skewing the histogram
    output_file <- paste0(id, "_histogram_of_read_counts_per_target_", count, ".pdf")
    median_value=median(data_v$values)
    #add median dashed line to histogram
    g <- ggplot(data = data_v, aes(x=values)) + geom_histogram(binwidth=50, color="black", fill="lightsteelblue3") + 
                                                geom_vline(aes(xintercept=median_value), color="red", linetype="dashed", linewidth=1.3) +
                                                annotate("text", median_value, y=10, color="red", label=paste("Median:", round(median_value)), vjust = 3.0, hjust=1.1, size = 4) +
                                                labs(title = paste0("Read counts per target, req. sample id: ", id), x = "Read counts per target", y = "Frequency") + 
                                                theme_minimal()
    #save plot to pdf file
    ggsave(output_file, plot = g)

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
files_path <- args[1] # to fetch list of ptc metrics for each sample
max_value <- args[2] # maximum value for histogram -> values larger than this value are set to the max specified
cat("Path to metrics files provided:", files_path, "\n")

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
    plot_histogram_read_counts_per_target(read.counts[[id]], c, id, max_value)
    c <- c+1
}

