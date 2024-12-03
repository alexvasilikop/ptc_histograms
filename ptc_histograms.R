library(ggplot2)
library(dplyr)
library(openxlsx)

############################################################################################
plot_histogram_read_counts_per_target <- function (x=vector(), count=numeric(), id=character(), max_value=numeric()) {
    #function to plot histogram of read counts per sample
    #validation of input vector
    if (length(x)==0) {
        stop("Input vector 'x' is empty. Please provide data.")
    }

    #print histogram for vector
    suppressWarnings(data_v <- data.frame( values = as.numeric(x) ))
    #remove NAs
    data_v <- data_v[!is.na(data_v$values),, drop=FALSE]
    #if values larger than max_value replace them with max_value
    data_v$values[data_v$values>as.numeric(max_value)] <- as.numeric(max_value) #to avoid extremely long tails on the histogram
    output_file <- paste0(id, "_histogram_of_read_counts_per_target_", count, ".pdf")
    median_value=median(data_v$values)
    #add median dashed line to histogram
    h <- ggplot(data = data_v, aes(x=values)) + geom_histogram(binwidth=100, color="black", fill="lightsteelblue3") + 
                                                geom_vline(aes(xintercept=median_value), color="red", linetype="dashed", linewidth=1.3) +
                                                annotate("text", median_value, y=10, color="red", label=paste("Median:", round(median_value)), vjust = 3.0, hjust=1.1, size = 4) +
                                                labs(title = paste0("Read counts per target, req. sample id: ", id), x = "Read counts per target", y = "Frequency") + 
                                                theme_minimal()
    #save plot to pdf file
    ggsave(output_file, plot = h)

    # Check if the file is created and is non-zero size
    if (file.info(output_file)$size == 0) {
        stop("Error: PDF file is empty, something went wrong with the plot creation.")
    }
}

plot_violinplot_across_samples <- function(x=data.frame(), selected_column=character()) {
    #plots violin plot of any provided column from df
    #data frame should have 2 columns one with normalized values for targets and one with the sample id the target belongs to

    if (!selected_column %in% colnames(x)){
        stop("Error in violinplot: requested column does not exist")
    }

    if (! "sample_id" %in% colnames(x)){
        stop("Error in violinplot: requested column 'sample_id' does not exist")
    }

    #transform for plotting
    new_df <- data.frame(sample_id=x$sample_id, values=as.numeric(x[[selected_column]]))
    new_df <- new_df[!is.na(new_df$values),]

    #create output file name
    output_file=paste0("violin_plots-", selected_column, "-across_samples.pdf")

    v <- ggplot(data=new_df, aes(x=values, y=sample_id)) + 
         geom_violin(scale="area", color="black", fill="lightsteelblue4") +
         geom_boxplot(width = 0.12, color = "black", outlier.colour = "red2", outlier.shape = 16, outlier.size = 1) +
         labs(title=paste("Violin plots of", selected_column, "across samples"), x=selected_column, y="Sample id") + 
         theme_minimal()

    ggsave(output_file, plot=v)

    if (file.info(output_file)$size == 0) {
        stop("Error: PDF file with violin plots is empty, something went wrong with the plot creation.")
    }
}

find_low_covered_regions <- function(df=data.frame(), low.threshold=numeric()) {
    #return df of low covered regions (chr, start, end, name) based on threshold of mean_coverage provided
    df$mean_coverage <- as.numeric(df$mean_coverage)
    new.df <- df[!is.na(df$mean_coverage),]
    df.low.covered <- filter(new.df, mean_coverage<as.numeric(low.threshold))
    #return data frame with chr, start , end and name of low covered region
    return(data.frame(chrom=df.low.covered$chrom, start=df.low.covered$start, end=df.low.covered$end, name=df.low.covered$name))
}

add_meancov_and_normcov <- function(df=data.frame(), filtered.df=data.frame()) {
    #create new filtered df with additional columns from matching rows in df (unfiltered)
    new.df.filtered <- filtered.df %>%
                       left_join(df, by=c("chrom", "start", "end", "name"))

    if (length(new.df.filtered$name) >= length(df$name)) {
        stop("New filtered set of regions in larger than unfiltered set!")
    }

    return(data.frame(chr=new.df.filtered$chrom, start=new.df.filtered$start, end=new.df.filtered$end, name=new.df.filtered$name, mean_coverage=new.df.filtered$mean_coverage, normalized_coverage=new.df.filtered$normalized_coverage))
}

add_exp_average_perbase_cov_depth <- function(x=data.frame(), read_length=numeric()) {
    #expected avg per base depth using region length, read length and read counts for region -> add to x as new column
    x <- x %>% mutate(
         avg_bp_cov = round((as.numeric(read_count)*as.numeric(read_length))/as.numeric(length)) 
         )
    return(x)
}

write_separate_samples_low_covered <- function(df=data.frame(), file_path=character()){
    #create workbook, add sheets and print df in in different excel sheets
    work.book <- createWorkbook()

    for (n in names(df)){
        addWorksheet(work.book, n)
        writeData(work.book, n, df[[n]])
    }

    saveWorkbook(work.book, file_path, overwrite=TRUE)
}

##########################################################################################################

#get cmd args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
    stop("Please provide correct number of arguments: 6")
}
#fetch path to files from cmd arguments
files_path  <- args[1] # to fetch list of ptc (per target coverage metrics) for each sample
max_value   <- args[2] # maximum value for histogram -> values larger than this value are set to the max specified
read_length <- args[3] # length of reads (e.g., 151)
t           <- args[4] # threshold for low coverage regions
outpath     <- args[5] # full output file path for common low regions (.csv)
excel_path  <- args[6] # full output file path for common low regions in each sample (xlsx)
cat("Path to metrics files provided:", files_path, "\n")

#List files and read files in a list of data frames
data.list <- list.files(files_path, pattern="*.metrics", full.names=T )
print(paste0("My files to read:", data.list))

#read files into list of dataframes
data.frames <- lapply(data.list, read.csv, sep="\t", header=T, row.names=NULL)


#generate list of req sample ids
#get list of lists after splitting with "/"
names.data.frames <- lapply(data.list, strsplit, "/")
#get last element after "/" --> file name
names.data.files <- lapply(names.data.frames, function(l){tail(tail(l[[1]],1),1)})
#using file names to extract first entry after splitting with "." -> this is the request sample id
sample.ids <- sapply(names.data.files, function (x){head(unlist(strsplit(x, "\\."))[[1]],1)})

####################################
#### Get read counts
#Get read counts in a list of vectors
read.counts <- lapply(data.frames, function(df) {df$read_count})
#add req sample ids as names
read.counts <- setNames(read.counts, sample.ids)
#plot histograms of read counts per target for each sample
c <- 1
for (id in names(read.counts)) {
    print(paste("Plotting histogram for sample ID:", id))
    plot_histogram_read_counts_per_target(read.counts[[id]], c, id, max_value)
    c <- c+1
}
#make empty boxplots
df_boxplot=data.frame(matrix(ncol=length(c("sample_id", "read_counts")),nrow=0))
colnames(df_boxplot) <- c("sample_id", "read_counts")

#plot_violinplots read counts across samples -> transform to 2 columns df with sample_id and read_counts
c <- 1
for (id in names(read.counts)) {
    v.names=rep(id, times=length(read.counts[[id]]))
    df_boxplot <- rbind(df_boxplot, data.frame(sample_id=v.names, read_counts=read.counts[[id]]))
    cat("Length of dataframe after adding ", c, " samples", length(df_boxplot$sample_id), "\n" )
    c <- c+1
}
#plot
plot_violinplot_across_samples(df_boxplot, "read_counts")


################
#### Get normalized coverages
#Get normalized coverages in a list of vectors (one vector for each sample)
norm.cov <- lapply(data.frames, function(df) {df$normalized_coverage})
#add req sample ids as names
norm.cov <- setNames(norm.cov, sample.ids)

#plot_violinplots_normalized_coverage_across_samples -> transform to 2 columns df with 1) sample_ids and 2) normalized_coverages
df_boxplot=data.frame(matrix(ncol=length(c("sample_id", "read_counts")),nrow=0))
colnames(df_boxplot) <- c("sample_id", "normalized_coverage")
#bind rows of each df to  df_boxplot
c <- 1
for (id in names(norm.cov)) {
    v.names=rep(id, times=length(norm.cov[[id]]))
    df_boxplot <- rbind(df_boxplot, data.frame(sample_id=v.names, normalized_coverage=norm.cov[[id]]))
    cat("Length of dataframe after adding ", c, " samples", length(df_boxplot$sample_id), "\n" )
    c <- c+1
}
#plot violin plots of normalized coverage
plot_violinplot_across_samples(df_boxplot, "normalized_coverage")

## uncomment if you want to add a column with expected coverage
#add_exp_average_perbase_cov_depth
#data.frames.with.names <- setNames(data.frames, sample.ids)
#data.frames_new <- lapply(data.frames.with.names, add_exp_average_perbase_cov_depth, read_length)

########################
#get low covered regions in df  -> list of dfs with chr, start, end, name of region
list.data.frames.with.names <- setNames(data.frames, sample.ids)
low.covered.regions <- lapply(list.data.frames.with.names, find_low_covered_regions, t)
#get common low covered
common.low.covered <- Reduce(function(df1, df2){merge(df1,df2)}, low.covered.regions)
#write common low covered in file (chr, start, end, name)
write.csv(file=outpath, common.low.covered, row.names=FALSE)


#write normalized and mean coverage per sample for low covered regions in an excel file (separate sheets for each sample)
list.filtered.data.frames.with.names <- lapply(list.data.frames.with.names, add_meancov_and_normcov, common.low.covered)
write_separate_samples_low_covered(list.filtered.data.frames.with.names, excel_path)