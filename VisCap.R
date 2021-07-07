# VisCap: Visualize normalized capture coverage
# Generates heatmap of exon coverage from a directory of sample_interval_summary files
# By Trevor Pugh, March 2012 - April 2013

###########
#Libraries#
###########

library("cluster")
library("gplots")
library("zoo")

##########
#Defaults#
##########

svn.revision                <- "$Id: VisCap.R 1372 2013-04-25 13:53:52Z tp908 $"
source('VisCap.cfg')
#lane_dir                    <- "\\\\rfanfs.research.partners.org\\NGS"
#out_dir                     <- "\\\\rfanfs.research.partners.org\\gigpad_clinical\\Facilities\\Laboratory of Molecular Medicine_4171303\\VisCap_Outputs"

#lane_dir                    <- "\\\\rfanfs.research.partners.org\\NGS"
#out_dir                     <- "\\\\rfanfs.research.partners.org\\gigpad_clinical\\Facilities\\Laboratory of Molecular Medicine_4171303\\VisCap_Outputs"
# interval_list_dir           <- "\\\\Sfa6\\lmm$\\DEVELOPMENT\\ACTIVE_DEVELOPMENT\\NEXT_GEN_COPY_NUMBER\\VisCap\\interval_lists"   
#explorer_file               <- "C:\\Windows\\explorer.exe"
#cov_file_pattern            <- ".target.cov.sample_interval_summary$"
#cov_field                   <- "_total_cvg"
#interval_file_pattern       <- ".interval_list$"
#ylimits                     <- c(-2, 2)
#iqr_multiplier              <- 3
#threshold.min_exons         <- 1
#iterative.calling.limit     <- 0        #Set to 0 to iterate until all failed samples are removed
#infer.batch.for.sub.out_dir <- TRUE     #Set to FALSE to prompt users for output directory
#clobber.output.directory    <- FALSE    #Set to FALSE to stop run when output directory already exists

#Setting a path to data skips prompts. Set to FALSE for deployment.
#dev_dir    <- "\\\\Sfa6\\lmm$\\DEVELOPMENT\\ACTIVE_DEVELOPMENT\\NEXT_GEN_COPY_NUMBER\\VisCap\\test_data"
dev_dir     <- FALSE


source("functions.R")

###########
#Arguments#
###########

#Argument collection and parsing
arguments <- commandArgs(trailingOnly = TRUE)

if(length(arguments) == 1) {  
    if(dev_dir != FALSE) {
        #Skips prompts if dev_dir is set
        lane_dir <- dev_dir
        out_dir <- dev_dir
    } else {
        #Collect input and output information from user

        #Input directory
        lane_dir <- choose.dir(caption = "Select a lane directory (e.g. L001):", default = lane_dir)
        if(is.na(lane_dir)) {
            try( winDialog.nonint(type="ok", "Run canceled. No input lane directory provided."), silent=TRUE)
            q(save="no")
        }
        
        #Output diretory: Attempt to derive batch information from file name, prompt user if unsuccessful
        file1 <- list.files(lane_dir, full.names=TRUE, pattern=cov_file_pattern, recursive=TRUE)[1]
        batch.regex <- "__(B*[0-9]+)"
        batch.match <- regexec(batch.regex, file1)[[1]]
        batch <- substring(file1, batch.match[2], batch.match[2] + attr(batch.match, "match.length")[2] - 1)
        if((infer.batch.for.sub.out_dir == FALSE) || is.na(batch)) {
            out_dir <- choose.dir(caption = "Select an output directory:", default = out_dir)
            batch   <- basename(out_dir)
        } else {
            out_dir <- paste(out_dir, batch, sep="/")
        }
        if(is.na(out_dir)) {
            try( winDialog.nonint(type="ok", "Run canceled. No output directory provided."), silent=TRUE)
            q(save="no")
        }

        #If output directory already exists, prompt user to overwrite
        if((clobber.output.directory == FALSE) & (file.exists(out_dir))) {
            overwrite <- try( winDialog.nonint(type="yesno", "Output directory already exists. Overwrite?"), silent=TRUE)
            if(overwrite == "NO") {
                shell(paste(explorer_file, out_dir, sep=" "), wait=FALSE)
                q(save="no")
            }
        }
		
		interval_list_dir <-arguments[1]
	}
} else if(length(arguments) == 3) {
    #Uses provided command line arguments
    lane_dir <- arguments[1]
    out_dir  <- arguments[2]
    viscap.cfg <- arguments[3]
    batch    <- basename(out_dir)
} else {
    #Usage statement
    try( winDialog.nonint(type="ok", "Usage: VisCap.R lane_directory output_directory interval_lists_directory"), silent=TRUE)
    q(save="no")
}

######
#Main#
######

# Read coverage tables
mat.all <- make_matrix_from_cov_files(lane_dir, cov_file_pattern, cov_field)
mat.all[which(mat.all==0)]=0.00001

# Read interval name files
interval_lookup <- load_interval_names(interval_list_dir, interval_file_pattern)

# Sort matrix by genome coordinates found in rownames
chroms <- c(1:22,"X","Y", "MT", "M")
chroms <- factor(chroms, levels=chroms, labels=chroms, ordered=TRUE)
coords <- matrix(unlist(strsplit(rownames(mat.all), ":|-")), ncol=3, byrow=TRUE, dimnames=list(rownames(mat.all)))
coords <- coords[order(match(coords[,1], chroms), as.numeric(coords[,2]), as.numeric(coords[,3])),]
mat <- mat.all[rownames(coords),]

#Iteratively run VisCap algorithm, removing bad samples after each run
if(iterative.calling.limit == 0) {
    iterative.calling.limit <- dim(mat)[2]
}

for(iteration in 1:iterative.calling.limit) {
    if(iterative.calling.limit == 1) {
        out_dir.iteration <- out_dir
    } else {
        out_dir.iteration <- paste(out_dir, "/", batch, "_run", iteration, sep="")
    }
    dir.create(out_dir.iteration, showWarnings = FALSE, recursive=TRUE)
    
    # Normalize exon coverage by exon
    nmat_badX <- log2(divide_by_batch_median(mat))
    nmat <- rescale_chrX(nmat_badX, ylimits, iqr_multiplier, out_dir.iteration)

    # Call cnvs then plot heatmaps by chromosome and per-sample exon coverage
    heatmap_by_chrom(nmat, "QC_cnv_heatmap", ylimits, out_dir.iteration)
    cnv_bplot_and_calls <- call_cnvs(nmat, ylimits, interval_lookup, threshold.min_exons, iqr_multiplier, threshold.cnv_log2_cutoffs, out_dir.iteration)
    exon_plot_per_sample(nmat, ylimits, interval_lookup, cnv_bplot_and_calls, out_dir.iteration)

    # Write out matrix of log2 ratios to file
    outfile <- paste(out_dir.iteration, "/", "log2_ratio_table", ".xls", sep="")
    gene_exon <- annotate_interval_names(rownames(nmat), interval_lookup)
    nmat.with.gene_exon <- cbind(gene_exon, nmat)
    write.table(nmat.with.gene_exon, outfile, , quote = FALSE, sep = "\t", col.names = NA)

    #Write out VisCap run information
    run_info_table <- matrix(ncol = 2, byrow=TRUE, data = c(
            "Date",                                         date(),
            "VisCap command",                               paste(commandArgs(), collapse=" "),
            "Subversion revision information",              svn.revision,
            "Batch directory",                              lane_dir,      
            "Coverage file pattern",                        cov_file_pattern,
            "Field within coverage file",                   cov_field,              
            "Output directory",                             out_dir.iteration,                
            "Interval name lookup files",                   interval_list_dir,
            "Interval name lookup file pattern",            interval_file_pattern,
            "Exons used for CNV detection",                 dim(nmat)[1],
            "Samples used for CNV detection",               dim(nmat)[2],
            "Samples not used for CNV detection",           paste(c(colnames(mat.all)[!(colnames(mat.all) %in% colnames(nmat))], ""), sep=","),
            "Plot y-axis limits",                           paste(ylimits, collapse=","),
            "Minimum consecutive exons to call CNV",        threshold.min_exons,
            "IQR multiplier used for boxplots",             iqr_multiplier,
            "Static log2 ratio thresholds to call CNVs",    paste(threshold.cnv_log2_cutoffs, collapse=","),
            "Iteration",                                    iteration
            ))
    run_info_outfile <- paste(out_dir.iteration, "/", "VisCap_run_info", ".xls", sep="")
    write.table(run_info_table, run_info_outfile, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    save.image(file=paste(out_dir.iteration, "/", "session", ".Rdata", sep=""))

    #Remove failed samples from matrix for next run
    passes <- cnv_bplot_and_calls[[1]]$names[cnv_bplot_and_calls[[1]]$qc == "PASS"]
    if(length(passes) == dim(mat)[2]) {
        break
    } else {
        #Restrict mat only to samples that pass boxplot qc
        mat <- mat[,passes]
    }
}

# Open output directory and quit R
if((dev_dir == FALSE) && (length(arguments) == 0)) {
    shell(paste(explorer_file, out_dir, sep=" "), wait=FALSE)
}
quit(save="no")
