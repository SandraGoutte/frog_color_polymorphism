################################################################################
###	Create a .bed file to have a custom sliding window to use in Pixy
################################################################################
## set working directory and bed file structure
setwd("/path/to/directory/")

## create a new .bed file
bed = matrix(ncol=3,nrow=1)
bed = as.data.frame(bed)
colnames(bed) = c("chrom", "chromStart","chromEnd")

## get sliding windows of 3000 bp with 1000 bp overlap
## set the parameters
window_size <- 3000
overlap <- 1000

## set the chromosome names
chroms=c('Scaffold_2__2_contigs__length_236349107',
        'Scaffold_1__1_contigs__length_194224858',
        'Scaffold_5__2_contigs__length_177367226',
        'Scaffold_3__2_contigs__length_170906183',
        'Scaffold_6__2_contigs__length_157086139',
        'Scaffold_7__2_contigs__length_128225229',
        'Scaffold_4__1_contigs__length_105873206',
        'Scaffold_10__2_contigs__length_97064778',
        'Scaffold_8__1_contigs__length_80766577',
        'Scaffold_9__1_contigs__length_67344957',
        'Scaffold_12__2_contigs__length_57092889',
        'Scaffold_11__1_contigs__length_50805362' )

## set the overall start and end SNP positions for each chromosome
start <- rep(1,12)
end <- c(236349107,194224858,177367226,170906183,157086139,128225229,105873206,
         97064778,80766577,67344957,57092889,50805362)

## create the windows
windows_matrix=list()
for(i in 1:12){
  # Calculate the starting positions of the windows
  starts <- seq(start[i], end[i] - window_size + 1, by = overlap)
  
  # Prepare a matrix to store start and end positions
  windows_matrix[[i]] <- matrix(NA, ncol = 3, nrow = length(starts))
  colnames(windows_matrix[[i]]) <- c("chromosome", "start", "end")
  
  # Populate the matrix with window coordinates
  windows_matrix[[i]][, "chromosome"] <- chroms[i]
  windows_matrix[[i]][, "start"] <- starts
  windows_matrix[[i]][, "end"] <- starts + window_size - 1
}

## check our matrix
head(windows_matrix)

## reformat
windows_matrix_genome <- ldply(windows_matrix, data.frame)
head(windows_matrix_genome)
tail(windows_matrix_genome)

## write the matrices to a BED file for each chromosome separately
for(i in 1:12){
  write.table(windows_matrix[[i]], 
              file = paste("sliding_windows_w3000_step1000_chr",i,
                           ".bed", sep=""), sep = "\t", col.names = FALSE,
              quote = FALSE)
}
## then open in Excel and remove the first column and scientific notations
## to avoid issues when using with Pixy
