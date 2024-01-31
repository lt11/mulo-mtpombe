## header ---------------------------------------------------------------------

### This script extracts the HIAF fields from the FORMAT
### data of a multisample vcf file. Then it plots the HIAF for
### all the samples.
### The correct temporal order of the samples is determined 
### from the vcf input file, where the samples IDs are correctly 
### sorted since the IDs are provided in the correct order in the wrapper.
### It works on the data produced (and automatically annotated) 
### by vardict.

rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(ggplot2)
library(here)

## function(s) ----------------------------------------------------------------

ExtractHIAF <- function(x) {
  ### It extracts the HIAF data from a vcf in a data-frame for each sample.
  ### All the other subfields of the FORMAT field are removed. Missing variants
  ### (e.g. "./.") are considered as HIAF = 0.
  ### 
  ### arguments:
  ### (1) a vcf in a data-frame
  ###
  ### returns:
  ### (1) the same vcf in a data-frame with a FORMAT comprising only the HIAF
  allFormats <- unlist(strsplit(x = x[1, 9], split = ":"))
  indHiaf <- which(allFormats == "HIAF")
  nCol <- ncol(x)
  for (indC in c(10:nCol)) {
    myCol <- x[, indC]
    ### to avoid a warning we convert the tags of missing variants to
    ### zeros, e.g. from "./.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:.:." to
    ### "0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0"
    indMissing <- grep(pattern = "[[:alnum:]]", x = myCol, value = F, invert = T)
    myCol[indMissing] <- paste(rep("0", length(allFormats)), collapse = ":")
    myHiaf <- as.numeric(sapply(strsplit(myCol, split = ":"), "[[", indHiaf))
    x[, indC] <- myHiaf
  }
  return(x)
}

INFOtoTYPE <- function(x) {
  ### It extract the data in the TYPE subfield in the INFO field. 
  ### All the other subfields of the INFO field are removed.
  ###
  ### arguments:
  ### (1) a vcf in a data-frame
  ###
  ### returns:
  ### (1) a data-frame where INFO has been replaced by TYPE
  allInfos <- unlist(strsplit(x = x[1, 8], split = ";"))
  indType <- grep("^TYPE=", allInfos)
  strType <- sub(pattern = "^.*TYPE=([^;]*);.*$",
                 replacement = "\\1",
                 x = x[[8]])
  x[, 8] <- strType
  return(x)
}

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(here())
### def off
# dirBase <- "/Users/Lorenzo/dev/dev-mulo-mtpombe"

### input and output folders
dirOutTop <- file.path(dirBase, "plots")
dir.create(dirOutTop, showWarnings = F)
dirOutMain <- file.path(dirOutTop, "af")
unlink(dirOutMain, recursive = T)
dir.create(dirOutMain)
dirIn <- file.path(dirBase, "var-calls", "mrg")
dirOutTab <- file.path(dirBase, "tabs", "af")
unlink(dirOutTab, recursive = T)
dir.create(dirOutTab, recursive = T, showWarnings = F)

### vcf file fixed column names (no samples)
hdVcfFix <- c("Chrom_id", "Pos_bp", "Var_id", "Ref_allele", "Alt_allele",
              "Qual_val", "Filter_tag", "Info_tags", "Format_def")

## clmnt ----------------------------------------------------------------------

### get all multisample files (one file has all the time-points)
allFiles <- list.files(dirIn, pattern = ".vcf.gz$", full.names = T)
for (indS in allFiles) {
  dfVcf <- read.table(file = indS)
  strVcf <- readLines(indS)
  allFields <- strsplit(grep(pattern = "^#CHROM", x = strVcf, value = T),
                        split = "\t")
  allFields <- unlist(allFields)
  nFields <- length(allFields)
  nSamples <- nFields - 9
  sampleFields <- allFields[c(10:nFields)]
  colnames(dfVcf)[1:9] <- hdVcfFix
  colnames(dfVcf)[10:nFields] <- sampleFields
  
  ### get replicate ID
  fileName <- sub("^([^\\.]*)(.*$)", "\\1", basename(indS))
  replicateID <- sub("^([^-]*)-(.*$)", "\\2", fileName)
  
  ### extract HIAF and variant type
  dfHiaf <- ExtractHIAF(dfVcf)
  dfHiafType <- INFOtoTYPE(dfHiaf)
  dfHiafType$Var_id <- paste0(dfHiaf$Chrom_id, "_", dfHiaf$Pos_bp)
  
  ### fake duplicates to test the next block
  ### dfHiafType <- rbind(dfHiafType[1:10, ],
  ###                     dfHiafType[10, ],
  ###                     dfHiafType[10,], 
  ###                     dfHiafType[11:nrow(dfHiaf), ])
  
  ### remove loci that appear more than once in the vcf, keeping the variant 
  ### showing the largest AF at any time-point
  dupVars <- dfHiafType$Var_id[duplicated(dfHiafType$Var_id)]
  indOut <- c()
  for (indD in unique(dupVars)) {
    indE <- which(dfHiafType$Var_id == indD)
    mtAfDup <- as.matrix(dfHiafType[indE, seq(10, nFields, 1)])
    indF <- which(mtAfDup == max(mtAfDup), arr.ind = T)[1]
    indOut <- c(indOut, indE[-indF])
  }
  if (length(indOut) != 0) {
    dfHiafType <- dfHiafType[-indOut, ]
  }
  ### make output folder and write the data to output
  dirOutTabSamp <- file.path(dirOutTab, replicateID)
  dir.create(dirOutTabSamp)
  dfHiafTypeOut <- dfHiafType[, c(1:5, 8, 10:nFields)]
  indTags <- which(colnames(dfHiafTypeOut) == "Info_tags")
  colnames(dfHiafTypeOut)[indTags] <- "Var_type"
  fileHiafTabOut <- file.path(dirOutTabSamp,
                              paste0(replicateID, "-hiaf-af.txt"))
  write.table(dfHiafTypeOut, file = fileHiafTabOut, append = F, quote = F,
              sep = "\t", row.names = F, col.names = T)
  
  ## prepare data-frame for plotting ------------------------------------------
  
  dfPlot <- data.frame()
  for (indR in 1:nrow(dfHiafType)) {
    
    dfPlot <- rbind(dfPlot,
                    cbind(rep(dfHiafType$Chrom_id[indR], nSamples),
                          rep(dfHiafType$Pos_bp[indR], nSamples),
                          rep(dfHiafType$Var_id[indR], nSamples),
                          rep(dfHiafType$Ref_allele[indR], nSamples),
                          rep(dfHiafType$Alt_allele[indR], nSamples),
                          rep(dfHiafType$Info_tags[indR], nSamples),
                          c(colnames(dfHiafType)[10:nFields]),
                          as.numeric(dfHiafType[indR, 10:nFields])))
  }
  colnames(dfPlot) <- c("Chrom_id", "Pos_bp", "Var_id", "Ref_allele",
                        "Alt_allele", "Var_type", "Sample_id", "Af_val")
  
  dfPlot$Var_id <- as.character(dfPlot$Var_id)
  dfPlot$Af_val <- as.numeric(dfPlot$Af_val)
  dfPlot$Sample_id <- factor(dfPlot$Sample_id,
                             levels = sampleFields)
  colourNamed <- c(c(SNV = "darkgreen"),
                   c(Deletion = "blue"),
                   c(Insertion = "red"),
                   c(Complex = "purple4"))
  dirOutRep <- file.path(dirOutMain, replicateID)
  dir.create(dirOutRep, showWarnings = F, recursive = T)
  
  ## plotting -----------------------------------------------------------------
  
  allPlotPath <- file.path(dirOutRep, "all-af.pdf")
  allPlot <- ggplot(data = dfPlot,
                    aes(x = Sample_id,
                        y = Af_val,
                        group = Var_id,
                        colour = Var_type)) +
    geom_point(alpha = 0.4) +
    geom_line(size = 1, alpha = 0.4) +
    scale_color_manual(name = "Variant type", values = colourNamed) +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
          # panel.grid.major = element_blank(), 
          # panel.grid.minor = element_blank(), 
          # panel.border = element_blank(), 
          panel.background = element_blank(),
          # title style
          # plot.title = element_text(size = 16,
          # face = "bold", hjust = 1), 
          # plot.subtitle = element_text(size = 12,
          # face = "bold", hjust = 1), 
          # axis style
          axis.title = element_text(size = 22), 
          axis.text.x = element_text(size = 14, angle = 45),
          axis.text.y = element_text(size = 22)) +
    coord_cartesian(ylim = c(.0, 1.)) +
    labs(title = "", subtitle = "",
         x = "Sample",
         y = "AF",
         size = 28)
  pdf(file = allPlotPath, width = 9, height = 9)
  print(allPlot)
  dev.off()
  
  for (indP in unique(dfPlot$Var_id)) {
    indB <- which(dfPlot$Var_id == indP)
    dfSingPlot <- dfPlot[indB, ]
    strPrefix <- sub("_", "-", indP)
    outName <- file.path(dirOutRep,
                         paste0(sub(":", "-", strPrefix),
                                ".pdf"))
    strTitle <- dfSingPlot$Var_type
    strSubTitle <- paste0(sub("_", ":", dfSingPlot$Var_id),
                          " - ",
                          dfSingPlot$Ref_allele[1], " / ",
                          dfSingPlot$Alt_allele[1])
    
    singPlot <- ggplot(data = dfSingPlot) +
      geom_line(aes(x = Sample_id, y = Af_val,
                    group = Var_id)) +
      theme_minimal() +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
            # panel.grid.major = element_blank(), 
            # panel.grid.minor = element_blank(), 
            # panel.border = element_blank(), 
            panel.background = element_blank(),
            # title style
            # plot.title = element_text(size = 16,
            # face = "bold", hjust = 1), 
            # plot.subtitle = element_text(size = 12,
            # face = "bold", hjust = 1), 
            # axis style
            axis.title = element_text(size = 22), 
            axis.text = element_text(size = 22)) +
      coord_cartesian(ylim = c(.0, 1.)) +
      labs(title = strTitle, subtitle = strSubTitle,
           x = "Sample",
           y = "AF",
           size = 28)
    pdf(file = outName, width = 9, height = 9)
    print(singPlot)
    dev.off()
  }
}
