library(metagear)
library(stringr)
library(parallel)
library(data.table)

# list WoS search results
searchFiles <- list.files(pattern = "search[0-9]*.txt")

parallel::cores <- detectCores()

searchList <- mclapply(searchFiles, read.delim, header = T, sep = "\t", quote = NULL, mc.cores = cores) # read data in parallel

searchData <- do.call("rbind", searchList) # merge into one df

# need authors, title, abstract, DOI, Journal and date of publication
litData <- searchData[, c(2, 10, 18, 29, 33, 34)]

colnames(litData) <- c("authors", "title", "journal", "DOI", "year", "abstract") # change column headings to sensible ones

# add reviewer column and review progress column
litData$reviewer <- "me"

litData$screened <- "not yet"

litData <- litData[unique(litData$DOI), ] # dereplicate by DOI
# write to file
write.csv(litData, "preScreened.csv", row.names = F)

# clear workspace
rm(list = ls())
gc()

abstract_screener("preScreened.csv", aReviewer = "me", reviewerColumnName = "reviewer", unscreenedColumnName = "screened", unscreenedValue = "not yet", abstractColumnName = "abstract", titleColumnName = "title")

##### MANUALLY SCREEN ABSTRACTS ######
## read search results df back in and split into yes/no/maybe
dat <- read.csv("preScreened.csv", header = T)

confirmedStudies <- dat[ dat$screened == "YES", ]
maybeStudies <- dat[ dat$screened == "MAYBE", ]

dir.create("pdfDwnlds") # create dir for pdf downloads

confirmedStudies$filename <- paste("study", 1:nrow(confirmedStudies), sep = "") # create column of filenames

confPFD <- PDFs_collect(confirmedStudies, DOIcolumn = "DOI", FileNamecolumn = "filename", directory = "./pdfDwnlds/", validatePDF = T, quiet = T, showSummary = T) # download pdfs, seems to hang at end!!

dir.create("maybePapers")

maybeStudies$filename <- paste("study", 1:nrow(maybeStudies), sep = "") # create column of filenames

maybePDFs <- PDFs_collect(maybeStudies, DOIcolumn = "DOI", FileNamecolumn = "filename", directory = "./maybePapers/", validatePDF = T, quiet = T, showSummary = T)

####################################### MANUAL SCREENING #######################################

#### UPDATE TO 2017-2019 papers ####
### manually remove unnecessary cols ####
searchFiles <- list.files(pattern = "*update.txt")

searchList <- lapply(searchFiles, fread) # read data in parallel

searchData <- rbindlist(searchList) # merge into one df

colnames(searchData) <- c(
  "authors", "title", "journal", "abstract", "year", "vol", "DOI")

# add reviewer column and review progress column
searchData$reviewer <- "me"

searchData$screened <- "not yet"

searchData <- searchData[!duplicated(DOI) & year != 2020, ] # dereplicate by DOI

prevData <- fread("allWos.csv")

searchData <- searchData[!DOI %in% prevData$DOI, ]

# write to file
write.csv(searchData, "preScreened_update.csv", row.names = F)

#### manual screening of data update
abstract_screener("preScreened_update.csv", aReviewer = "me", reviewerColumnName = "reviewer", unscreenedColumnName = "screened", unscreenedValue = "not yet", abstractColumnName = "abstract", titleColumnName = "title")

# write studies to use to file, and add to extractedData file
updateData <- fread("preScreened_update.csv")
finalData <- updateData[screened %in% c("YES", "maybe")]

fwrite(finalData, "accept_studies_update.txt")
