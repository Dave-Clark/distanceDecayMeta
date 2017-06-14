library(metagear)

# list WoS search results
searchFiles <- list.files(pattern = "set[0-9].txt")

cores <- detectCores()

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
library(fossil)

s1Coords <- data.frame(long = c(57.3, 58.6, 59.1, 60.03, 61.5, 60.5, 57.5, 56.5, 55.75, 55.01, 54.99, 53.99, 54, 54), lat = c(62.67, 63.01, 63.24, 63.25, 63, 61, 60.75, 60.99, 60, 60.75, 61.75, 60, 60.99, 62.25))

s1Coords <- s1Coords * -1

s1Dist <- earth.dist(s1Coords)
