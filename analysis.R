library(metagear)
library(data.table)
library(ggplot2)
library(cowplot)

# detect web of science search results
dataFiles <- list.files(pattern = "search[0-9]")

# read in data
wosData <- lapply(dataFiles, read.delim, header = T)

# check df dims
# lapply(wosData, dim)

# concatenate into one df
allWos <- do.call("rbind", wosData)

# dereplicate search results by title
allWos <- allWos[!duplicated(allWos$TI), ]

# rename other cols for compliance with metagear
colnames(allWos) <- c("authors", "title", "journal", "abstract", "year", "volume", "DOI")

# add reviewer and screened cols
allWos$screened <- "NO"
allWos$reviewer <- "Dave"

# write data to file
write.csv(allWos, "allWos.csv", quote = F, row.names = F)

# launch abstract screener to manually curate abstracts
abstract_screener("allWos.csv", aReviewer = "Dave",
  reviewerColumnName = "reviewer", unscreenedColumnName = "screened",
  unscreenedValue = "not_screened", abstractColumnName = "abstract",
  titleColumnName = "title")

# prelim results
dat <- read.csv("allWos.csv")
accept <- dat[!dat$screened == "NO", ]
write.table(accept, "accept_studies.txt", row.names = F, sep = "\t", quote = F)

# read in extracted data
acceptDat <- fread("extractedData.csv")

# create directory for figures
dir.create("graphics/")

# get n studies and data points by year
hitsYear <- acceptDat[order(year),
  .(studies = length(unique(title)), dataPoints = .N), by = year]

# melt to long format
hitsYear <- melt(hitsYear, id = "year")

# plot studies + data points by year
hitYear <- ggplot(hitsYear,
    aes(x = year, y = cumsum(value), group = variable)) +
  geom_line(aes(linetype = variable), size = 1.2, color = "lightgrey") +
  geom_point(aes(shape = variable), size = 4) +
  theme_bw() +
  scale_shape(labels = c("studies", "data points")) +
  scale_linetype(labels = c("studies", "data points")) +
  labs(x = "Publication year", y = "Cumulative total") +
  scale_x_continuous(breaks = seq(2005, 2017, 2), labels = seq(2005, 2017, 2)) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = c(0.15, 0.9),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"))

ggsave("graphics/data_year.pdf", hitYear, device = "pdf")

# get number of studies per journal
hitsJourn <- acceptDat[, .(nStudies = length(unique(title))), by = journal]

# order by n studies
hitsJourn <- hitsJourn[order(nStudies)]

# plot studies per journal (top 15)
journPlot <- ggplot(tail(hitsJourn, 15),
    aes(x = factor(journal, levels = journal), y = nStudies)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Journal", y = "Number of studies") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 15.5)) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()) +
  coord_flip()

ggsave("graphics/journal_studies.pdf", journPlot, device = "pdf", width = 10)

# categorise methods into fingerprinting, HTS, or other
acceptDat[, technique := "fingerprinting"]
acceptDat[method %in% c("illumina", "pyrosequencing", "Pac-Bio", "Ion Torrent"),
  technique := "HTS"]
acceptDat[method %in% c("morphology", "sanger"), technique := "Other"]

# add "signficance" column
acceptDat[, signficant := ifelse(
  pValue <= 0.05, "signficant", "non-signficant")]

# test whether HTS produces significantly higher mantel coefs than other methods
anova1 <- aov(mantelR ~ technique, acceptDat)
summary(anova1)

# same test but only for "significant" mantel coefs
anova2 <- aov(mantelR ~ technique, acceptDat[pValue <= 0.05])
summary(anova2) # Approaching significance

# Tukey test to explore results
TukeyHSD(anova2) # HTS almost different to fingerprinting

# bw plot to show vals
techPlot <- ggplot(acceptDat, aes(x = technique, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Method", y = expression(R[Mantel])) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

sigTech <- ggplot(acceptDat[pValue <= 0.05], aes(x = technique, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Method", y = expression(R[Mantel])) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

bwPlots <- plot_grid(techPlot, sigTech, labels = "AUTO", label_size = 16,
  align = "hv")

ggsave("graphics/tech_mantel.pdf", bwPlots, device = "pdf", width = 10,
  height = 5)

# test whether mantel R is related to sample depth
lm1 <- lm(mantelR ~ log(seqDepth), acceptDat)
summary(lm1)

# repeat test for significant mantel vals
lm2 <- lm(mantelR ~ log(seqDepth), acceptDat[pValue <= 0.05])
summary(lm2)

# plot mantel r as function of seq depth
seqDepth <- ggplot(acceptDat, aes(x = log(seqDepth), y = mantelR)) +
  geom_point(col = "grey", size = 3, alpha = 0.7) +
  stat_smooth(method = "lm", se = F, col = "black", size = 1.1) +
  scale_x_continuous(breaks = seq(2, 16, 2), labels = seq(2, 16, 2)) +
  labs(x = "log(sampling depth)", y = expression(R[Mantel])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

# test whether indexes produce different results
# exclude those with 3 or fewer occurences
# perform anova test
anova3 <- aov(mantelR ~ simIndex, acceptDat[, if(.N > 3) .SD, by = simIndex])
summary(anova3) # sig differences found

# perform Tukey test
TukeyHSD(anova3)

# bw plot of mantel R vs sim indexes
simPlot <- ggplot(acceptDat[, if(.N > 3) .SD, by = simIndex],
    aes(x = simIndex, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(
    labels = c(expression(beta~MNTD), "Bray-Curtis", "Euclidean", "Hellinger",
      "Jaccard", "Raup-Crick", "S\u00F8rensen", "Unifrac")) +
  labs(x = "Similarity index", y = expression(R[Mantel])) +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("graphics/sim_index.pdf", simPlot, device = "pdf", width = 5)

# divide into different index types "phylogenetic", "abundance" or "binary"
# phylo = unifrac, betamntd betampd, rao
# abund = bray. hornmorista, euclid, hellinger, theta
# binary = jaccard, Raup-Crick, sorensen, simpson, beta sim

# now test biological hypotheses
