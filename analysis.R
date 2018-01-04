library(metagear)

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

# for analysis
library(data.table)
library(ggplot2)
library(cowplot)
library(plotrix)

# read in extracted data
acceptDat <- fread("extractedData.csv")

# summary stats
acceptDat[, summary(mantelR)]
acceptDat[, std.error(mantelR)]

# create directory for figures
dir.create("graphics/")

# get n studies and data points by year
hitsYear <- acceptDat[order(year),
  .(studies = length(unique(title)), dataPoints = .N), by = year]

# calculate cumulative totals
hitsYear[, dataPoints := cumsum(dataPoints)]
hitsYear[, studies := cumsum(studies)]

# melt to long format
hitsYear <- melt(hitsYear, id = "year")

# plot studies + data points by year
hitYear <- ggplot(hitsYear, aes(x = year, y = value, group = variable)) +
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
    legend.position = c(0.16, 0.9),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"))

ggsave("graphics/data_year.pdf", hitYear, device = "pdf", width = 7, height = 7)

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
acceptDat[, resolution := "Intermediate"]
acceptDat[method %in%
  c("sanger", "illumina", "pyrosequencing", "Pac-Bio", "Ion Torrent"),
  resolution := "High"]
acceptDat[method %in% c("morphology"), resolution := "Low"]

# test whether HTS produces significantly higher mantel coefs than other methods
aov1 <- aov(mantelR ~ resolution, acceptDat)
summary(aov1)

# same test but only for "significant" mantel coefs
aov2 <- aov(mantelR ~ resolution, acceptDat[pValue <= 0.05])
summary(aov2) # Approaching significance

# bw plot to show vals
techPlot <- ggplot(acceptDat, aes(x = resolution, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  scale_x_discrete(labels = paste0(unique(acceptDat$resolution), "\n", " (n = ",
    acceptDat[, .N, by = resolution]$N, ")")) +
  theme_bw() +
  labs(y = expression(R[Mantel])) +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    aspect.ratio = 1)

sigTech <- ggplot(acceptDat[pValue <= 0.05], aes(x = resolution, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  scale_x_discrete(labels = paste0(unique(acceptDat$resolution), "\n", " (n = ",
    acceptDat[pValue <= 0.05, .N, by = resolution]$N, ")")) +
  theme_bw() +
  labs(y = expression(R[Mantel])) +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    aspect.ratio = 1)

bwPlots <- plot_grid(techPlot, sigTech, labels = "AUTO", label_size = 16,
  align = "hv")

ggsave("graphics/tech_mantel.pdf", bwPlots, device = "pdf", width = 8.5,
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
  labs(x = "log(Coverage)", y = expression(R[Mantel])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("graphics/seq_depth.pdf", seqDepth, device = "pdf", height = 5,
  width = 5)

# does adding sample effort improve model
lm3 <- lm(mantelR ~ log(nSamples) + log(seqDepth), acceptDat)
summary(lm3)

# compare to seqDepth alone
AIC(lm1)
AIC(lm3)

# effect of sample effort
lm4 <- lm(mantelR ~ log(nSamples), acceptDat)
summary(lm4) # NS

# no real improvement
lm5 <- lm(mantelR ~ log(nSamples) * log(seqDepth), acceptDat)

# test correlation between seq depth and sample effort
cor.test(log(acceptDat$nSamples), log(acceptDat$seqDepth))

# test whether indexes produce different results
# exclude those with 3 or fewer occurences
# perform anova test
aov5 <- aov(mantelR ~ simIndex, acceptDat)
summary(aov5) # sig differences found

# perform tukey test
aov5Tukey <- setDT(as.data.frame(TukeyHSD(aov5)$simIndex), keep.rownames = T)
names(aov5Tukey)[ncol(aov5Tukey)] <- "p_adj"

# get comparisons that are significantly different
aov5Tukey[p_adj <= 0.05]

# divide into different index types "phylogenetic", "abundance" or "binary"
# phylo = unifrac, betamntd betampd, rao
# abund = bray. hornmorista, euclid, hellinger, theta
# binary = jaccard, Raup-Crick, sorensen, simpson, beta sim
acceptDat[, indType := "Binary"]
acceptDat[
  simIndex %in% c("bray", "Morisita-horn", "euclidean", "hellinger", "theta"), indType := "Abundance"]
acceptDat[simIndex %in% c("unifrac", "betaMNTD", "betaMPD", "rao"),
  indType := "Phylogenetic"]

# test different index types
aov6 <- aov(mantelR ~ indType, acceptDat)
TukeyHSD(aov6)

# subet indices with more than 3 coefficients
simData <- acceptDat[, if(.N > 3) .SD, by = simIndex]

# reorder the factor levels
simData$simIndex <- factor(simData$simIndex,
  levels = c("bray", "euclidean", "hellinger", "jaccard", "Raup-Crick",
    "sorensen", "betaMNTD", "unifrac"),
  labels = c("Bray-Curtis", "Euclidean", "Hellinger", "Jaccard", "Raup-Crick",
    "S\U00F8rensen", "\U03B2-MNTD", "Unifrac"))

# bw plot of mantel R vs sim indexes
simPlot <- ggplot(simData, aes(x = simIndex, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  scale_x_discrete(labels = paste0(levels(simData$simIndex), "\n", "(n = ",
    simData[, .N, by = simIndex][order(simIndex)]$N, ")")) +
  geom_vline(xintercept = c(3.5, 6.5), col = "grey") +
  annotate("text", x = c(3.4, 6.4, 8.5), y = 0.95, hjust = 1, size = 4,
    label = c("Abundance", "Binary", "Phylogenetic")) +
  theme_bw() +
  labs(y = expression(R[Mantel])) +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = 18),
    axis.title.x = element_blank(),
    strip.text.x = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

# reorder factor levels to match
acceptDat$indType <- factor(acceptDat$indType,
  levels = c("Abundance", "Binary", "Phylogenetic"))

# plot different index types
indexType <- ggplot(acceptDat, aes(x = indType, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  scale_x_discrete(labels = paste0(levels(acceptDat$indType), "\n", " (n = ",
    acceptDat[, .N, by = indType][order(indType)]$N, ")")) +
  theme_bw() +
  labs(y = expression(R[Mantel])) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

# combine index and indexType plots into panel
indexPlots <- plot_grid(indexType, simPlot, labels = "AUTO", label_size = 16,
  align = "hv", rel_widths = c(1, 0.6))

# write panel plot to file
ggsave("graphics/index_panel.pdf", indexPlots, device = cairo_pdf, height = 4,
  width = 8)

### biological factors ###
# test for taxa influence
aov7 <- aov(mantelR ~ taxa, acceptDat)
summary(aov7)

# test taxa influence for only significant results
aov8 <- aov(mantelR ~ taxa, acceptDat[pValue <= 0.05])

# Tukey test to find which groups different
aov8Tukey <- setDT(as.data.frame(TukeyHSD(aov8)$taxa), keep.rownames = T)
names(aov8Tukey)[ncol(aov8Tukey)] <- "p_adj"

# get significant diff groups
aov8Tukey[p_adj <= 0.05]

# remove bac-fungi coefs and re-test
aov9 <- aov(mantelR ~ taxa, acceptDat[pValue <= 0.05 & taxa != "bac_fungi"])
summary(aov9)

# test for biome effect
aov10 <- aov(mantelR ~ biome, acceptDat[, if(.N > 3) .SD, by = biome])
summary(aov10)

# conduct tukey test
aov10Tukey <- setDT(as.data.frame(TukeyHSD(aov10)$biome), keep.rownames = T)
names(aov10Tukey)[ncol(aov10Tukey)] <- "p_adj"

# get sig different groups
aov10Tukey[p_adj <= 0.05]

# test for different env medium
aov11 <- aov(mantelR ~ medium, acceptDat)
summary(aov11)

medPlot <- ggplot(acceptDat[!is.na(medium), ], aes(x = medium, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(labels = paste0(
    unique(acceptDat[!is.na(medium), medium]), "\n", " (n = ",
    acceptDat[!is.na(medium), .N, by = medium]$N, ")")) +
  labs(y = expression(R[Mantel]), x = "Micro-habitat") +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("graphics/medium.pdf", medPlot, device = "pdf", height = 5, width = 5)

# Tukey test to find sig differences
aov11Tukey <- TukeyHSD(aov11)

# test for studies that include more than one focal taxa
multStudies <- dat[, if(length(unique(taxa)) > 1) .SD,
  .SDcols = c("taxa", "mantelR", "title"), by = title]

summary(aov(mantelR ~ taxa, multStudies))

# test for relationship with scale
lm6 <- lm(mantelR ~ log(scale), acceptDat)
summary(lm6)

#plot scale against mantelr
scalePlot <- ggplot(acceptDat, aes(x = log(scale), y = mantelR)) +
  geom_point(col = "grey", size = 3, alpha = 0.7) +
  stat_smooth(method = "lm", se = F, col = "black", size = 1.1) +
  labs(x = "log(Scale (km))", y = expression(R[Mantel])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("graphics/scale.pdf", scalePlot, device = "pdf", width = 5, height = 5)

# test if scale and sample effort are correlated
acceptDat[, cor.test(log(nSamples), log(seqDepth))]

# add sampling effort to scale
lm7 <- lm(mantelR ~ log(nSamples) + log(scale), acceptDat)
summary(lm7)
