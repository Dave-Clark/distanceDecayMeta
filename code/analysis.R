library(metagear)

# detect web of science search results
dataFiles <- list.files(pattern = "search[0-9]")

# read in data
wosData <- lapply(dataFiles, read.delim, header = T)

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
library(patchwork)
remotes::install_github("erocoar/ggparl")
library(ggparl)
library(plotrix)
library(DescTools)

# read in extracted data
acceptDat <- fread("extractedData.csv")

# z transform correlation coefficients
acceptDat[, mantelZ := FisherZ(mantelR)]

# summary stats
acceptDat[, summary(mantelR)]
acceptDat[, std.error(mantelR)]

# get n studies and data points by year
hitsYear <- acceptDat[order(year),
  .(studies = length(unique(title)), dataPoints = .N), by = year]

# calculate cumulative totals
hitsYear[, ":="(dataPoints = cumsum(dataPoints), studies = cumsum(studies))]

# melt to long format
hitsYear <- melt(hitsYear, id = "year")

# plot studies + data points by year
hitYear <- ggplot(hitsYear, aes(x = year, y = value, group = variable)) +
  geom_line(aes(linetype = variable), size = 1.2, color = "grey") +
  geom_point(aes(shape = variable), size = 4, col = "grey", alpha = 0.7) +
  theme_bw() +
  scale_shape(labels = c("studies", "data points")) +
  scale_linetype(labels = c("studies", "data points")) +
  labs(x = "Publication year", y = "Cumulative total",
    title = "B") +
  scale_x_continuous(breaks = seq(2005, 2017, 2), labels = seq(2005, 2017, 2)) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title.position = "panel",
    plot.title = element_text(size = 20),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = c(0.3, 0.9),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"))

# get number of studies per journal
hitsJourn <- acceptDat[, .N, by = journal][order(N), ]

# plot studies per journal (top 15)
journPlot <- ggplot(hitsJourn[N >= 4],
    aes(y = factor(journal, levels = journal), x = N)) +
  geom_bar(stat = "identity", fill = "grey", col = "grey") +
  theme_bw() +
  labs(x = "Number of distance-decay relationships", title = "A") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 51)) +
  theme(axis.text = element_text(size = 16),
    plot.title.position = "panel",
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

metaPanel <- journPlot + hitYear

ggsave("../graphics/Figure_1.pdf", metaPanel, height = 5.5, width = 14,
  device = "pdf")

# categorise methods into fingerprinting, HTS, or other
acceptDat[, resolution := "Intermediate"]
acceptDat[method %in%
  c("sanger", "illumina", "pyrosequencing", "Pac-Bio", "Ion Torrent"),
  resolution := "High"]
acceptDat[method %in% c("morphology"), resolution := "Low"]

####################### analyses of context related variales ###################
taxon <- aov(mantelZ ~ taxa, acceptDat)
summary(taxon)

# test for biome effect
biome <- aov(mantelZ ~ biome, acceptDat[, if(.N > 3) .SD, by = biome])
summary(biome)

# conduct tukey test
biomeTukey <- setDT(as.data.frame(TukeyHSD(biome)$biome), keep.rownames = T)
names(biomeTukey)[ncol(biomeTukey)] <- "p_adj"

# get maximum p val for all sponge comparisons
biomeTukey[grep("sponge", rn), max(p_adj)]

biomeTukey[grep("grassland", rn), ][order(p_adj), ]

# order factor levels from highest to lowest mean
acceptDat[, biome := factor(biome,
  levels = acceptDat[,
    mean(mantelZ), by = biome][order(V1, decreasing=T), biome])]

material <- aov(mantelZ ~ medium, acceptDat[, if(.N > 2) .SD, by = biome])
materialTukey <- as.data.table(TukeyHSD(material)$medium, keep.rownames = T)
names(materialTukey)[ncol(materialTukey)] <- "p_adj"

# order factor levels from highest to lowest mean
acceptDat[, medium := factor(medium,
  levels = acceptDat[,
    mean(mantelZ), by = medium][order(V1, decreasing=T), medium])]

biomePlot <- ggplot(acceptDat[, if(.N > 3) .SD, by = biome],
    aes(x = biome, y = mantelZ)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  geom_boxjitter(width = 0.6, outlier.shape = NA, jitter.alpha = 0.5) +
  labs(x = "Biome", y = expression(Effect~size~(italic(Z[r]))), title = "A") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

habitatPlot <- ggplot(acceptDat[, if(.N > 2) .SD, by = medium],
    aes(x = medium, y = mantelZ)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  geom_boxjitter(width = 0.6, outlier.shape = NA, jitter.alpha = 0.5) +
  labs(x = "Habitat type", y = expression(Effect~size~(italic(Z[r]))),
    title = "B") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(size = 20),
    panel.grid = element_blank())

envPanel <- biomePlot + habitatPlot

ggsave("../graphics/Figure_2.pdf", envPanel, height = 5, width = 9,
  device = "pdf")

scaleLm <- lm(mantelZ ~ log10(scale), acceptDat)
summary(scaleLm)

scalePlot <- ggplot(acceptDat, aes(x = scale, y = mantelZ)) +
  geom_point(size = 3, alpha = 0.8, colour = "grey") +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) +
  stat_smooth(method = "lm", se = F, colour = "black") +
  labs(x = "Spatial extent (km)", y = expression(Effect~size~(italic(Z[r])))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

ggsave("../graphics/Figure_3.pdf", scalePlot, height = 4, width = 4.5,
  device = "pdf")

# test for correlation between scale and sampling effort
acceptDat[, cor.test(nSamples, scale)]

########################### Methodological differences #########################
resolution <- aov(mantelZ ~ resolution, acceptDat)
summary(resolution)
TukeyHSD(resolution)

coverage <- lm(mantelZ ~ log10(seqDepth), acceptDat)
summary(coverage)

# relabel factor levels for plot
plotMethods <- unique(acceptDat[! method %in% c("DGGE", "TRFLP"),
  if(.N > 10) .SD, by = method][, method])

acceptDat[, methodLab := ifelse(method %in% plotMethods, method, "Other")]
acceptDat[, methodLab := factor(methodLab,
  levels = c("illumina", "pyrosequencing", "sanger", "morphology", "Other"),
  labels = c("Illumina MiSeq/HiSeq", "454 Pyrosequencing", "Sanger sequencing",
    "Morphology/microscopy", "Other"))]

depthPlot <- ggplot(acceptDat,
    aes(x = seqDepth, y = mantelZ, col = methodLab)) +
  geom_point(size = 3, alpha = 0.6) +
  stat_smooth(method = "lm", se = F, colour = "black") +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000)) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Community coverage\n(sequences/individuals per sample)",
    y = expression(Effect~size~(italic(Z[r])))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../graphics/Figure_4.pdf", depthPlot, height = 4.5, width = 7.5,
  device = "pdf")

# test whether HTS produces significantly higher mantel coefs than other methods
aov1 <- aov(mantelZ ~ resolution, acceptDat)
summary(aov1)

# same test but only for "significant" mantel coefs
aov2 <- aov(mantelR ~ resolution, acceptDat[pValue <= 0.05])
summary(aov2) # Approaching significance

# does adding sample effort improve model
lm3 <- lm(mantelR ~ log(nSamples) + log(seqDepth), acceptDat)
summary(lm3)

# effect of sample effort
lm4 <- lm(mantelR ~ log(nSamples), acceptDat)
summary(lm4) # NS

# no real improvement
lm5 <- lm(mantelR ~ log(nSamples) * log(seqDepth), acceptDat)

# test ation between seq depth and sample effort
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

indType <- ggplot(simData, aes(x = indType, y = mantelZ)) +
  geom_boxjitter(width = 0.7, outlier.shape = NA) +
  labs(x = "Similarity index type", y = expression(Effect~size~(italic(Z[r]))),
    title = "A") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45,  hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20),
    plot.title.position = "plot",
    panel.grid = element_blank())

indIdent <- ggplot(simData, aes(x = simIndex, y = mantelZ)) +
  geom_boxjitter(width = 0.6, outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~indType, scales = "free_x") +
  labs(x = "Similarity index", y = expression(Effect~size~(italic(Z[r]))),
    title = "B") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45,  hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20),
    plot.title.position = "plot",
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

indexPanel <- indType + indIdent + plot_layout(widths = c(0.35, 1))

ggsave("../graphics/Figure_5.pdf", indexPanel, height = 4.5, width = 10,
  device = cairo_pdf)
