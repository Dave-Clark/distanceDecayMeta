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
# re-look at biome classifications
# weighted vs unweighted unifrac


library(data.table)
library(ggplot2)
library(patchwork)
# remotes::install_github("erocoar/ggparl")
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
  labs(x = "Publication year", y = "Cumulative total") +
  scale_x_continuous(breaks = seq(2005, 2019, 2), labels = seq(2005, 2019, 2)) +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = c(0.3, 0.9),
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"))

# get number of studies per journal
hitsJourn <- acceptDat[, .N, by = journal][order(N), ]

# plot studies per journal (top 15)
journPlot <- ggplot(hitsJourn[N >= 5],
    aes(y = factor(journal, levels = journal), x = N)) +
  geom_bar(stat = "identity", fill = "grey", col = "grey") +
  theme_bw() +
  labs(x = "Number of distance-decay relationships") +
  theme(axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

metaPanel <- journPlot + hitYear + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20))

ggsave("../graphics/Figure_1.pdf", metaPanel, height = 5.5, width = 14,
  device = "pdf")

# categorise methods into fingerprinting, HTS, or other
acceptDat[, resolution := "Intermediate"]
acceptDat[method %in%
  c("sanger", "illumina", "pyrosequencing", "Pac-Bio", "Ion Torrent"),
  resolution := "High"]
acceptDat[method %in% c("morphology"), resolution := "Low"]

####################### analyses of context related variales ###################
taxon <- aov(mantelZ ~ taxa, acceptDat[, if(.N > 3) .SD, by = taxa])
summary(taxon)

acceptDat[, taxa := factor(taxa,
  levels = acceptDat[, .N, by = taxa][order(N, decreasing = T), taxa],
  labels = c("Bacteria", "Fungi", expression(mu*"-Eukarya"), "Archaea", "Bacteria/\nArchaea",  "Bacteria/\nFungi", "Bacteria/\nEukarya", "All"))]

taxonFig <- ggplot(acceptDat, aes(x = taxa, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  geom_boxjitter(jitter.alpha = 0.5, outlier.shape = NA) +
  scale_x_discrete(labels = parse(text = levels(acceptDat$taxa))) +
  labs(x = "Taxon", y = expression(Mantel[italic(r)])) +
  theme_bw() +
  theme(axis.text.x = element_text(
      size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

ggsave("../graphics/Figire_S1.pdf", taxonFig, height = 5, width = 4,
  device = "pdf")

# test for biome effect
biome <- aov(mantelZ ~ biome, acceptDat[, if(.N > 3) .SD, by = biome])
summary(biome)

# conduct tukey test
biomeTukey <- setDT(as.data.frame(TukeyHSD(biome)$biome), keep.rownames = T)
names(biomeTukey)[ncol(biomeTukey)] <- "p_adj"

# order factor levels from highest to lowest mean
acceptDat[, biome := factor(biome,
  levels = acceptDat[,
    mean(mantelR), by = biome][order(V1, decreasing=T), biome])]

# make dummy variable to remove interactions with small sample sizes
acceptDat[, biome_medium := paste(biome, medium, sep = "_")]
material <- aov(mantelZ ~ biome * medium,
  acceptDat[, if(.N > 3) .SD, by = biome_medium])

materialTukey <- setDT(
  as.data.frame(TukeyHSD(material)[[2]]), keep.rownames = T)
biome_mediumTukey <- setDT(
  as.data.frame(TukeyHSD(material)[[3]]), keep.rownames = T)
names(biome_mediumTukey)[ncol(biome_mediumTukey)] <- "p_adj"
biome_mediumTukey <- biome_mediumTukey[!is.na(p_adj)]

envPlot <- ggplot(acceptDat[, if(.N > 3) .SD, by = biome_medium],
    aes(x = biome, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  geom_boxjitter(jitter.alpha = 0.5, outlier.shape = NA) +
  labs(x = "", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text.x = element_text(
      size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

medPlot <- ggplot(acceptDat[, if(.N > 3) .SD, by = biome_medium],
    aes(x = medium, y = mantelR, col = medium)) +
  geom_boxjitter(jitter.alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ biome, scales = "free_x") +
  labs(x = "", y = expression(Mantel[r]), col = "Habitat") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 18),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14))

envPanel <- envPlot + medPlot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20))

ggsave("../graphics/Figure_2.pdf", envPanel, height = 6, width = 12,
  device = "pdf")

scaleLm <- lm(mantelZ ~ log10(scale), acceptDat)
summary(scaleLm)

scalePlot <- ggplot(acceptDat, aes(x = scale, y = mantelR)) +
  geom_point(size = 3, alpha = 0.8, colour = "grey") +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) +
  stat_smooth(method = "lm", se = T, colour = "black") +
  labs(x = "Spatial extent (km)", y = expression(Mantel[r])) +
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

method <- aov(mantelZ ~ method, acceptDat[, if(.N > 4) .SD, by = method])

# relabel factor levels for plot
plotMethods <- unique(acceptDat[, if(.N > 4) .SD, by = method][, method])

acceptDat[, methodLab := ifelse(method %in% plotMethods, method, "Other")]
acceptDat[, methodLab := factor(methodLab,
  levels = c("illumina", "pyrosequencing", "TRFLP", "ARISA", "DGGE", "sanger", "Other", "morphology", "plfa"),
  labels = c("Illumina MiSeq/HiSeq", "454 Pyrosequencing", "TRFLP", "ARISA", "DGGE", "Sanger", "Other", "Morphology/\nmicroscopy", "PLFA"))]

methodPlot <- ggplot(acceptDat[, if(.N > 4) .SD, by = method],
  aes(x = methodLab, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  geom_boxjitter(jitter.alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ resolution, scales = "free_x") +
  labs(x = "Molecular method", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text.x = element_text(
      size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 14))

ggsave("../graphics/Figure_S2.pdf", methodPlot, height = 5, width = 8,
  device = "pdf")

acceptDat[, seqDepth := as.numeric(seqDepth)]
coverage <- lm(mantelZ ~ log10(seqDepth), acceptDat)
summary(coverage)

sampleDepth <- lm(mantelZ ~ log10(nSamples), acceptDat)
summary(sampleDepth)


depthPlot <- ggplot(acceptDat[!is.na(seqDepth)],
    aes(x = seqDepth, y = mantelR, col = methodLab)) +
  geom_point(size = 3, alpha = 0.6) +
  stat_smooth(method = "lm", se = T, colour = "black") +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Community coverage\n(sequences/individuals per sample)",
    y = expression(Mantel[r])) +
  scale_y_continuous(breaks = seq(-0.25, .75, 0.25)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank())

samplePlot <- ggplot(acceptDat,
    aes(x = nSamples, y = mantelR, col = methodLab)) +
  geom_point(size = 3, alpha = 0.5) +
  stat_smooth(method = "lm", se = T, colour = "black", linetype = 2) +
  scale_x_log10(breaks = c(10, 100, 1000)) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Number of samples",
    y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank())

samplePanel <- depthPlot + samplePlot + plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20))

ggsave("../graphics/Figure_4.pdf", samplePanel, height = 8, width = 7,
  device = "pdf")


# test whether indexes produce different results
# exclude those with 3 or fewer occurences
# perform anova test
aov5 <- aov(mantelR ~ simIndex, acceptDat[, if(.N > 3) .SD, by = simIndex])
summary(aov5) # sig differences found

# perform tukey test
aov5Tukey <- setDT(as.data.frame(TukeyHSD(aov5)$simIndex), keep.rownames = T)
setnames(aov5Tukey, old = "p adj", new = "p_adj")

# get comparisons that are significantly different
aov5Tukey[p_adj <= 0.05]

# divide into different index types "phylogenetic", "abundance" or "binary"
# phylo = unifrac, betamntd betampd, rao
# abund = bray. hornmorista, euclid, hellinger, theta
# binary = jaccard, Raup-Crick, sorensen, simpson, beta sim
acceptDat[, indType := "Binary"]
acceptDat[
  simIndex %in% c("bray", "Morisita-horn", "Horn-morisita", "euclidean", "hellinger", "theta", "beta_bray", "canberra", "nes_bray"), indType := "Abundance"]
acceptDat[simIndex %in% c("w_unifrac", "u_unifrac", "betaMNTD", "betaMPD", "rao", "Mash"),
  indType := "Phylogenetic"]

indexAov <- aov(mantelZ ~ indType/simIndex,
  acceptDat[, if(.N > 3) .SD, by = simIndex])
summary(indexAov)

indexTukey <- TukeyHSD(indexAov)
indexTukey <- setDT(as.data.frame(indexTukey[[2]]), keep.rownames = T)
setnames(indexTukey, old = "p adj", new = "p_adj")
indexTukey <- indexTukey[!is.na(p_adj)]


# subet indices with more than 3 coefficients
simData <- acceptDat[, if(.N > 3) .SD, by = simIndex]

# reorder the factor levels
simData$simIndex <- factor(simData$simIndex,
  levels = c("bray", "euclidean", "hellinger", "jaccard", "Raup-Crick",
    "sorensen", "beta_sim", "betaMNTD", "u_unifrac", "w_unifrac"),
  labels = c("Bray-Curtis", "Euclidean", "Hellinger", "Jaccard", "Raup-Crick",
    "S\U00F8rensen", "\U03B2-sim", "\U03B2-MNTD", "unweighted Unifrac", "weighted Unifrac"))

indType <- ggplot(simData, aes(x = indType, y = mantelR)) +
  geom_boxjitter(width = 0.7, outlier.shape = NA, jitter.alpha = 0.5) +
  labs(x = "Similarity index type", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45,  hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

indIdent <- ggplot(simData, aes(x = simIndex, y = mantelR)) +
  geom_boxjitter(width = 0.6, outlier.shape = NA, jitter.alpha = 0.5) +
  facet_wrap(~indType, scales = "free_x") +
  labs(x = "Similarity index", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45,  hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

indexPanel <- indType + indIdent + plot_layout(widths = c(0.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20))

ggsave("../graphics/Figure_5.pdf", indexPanel, height = 5, width = 10.5,
  device = cairo_pdf)

######################################### Model comparison #####################
allVars <- c("method", "nSamples", "seqDepth", "simIndex", "taxa", "medium",
  "biome", "scale", "resolution", "indType")

subData <- acceptDat[complete.cases(acceptDat[, .SD, .SDcols = allVars]), ]

ecoModel <- lm(mantelZ ~ log10(scale) + taxa + biome + medium, subData)
methModel <- lm(mantelZ ~ method + log(nSamples) + log10(seqDepth) + simIndex, subData)
null <- lm(mantelZ ~ 1, subData)

AIC(ecoModel, methModel)
summary(ecoModel)
summary(methModel)

anova(null, ecoModel)
anova(methModel, null)
