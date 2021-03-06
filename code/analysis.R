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
# re-look at environment classifications
# weighted vs unweighted unifrac


library(data.table)
library(lme4)
library(ggplot2)
library(patchwork)
# remotes::install_github("erocoar/ggparl")
library(ggparl)
library(paletteer)
library(plotrix)
library(DescTools)

# read in extracted data
acceptDat <- fread("extractedData.csv")

# z transform correlation coefficients
acceptDat[, mantelZ := FisherZ(mantelR)]

# rough funnel plot
funPlot <- ggplot(acceptDat, aes(x = mantelR, y = nSamples)) +
  geom_vline(xintercept = mean(acceptDat$mantelR), linetype = 2) +
  geom_point(alpha = 0.7) +
  labs(x = expression(Mantel[italic(r)]), y = "Sample size") +
  spatialExtent_x_continuous(limits = c(-1, 1)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

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
  spatialExtent_shape(labels = c("studies", "data points")) +
  spatialExtent_linetype(labels = c("studies", "data points")) +
  labs(x = "Publication year", y = "Cumulative total") +
  spatialExtent_x_continuous(breaks = seq(2005, 2019, 2), labels = seq(2005, 2019, 2)) +
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

hitsJourn[journal == "Science of the total environment",
  journal := "Science of The Total Environment"]

hitsJourn[journal == "ISME",
  journal := "ISME J"]

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

ggsave("../graphics/Figure_S1.pdf", metaPanel, height = 5.5, width = 14,
  device = "pdf")

# categorise methods into fingerprinting, HTS, or other
acceptDat[, resolution := "Intermediate"]
acceptDat[method %in%
  c("sanger", "illumina", "pyrosequencing", "Pac-Bio", "Ion Torrent"),
  resolution := "High"]
acceptDat[method %in% c("morphology"), resolution := "Low"]

####################### analyses of context related variales ###################
taxonLm <- lm(mantelZ ~ taxa, acceptDat[, if(.N > 3) .SD, by = taxa])

# re run as aov for post hoc analysis
taxon <- aov(mantelZ ~ taxa, acceptDat[, if(.N > 3) .SD, by = taxa])
summary(taxon)

acceptDat[, taxa := factor(taxa,
  levels = acceptDat[, .N, by = taxa][order(N, decreasing = T), taxa],
  labels = c("Bacteria", "Fungi", " Other Eukarya", "Archaea", "Bacteria/\nArchaea",  "Bacteria/\nFungi", "Bacteria/\nEukarya", "All"))]

taxonFig <- ggplot(acceptDat, aes(x = taxa, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey") +
  geom_boxjitter(jitter.alpha = 0.5, outlier.shape = NA) +
  labs(x = "Taxon", y = expression(Mantel[italic(r)])) +
  theme_bw() +
  theme(axis.text.x = element_text(
      size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

ggsave("../graphics/Figure_1.pdf", taxonFig, height = 5, width = 5,
  device = "pdf")

# test for environment effect
environmentLm <- lm(mantelZ ~ environment, acceptDat[, if(.N > 3) .SD, by = environment])
summary(environment)

# run as anova for post hoc analysis
environment <- aov(mantelZ ~ environment, acceptDat[, if(.N > 3) .SD, by = environment])

# conduct tukey test
environmentTukey <- setDT(as.data.frame(TukeyHSD(environment)$environment), keep.rownames = T)
names(environmentTukey)[ncol(environmentTukey)] <- "p_adj"

# order factor levels from highest to lowest mean
acceptDat[, environment := factor(environment,
  levels = acceptDat[,
    mean(mantelR), by = environment][order(V1, decreasing=T), environment])]

# make dummy variable to remove interactions with small sample sizes
acceptDat[, environment_habitat := paste(environment, habitat, sep = "_")]

# run as linear model
materialLm <- lm(mantelZ ~ environment * habitat,
  acceptDat[, if(.N > 3) .SD, by = environment_habitat])

# run as anova for post-hoc analysis
material <- aov(mantelZ ~ environment * habitat,
  acceptDat[, if(.N > 3) .SD, by = environment_habitat])

materialTukey <- setDT(
  as.data.frame(TukeyHSD(material)[[2]]), keep.rownames = T)
environment_habitatTukey <- setDT(
  as.data.frame(TukeyHSD(material)[[3]]), keep.rownames = T)
names(environment_habitatTukey)[ncol(environment_habitatTukey)] <- "p_adj"
environment_habitatTukey <- environment_habitatTukey[!is.na(p_adj)]

# create custom simpsons palette
simpsonFullPal <- paletteer_d("ggsci::springfield_simpsons")
simpsonPal <- simpsonFullPal[c(1:3, 5, 9, 11, 14:16)]

envPlot <- ggplot(acceptDat[, if(.N > 3) .SD, by = environment_habitat],
    aes(x = environment, y = mantelR)) +
  geom_hline(yintercept = 0, linetype = 2, col = "grey") +
  geom_boxjitter(jitter.alpha = 0.5, outlier.shape = NA) +
  labs(x = "", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text.x = element_text(
      size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

medPlot <- ggplot(acceptDat[, if(.N > 3) .SD, by = environment_habitat],
    aes(x = habitat, y = mantelR, col = habitat)) +
  geom_boxjitter(jitter.alpha = 0.7, outlier.shape = NA, alpha = 1) +
  facet_wrap(~ environment, spatialExtents = "free_x") +
  spatialExtent_color_manual(values = simpsonPal) +
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

lakeLm <- lm(mantelZ ~ within_lake, acceptDat[!is.na(within_lake)])
anova(lakeLm)

acceptDat[, within_lake := fcase(within_lake == "within", "Within lake",
  within_lake == "across", "Across lakes")]

# test comparison of within-lake and between-lake studies
lakePlot <- ggplot(acceptDat[!is.na(within_lake)],
  aes(x = within_lake, y = mantelR)) +
  geom_boxjitter(jitter.alpha = 0.7, outlier.shape = NA, alpha = 1) +
  labs(x = "", y = expression(Mantel[r]), col = "Habitat") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

ggsave("../graphics/Figure_S2.pdf", lakePlot, height = 4, width = 4,
  device = "pdf")

spatialExtentLm <- lm(mantelZ ~ log10(spatialExtent), acceptDat)
summary(spatialExtentLm)

#### NEED RANDOM INTERCEPT FOR STUDY EFFECT ####
spatialExtentLmer <- lmer(mantelZ ~ (1 + log10(spatialExtent) | title) + log10(spatialExtent),
  acceptDat, REML = F)
spatialExtentLmer2 <- lmer(mantelZ ~ (1 | title) + log10(spatialExtent),
  acceptDat, REML = F)
spatialExtentNull <- lmer(mantelZ ~ (1 | title), acceptDat, REML = F)
AIC(spatialExtentLmer2, spatialExtentLmer) # use simpler random intercept only model
#
MuMIn::r.squaredGLMM(spatialExtentLmer2)

# get p value for spatialExtent
coefs <- data.frame(coef(summary(spatialExtentLmer2)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))

# using mantelR instead of Z
spatialExtentLmer2R <- lmer(mantelR ~ (1 | title) + log10(spatialExtent),
  acceptDat, REML = F)

spatialExtentPred <- data.table(spatialExtent = seq(
  min(acceptDat$spatialExtent, na.rm = T),
  max(acceptDat$spatialExtent, na.rm = T),
  100))

spatialExtentPred$prediction <- predict(spatialExtentLmer2R, newdata = spatialExtentPred, re.form = NA)

spatialExtentPlot <- ggplot(acceptDat,
    aes(x = spatialExtent * 1000, y = mantelR)) +
  geom_point(size = 3, alpha = 0.8, col = "grey") +
  geom_line(data = spatialExtentPred, aes(x = spatialExtent * 1000, y = prediction),
    linetype = 2, size = 1.1) +
  spatialExtent_x_log10(
    breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  labs(x = "Spatial extent (m)", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

ggsave("../graphics/Figure_3.pdf", spatialExtentPlot, height = 4, width = 5,
  device = "pdf")

# test for correlation between spatialExtent and sampling effort
acceptDat[, cor.test(nSamples, spatialExtent)]

# plot spatialExtent against environment and Habitat
spatialExtentEnv <- ggplot(acceptDat, aes(x = environment, y = spatialExtent * 1000)) +
  geom_boxplot() +
  spatialExtent_y_log10() +
  labs(x = "Environment", y = "Spatial extent (m)") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

spatialExtentHabitat <- ggplot(acceptDat[!is.na(habitat)],
    aes(x = habitat, y = spatialExtent * 1000)) +
  geom_boxplot() +
  spatialExtent_y_log10() +
  labs(x = "Habitat", y = "Spatial extent (m)") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    panel.grid = element_blank())

spatialExtentPanel <- spatialExtentEnv + spatialExtentHabitat +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20))

ggsave("../graphics/Figure_S3.pdf", spatialExtentPanel, height = 10, width = 6,
  device = "pdf")

########################### Methodological differences #########################
resolutionLm <- lm(mantelZ ~ resolution, acceptDat)

resolution <- aov(mantelZ ~ resolution, acceptDat)
summary(resolution)
TukeyHSD(resolution)

methodLm <- lm(mantelZ ~ method, acceptDat[, if(.N > 4) .SD, by = method])
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
  facet_wrap(~ resolution, spatialExtents = "free_x") +
  labs(x = "Molecular method", y = expression(Mantel[r])) +
  theme_bw() +
  theme(axis.text.x = element_text(
      size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 14))

ggsave("../graphics/Figure_S4.pdf", methodPlot, height = 5, width = 8,
  device = "pdf")

acceptDat[, seqDepth := as.numeric(seqDepth)]
coverage <- lmer(mantelZ ~ (1 | title) + log10(seqDepth), acceptDat)
summary(coverage)

# get p value for coverage
coverageCoefs <- data.frame(coef(summary(coverage)))
# use normal distribution to approximate p-value
coverageCoefs$p.z <- 2 * (1 - pnorm(abs(coverageCoefs$t.value)))

EnvStats::rosnerTest(na.omit(log10(acceptDat$seqDepth)), k = 2)

coverage2 <- lmer(mantelZ ~ (1 | title) + log10(seqDepth),
  acceptDat[log10(seqDepth) < 7.5, ])
# get p value for coverage
coverageCoefs2 <- data.frame(coef(summary(coverage2)))
# use normal distribution to approximate p-value
coverageCoefs2$p.z <- 2 * (1 - pnorm(abs(coverageCoefs2$t.value)))

# generate predictions
coveragePred <- data.table(seqDepth = seq(
  min(na.omit(acceptDat$seqDepth)), max(na.omit(acceptDat$seqDepth)), 1000))

coveragePred[, prediction := predict(coverage, newdata = coveragePred,
  re.form = NA)]

sampleDepth <- lmer(mantelZ ~ (1 | title) + log10(nSamples), acceptDat)
summary(sampleDepth)

# get p value for sample depth
depthCoefs <- data.frame(coef(summary(sampleDepth)))
# use normal distribution to approximate p-value
depthCoefs$p.z <- 2 * (1 - pnorm(abs(depthCoefs$t.value)))

depthPlot <- ggplot(data = acceptDat[!is.na(seqDepth)],
    aes(x = seqDepth, y = mantelR, col = methodLab)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(data = coveragePred, aes(x = seqDepth, y = prediction),
    linetype = 1, size = 1.1, col  = "black") +
  spatialExtent_x_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  spatialExtent_color_manual(values = simpsonPal) +
  labs(x = "Community coverage\n(sequences/individuals per sample)",
    y = expression(Mantel[r])) +
  spatialExtent_y_continuous(breaks = seq(-0.25, .75, 0.25)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank())

samplePlot <- ggplot(acceptDat,
    aes(x = nSamples, y = mantelR, col = methodLab)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_smooth(method = "lm", se = T, colour = "black", linetype = 2) +
  spatialExtent_x_log10(breaks = c(10, 100, 1000)) +
  spatialExtent_color_manual(values = simpsonPal) +
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


# test whether indices produce different results
# exclude those with 3 or fewer occurences
# perform anova test
simIndexLm <- lm(mantelR ~ simIndex, acceptDat[, if(.N > 3) .SD, by = simIndex])

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
  facet_wrap(~indType, spatialExtents = "free_x") +
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

# Reviewers suggested analysing differences between correlation types used
# e.g. Pearson vs Spearman
corType <- lm(mantelZ ~ correlation, data = acceptDat)
summary(corType)

######################################### Model comparison #####################
allVars <- c("method", "nSamples", "seqDepth", "simIndex", "taxa", "habitat",
  "environment", "spatialExtent", "resolution", "indType")

subData <- acceptDat[complete.cases(acceptDat[, .SD, .SDcols = allVars]), ]

ecoModel <- lm(mantelZ ~ log10(spatialExtent) + taxa + environment + habitat, subData)
methModel <- lm(mantelZ ~ method + log(nSamples) + log10(seqDepth) + simIndex, subData)
null <- lm(mantelZ ~ 1, subData)

AIC(ecoModel, methModel)
summary(ecoModel)
summary(methModel)

anova(null, ecoModel)
anova(methModel, null)
