install.packages("dplyr")
install.packages("tidyverse")
library(dplyr)
library(tidyverse)
if (!require("ßBiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(tibble)
BiocManager::install("ShortRead")
BiocManager::install("dada2")
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
library(Biostrings)
library(DECIPHER)
library(dada2)
library(ShortRead)
library(phyloseq)
library(ggplot2)
library(DESeq2)
devtools::install_github("gmteunisse/fantaxtic")
devtools::install_github("kassambara/rstatix")
library(fantaxtic)
library(edgeR)
library(stringr)
library(ggpubr)
library(DESeq2)


seqtab.nochim <- readRDS("seqtab_nochim.rds")
taxonomy.test <- read.csv("taxonomy_test.csv")

seqtab_df <- data.frame(
  full_seq = colnames(seqtab.nochim),
  prefix = substr(colnames(seqtab.nochim), 1, 50),
  stringsAsFactors = FALSE
)

rownames(taxonomy.test) <- taxonomy.test$X

taxa_df <- data.frame(
  full_seq = rownames(taxonomy.test),
  prefix = substr(rownames(taxonomy.test), 1, 50),
  stringsAsFactors = FALSE
)

common_prefixes <- intersect(seqtab_df$prefix, taxa_df$prefix)

seqtab_mapped <- seqtab_df[seqtab_df$prefix %in% common_prefixes, ]
seqtab_mapped <- seqtab_mapped[!duplicated(seqtab_mapped$prefix), ]

taxa_mapped <- taxa_df[taxa_df$prefix %in% common_prefixes,]
taxa_mapped <- taxa_mapped[!duplicated(taxa_mapped$prefix),]

seqtab_matched <- seqtab.nochim[, seqtab_mapped$full_seq]
taxa_matched <- taxonomy.test[taxa_mapped$full_seq,]

colnames(seqtab_matched) <- seqtab_mapped$prefix
rownames(taxa_matched) <- taxa_mapped$prefix
taxa_matched <- as.matrix(taxa_matched)


tax <- t(seqtab.nochim) %>%
  as.data.frame() %>%
  rownames_to_column("X") %>% inner_join(taxonomy.test, by = "X")

dna <- DNAStringSet(getSequences(seqtab.nochim))

samples.out <- rownames(seqtab.nochim)
location <- sapply(strsplit(samples.out, "-"), `[`, 1)
health <- sapply(strsplit(samples.out, "-"), `[`, 2)
fhl.num <- sapply(strsplit(samples.out, "[[:digit:]]-"), `[`, 2)
fhl.num <- sapply(strsplit(fhl.num, "_R1.fastq.gz"), `[`, 1)


sample.info <- data.frame(Location = location, Health = health, FHL = fhl.num)
rownames(sample.info) <- samples.out
sample.info <- sample.info %>% filter(stringr::str_detect(Health, "Dying|Dead|Healthy"))


samples.2health <- sample.info %>% 
  mutate(Health = case_when(
    Health == "Dying" ~ "Unhealthy",
    Health == "Dead" ~ "Unhealthy",
    TRUE ~ "Healthy"
))

ps <- phyloseq(
  otu_table(seqtab_matched, taxa_are_rows = FALSE),
  sample_data(sample.info),
  tax_table(taxa_matched)
)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
theme_set(theme_bw())

otu.table <- as.data.frame(otu_table(ps))
#taxa_matched$Count <- colSums(otu.table)

taxa_matched <- as.matrix(taxa_matched)
#taxa_matched$Count <- as.character(taxa_matched$Count)

ps.count <- phyloseq(
  otu_table(seqtab_matched, taxa_are_rows = FALSE),
  sample_data(sample.info),
  tax_table(taxa_matched)
)

ps.2health <- phyloseq(
  otu_table(seqtab_matched, taxa_are_rows = FALSE),
  sample_data(samples.2health),
  tax_table(taxa_matched)
)

tax.table.count <- psmelt(ps.2health)
write.csv(tax.table.count, "tax_table_count.csv")

topbyhealth <- tax.table.count %>%
  filter(Abundance > 0) %>%
  group_by(Health) %>%
  summarize(Taxa = list(unique(Genus)))

topbyhealth.species <- tax.table.count %>%
  filter(Abundance > 0) %>%
  group_by(Health) %>%
  summarize(Taxa = list(unique(Species)))

tax.health <- nested_top_taxa(ps_obj = ps.2health,
                              top_tax_level = "Genus",
                              nested_tax_level = "Species",
                              n_top_taxa = 20,
                              n_nested_taxa = 20,
                              include_na_taxa = F,
                              grouping = "Health",
                              by_proportion = T)

top.taxa <- tax.health$top_taxa

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "Health", title = "Bray NMDS") 

ps.sub <- subset_samples(ps, !is.na(Health))
expt <- prune_taxa(names(sort(taxa_sums(ps.sub), TRUE)[1:50]), ps.sub)
ord <- ordinate(expt, formula = ~Health, "NMDS", "bray")
ordplot <- plot_ordination(expt, ord, "samples", color = "Health", shape = "Health")
ordplot + stat_ellipse(geom = "polygon", type = "norm", linetype = 2, alpha = 0.2, aes(fill=Health)) +
  stat_ellipse(type = "t", level = 0.95) + theme_bw()

#NMDS overall taxa
ps.2sub <- subset_samples(ps.2health, !is.na(Health))
expt2 <- prune_taxa(names(sort(taxa_sums(ps.2sub), TRUE)[1:50]), ps.2sub)
ord2 <- ordinate(expt2, formula = ~Health, "NMDS", "bray")
ordplot2 <- plot_ordination(expt2, ord2, "samples", color = "Health", shape = "Health")
ordplot2 + stat_ellipse(geom = "polygon", type = "norm", linetype = 2, alpha = 0.2, aes(fill = Health)) +
  stat_ellipse(type = "t", level = 0.95) +theme_bw() + labs(title = "NMDS Ordination by Sample Health") + theme(plot.title = element_text(hjust = 0.5))

#NMDS at species level
ps.species <- tax_glom(ps.2health, taxrank = "Species")
ps.speciesrel <- transform_sample_counts(ps.species, function(otu) otu/sum(otu))
species.ord <- ordinate(ps.speciesrel, method = "NMDS", distance = "bray")
speciesordplot <- plot_ordination(ps.speciesrel, species.ord, "samples", color = "Health", shape = "Health")
speciesordplot + stat_ellipse(geom = "polygon", type = "norm", linetype = 1, alpha = 0.1, aes(fill = Health)) + theme_bw() +
  labs(title = "NMDS Ordination at Species Level") + theme(plot.title = element_text(hjust = 0.5))

speciesordplot <- plot_ordination(ps.speciesrel, species.ord, "samples", color = "Location", shape = "Location")
speciesordplot + stat_ellipse(geom = "polygon", type = "norm", linetype = 1, alpha = 0.1, aes(fill = Location)) + theme_bw() +
  labs(title = "NMDS Ordination at Species Level") + theme(plot.title = element_text(hjust = 0.5))

table1 <- summarise()

#top 20
ps.top40 <- subset_taxa(ps.2health, !is.na(Species) & !is.na(Phylum) & !is.na(Genus) & !is.na(Class) & !is.na(Family) & !is.na(Order) & !is.na(Kingdom))
ps.top40 <- subset_taxa(ps.2top20, !apply(tax_table(ps.2top20), 1, function(x) any(grepl("Incertae_sedis", x, ignore.case = TRUE))))
top40 <- names(sort(taxa_sums(ps.2top20), decreasing = TRUE)) [1:40]
ps.top40 <- transform_sample_counts(ps.2top20, function(otu) otu/sum(otu))
?transform_sample_counts()
ps.top40 <- prune_taxa(top20.2, ps.2top20)
top40df <- psmelt(ps.top40) 
total.abundance <- top40df %>%
  group_by(Health, Genus, Species) %>%
  summarize(TotalAbundance = sum(Abundance))

ggplot(top40df, aes(x = Species, y = Abundance, fill = Genus)) + geom_bar(stat = "identity", position = "stack") +
  geom_text(data = total.abundance, aes(x = Species, y = TotalAbundance, label = round(TotalAbundance, digits = 2), fill = NULL), hjust = 0, size = 3) +
  theme_bw() + coord_flip() + facet_wrap(~Health)

top20 <- names(sort(taxa_sums(ps), decreasing = TRUE)) [1:20]
ps.top20 <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ps.top20 <- prune_taxa(top20, ps.top20)
sample_variables(ps.top20)
plot_bar(ps.top20, x = "Health", fill = "Phylum") + facet_wrap(~Health, scales = "free_x")  

moreunhealthy <- names(sort(taxa_sums(ps.2top20), decreasing = TRUE))


##Get OTUs that were more present in unhealthy than healthy samples
health2.speciesnoNA <- subset_taxa(ps.2health, !is.na(Species))
table.nospeciesNA <- psmelt(health2.speciesnoNA)
table1 <- table.nospeciesNA %>% filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(n=n())
table1

health2.noNA <- subset_taxa(ps.2health, !is.na(Species) & !is.na(Phylum) & !is.na(Genus) & !is.na(Class) & !is.na(Family) & !is.na(Order) & !is.na(Kingdom))
health2.noNA <- subset_taxa(health2.noNA, !apply(tax_table(health2.noNA), 1, function(x) any(grepl("Incertae_sedis", x, ignore.case = TRUE))))
health2df <- psmelt(health2.noNA)
nrow(health2df)

table1 <- health2df %>%
  group_by(Location) %>%
  summarise(n=n_distinct(Sample))
table1

names.healthy <- sample_names(subset_samples(health2.noNA, Health == "Healthy"))
names.unhealthy <- sample_names(subset_samples(health2.noNA, Health == "Unhealthy"))

otu.matrix <- as(otu_table(health2.noNA), "matrix")
if (taxa_are_rows(health2.noNA) == FALSE) {otu.matrix <- t(otu.matrix)}
sum.healthy <- rowSums(otu.matrix[, names.healthy, drop = FALSE])
sum.unhealthy <- rowSums(otu.matrix[, names.unhealthy, drop = FALSE])
moreunhealthy <- sum.unhealthy > sum.healthy
ps.moreunhealthy <- prune_taxa(moreunhealthy, health2.noNA)
moreunhealthytable <- psmelt(ps.moreunhealthy)


plot_bar(ps.moreunhealthy, x = "Genus", fill = "Species") + facet_wrap(~Health, scales = "free_x") + coord_flip() + theme_bw()
top20unhealthy <- names(sort(taxa_sums(ps.moreunhealthy), decreasing = TRUE)) [1:40]
ps.unhealthy40 <- prune_taxa(top20unhealthy, ps.moreunhealthy)
plot_bar(ps.unhealthy40, x = "Genus", fill = "Species") + facet_wrap(~Health) + coord_flip() + theme_bw()

#Get unique OTU between health groups 
ps.healthy <- subset_samples(ps, Health == "Healthy")
ps.unhealthy <- subset_samples(ps, Health %in% c("Dead", "Dying"))

taxa.healthy <- taxa_names(prune_taxa(taxa_sums(ps.healthy) > 0, ps.healthy))
taxa.unhealthy <- taxa_names(prune_taxa(taxa_sums(ps.unhealthy) > 0, ps.unhealthy))

unique.healthy <- setdiff(taxa.healthy, taxa.unhealthy)
unique.unhealthy <- setdiff(taxa.unhealthy, taxa.healthy)

tax.table <- as.data.frame(tax_table(ps))
seq.table <- as.data.frame(otu_table(ps))

tax2.table <- as.data.frame(tax_table(ps.2health))

unique.healthy.table <- tax.table[unique.healthy, ]
unique.unhealthy.table <- tax.table[unique.unhealthy, ]

unique.healthy.table <- as.data.frame(unique.healthy.table)
unique.unhealthy.table <- as.data.frame(unique.unhealthy.table)

otu.healthy <- as.data.frame(otu_table(ps.healthy))
otu.unhealthy <- as.data.frame(otu_table(ps.unhealthy))

tax.noNA <- tax %>%
  drop_na("Species", "Phylum", "Genus", "Class", "Family", "Order", "Kingdom") %>%
  filter(!if_any(c(Species, Phylum, Genus, Class, Family, Order, Kingdom), ~ grepl("incertae_sedis", ignore.case = TRUE, .)))


length(unique(tax.noNA$X))
otu.noNA <- as.data.frame(otu_table(health2.noNA))
samp.names <- colnames(otu.noNA[1:246])

otu.samples <- otu.noNA %>%
  mutate(across(c(samp.names), ~ifelse(.x !=0,1,.x)
  ))
otu.samples[nrow(otu.samples) + 1, 1:246] <- colSums(otu.samples[,1:246], na.rm=TRUE)
rownames(otu.samples)[97] <- "Total"

#Phytophthora custom BLAST
blastoutre <- read.csv("blastoutre.csv")
blastoutre <- blastoutre %>% mutate(qseqid = if_else(
  stringr::str_detect(qseqid, ";size"),
  stringr::str_replace(qseqid, ";size", ",size"),
  qseqid))

blastoutsplit <- blastoutre %>%
  separate_wider_delim(qseqid, ";", names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"), too_few = "align_start")

##DESeq2 - differential abundance
health2.noNA <- subset_taxa(ps.2health, !is.na(Species) & !is.na(Phylum) & !is.na(Genus) & !is.na(Class) & !is.na(Family) & !is.na(Order) & !is.na(Kingdom))
health2.noNA <- subset_taxa(health2.noNA, !Genus %in% c("gen_incertae_sedis"))
diagdds <- phyloseq_to_deseq2(health2.noNA, ~Health)
diagdds <- estimateSizeFactors(diagdds, type = "poscounts")
diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")

res <- results(diagdds, cooksCutoff = FALSE)
alpha <- 0.05
sigtab <- res[which(res$padj < alpha),]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(health2.noNA)[rownames(sigtab),],"matrix"))
head(sigtab)
dim(sigtab)

x <- tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max (x))
x <- sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels = names(x))
x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab$Genus <- factor(as.character(sigtab$Genus), levels = names(x))
ggplot(sigtab, aes(x = Genus, y=log2FoldChange, color = Phylum)) + geom_point(size = 6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) + theme_bw()

#alpha diversity
health2prune <- prune_species(speciesSums(ps.2health) > 0, ps.2health)
theme_set(theme_bw())
pal <- "Set1"
plot_richness(health2prune, measures = "Shannon", color = "Location") + facet_wrap(~Health, scales = "free_x")
plot_richness(health2prune, x = "Health", measures = "Shannon", color = "Location") + facet_wrap(~Health, scales = "free_x")
