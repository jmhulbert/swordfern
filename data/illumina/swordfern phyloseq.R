library(dplyr)
library(tidyverse)
if (!require("ßBiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(tibble)
BiocManager::install("ShortRead")
BiocManager::install("dada2")
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
BiocManager::install("rBLAST")
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
library(Biostrings)
library(DECIPHER)
library(dada2)
library(ShortRead)
library(phyloseq)
library(rBLAST)
library(ggplot2)
library(DESeq2)
devtools::install_github("gmteunisse/fantaxtic")
library(fantaxtic)
library(edgeR)
library(stringr)
library(ggpubr)


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

topbyhealth <- tax.table.count %>%
  filter(Abundance > 0) %>%
  group_by(Health) %>%
  summarize(Taxa = list(unique(Genus)))

topbyhealth.species <- tax.table.count %>%
  filter(Abundance > 0) %>%
  group_by(Health) %>%
  summarize(Taxa = list(unique(Species)))

topbyhealth$Taxa[[1]]
topbyhealth$Taxa[[2]]
topbyhealt.species$Taxa[[1]]

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
f
ps.sub <- subset_samples(ps, !is.na(Health))
expt <- prune_taxa(names(sort(taxa_sums(ps.sub), TRUE)[1:50]), ps.sub)
ord <- ordinate(expt, formula = ~Health, "NMDS", "bray")
ordplot <- plot_ordination(expt, ord, "samples", color = "Health", shape = "Health")
ordplot + stat_ellipse(geom = "polygon", type = "norm", linetype = 2, alpha = 0.2, aes(fill=Health)) +
  stat_ellipse(type = "t", level = 0.95) + theme_bw()

#explore this but species lvl
ps.2sub <- subset_samples(ps.2health, !is.na(Health))
expt2 <- prune_taxa(names(sort(taxa_sums(ps.2sub), TRUE)[1:50]), ps.2sub)
ord2 <- ordinate(expt2, formula = ~Health, "NMDS", "bray")
ordplot2 <- plot_ordination(expt2, ord2, "samples", color = "Health", shape = "Health")
ordplot2 + stat_ellipse(geom = "polygon", type = "norm", linetype = 2, alpha = 0.2, aes(fill = Health)) +
  stat_ellipse(type = "t", level = 0.95) +theme_bw()

#NMDS at species level
ps.species <- tax_glom(ps.2health, taxrank = "Species")
ps.speciesrel <- transform_sample_counts(ps.species, function(otu) otu/sum(otu))
species.ord <- ordinate(ps.speciesrel, method = "NMDS", distance = "bray")
speciesordplot <- plot_ordination(ps.speciesrel, species.ord, "samples", color = "Health", shape = "Health")
speciesordplot + stat_ellipse(geom = "polygon", type = "norm", linetype = 2, alpha = 0.2, aes(fill = Health)) + theme_bw() + stat_ellipse(type = "t", level = 0.95)


#top 20
ps.2top20 <- subset_taxa(ps.2health, !is.na(Species) & !is.na(Phylum) & !is.na(Genus) & !is.na(Class) & !is.na(Family) & !is.na(Order) & !is.na(Kingdom))
ps.2top20 <- subset_taxa(ps.2top20, !Genus %in% c("gen_incertae_sedis"))
top20.2 <- names(sort(taxa_sums(ps.2top20), decreasing = TRUE)) [1:40]
top20.2
?taxa_sums
ps.2top20 <- transform_sample_counts(ps.2top20, function(otu) otu/sum(otu))
?transform_sample_counts()
ps.2top20 <- prune_taxa(top20.2, ps.2top20)
plot_bar(ps.2top20, x = "Genus", fill = "Species") + facet_wrap(~Health, scales = "free_x") + coord_flip() + theme_bw()


top20 <- names(sort(taxa_sums(ps), decreasing = TRUE)) [1:20]
ps.top20 <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ps.top20 <- prune_taxa(top20, ps.top20)
sample_variables(ps.top20)
plot_bar(ps.top20, x = "Health", fill = "Phylum") + facet_wrap(~Health, scales = "free_x")  


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

colorfunc <- function(palname = "Set1", ...){
  scale_fill_brewer(palatte = palname, ...)
}
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
<<<<<<< HEAD
pal <- "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}
plot_richness(health2prune, measures = "Shannon") + facet_wrap(~Health, scales = "free_x") 



=======
plot_richness(health2prune, x = "Health", measures = "Shannon") + facet_wrap(~Health, scales = "free_x")
>>>>>>> d0f2b920b78c74699d6535d10d58ee99855171da
