#### get altered pathway in FOXA1 mutants under DHT treatment
#### library(edgeR);
library(clusterProfiler);
library(gProfileR);
setwd("/.mounts/labs/cpcgene/private/projects/RNA-Seq_PipelineDevelopment/analysis/star/star-circ/lsd1/mutant/analysis");
#### only one replicate for each condition, use fold change 
rpkm.gene <- read.table('rnaseq/count/2018-11-23_STAR_FPKM_TMM_stranded.txt', row.names = 1, as.is = TRUE);
rpkm.gene <- rpkm.gene[, grep('DHT|Symbol', colnames(rpkm.gene))];
diff.gene <- data.frame(gene = rpkm.gene$Symbol, D226G_DHT = rowMeans(rpkm.gene[, 1:2]), 
	M253K_DHT = rowMeans(rpkm.gene[, 3:4]), WT_DHT = rowMeans(rpkm.gene[, 5:6]));
diff.gene$lfc_D226GvsWT <- apply(diff.gene[, -1] + 1, 1, function(x) log2(x[1]/x[3]));
diff.gene$lfc_M253KvsWT <- apply(diff.gene[, -1] + 1, 1, function(x) log2(x[2]/x[3]));
diff.gene$mean <- rowMeans(diff.gene[, 2:4]);
### tried with different cut-offs
##### mean>1&abs(lfc)>0.58; mean>0.1&abs(lfc)>0.58;
d226g_up <- diff.gene[diff.gene$mean>1&diff.gene$lfc_D226GvsWT>0, ]$gene;
d226g_down <- diff.gene[diff.gene$mean>1&diff.gene$lfc_D226GvsWT<(-0), ]$gene;
m253k_up <- diff.gene[diff.gene$mean>1&diff.gene$lfc_M253KvsWT>0, ]$gene;
m253k_down <- diff.gene[diff.gene$mean>1&diff.gene$lfc_M253KvsWT<(-0), ]$gene;
####perform pathway enrichment analysis
ego_d226g_up <-  gprofiler(query = as.vector(d226g_up),
    organism = "hsapiens", significant = TRUE,
    correction_method = "fdr",
    src_filter = c("GO", 'REAC'),
    underrep = FALSE,
    exclude_iea = FALSE,
    custom_bg = as.vector(diff.gene$gene)
    );
###
ego_d226g_down <-  gprofiler(query = as.vector(d226g_down),
    organism = "hsapiens", significant = TRUE,
    correction_method = "fdr",
    src_filter = c("GO", 'REAC'),
    underrep = FALSE,
    exclude_iea = FALSE,
    custom_bg = as.vector(diff.gene$gene)
    );
###
ego_m253k_up <-  gprofiler(query = as.vector(m253k_up),
    organism = "hsapiens", significant = TRUE,
    correction_method = "fdr",
    src_filter = c("GO", 'REAC'),
    underrep = FALSE,
    exclude_iea = FALSE,
    custom_bg = as.vector(diff.gene$gene)
    );
##
ego_m253k_down <-  gprofiler(query = as.vector(m253k_down),
    organism = "hsapiens", significant = TRUE,
    correction_method = "fdr",
    src_filter = c("GO", 'REAC'),
    underrep = FALSE,
    exclude_iea = FALSE,
    custom_bg = as.vector(diff.gene$gene)
    );
save(list = ls(), file = paste0(Sys.Date(), '_foxa1_mutant_pathway.rda'));
####
ego <- ego_d226g_up;
ego <- ego[ego$term.size>=10&ego$term.size<=500, ];
ego <- ego[order(ego$p.value), ];
write.csv(ego, paste0(Sys.Date(), '_foxa1_mutant_pathway_d226g_up.csv'));
ego <- ego_d226g_down;
ego <- ego[ego$term.size>=10&ego$term.size<=500, ];
ego <- ego[order(ego$p.value), ];
write.csv(ego, paste0(Sys.Date(), '_foxa1_mutant_pathway_d226g_down.csv'));
ego <- ego_m253k_up;
ego <- ego[ego$term.size>=10&ego$term.size<=500, ];
ego <- ego[order(ego$p.value), ];
write.csv(ego, paste0(Sys.Date(), '_foxa1_mutant_pathway_m253k_up.csv'));
ego <- ego_m253k_down;
ego <- ego[ego$term.size>=10&ego$term.size<=500, ];
ego <- ego[order(ego$p.value), ];
write.csv(ego, paste0(Sys.Date(), '_foxa1_mutant_pathway_m253k_down.csv'));
###
gene.rea <- read.gmt('~/circRNA/star-circ/circRNA_landscape/rnaseq_landscape/data/ref/gp_ref/hsapiens.REAC.NAME.gmt');
diff.gene$tgfb <- ifelse(diff.gene$gene%in%gene.rea[gene.rea$ont=='REAC:170834', ]$gene, 1, 0);
write.csv(diff.gene, paste0(Sys.Date(), '_foxa1_mutant_rpkm_diff.csv'))