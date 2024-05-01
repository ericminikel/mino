
# BiocManager::install("DESeq2")

library(AnVIL)
library(tidyverse)
library(janitor)
library(DESeq2)

setwd('R/')

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35

percent = function(x, format='f', digits=0, signed=F) { paste0(ifelse(signed & x > 0, '+', ''), formatC(x*100,format=format,digits=digits),'%') }

clipdist = function(x, minx, maxx) pmin(pmax(x,minx),maxx)


# gsutil_cp('gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/meta.tsv', 'R/output/')

# gsutil_cp('gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/sorce-nuvolone-2020-table-s2.txt', 'R/output/')

sn20 = read_tsv('R/output/sorce-nuvolone-2020-table-s2.txt') %>% clean_names()
meta = read_tsv('R/output/meta.tsv')
qmat_raw = read_tsv('R/output/qmat.tsv')


# quick look to convince yourself - changes may only be seen in animals actively on drug,
# washout is very fast
qmat_raw %>% filter(symbol=='Aif1') %>% pivot_longer(2:16) %>% inner_join(meta, by=c('name'='sample_id')) %>% group_by(treatment, hours) %>% summarize(.groups='keep', mean=mean(value))

meta$keep = meta$treatment=='saline' | meta$hours==312

qmat = as.matrix(qmat_raw[,which(colnames(qmat_raw) %in% meta$sample_id[meta$keep])])
rownames(qmat) = qmat_raw$symbol
mode(qmat) = 'integer'

qmat = qmat[rowMeans(qmat) >= 1,]

pca = prcomp(qmat, scale.=T)
pcvar = pca$sdev^2
var_exp = pcvar / sum(pcvar)
pc = pca$rotation
plot(pc[,1], pc[,2], pch=20)
text(pc[,1], pc[,2], labels=meta$sample_id, pos=2)

meta$pc1 = pc[,1]

meta = meta[meta$keep,]

# confirm that order matches
stopifnot(all(meta$sample_id == colnames(qmat)))

dds = DESeqDataSetFromMatrix(countData = qmat,
                             colData = meta,
                             design = ~ sex + treatment)
dds = DESeq(dds)
results_names = resultsNames(dds) # lists the coefficients
res = results(dds, contrast=c("treatment","minocycline","saline"))


res$padj[is.na(res$padj)] = 1

res %>%
  as_tibble(rownames='gene') %>%
  arrange(padj) -> restbl


enr_meta = tibble(code=c('_NE','_NP','AS','EC','MG','N','OL'),
                  name=c('not enriched','not assigned','astrocytes','endothelial','microglia','neurons','oligodendrocytes'),
                  colors=c('#7C7C7C','#7C7C7C','#CC6576','#49A999','#882354','#332488','#88CDEE'))


dge_contrast = function(dds_obj, condition1, condition2) {
  results(dds_obj, contrast=c("treatment",condition1,condition2)) %>% as_tibble(rownames='gene') %>% arrange(padj) %>% clean_names() -> return_value
  return_value$enr = sn20$enrichment[match(return_value$gene, sn20$gene_id)]
  return_value$enr[is.na(return_value$enr)] = '_NP'
  return_value$padj[is.na(return_value$padj)] = 1
  return (return_value)
}

volcano = function(res_obj, save=T, name=deparse(substitute(res_obj))) {
  res_obj$color = enr_meta$colors[match(res_obj$enr, enr_meta$code)]
  res_obj$color_disp = alpha(res_obj$color, ifelse(res_obj$padj < 0.05, 1, ci_alpha))
  res_obj$cex = ifelse(res_obj$padj < 0.05, 1, 0.5)
  
  
  xlims = c(-3, 3)
  ylims = c(0, max(5,max(-log10(res_obj$padj),na.rm=T)*1.1))
  par(mar=c(4,4,3,3))
  plot(x=clipdist(res_obj$log2fold_change, xlims[1], xlims[2]), y=-log10(res_obj$padj), pch=20, xlim=xlims, ylim=ylims,
       col=res_obj$color_disp, cex=res_obj$cex,
       xlab = 'L2FC', ylab = '-log10(padj)', xaxs='i', yaxs='i', axes=F)
  mtext(side=3, line=0.5, text='', cex=0.8, font=2)
  axis(side=1, at=-5:5)
  abline(v=0, lwd=.125)
  abline(h=-log10(0.05), lwd=.125)
  axis(side=2, las=2)
  axis(side=2, at=ylims, labels=NA, lwd=1, lwd.ticks=0)
  res_obj %>%
    filter(padj < 0.05 & base_mean > 3) %>%
    mutate(pos=ifelse(log2fold_change < 0, 2, 4)) %>%
    arrange(padj) %>%
    head(20) -> to_label
  if (nrow(to_label) > 0) {
    par(xpd=T)
    text(x=clipdist(to_label$log2fold_change, xlims[1], xlims[2]), y=-log10(to_label$padj), labels=to_label$gene, font=3, pos=to_label$pos, cex=0.75)
    par(xpd=F)
  }
}

handle_contrasts = function(group1, group2) {
  name = paste(group1,group2,sep='_')
  res = dge_contrast(dds, group1, group2)
  assign(name, res, envir = parent.frame())
  resx=300
  png(paste0('display_items/',name,'.png'), width=3.5*resx, height=3.5*resx, res=resx)
  volcano(res, name=name)
  dev.off()
  write_tsv(res, paste0('output/res_',group1,'_vs_',group2,'.tsv'))
}

handle_contrasts('minocycline', 'saline')

minocycline_saline$sig = minocycline_saline$padj < 0.05
minocycline_saline$sigdown = minocycline_saline$padj < 0.05 & minocycline_saline$log2fold_change < 0
minocycline_saline %>%
  group_by(enr) %>%
  summarize(.groups='keep',
            n_sig = sum(padj < 0.05),
            p_sig = mean(padj < 0.05))

table(minocycline_saline[,c('sig','enr')])
minocycline_saline$mg = minocycline_saline$enr=='MG'
chisq.test(table(minocycline_saline[,c('sig','enr')]))
chisq.test(table(minocycline_saline[,c('sigdown','mg')]))

t.test(minocycline_saline$log2fold_change[minocycline_saline$enr=='MG' & minocycline_saline$padj < 0.05],
       minocycline_saline$log2fold_change[minocycline_saline$enr!='MG' & minocycline_saline$padj < 0.05])

t.test(minocycline_saline$log2fold_change[minocycline_saline$enr=='MG' & minocycline_saline$padj < 0.05],
       minocycline_saline$log2fold_change[minocycline_saline$enr!='MG' & minocycline_saline$padj < 0.05])



t.test(minocycline_saline$log2fold_change[minocycline_saline$enr=='EC' & minocycline_saline$padj < 0.05],
       minocycline_saline$log2fold_change[minocycline_saline$enr!='EC' & minocycline_saline$padj < 0.05])

write_tsv(minocycline_saline, 'output/minocycline_saline.tsv')

gsutil_cp( 'src/merge_count_matrix.R', 'gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/src/merge_count_matrix.R')
gsutil_cp( 'src/mino_dge.R', 'gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/src/mino_dge.R')
gsutil_cp( 'output/qmat.tsv', 'gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/output/qmat.tsv')
gsutil_cp( 'output/minocycline_saline.tsv', 'gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/output/minocycline_saline.tsv')
gsutil_cp( 'display_items/minocycline_saline.png', 'gs://fc-3effa716-6fbc-4071-a224-96c78e5129c7/output/minocycline_saline.png')
