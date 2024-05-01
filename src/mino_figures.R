overall_start_time = Sys.time()
tell_user = function(x) { cat(file=stderr(), x); flush.console() }

tell_user('Loading required packages...')

options(stringsAsFactors=F)
setwd('~/d/sci/src/mino/')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
suppressMessages(library(robustrank))
suppressMessages(library(magick))
#BiocManager::install('EmpiricalBrownsMethod')
suppressMessages(library(EmpiricalBrownsMethod))

# for debugging - can remove it and all references to it later
stderrmsg = function(msg, verbose=T) {
  if (verbose) {
    cat(file=stderr(),msg)
    flush.console()
  }
}

tell_user('done.\nLoading functions...')

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}


alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot
pt_alpha = 0.7 # degree of transparency for shading points in plot

rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}

p_to_symbol = function(p) {
  case_when(p < 0.001 ~ '***',
            p < 0.01 ~ '**',
            p < 0.05 ~ '*',
            TRUE ~ '')
}


clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

# one fewer parameter - do clipdist symmetrically around 0
clipsym = function(x, absx) {
  return (pmin(absx, pmax(-absx, x)))
}


bin_age = function(age) {
  case_when(age < 20 ~ '<20',
            age >= 20 & age < 80 ~ paste0(floor(age/5)*5,'-',floor(age/5)*5+4),
            age >= 80 ~ '≥80')
}

unbin_age = function(age_bin) {
  case_when(age_bin == '<20' ~ 18,
            age_bin == '≥80' ~ 82,
            TRUE ~ suppressWarnings((as.numeric(gsub('-.*','',age_bin)) + as.numeric(gsub('.*-','',age_bin)))/2))
}

makecano = function(peptides, 
                    meta, 
                    contrast,
                    dataset_id=NULL, 
                    clipwidth=4, 
                    ymax=NA, 
                    verbose=T,
                    saveimgs=F) {
  
  # peptides can be either a path or a tibble. if path, read in
  if (length(class(peptides))==1) {
    if (class(peptides)=='character') {
      peptides_path = peptides
      peptides = read.xlsx(peptides, startRow = 3) %>% clean_names() %>% as_tibble()
      # if dataset_id is not provided, try to extract it from the peptides path
      if (is.null(dataset_id)) {
        dataset_id = substr(peptides_path,6,9)
      }
    }
  }
  stderrmsg(paste0('Processing dataset ',dataset_id,'\n'),verbose)
  
  
  # meta should have exactly these two columns
  stopifnot(identical(colnames(meta), c('tmt_colname','grp')))
  
  
  # some datasets have an identical peptide that appears multiple times in the LC gradient
  # with different scan_number values. samples will generally be high in one scan and
  # low in the other, and if this is confounded with sample group, this leads the different rows 
  # to have opposite direction of effect but both have significant P values. when combined with
  # EBM, this yields extremely high P values with virtually zero effect sizes. to remove
  # this artifact, group at the level of peptide sequence before running EBM.
  peptides %>%
    select(-charge, -scan_number) %>%
    group_by(gene_symbol, peptide_sequence) %>%
    summarize(.groups='keep', across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
    ungroup() -> peptides_unique
  
  # contrast should indicate which groups to compare; groups to be collapsed get a "+" between them
  if (nchar(contrast)==5) {
    comparison = strsplit(contrast,'/')[[1]]
    tmt_colnames = meta$tmt_colname[meta$grp %in% comparison]
  } else {
    comparison = strsplit(contrast,'/')[[1]]
    tmt_colnames = meta$tmt_colname
    latter = comparison[2]
    latters = strsplit(latter,'\\+')[[1]]
    meta$grp = ifelse(meta$grp %in% latters, latter, meta$grp)
  }
  
  stderrmsg(paste0('Groups: ',comparison[1],' N = ',sum(meta$grp==comparison[1]),', ',
                   comparison[2],' N = ',sum(meta$grp==comparison[2]),'\n'),verbose)
  
  stderrmsg(paste0('Producing individual peptide volcano\n'),verbose)
  
  # individual peptide-based volcano  
  peptides_unique %>%
    select(gene=gene_symbol, peptide_sequence, all_of(tmt_colnames)) %>%
    pivot_longer(cols = all_of(tmt_colnames)) %>%
    inner_join(meta, by=c('name'='tmt_colname')) %>%
    filter(grp %in% comparison) %>%
    group_by(gene, peptide_sequence, grp) %>%
    summarize(.groups='keep', values = list(value)) %>%
    spread(grp, values) %>%
    group_by(gene, peptide_sequence) %>%
    summarize(.groups='keep',
              p_value = t.test(unlist(get(comparison[1])), unlist(get(comparison[2])))$p.value,
              l2fc = log2(mean(unlist(get(comparison[1])) / mean(unlist(get(comparison[2])))))) -> pepcano
  
  stderrmsg(paste0('Collapsing by gene symbol\n'),verbose)
  
  peptides_unique %>%
    select(gene=gene_symbol, peptide_sequence, all_of(tmt_colnames)) %>%
    inner_join(pepcano, by=c('gene','peptide_sequence')) %>%
    group_by(gene) %>%
    reframe(.groups='keep',
              mat = as.matrix(across(all_of(tmt_colnames))),
              pvals = as.numeric(p_value)) %>%
    group_by(gene) %>%
    summarize(.groups='keep',
              p_ebm = empiricalBrownsMethod(data_matrix=mat, p_values=pvals, extra_info=F)) -> ebmcano
  
  peptides_unique %>%
    select(gene=gene_symbol, peptide_sequence, all_of(tmt_colnames)) %>%
    pivot_longer(cols = all_of(tmt_colnames)) %>%
    inner_join(meta, by=c('name'='tmt_colname')) %>%
    filter(grp %in% comparison) %>%
    group_by(gene, peptide_sequence, grp) %>%
    summarize(.groups='keep', values = list(value)) %>%
    spread(grp, values) %>%
    group_by(gene, peptide_sequence) %>%
    summarize(.groups='keep',
              tot = mean(c(unlist(get(comparison[1])), unlist(get(comparison[2])))),
              l2fc = log2(mean(unlist(get(comparison[1])) / mean(unlist(get(comparison[2])))))) %>%
    group_by(gene) %>%
    summarize(.groups='keep',
              n_pep = n(),
              abun = sum(tot),
              l2fc = mean(l2fc)) -> pepstats
  
  corrected_threshold = 0.05 / nrow(pepstats)
  ebmcano %>%
    inner_join(pepstats, by='gene') %>%
    mutate(color = case_when(p_ebm >= 0.05 ~ '#C9C9C9',
                             p_ebm < 0.05 & p_ebm >= corrected_threshold & l2fc < 0 ~ '#FFCCCC',
                             p_ebm < 0.05 & p_ebm >= corrected_threshold & l2fc > 0 ~ '#CCCCFF',
                             p_ebm < corrected_threshold & l2fc < 0 ~ '#BB1111',
                             p_ebm < corrected_threshold & l2fc > 0 ~ '#1111BB')) %>%
    ungroup() %>%
    arrange(p_ebm) -> allcano
  
  stderrmsg(paste0('Making QQ and volcano plots\n'),verbose)
  
  if(saveimgs) {
    resx=300
    png(paste0('display_items/tmt_',dataset_id,'_',paste0(comparison,collapse='_vs_'),'_qq.png'), width=3.5*resx, height=4.0*resx, res=resx)
  n = nrow(allcano)
  observed = sort(allcano$p_ebm) # actual p values
  expected = seq(1/n,1,1/n) # under null, p values are uniformly distributed between 0 and 1
  o = -log10(observed) # log transform the p values
  e = -log10(expected) # log transform the expected p values
  plot(e,o,pch=19,xlab="expected -log10(p)",ylab="observed -log10(p)",main="qq plot")
  abline(a=0,b=1,col="red") # add a line representing the null
  silence_message = dev.off()
  
  resx=300
  png(paste0('display_items/tmt_',dataset_id,'_',paste0(comparison,collapse='_vs_'),'_volcano.png'), width=3.5*resx, height=4.0*resx, res=resx)
  par(mar=c(3,3,1,3))
  xlims = clipwidth*c(-1.05,1.05)
  ylims = c(0, ifelse(is.na(ymax), max(-log10(allcano$p_ebm))*1.05, ymax))
  plot(x=NA, y=NA, axes=F, ann=F, xaxs='i', yaxs='i', xlim=xlims, ylim=ylims)
  points(x=clipsym(allcano$l2fc,clipwidth), 
         y=-log10(allcano$p_ebm), 
         pch=20, 
         col=allcano$color)
  axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
  axis(side=1, at=floor(min(xlims)):ceiling(max(xlims)), tck=-0.025, labels=NA)
  axis(side=1, line=-1, at=floor(min(xlims)):ceiling(max(xlims)), lwd=0, cex.axis=0.7)
  mtext(side=1, line=1, text='log2(fold change)')
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  axis(side=2, at=floor(min(ylims)):ceiling(max(ylims)), tck=-0.025, labels=NA)
  axis(side=2, at=seq(floor(min(ylims)),ceiling(max(ylims)),5), tck=-0.05, cex.axis=0.7, las=2)
  mtext(side=2, line=2, text='-log10(P)')
  abline(h=-log10(c(0.05, corrected_threshold)), lwd=0.25)
  mtext(side=4, at=-log10(c(0.05, corrected_threshold)), line=0.25, las=2,
        text=c('nominal\nP < 0.05', 'Bonferroni\nP < 0.05'), cex=0.5)
  par(xpd=T)
  # label those that are P < 0.1 after Bonferroni, AND
  # (are in the top 10 either by P value OR by absolute L2FC)
  callouts = allcano$p_ebm < 0.1 / nrow(allcano) &
    (rank(allcano$p_ebm) <= 10 | rank(-abs(allcano$l2fc)) <= 10) | allcano$gene %in% 'Nefl'
  text(x=clipsym(allcano$l2fc[callouts],clipwidth), 
       y=-log10(allcano$p_ebm[callouts]), 
       labels=allcano$gene[callouts], 
       pos=ifelse(allcano$l2fc[callouts] < 0 & allcano$l2fc[callouts] > -clipwidth, 2, 4), 
       font=3, cex=0.6, family='mono')
  par(xpd=F)
  silence_message = dev.off()
  }

  
  return (allcano)
}


dviz = function(tbl,
                xlims,
                ylims,
                xvar,
                yvar,
                colorvar=NULL,
                pchvar=NULL,
                pchbg='#FFFFFF',
                pchcex=1,
                xcols=character(0), # columns that group with x TBD
                xats=xbigs,
                xbigs,
                xlwds=1,
                xbiglabs=xbigs,
                xaxcex=1,
                xlabcex=1,
                xlabline=1.5,
                yats=ybigs,
                ybigs,
                ylwds=1,
                ybiglabs=ybigs,
                yaxcex=1,
                ylabcex=1,
                ylabline=2.25,
                log,
                mar=c(3,3,1,1),
                jitamt=0.1,
                randseed=1,
                boxwidth=NA,
                barwidth=NA,
                polygon=NA,
                arrowlength=0.05,
                test=NA,
                control_group=NA,
                xlab='',
                ylab='',
                legtext=NULL,
                legcol='#000000',
                legtextcol=legcol,
                leglty=1,
                legpch=20,
                leglwd=1,
                legcex=1,
                crosshairs=F
) {
  
  if (is.null(pchvar)) {
    pchvar='pch'
    tbl$pch = 19
  }
  
  if (is.null(colorvar)) {
    colorvar='color'
    tbl$color = '#000000'
  }
  
  tbl %>%
    mutate(x=!!as.name(xvar), y=!!as.name(yvar), color=!!as.name(colorvar), pch=!!as.name(pchvar)) %>%
    select(x, y, color, pch, all_of(xcols)) -> tbl
  
  if (!crosshairs) {
    xcols = c('x', xcols)
  }
  
  tbl %>%
    group_by(color, across(all_of(xcols))) %>%
    summarize(.groups='keep',
              n = n(),
              mean = mean(y),
              l95 = mean(y) - 1.96 * sd(y) / sqrt(n()),
              u95 = mean(y) + 1.96 * sd(y) / sqrt(n()),
              median = median(y),
              q25 = quantile(y, .25),
              q75 = quantile(y, .75),
              cv = sd(y) / mean(y),
              x_mean = mean(x),
              x_l95  = mean(x) - 1.96 * sd(x) / sqrt(n()),
              x_u95  = mean(x) + 1.96 * sd(x) / sqrt(n())) %>%
    ungroup() %>%
    mutate(l95 = ifelse(l95 < 0 & log %in% c('xy','y'), min(ylims), l95)) %>%
    #mutate(x = case_when(crosshairs ~ x_mean, 
    #                     TRUE ~ x)) %>%
    arrange(x) -> tbl_smry
  
  if (crosshairs) {
    tbl_smry$x = tbl_smry$x_mean
  }
  
  par(mar=mar)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log=log)
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  if (!is.null(xats)) {
    axis(side=1, at=xats, tck=-0.025, lwd.ticks=xlwds, labels=NA)
  }
  if (!is.null(xbigs)) {
    axis(side=1, at=xbigs, tck=-0.05, lwd.ticks=xlwds, labels=NA)
    axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbiglabs, cex.axis=xaxcex)
  }
  mtext(side=1, line=xlabline, text=xlab, cex=xlabcex)
  if (!is.null(yats)) {
    axis(side=2, at=yats, tck=-0.025, labels=NA)
  }
  if (!is.null(ybigs)) {
    axis(side=2, at=ybigs, tck=-0.05, labels=NA)
    axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.3, labels=ybiglabs, cex.axis=yaxcex)
  }
  mtext(side=2, line=ylabline, text=ylab, cex=ylabcex)
  if (crosshairs) {
    jitamt = 0
  }
  set.seed(randseed)
  if (jitamt==0) { # ?jitter won't take 0 for an answer, so work around it
    jittered_x = tbl$x
  } else {
    jittered_x = jitter(tbl$x, amount=jitamt)
  }
  points(x=jittered_x, y=tbl$y, col=alpha(tbl$color,ci_alpha), pch=tbl$pch, bg=pchbg)
  
  if (crosshairs) {
    segments(x0=tbl_smry$x_l95, x1=tbl_smry$x_u95, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    segments(x0=tbl_smry$x_mean, y0=tbl_smry$l95, y1=tbl_smry$u95, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(barwidth)) {
    segments(x0=tbl_smry$x-barwidth, x1=tbl_smry$x+barwidth, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    arrows(x0=tbl_smry$x, y0=tbl_smry$l95, y1=tbl_smry$u95, code=3, angle=90, length=arrowlength, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(boxwidth)) {
    rect(xleft=tbl_smry$x-boxwidth, xright=tbl_smry$x+boxwidth, ybottom=tbl_smry$q25, ytop=tbl_smry$q75, border=tbl_smry$color, lwd=1.5, col=NA)
    segments(x0=tbl_smry$x-boxwidth, x1=tbl_smry$x+boxwidth, y0=tbl_smry$median, col=tbl_smry$color, lwd=1.5)
  }
  
  if (!is.na(polygon)) {
    for (clr in unique(tbl_smry$color)) {
      if (polygon=='iqr') {
        subs = subset(tbl_smry, color==clr & !is.na(q25) & !is.na(q75))
        points( x=  subs$x, y=  subs$median, type='l', lwd=2, col=subs$color)
        polygon(x=c(subs$x, rev(subs$x)), y=c(subs$q25, rev(subs$q75)), col=alpha(subs$color, ci_alpha), border=NA)
      } else if (polygon=='ci') {
        subs = subset(tbl_smry, color==clr & !is.na(l95) & !is.na(u95))
        points( x=  subs$x, y=  subs$mean, type='l', lwd=2, col=subs$color)
        polygon(x=c(subs$x, rev(subs$x)), y=c(subs$l95, rev(subs$u95)), col=alpha(subs$color, ci_alpha), border=NA)
      }
    }
  }
  
  if (!is.na(test)) {
    testfun = get(test) # e.g. ks.test
    control_color = control_group
    tbl_smry$p = as.numeric(NA)
    for (i in 1:nrow(tbl_smry)) {
      this_x = tbl_smry$x[i]
      this_rows = tbl$x == this_x & tbl$color == tbl_smry$color[i]
      ctrl_rows = tbl$x == this_x & tbl$color == control_group
      test_obj = suppressWarnings(testfun(tbl$y[this_rows],  tbl$y[ctrl_rows]))
      tbl_smry$p[i] = test_obj$p.value
    }
    tbl_smry$p_symb = ''
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.05] = '*'
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.01] = '**'
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.001] = '***'
    text(x=tbl_smry$x[tbl_smry$color != control_color], y=max(ylims)*.95, labels=tbl_smry$p_symb[tbl_smry$color != control_color])
  }
  
  if (!is.null(legtext)) {
    par(xpd=T)
    legend(x=max(xlims),y=max(ylims),legtext,col=legcol,text.col=legtextcol,pch=legpch,lwd=leglwd, cex=legcex,lty=leglty, bty='n')
    par(xpd=F)
  }
  
  return(tbl_smry)
  
}






##############
# OUTPUT STREAMS
##############



tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/table-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}









#####
# FIGURE 1
#####

tell_user('done.\nCreating Figure 1...')

resx=300
png('display_items/figure-1.png',width=6.5*resx,height=4.5*resx,res=resx)

par(mfrow=c(1,2), mar=c(4,5,2.5,1))
panel = 1

tmt_pep = read_tsv('data/tmt/tmt_1622_peptides.tsv', col_types=cols())
tmt_meta = read_tsv('data/tmt/tmt_1622_meta.tsv', col_types=cols())

tmt_pep %>%
  select(scan_number, gene_symbol, peptide_sequence, protein_id, a:j) %>%
  pivot_longer(cols=a:j) %>%
  inner_join(tmt_meta, by=c('name'='sample')) -> pep

pep %>%
  group_by(scan_number, gene_symbol, peptide_sequence, protein_id, lp_order, description, color) %>%
  summarize(.groups='keep', 
            meanval = mean(value),
            techrep_cv = sd(value)/mean(value)) %>%
  ungroup() -> pep_means

pep_means %>%
  group_by(scan_number, gene_symbol, peptide_sequence, protein_id) %>%
  summarize(.groups='keep', 
            mean_abundance = mean(meanval),
            mean_techrep_cv = mean(techrep_cv)) %>%
  ungroup() -> pep_stats

# sum(pep_stats$mean_techrep_cv < 0.25)

pep_means %>%
  inner_join(pep_stats, by=c('scan_number', 'gene_symbol', 'peptide_sequence', 'protein_id')) %>%
  filter(mean_techrep_cv < 0.25) %>%
  group_by(gene_symbol) %>%
  summarize(.groups='keep', 
            n_peptides = n(),
            baseline_cv = sd(meanval[lp_order %in% 1:3]) / mean(meanval[lp_order %in% 1:3]),
            baseline_abundance = mean(meanval[lp_order %in% 1:3]),
            ratio1 = mean(meanval[lp_order==1]/mean(meanval[lp_order %in% 1:3])),
            ratio2 = mean(meanval[lp_order==2]/mean(meanval[lp_order %in% 1:3])),
            ratio3 = mean(meanval[lp_order==3]/mean(meanval[lp_order %in% 1:3])),
            ratio4 = mean(meanval[lp_order==4]/mean(meanval[lp_order %in% 1:3])),
            ratio5 = mean(meanval[lp_order==5]/mean(meanval[lp_order %in% 1:3]))) %>%
  ungroup() %>% 
  filter(baseline_cv <= 1) -> prot_ratios



xlims = c(0, 1.0)
ylims = c(.07, 15)
yats = rep(1:9, 3) * 10^(rep(-1:1, each=9))
ybigs = 10^(-1:1)
xats = seq(0, 1.5, .1)
xbigs = seq(0, 1.5, .5)

plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, labels=c('0%','50%','100%','150%'), tck=-0.05, lwd=0, line=-0.5, cex.axis=.7)
mtext(side=1, line=1.5, text='baseline %CV')
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=c('0.1', '1', '10'), tck=-0.05, lwd=0, line=-0.25, las=2)
mtext(side=2, line=2.25, text='fold change from baseline\nat LP 4 days post-minocycline')
points(prot_ratios$baseline_cv, prot_ratios$ratio4, col=alpha('#E38217',pt_alpha), pch=20)
callouts = tibble(gene_symbol=c('NEFL','NEFM','TREM2','AIF1','CRP','VIM'),
                  pos=c(2,4,3,1,2,2))
callouts %>%
  inner_join(prot_ratios, by='gene_symbol') -> callout_plot
points(x=callout_plot$baseline_cv, y=callout_plot$ratio4, pch=1, col='#000000')
text(x=callout_plot$baseline_cv, y=callout_plot$ratio4, pos=callout_plot$pos, labels=callout_plot$gene_symbol, font=3, cex=.8)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(pep_stats,   'Peptide metrics for TMT in CSF from index case.')
write_supp_table(prot_ratios, 'Protein abundance quantified by TMT in serial CSF samples from index case.')


par(mar=c(4,3,2.5,1))

prot_ratios %>%
  select(gene_symbol, ratio1:ratio5) %>%
  pivot_longer(cols=ratio1:ratio5) %>%
  mutate(lp_order = as.integer(gsub('ratio','',name)),
         ratio = value) %>%
  select(-value,-name) -> time_series

xlims = c(1, 5)
ylims = c(.07, 15)
yats = rep(1:9, 3) * 10^(rep(-1:1, each=9))
ybigs = 10^(-1:1)
xats = NULL
xbigs = 1:5

plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=1, at=xbigs, tck=-0.025, labels=NA)
axis(side=1, at=tmt_meta$lp_order, labels=gsub(' ','\n',tmt_meta$description), padj=1, lwd=0, line=-1, cex.axis=.7)
mtext(side=1, line=2.5, text='LP order')
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=c('0.1', '1', '10'), tck=-0.05, lwd=0, line=-0.25, las=2)
mtext(side=2, line=2.25, text='fold change from baseline')
for (gene in unique(time_series$gene_symbol)) {
  subs = subset(time_series, gene_symbol==gene)
  points(subs$lp_order, subs$ratio, type='l', col='#E38217', lwd=.75)
}

par(xpd=T)
callouts = tibble(gene_symbol=c('NEFL','NEFM','TREM2','AIF1','VIM','CRP','CRP'),
                  lp_order = c(4,4,4,4,4,5,2),
                  pos=c(2,2,3,1,1,1,3))
callouts %>%
  inner_join(time_series, by=c('gene_symbol','lp_order')) -> callout_plot
points(x=callout_plot$lp_order, y=callout_plot$ratio, pch=1, col='#000000', font=3)
text(x=callout_plot$lp_order, y=callout_plot$ratio, pos=callout_plot$pos, labels=callout_plot$gene_symbol, font=3, cex=.8)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

unnecessary_message_end_figure_1 =dev.off()






#####
# FIGURE 2
#####

tell_user('done.\nCreating Figure 2...')

resx=300
png('display_items/figure-2.png',width=6.5*resx,height=3.5*resx,res=resx)

par(mfrow=c(1,2))
panel = 1

min_days_on = 0
max_included_days_off = 7
max_days_assumed_on = 30

patient = read_tsv('data/bwh/deid_patient_data.tsv', col_types=cols()) %>%
  rename(age_bin = age) %>%
  mutate(age = unbin_age(age_bin))
nfl = read_tsv('data/bwh/deid_nfl_data.tsv', col_types=cols())

patient %>%
  group_by(query_type, drug) %>%
  summarize(.groups='keep', 
            n_patients = length(unique(subject_id))) %>%
  ungroup() %>%
  mutate(drug = replace_na(drug, 'other/unknown')) %>%
  arrange(desc(query_type), desc(n_patients)) -> query_smry

write_supp_table(query_smry, 'Number of patients by query type and drug prescribed.')

patient %>%
  filter(query_type == 'derms') %>%
  mutate(use_med_end = case_when(!is.na(med_end) ~ med_end,
                                 TRUE ~ med_start + max_days_assumed_on)) %>%
  inner_join(nfl, by='subject_id', relationship='many-to-many') %>%
  filter(sample_date > med_start + min_days_on & (sample_date <= use_med_end + max_included_days_off)) %>%
  mutate(days_since_start = sample_date - med_start) %>%
  select(subject_id, query_type, sex, age_bin, age, drug, route, daily_dose, nfl, days_since_start) -> derms_nfl

patient %>%
  filter(query_type == 'wells') %>%
  inner_join(nfl, by='subject_id', relationship='many-to-many') %>%
  mutate(days_since_start = as.integer(NA)) %>%
  select(subject_id, query_type, sex, age_bin, age, drug, route, daily_dose, nfl, days_since_start) -> wells_nfl

rbind(derms_nfl, wells_nfl) -> derms_wells_nfl

derms_wells_nfl %>%
  filter(drug != 'azithromycin') %>%
  filter(route=='Oral' | drug=='clindamycin' | query_type == 'wells') %>%
  filter(days_since_start >= min_days_on | query_type=='wells') %>%
  group_by(subject_id, sex, age_bin, age, drug) %>%
  summarize(.groups='keep', nfl = mean(nfl)) %>%
  ungroup() -> bwh_dw_indiv_means

write_supp_table(bwh_dw_indiv_means %>% select(-age), 'Mean plasma NfL concentrations for well visit and dermatology patients.')

bwh_dw_indiv_means$drug = factor(bwh_dw_indiv_means$drug, levels=c('well visit control','clindamycin','isotretinoin','doxycycline','minocycline'))
bwh_dw_indiv_means %>%
  arrange(drug) -> bwh_dw_indiv_means

m = lm(log(nfl) ~ drug + age, data=bwh_dw_indiv_means)
dw_coefs = summary(m)$coefficients  %>% 
  as_tibble(rownames='parameter') %>%
  clean_names() 
write_supp_table(dw_coefs, 'Log-linear model coefficients for well visit and dermatology patients.')

meta = read_tsv('data/bwh/bwh_meta.tsv', col_types=cols())

bwh_dw_indiv_means$color = meta$color[match(bwh_dw_indiv_means$drug, meta$drug)]


par(mar=c(4,4,3,1))
xlims = c(18, 82)
xbigs = c(18, 3:9*10)
xats = c(18,seq(22.5,78.5,5),82)
xatlabs = bin_age(xats)
ylims = c(0,160)
ybigs = seq(0,175,25)
yats = seq(0,175,5)
pt_alpha = .7
plot(NA, NA, xlim=xlims, ylim=ylims, yaxs='i', axes=F, ann=F)
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xats, lwd=0, labels=xatlabs, las=2, cex.axis=0.8, line=-0.5)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, lwd=0, line=-0.25, las=2)
par(xpd=T)
points(x=bwh_dw_indiv_means$age, y=bwh_dw_indiv_means$nfl, col=alpha(bwh_dw_indiv_means$color,pt_alpha), pch=19)
par(xpd=F)
x = seq(18,80,1)
for (this_drug in meta$drug[meta$drug != 'minocycline']) {
  this_color = meta$color[meta$drug==this_drug]
  y = exp(predict.lm(m, newdata=list(age=x, drug=rep(this_drug,length(x)))))
  points(x, y, col=this_color, type='l', lwd=1.5)
}
mtext(side=1, line=2.5, text='age')
mtext(side=2, line=2.5, text='plasma NfL (pg/mL)')
par(xpd=T)
legend(x=min(xlims),y=max(ylims)*1.1,meta$disp,col=alpha(meta$color,pt_alpha),pch=19,cex=.7,bty='n')
par(xpd=F)

mtext(side=3, line=1, text='dermatology and well visits')
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


write_stats("Linear model predicted baseline NfL at age 18: ",formatC(y[x==18],format='g',digits=3))

patient %>%
  filter(query_type=='inpt') %>% 
  filter(drug=='minocycline') %>% # some of these patients also took other drugs; we only care about their minocycline usage dates
  mutate(use_med_end = case_when(!is.na(med_end) ~ med_end,
                                 TRUE ~ med_start + max_days_assumed_on)) %>%
  inner_join(nfl, by='subject_id', relationship='many-to-many') %>%
  mutate(on = sample_date > med_start & (sample_date <= use_med_end + max_included_days_off)) %>%
  mutate(days_since_start = sample_date - med_start) %>%
  mutate(days_since_end = sample_date - use_med_end) %>%
  select(subject_id, query_type, sex, age_bin, drug, route, daily_dose, nfl, on, sample_date, days_since_start, days_since_end) %>%
  group_by(subject_id, sample_date) %>%
  arrange(desc(on), days_since_start) %>% # there will be many-to-many match onto both current and expired medication orders; desc(on) puts TRUE first, prioritizing active orders; asc(days_since_start) prioritizes most recent course of drug
  slice(1) -> inpt_nfl


inpt_details = read_tsv('data/bwh/deid_inpatient_details.tsv', col_types=cols()) %>%
  filter(subject_id %in% inpt_nfl$subject_id) %>%
  select(subject_id, summary_of_most_severe_health_problems, indication_for_minocycline, day_of_death)

write_supp_table(inpt_details, 'Health problems and indications for minocycline in hospitalized patients.')

inpt_nfl %>%
  mutate(start = sample_date - days_since_start) %>%
  group_by(subject_id) %>%
  summarize(.groups='keep', new_zero = min(start)-1) %>%
  ungroup() -> new_zeroes

inpt_nfl$new_zero = new_zeroes$new_zero[match(inpt_nfl$subject_id, new_zeroes$subject_id)]
inpt_nfl$plot_day = inpt_nfl$sample_date - inpt_nfl$new_zero

on_color = meta$color[meta$drug=='minocycline']
off_color = '#787878'
xlims = range(inpt_nfl$plot_day) + c(0,5)
xbigs = seq(0,180,30)
xats = seq(0,180,10)
ylims = c(0, 3000)
ybigs = seq(0,3000,500)
yats = seq(0,3000,100)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='day')
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=formatC(ybigs, big.mark=','), tck=-0.05, las=2, lwd=0, line=-0.25)
mtext(side=2, line=3.0, text='plasma NfL (pg/mL)')
for (this_subj in unique(inpt_nfl$subject_id)) {
  subs = subset(inpt_nfl, subject_id==this_subj)
  points(subs$plot_day, subs$nfl, pch=19, col=alpha(ifelse(subs$on, on_color, off_color), pt_alpha))
  points(subs$plot_day, subs$nfl, type='l', lwd=0.5, col='#000000')
}
par(xpd=T)
legend(x=mean(xlims),y=max(ylims)*1.1,c('on minocycline','off minocycline'),col=alpha(c(on_color, off_color),pt_alpha),pch=19,cex=.7,bty='n')
par(xpd=F)

mtext(side=3, line=1, text='inpatients')

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write_supp_table(inpt_nfl, 'Serial inpatient NfL measurements.')

# reviewed partially paired methods; t test is less appropriate because NfL is non-normal; decided on partially paired wilcoxon
# https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.865.734&rep=rep1&type=pdf
# https://stats.stackexchange.com/questions/25941/t-test-for-partially-paired-and-partially-unpaired-data
# https://rdrr.io/cran/robustrank/man/pm.wilcox.test.html

inpt_nfl %>%
  group_by(subject_id, on) %>%
  summarize(.groups='keep', mean_nfl = mean(nfl)) %>%
  ungroup() -> inpt_onoff_means

inpt_onoff_means %>%
  mutate(mino = ifelse(on, 'on', 'off')) %>%
  select(subject_id, mino, mean_nfl) %>%
  pivot_wider(names_from=mino, values_from=mean_nfl) -> wider

pm_wilcoxon = pm.wilcox.test(wider$on, wider$off, alternative='two.sided')

write_stats('Partially paired Wilcoxon for inpatients on/off minocycline: P = ',formatC(pm_wilcoxon$p.value, format='g', digits=2))

# https://stackoverflow.com/a/17042397/3806692
u = par("usr")
v = c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v = c( v[1] + (v[2]-v[1])*.75, v[2], v[3] + (v[4]-v[3])*.25, v[3] + (v[4]-v[3])*.75 )
par( fig=v, new=TRUE, mar=c(0,0,0,0) )

xlims = c(0.5, 2.5)
ylims = c(0,3000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
points(x=rep(1, nrow(wider)), y=wider$on, col=alpha(on_color,pt_alpha),  pch=19)
points(x=rep(2, nrow(wider)), y=wider$off, col=alpha(off_color,pt_alpha),  pch=19)
segments(x0=1, x1=2, y0=wider$on, y1=wider$off, col='#000000', lwd=.5)
axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
mtext(side=3, line=0.1, text='means', cex=0.8)
mtext(side=1, at=c(1,2), text=c('on','off'), col=c(on_color, off_color))
box()

write_supp_table(wider, 'Inpatient NfL averaged by on/off minocycline status for each individual.')


unnecessary_message_end_figure_2 = dev.off()











#####
# FIGURE S1
#####

tell_user('done.\nCreating Figure S1...')

resx=300
png('display_items/figure-s1.png',width=6.5*resx,height=4.5*resx,res=resx)


layout_matrix = matrix(c(1:20,rep(21,5)), nrow=5, byrow=T)
layout(layout_matrix, heights=c(1,1,1,1,.5))

par(mar=c(2,3,2,1))
for (subj in unique(inpt_nfl$subject_id)) {
  subs = subset(inpt_nfl, subject_id==subj)
  anno = subset(patient, subject_id==subj & drug=='minocycline')
  anno$color = case_when(anno$route=='Oral' ~ '#FF000015',
                         anno$route=='Intravenous' ~ '#FF000088')
  hist = subset(inpt_details, subject_id==subj)
  
  xlims = c(0, max(c(subs$sample_date, anno$med_end, hist$day_of_death), na.rm=T)+5)
  if (max(xlims) > 100) {
    xats = seq(0,500,7)
    xbigs = seq(0,500,28)
  } else {
    xats = 0:100
    xbigs = seq(0,100,7)
  }
  ylims = c(0, max(c(subs$nfl*1.1,100), na.rm=T))
  if (max(ylims) > 1000) {
    yats = 0:30*100
    ybigs = seq(0,3000,500)
  } else {
    yats = 0:100*10
    ybigs = seq(0,1000,100)
  }
  
  plot(NA, NA, xlim = xlims, ylim = ylims, yaxt='n', xaxt='n', yaxs='i', xaxs='i', axes=F, ann=F)
  points(x=subs$sample_date, y=subs$nfl, pch=20, cex=.6)
  points(x=subs$sample_date, y=subs$nfl, type='l', lwd=0.5)
  axis(side=1, labels=NA, tck=-0.025, at=xats)
  axis(side=1, tck=-0.05, at=xbigs, labels=NA)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-.75)
  mtext(side=1, line=1.0, text='day', cex=0.6)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-.75, labels=ybigs, cex.axis=0.7)
  mtext(side=2, line=1.55, text='NfL', cex=0.6)
  
  if(!is.na(hist$day_of_death[hist$subject_id==subj])) {
    death_day = as.integer(hist$day_of_death[hist$subject_id==subj])
    abline(v=death_day, col='#3D59AB', lwd=2)
  } else {
    death_day = Inf
  }
  
  rect(xleft=anno$med_start, xright=pmin(anno$med_end, min(c(max(xlims), death_day), na.rm=T)), ybottom=rep(0.1, nrow(anno)), ytop=rep(1e4, nrow(anno)), col=anno$color, border=NA)
  mtext(side=3, line=0.125, text=paste(subj, anno$sex[1], anno$age_bin[1]), cex=0.7)
}

par(mar=c(0,0,0,0))
plot(NA, NA, xlim=0:1, ylim=0:1, ann=F, axes=F, xaxs='i', yaxs='i')
legend('center', 
       c('oral minocycline','iv minocycline','plasma NfL reading','patient death'),
       pt.bg=c('#FF000015','#FF000088','#000000','#3D59AB'),
       col=c('#000000','#000000','#000000','#3D59AB'),
       pch=c(22,22,20,124),
       pt.lwd=c(.5,.5,1,0),
       lwd=c(0,0,1,0),
       pt.cex=2,
       horiz=T, bty='n')

unnecessary_message_end_figure_2 = dev.off()









#####
# FIGURE 3
#####

tell_user('done.\nCreating Figure 3...')

resx=300
png('display_items/figure-3.png',width=6.5*resx,height=3.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,
                         1,4,5), byrow=T, nrow=2)
layout(layout_matrix, widths=c(4,1,1))

panel = 1


main = read_tsv('data/animals/CMR-2285-main.tsv', col_types=cols())
ella = read_tsv('data/ella/CMR-2285-ella.tsv', col_types=cols())
sex_meta = read_tsv('data/animals/sex_meta.tsv', col_types=cols())


ella %>%
  inner_join(main, by=c('animal')) %>%
  mutate(dob = as.Date(dob)) %>%
  mutate(days = as.integer(date - dob)) %>%
  mutate(weeks = floor(days/7)) %>%
  filter(!is.na(nfl)) %>%
  inner_join(sex_meta, by='sex') %>%
  select(animal, sex, color, pch, weeks, days, nfl) %>%
  arrange(animal, days) -> mashup

mashup %>%
  select(animal, sex, weeks, nfl) -> mashup_out
write_supp_table(mashup_out, 'Longitudinal NfL measurements in naive mice.')

sex_offset = .125
mashup %>%
  filter(!is.na(nfl)) %>%
  group_by(weeks, sex, color) %>%
  summarize(.groups='keep',
            n=n(),
            mean=mean(nfl),
            sd=sd(nfl),
            median=median(nfl),
            q25 = quantile(nfl, .25),
            q75 = quantile(nfl, .75)) %>%
  ungroup() %>%
  mutate(x = weeks + sex_offset * ifelse(sex=='M',-1,1)) -> long_smry


long_smry %>%
  select(-color, -x) -> long_smry_out
write_supp_table(long_smry_out, 'Summary of longitudinal NfL measurements in naive mice.')

par(mar=c(3,5,3,1))
xlims = c(5,24)
ylims = c(10,1500)
xats = 6:24
xbigs = seq(6,24,6)
yats = rep(1:9, 4) * 10^rep(0:3,each=9)
ybigs = c(1, 10, 100, 1000)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=1, at=xlims, lwd.ticks = 0, labels=NA)
axis(side=1, at=xats, tck=-0.0125, labels=NA)
axis(side=1, at=xbigs, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, lwd=0, line=-0.25)
mtext(side=1, line=1.5, text='age (weeks)')
axis(side=2, at=yats, tck=-0.0125, labels=NA)
axis(side=2, at=ybigs, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, las=2, lwd=0)
mtext(side=2, line=3, text='plasma NfL (pg/mL)')
points(mashup$weeks, mashup$nfl, col=alpha(mashup$color,ci_alpha), pch=mashup$pch)
for (anml in unique(mashup$animal)) {
  subs = subset(mashup, animal==anml)
  points(subs$weeks, subs$nfl, type='l', lwd=.125, col=subs$color)
}
boxwidth = 0.1
rect(xleft=long_smry$x-boxwidth, xright=long_smry$x+boxwidth, ybottom=long_smry$q25, ytop=long_smry$q75, border=long_smry$color, lwd=1.5, col=NA)
segments(x0=long_smry$x-boxwidth, x1=long_smry$x+boxwidth, y0=long_smry$median, col=long_smry$color, lwd=1.5)

x = seq(6,24,.1)
loess_f = loess(nfl ~ weeks, data=subset(mashup,sex=='F'), span = 1.5)
y_f = predict(loess_f, x)
points(x, y_f, type='l', lwd=5, col=alpha(sex_meta$color[sex_meta$sex=='F'],.7))

loess_m = loess(nfl ~ weeks, data=subset(mashup,sex=='M'), span = 1.5)
y_m = predict(loess_m, x)
points(x, y_m, type='l', lwd=5, col=alpha(sex_meta$color[sex_meta$sex=='M'],.7))

par(xpd=T)
legend(x=20,y=max(ylims),sex_meta$sex,pch=sex_meta$pch,col=sex_meta$color,bty='n')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



mashup %>%
  group_by(sex, color) %>%
  summarize(.groups='keep',
            mean = mean(nfl),
            mean_cv = sd(nfl)/mean(nfl)) %>%
  ungroup() %>%
  mutate(x = case_when(sex=='F' ~ 2, sex=='M' ~ 1)) -> overall_sex_smry

mashup %>%
  group_by(animal, sex, color) %>%
  summarize(.groups='keep',
            mean = mean(nfl),
            cv = sd(nfl)/mean(nfl)) %>%
  ungroup() %>%
  group_by(sex, color) %>%
  summarize(.groups='keep',
            mean_mean = mean(mean),
            mean_between_cv = sd(mean)/mean(mean),
            mean_within_cv = mean(cv)) %>%
  ungroup() %>%
  mutate(x = case_when(sex=='F' ~ 2, sex=='M' ~ 1)) -> within_between_sex_smry

overall_sex_smry %>%
  select(-x) %>%
  inner_join(within_between_sex_smry, by=c('sex','color')) %>%
  select(-mean_mean, -x, -color) %>%
  rename(mean_overall_cv = mean_cv) -> sex_smry_out

write_supp_table(sex_smry_out, 'Longitudinal NfL in naive mice summarized by sex.')

par(mar=c(2,5,3,1))
xlims = range(overall_sex_smry$x) + c(-0.5, 0.5)
ylims = c(0,200)
yats = 0:20*10
ybigs = c(0,100,200)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks = 0, labels=NA)
mtext(side=1, line=0.25, at=overall_sex_smry$x, text=overall_sex_smry$sex)
axis(side=2, at=yats, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=NA, tck=-0.2)
axis(side=2, at=ybigs, labels=ybigs, lwd=0, line=-0.5, las=2)
mtext(side=2, line=3, text='mean NfL (pg/mL)', cex=0.9)
barwidth = 0.3
rect(xleft=overall_sex_smry$x-barwidth, xright=overall_sex_smry$x+barwidth, ybottom=rep(0,nrow(overall_sex_smry)), ytop=overall_sex_smry$mean, col=overall_sex_smry$color, border=NA)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

xlims = range(overall_sex_smry$x) + c(-0.5, 0.5)
ylims = c(0,1.2)
yats = 0:20/10
ybigs = c(0,.5,1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks = 0, labels=NA)
mtext(side=1, line=0.25, at=overall_sex_smry$x, text=overall_sex_smry$sex)
axis(side=2, at=yats, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=NA, tck=-0.2)
axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, line=-0.5, las=2)
mtext(side=2, line=3, text='all-sample CV (%)', cex=0.9)
barwidth = 0.3
rect(xleft=overall_sex_smry$x-barwidth, xright=overall_sex_smry$x+barwidth, ybottom=rep(0,nrow(overall_sex_smry)), ytop=overall_sex_smry$mean_cv, col=overall_sex_smry$color, border=NA)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



xlims = range(within_between_sex_smry$x) + c(-0.5, 0.5)
ylims = c(0,1.2)
yats = 0:20/10
ybigs = c(0,.5,1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks = 0, labels=NA)
mtext(side=1, line=0.25, at=within_between_sex_smry$x, text=within_between_sex_smry$sex)
axis(side=2, at=yats, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=NA, tck=-0.2)
axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, line=-0.5, las=2)
mtext(side=2, line=3, text='between-animal CV (%)', cex=0.9)
barwidth = 0.3
rect(xleft=within_between_sex_smry$x-barwidth, xright=within_between_sex_smry$x+barwidth, ybottom=rep(0,nrow(within_between_sex_smry)), ytop=within_between_sex_smry$mean_between_cv, col=within_between_sex_smry$color, border=NA)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


xlims = range(within_between_sex_smry$x) + c(-0.5, 0.5)
ylims = c(0,1.2)
yats = 0:20/10
ybigs = c(0,.5,1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks = 0, labels=NA)
mtext(side=1, line=0.25, at=within_between_sex_smry$x, text=within_between_sex_smry$sex)
axis(side=2, at=yats, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, labels=NA, tck=-0.2)
axis(side=2, at=ybigs, labels=percent(ybigs), lwd=0, line=-0.5, las=2)
mtext(side=2, line=3, text='within-animal CV (%)', cex=0.9)
barwidth = 0.3
rect(xleft=within_between_sex_smry$x-barwidth, xright=within_between_sex_smry$x+barwidth, ybottom=rep(0,nrow(within_between_sex_smry)), ytop=within_between_sex_smry$mean_within_cv, col=within_between_sex_smry$color, border=NA)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

unnecessary_message_end_figure_3 =dev.off()





######
# FIGURE 4
######

tell_user('done.\nCreating Figure 4...')

resx=300
png('display_items/figure-4.png', width=6.5*resx, height=8*resx, res=resx)

layout_matrix = matrix(c(1,2,
                         3,4,
                         5,6), nrow=3, byrow=T)
# layout_matrix = matrix(c(1,1,2,3,
#                          4,5,5,6), nrow=2, byrow=T)
layout(layout_matrix)
panel = 1


# Minocycline pilot study in mice EF-0047-DA-ME
meta = read_tsv('data/animals/EF-0047-DA-ME-meta.tsv', col_types=cols())
main = read_tsv('data/animals/EF-0047-DA-ME-main.tsv', col_types=cols())
ella = read_tsv('data/ella/EF-0047-DA-ME-ella.tsv', col_types=cols()) 

predose_color = '#A9A9E2'

ella %>%
  inner_join(main, by='animal') %>%
  inner_join(meta, by='group') %>%
  select(animal, group, drug, day, nfl, color) %>%
  filter(!is.na(nfl)) %>%
  mutate(pch = ifelse(day > 0 & day < 14, 19, 1),
         ptcol = ifelse(day == 0, predose_color, color),
         grouped_group = ifelse(day==0, 'predose', drug),
         grouped_color = ifelse(day==0, predose_color, color)) %>%
  mutate(grouped_disp = ifelse(day %in% c(16,20), '16-20', as.character(day))) %>%
  mutate(grouped_day = ifelse(day %in% c(16,20), 20, day)) -> mp_nfl

mp_nfl %>%
  select(animal, group, drug, day, nfl) -> mp_nfl_out
write_supp_table(mp_nfl_out, 'Individual NfL values in mouse treatment study 1')

par(mar=c(3,4,3,1))
xlims=c(-3, 45)
ylims=c(20,3000)
yats=rep(1:9,4) * rep(10^(1:4),each=9)
ybigs=10^(1:4)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.3, labels=ybigs, cex.axis=0.8)
mtext(side=2, line=2.5, text='plasma NfL (pg/mL)')
par(xpd=T)
points(mp_nfl$grouped_day, mp_nfl$nfl, col=alpha(mp_nfl$ptcol, ci_alpha), pch=mp_nfl$pch)
par(xpd=F)
for (anml in unique(mp_nfl$animal)) {
  subs = subset(mp_nfl, animal==anml)
  points(subs$grouped_day, subs$nfl, col=subs$color, type='l', lwd=0.25)
}


mp_nfl %>%
  group_by(grouped_day) %>%
  summarize(.groups='keep',
            n = n(),
            ratio = mean(nfl[drug=='minocycline'])/mean(nfl[drug=='none']),
            wilcoxon_p = suppressWarnings(wilcox.test(nfl[drug=='minocycline'], nfl[drug=='none'])$p.value)) %>%
  ungroup() %>%
  mutate(psymb = p_to_symbol(wilcoxon_p)) -> mp_wilcoxons
  
mp_nfl %>%
  group_by(grouped_day, grouped_disp, grouped_group, grouped_color) %>%
  summarize(.groups='keep', 
            n=n(),
            cv = sd(nfl) / mean(nfl),
            median = median(nfl),
            q25 = quantile(nfl, .25),
            q75 = quantile(nfl, .75),
            mean = mean(nfl),
            l95 = lower(nfl),
            u95 = upper(nfl)) %>%
  ungroup() -> mp_smry
boxwidth = 2
rect(xleft=mp_smry$grouped_day-boxwidth, xright=mp_smry$grouped_day+boxwidth, ybottom=mp_smry$q25, ytop=mp_smry$q75, border=mp_smry$grouped_color, lwd=1.5, col=NA)
segments(x0=mp_smry$grouped_day-boxwidth, x1=mp_smry$grouped_day+boxwidth, y0=mp_smry$median, col=mp_smry$grouped_color, lwd=1.5)
mtext(side=3, line=0, at=mp_wilcoxons$grouped_day, text=mp_wilcoxons$psymb, cex=1.0)

washout_color = '#FCE412'
mp_smry %>%
  distinct(grouped_day, grouped_disp) -> xleg
ybot = 8
ytop = 12
par(xpd=T)
rect(xleft=0.5, xright = 15, ybottom = ybot, ytop=ytop, col=alpha(meta$color[meta$drug=='minocycline'], 0.3), border=NA)
text(x=(0.5+15)/2, y=sqrt(ybot*ytop), labels='treatment', col=meta$color[meta$drug=='minocycline'], cex=0.8)
rect(xleft=15, xright = 42, ybottom = ybot, ytop=ytop, col=alpha(washout_color, 0.3), border=NA)
text(x=(42+15)/2, y=sqrt(ybot*ytop), labels='washout', col='black', cex=0.8)
par(xpd=F)

axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, adj=-0.17, text='day:', cex=0.8)
mtext(side=1, at=xleg$grouped_day, text=xleg$grouped_disp, cex=0.8)
#axis(side=3, line=0.5, at=c(0.5, 13.5), tck=0.05, col=meta$color[meta$drug=='minocycline'], lwd=2, labels=NA)
#mtext(side=3, line=0.5, at=7.0, text='treatment', col=meta$color[meta$drug=='minocycline'])

mp_nfl %>%
  group_by(grouped_group, grouped_color, pch, ptcol) %>%
  summarize(.groups='keep', minday = min(grouped_day)) %>%
  ungroup() %>%
  arrange(minday, desc(grouped_group)) %>%
  mutate(disp = case_when(minday==20 ~ paste0(grouped_group, ' washout'),
                          TRUE ~ grouped_group)) -> leg

par(xpd=T)
legend(x=28, y=5000, leg$disp, col=leg$ptcol, pch=leg$pch, bty='n', cex=0.8)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

mp_smry %>%
  select(grouped_day, grouped_group, n, cv, median, q25, q75, mean, l95, u95) -> mp_smry_out
write_supp_table(mp_smry_out, 'Summarized NfL values in mouse treatment study 1')
write_supp_table(mp_wilcoxons, 'Wilcoxon tests on NfL values in mouse treatment study 1')


# this is for Figure S2 but placing it here so that supp tables are in order:
weights = read_tsv('data/animals/EF-0047-DA-ME-weights.tsv', col_types=cols())
weights %>%
  filter(day==1) %>%
  rename(baseline=weight_g) %>% 
  select(animal, baseline) -> baselines
weights %>%
  inner_join(baselines, by='animal') %>%
  mutate(delta = (weight_g - baseline)/baseline) -> weights_rel_study1
weights_rel_study1 %>%
  group_by(day) %>%
  summarize(.groups='keep', 
            n_animals = n(),
            mean = mean(delta),
            l95 = lower(delta),
            u95 = upper(delta)) %>%
  ungroup() %>%
  filter(!is.na(l95)) -> weights_rel_avg_study1
weights_rel_study1 %>%
  select(animal, day, weight_g, delta) %>%
  arrange(animal, day) -> weights_rel_out
write_supp_table(weights_rel_out, 'Individual treated animal weights in mouse treatment study 1')
write_supp_table(weights_rel_avg_study1, 'Summary of treated animal weights in mouse treatment study 1')


mg_genes = read_tsv('data/tmt/microglia_genes.tsv', col_types=cols())
mg_genes$mg_enriched = mg_genes$microglia_p_value < 0.05 & mg_genes$microglia_average > 0.1 & mg_genes$microglia_log2_fold_change > 2

deg = read_tsv('data/rnaseq/mp_deseq2_dge.tsv', col_types=cols()) %>%
  left_join(select(mg_genes, feature_name, mg_enriched), by=c('gene'='feature_name')) %>%
  mutate(mg_enriched = replace_na(mg_enriched, F)) %>%
  mutate(color = case_when(pvalue >= 0.05 ~ '#C9C9C999',
                           pvalue < 0.05 & padj >= 0.05 & l2fc < 0 ~ '#FFCCCC99',
                           pvalue < 0.05 & padj >= 0.05 & l2fc > 0 ~ '#CCCCFF99',
                           padj < 0.05 & l2fc < 0 ~ '#BB111199',
                           padj < 0.05 & l2fc > 0 ~ '#1111BB99')) %>%
  arrange(padj)
deg %>%
  select(gene, base_mean, l2fc, lfc_se, stat, pvalue, padj, mg_enriched) %>%
  arrange(padj) -> deg_out
write_supp_table(deg_out, 'RNA-seq differential gene expression analysis from mouse treamtent study 1')










# Dose-response mino & iso study EF-0048-DA-ME
meta = read_tsv('data/animals/EF-0048-DA-ME-meta.tsv', col_types=cols())
main = read_tsv('data/animals/EF-0048-DA-ME-main.tsv', col_types=cols())
ella = read_tsv('data/ella/EF-0048-DA-ME-ella.tsv', col_types=cols()) 


ella %>%
  mutate(animal = as.integer(gsub('-','',substr(sample_name, 4,5)))) %>%
  inner_join(main, by='animal') %>%
  inner_join(meta, by='drug') %>%
  select(animal, group, drug, dose_mgkg, nfl, color) -> dr_nfl

dr_nfl %>%
  select(animal, group, drug, dose_mgkg, nfl) -> dr_nfl_out
write_supp_table(dr_nfl_out, 'Individual NfL values in mouse treatment study 2')

main %>%
  inner_join(meta, by='drug') %>%
  distinct(group, drug, dose_mgkg, color) -> xinfo

xinfo %>%
  group_by(drug, color) %>%
  summarize(.groups='keep',
            xmin = min(group),
            xmax = max(group),
            xmid = mean(group)) %>%
  ungroup() -> xtranches

ella %>%
  mutate(animal = as.integer(gsub('-','',substr(sample_name, 4,5)))) %>%
  inner_join(main, by='animal') %>%
  inner_join(meta, by='drug') %>%
  select(animal, group, drug, dose_mgkg, nfl, color) -> dr_nfl

control_values = dr_nfl$nfl[dr_nfl$drug=='saline']

dr_nfl %>%
  group_by(group, drug, dose_mgkg) %>%
  summarize(.groups='keep', wilcoxon_p = suppressWarnings(wilcox.test(control_values, nfl)$p.value)) %>%
  ungroup() -> wilcoxons

dr_smry = dviz(dr_nfl, xlims=c(0.5, 9.5), ylims=c(20,3000),
               xvar='group', yvar='nfl', colorvar='color', xcols=c('drug','dose_mgkg'),
               xats=NULL, xbigs=NULL, yats=rep(1:9,4) * rep(10^(1:4),each=9), ybigs=10^(1:4),
               xlab = '', ylab='plasma NfL (pg/mL)', xlwds=0, yaxcex=0.8,
               log='y', boxwidth=0.33, mar=c(3,4,3,1))
mtext(side=1, line=-0.5, at=xinfo$group, text=paste0(xinfo$dose_mgkg,'\n'), padj=1, cex=0.6)
mtext(side=1, line=-0.5, adj=-0.17, padj=1, text='mg/kg:\n', cex=0.6)
tranche_line = 1.25
overhang = .4
for (i in 1:nrow(xtranches)) {
  axis(side=1, line=tranche_line, tck=0.025, at=unlist(xtranches[i,c('xmin','xmax')]) + overhang*c(-1,1), labels=NA)
  mtext(side=1, line=tranche_line, at=xtranches$xmid[i], text=xtranches$drug[i], col=xtranches$color[i], cex=0.8)
}

dr_smry$wilcoxon_p = wilcoxons$wilcoxon_p[match(dr_smry$x, wilcoxons$group)]
dr_smry$psymb = p_to_symbol(dr_smry$wilcoxon_p)
mtext(side=3, line=0, at=dr_smry$x, text=dr_smry$psymb, cex=1.0)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

dr_smry %>%
  select(x, drug, dose_mgkg, n, cv, median, q25, q75, mean, l95, u95, wilcoxon_p) -> dr_smry_out
write_supp_table(dr_smry_out, 'Summary of NfL values in mouse treatment study 2')



meta = read_tsv('data/animals/CMR-1982-meta.tsv', col_types=cols())
main = read_tsv('data/animals/CMR-1982-main.tsv', col_types=cols()) %>%
  mutate(animal = as.character(animal))
ella = read_tsv('data/ella/CMR-1982-ella.tsv', col_types=cols()) %>%
  mutate(animal = as.character(animal))

ella %>%
  inner_join(main, by='animal') %>%
  inner_join(meta, by=c('tx','day')) %>%
  select(animal, cohort, sex, tx, day, nfl, x, color, pch) -> cmr_1982_nfl

cmr_1982_nfl %>%
  rename(drug=tx) %>%
  select(animal, sex, drug, day, nfl) %>%
  arrange(day, desc(drug), animal) -> cmr1982_nfl_out
write_supp_table(cmr1982_nfl_out, 'Individual NfL values in mouse treatment study 3')


cmr_1982_smry = dviz(cmr_1982_nfl, xlims=c(-2,7), ylims=c(20,3000),
                     xvar='x', yvar='nfl', colorvar='color', pchvar='pch', xcols=c('tx','day','cohort'),
                     xats=NULL, xbigs=NULL, yats=rep(1:9,4) * rep(10^(1:4),each=9), ybigs=10^(1:4),
                     xlab = '', ylab='plasma NfL (pg/mL)', xlwds=0, yaxcex=0.8,
                     log='y', boxwidth=0.33, mar=c(3,4,3,1))
meta %>% distinct(x, disp, color) -> xleg
mtext(side=1, line=0.25, at=xleg$x, text=xleg$disp, col=xleg$color, cex=0.8)
for (anml in unique(cmr_1982_nfl$animal)) {
  subs = subset(cmr_1982_nfl, animal==anml)
  points(subs$x, subs$nfl, col=subs$color, type='l', lwd=0.25)
}


meta %>%
  filter(day != 0) %>%
  group_by(day) %>%
  summarize(.groups='keep', xmin=min(x), xmax=max(x)) %>%
  ungroup() -> xtranches
# mtext(side=1, line=1.0, at=xtranches$day, text=xtranches$day, col='#000000', cex=0.8)
# mtext(side=1, line=1.5, adj=-0.17,  text='day:', cex=0.8)
tranche_line = 1.5
overhang = .8
for (i in 1:nrow(xtranches)) {
  axis(side=1, line=tranche_line, tck=0.025, at=unlist(xtranches[i,c('xmin','xmax')]) + overhang*c(-1,1), labels=NA)
  mtext(side=1, line=tranche_line, at=xtranches$day[i], text=paste0('day ',xtranches$day[i]), col='#000000', cex=0.8)
}

cmr_1982_nfl %>%
  group_by(day) %>%
  summarize(.groups='keep', 
            wilcoxon_p = suppressWarnings(wilcox.test(nfl[tx=='mino'], nfl[tx=='none'])$p.value)) %>%
  ungroup() -> wilcoxons
cmr_1982_smry$wilcoxon_p = wilcoxons$wilcoxon_p[match(cmr_1982_smry$day, wilcoxons$day)]
cmr_1982_smry$wilcoxon_p[cmr_1982_smry$tx=='none'] = 1.0
cmr_1982_smry$psymb = p_to_symbol(cmr_1982_smry$wilcoxon_p)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

cmr_1982_smry %>%
  rename(drug=tx) %>%
  select(drug, day, cohort, n, median, q25, q75, mean, l95, u95, cv, wilcoxon_p) %>%
  arrange(day, desc(drug)) -> cmr_1982_smry_out
write_supp_table(cmr_1982_smry_out, 'Summary of NfL values in mouse treatment study 3')


meta = read_tsv('data/animals/CMR-2535-meta.tsv', col_types=cols())
main = read_tsv('data/animals/CMR-2535-main.tsv', col_types=cols()) %>%
  mutate(sex=ifelse(sex,'M','F'))
ella = read_tsv('data/animals/CMR-2535-ella.tsv', col_types=cols())

main %>%
  inner_join(ella, by=c('animal')) %>%
  select(animal, tx, day, nfl, gn_rs, sex, cohort) %>%
  inner_join(select(meta, -day), by=c('cohort','tx')) -> cmr2535

cmr2535 %>%
  rename(drug=tx) %>%
  select(animal, sex, day, drug, nfl) -> cmr2535_out

write_supp_table(cmr2535_out, 'Individual NfL values in mouse treatment study 4')

cmr2535 %>%
  filter(day==5) -> cmr2535_day5
meta %>%
  filter(day==5) -> xleg

cmr2535_day5_smry = dviz(cmr2535_day5, xlims=range(xleg$x) + c(-0.5, 0.5), ylims=c(20,3000), log='y',
                         xvar='x', yvar='nfl', colorvar='color', xcols=c('cohort','day','tx'),
                         xbigs=NULL, yats=rep(1:9,4) * rep(10^(1:4),each=9), ybigs=10^(1:4),
                         xlab='', ylab='plasma NfL (pg/mL)', xaxcex=0.8, ylabcex=1, xlwds=0, yaxcex=0.8, 
                         boxwidth=.33, mar=c(3,3,3,1))
mtext(side=1, line=0.25, at=xleg$x, text=xleg$tx, col=xleg$color)
wilcoxon_p = suppressWarnings(wilcox.test(cmr2535_day5$nfl[cmr2535_day5$tx=='mino'],cmr2535_day5$nfl[cmr2535_day5$tx=='mock'])$p.value)
mtext(side=3, at=xleg$x[xleg$tx=='mino'], line=0, text=p_to_symbol(wilcoxon_p))
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


cmr2535 %>%
  group_by(tx, day) %>%
  summarize(.groups='keep', 
            n=n(), 
            mean_nfl = mean(nfl), nfl_l95=lower(nfl), nfl_u95=upper(nfl),
            median_nfl = median(nfl), nfl_q25=quantile(nfl,.25), nfl_q75=quantile(nfl,.75)) %>%
  ungroup() -> cmr2535_smry


cmr2535 %>%
  group_by(day) %>%
  summarize(.groups='keep', 
            wilcoxon_p = suppressWarnings(wilcox.test(nfl[tx=='mino'], nfl[tx=='mock'])$p.value)) %>%
  ungroup() -> wilcoxons_versus_mock
cmr2535_smry$wilcoxon_p_versus_mock = wilcoxons_versus_mock$wilcoxon_p[match(cmr2535_smry$day, wilcoxons_versus_mock$day)]
cmr2535_smry$wilcoxon_p_versus_mock[cmr2535_smry$tx=='mock'] = 1.0
cmr2535 %>%
  group_by(tx) %>%
  summarize(.groups='keep', 
            wilcoxon_p = suppressWarnings(wilcox.test(nfl[day==10], nfl[day==24])$p.value)) %>%
  ungroup() -> wilcoxons_washout
cmr2535_smry$wilcoxon_p_washout = wilcoxons_washout$wilcoxon_p[match(cmr2535_smry$tx, wilcoxons_washout$tx)]
cmr2535_smry$wilcoxon_p_washout[cmr2535_smry$day != 24] = 1.0


write_supp_table(cmr2535_smry, 'Summary of NfL values in mouse treatment study 4')

write_stats('Difference between mock vs. minocycline NfL in mouse treatment study 4, day 5 timepoint: ',formatC(wilcoxon_p,format='g',digits=2))

cmr2535_weights = read_tsv('data/animals/CMR-2535-weights.tsv', col_types=cols()) 

cmr2535_weights %>%
  inner_join(main, by='animal') -> cmr2535_weights_with_meta

cmr2535_weights_with_meta %>%
  select(animal, sex, tx, bleed_day, weight_day=day, weight_g = weight, delta) -> cmr2535_weights_out

write_supp_table(cmr2535_weights_out, 'Individual weights in mouse treatment study 4')


cmr2535_weights_with_meta %>%
  group_by(tx, day) %>%
  summarize(.groups='keep',
            n = n(),
            mean_change = mean(delta),
            l95 = lower(delta),
            u95 = upper(delta)) %>%
  ungroup() -> cmr2535_weight_smry

write_supp_table(cmr2535_weight_smry, 'Summary of weights in mouse treatment study 4')

mock_mino_day5_weight_tobj = t.test(cmr2535_weights_with_meta$delta[cmr2535_weights_with_meta$day==5 & cmr2535_weights_with_meta$tx=='mino'],
                                    cmr2535_weights_with_meta$delta[cmr2535_weights_with_meta$day==5 & cmr2535_weights_with_meta$tx=='mock'])
  
write_stats('Mean weight change by day 5 in mouse treatment study 4: mock ',
            percent(mean(cmr2535_weights_with_meta$delta[cmr2535_weights_with_meta$day==5 & cmr2535_weights_with_meta$tx=='mock']), signed=T, digits=1),
            ' vs. mino ',
            percent(mean(cmr2535_weights_with_meta$delta[cmr2535_weights_with_meta$day==5 & cmr2535_weights_with_meta$tx=='mino']), signed=T, digits=1),
            ', P = ',formatC(mock_mino_day5_weight_tobj$p.value, format='e', digits=1))






tmt = read_tsv('data/tmt/tmt_1997_peptides.tsv', col_types=cols())
main %>%
  mutate(tmt_colname = paste0('x',gsub('\\.','_',animal))) %>%
  filter(tmt_colname %in% colnames(tmt)) %>%
  rename(grp=tx) %>%
  select(tmt_colname, grp) -> tmt_meta
allcano = makecano(tmt, tmt_meta, 'mino/mock', dataset_id='1997', verbose=F, saveimgs=F)
dataset_id = 1997
comparison = c('mino','mock')


# 
# 
# function(peptides, 
#          meta, 
#          contrast,
#          dataset_id=NULL, 
#          clipwidth=4, 
#          ymax=NA, 
#          verbose=T,
#          saveimgs=F) {



# QQ plot
# par(mar=c(3,3,1,1))
# n = nrow(allcano)
# observed = sort(allcano$p_ebm) # actual p values
# expected = seq(1/n,1,1/n) # under null, p values are uniformly distributed between 0 and 1
# o = -log10(observed) # log transform the p values
# e = -log10(expected) # log transform the expected p values
# plot(NA, NA, xlim=c(0,4), ylim=c(0,11),xlab="expected -log10(p)",ylab="observed -log10(p)", axes=F, ann=F, xaxs='i', yaxs='i')
# par(xpd=T)
# points(e,o,pch=19)
# par(xpd=F)
# abline(a=0,b=1,col="red") # add a line representing the null
# axis(side=1, at=0:20)
# axis(side=2, las=2, at=0:20)


allcano$mg = allcano$gene %in% mg_genes$feature_name[mg_genes$mg_enriched]
clipwidth = 3
ymax = 10
corrected_threshold = 0.05 / nrow(allcano)

allcano %>%
  select(gene, p_ebm, n_pep, abun, l2fc, mg) -> allcano_out
write_supp_table(allcano_out, 'TMT differential protein expression in mouse brain after 5 days minocycline in mouse treatment study 4.')

write_supp_table(mg_genes, 'Microglial enrichment data in mouse cortex from Mortberg & Gentile 2023')


# resx=300
# png(paste0('display_items/tmt_mino_cmr_2535.png'), width=3.5*resx, height=4.0*resx, res=resx)

par(mar=c(3,4,3,4))
xlims = clipwidth*c(-1.05,1.05)
ylims = c(0, ifelse(is.na(ymax), max(-log10(allcano$p_ebm))*1.05, ymax))
plot(x=NA, y=NA, axes=F, ann=F, xaxs='i', yaxs='i', xlim=xlims, ylim=ylims)
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=floor(min(xlims)):ceiling(max(xlims)), tck=-0.025, labels=NA)
axis(side=1, line=-0.75, at=floor(min(xlims)):ceiling(max(xlims)), lwd=0, cex.axis=0.8)
mtext(side=1, line=1.4, text='log2(fold change)')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
axis(side=2, at=floor(min(ylims)):ceiling(max(ylims)), tck=-0.025, labels=NA)
axis(side=2, at=seq(floor(min(ylims)),ceiling(max(ylims)),5), tck=-0.05, cex.axis=0.8, las=2)
mtext(side=2, line=2, text='-log10(P)')
abline(h=-log10(c(0.05, corrected_threshold)), lwd=0.25)
mtext(side=4, at=-log10(c(0.05, corrected_threshold)), line=0.25, las=2,
      text=c('nominal\nP < 0.05', 'Bonferroni\nP < 0.05'), cex=0.5)
par(xpd=T)
points(x=clipsym(allcano$l2fc,clipwidth), 
       y=-log10(allcano$p_ebm), 
       pch=20, 
       col=allcano$color)
callouts = rank(allcano$p_ebm) < 5 | (allcano$p_ebm < corrected_threshold & abs(allcano$l2fc) > .3) | (allcano$p_ebm < 0.05 & abs(allcano$l2fc) > 2)
allcano$callout_pos = ifelse(allcano$l2fc < 0 & allcano$l2fc > -clipwidth, 2, 4)
allcano$label_x_offset = 0
allcano$label_y_offset = 0
allcano$label_y_offset[allcano$gene=='Sult1a1'] = -0.35
allcano$label_x_offset[allcano$gene=='Sult1a1'] = -0.2
allcano$label_y_offset[allcano$gene=='Slc12a5'] = 0.3
allcano$label_x_offset[allcano$gene=='Slc12a5'] = -0.3
allcano$label_x_offset[allcano$gene=='Itih4'] = -0.1
allcano$label_y_offset[allcano$gene=='Itih4'] = 0.05
points(x=clipsym(allcano$l2fc[allcano$mg],clipwidth), 
       y=-log10(allcano$p_ebm[allcano$mg]), pch=1)
text(x=clipsym(allcano$l2fc[callouts],clipwidth)  + allcano$label_x_offset[callouts], 
     y=-log10(allcano$p_ebm[callouts]) + allcano$label_y_offset[callouts], 
     labels=allcano$gene[callouts], 
     pos=allcano$callout_pos[callouts], 
     font=3, cex=0.8, family='mono')
legend(x=-3, y=4, legend='microglial\ngenes', pch=1, pt.cex=1, cex=0.9, bty='n')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


# tmt_allcano %>% filter(p_ebm < 0.05) %>% arrange(desc(l2fc)) %>% View()
# 
# unnecessary_message_end_figure_4 = dev.off()
# 
# tmt_meta$order = 1:18
# tmt %>%
#   select(charge, gene_symbol, scan_number, peptide_sequence, x116322_1:x116329_2) %>%
#   pivot_longer(x116322_1:x116329_2) %>%
#   inner_join(tmt_meta, by =c('name'='tmt_colname')) %>%
#   mutate(x = ifelse(grp=='mock',1,2)) %>%
#   arrange(gene_symbol, scan_number, peptide_sequence, order) -> tmt_long
# 
# tmt_long %>% filter(gene_symbol=='Psd3') %>% View()
# tmt_long %>% filter(gene_symbol=='Psd3') %>% group_by(charge, gene_symbol, scan_number, peptide_sequence, grp) %>% summarize(.groups='keep', mean=mean(value))
# tmt_long %>% filter(gene_symbol=='Ero1a') %>% group_by(charge, gene_symbol, scan_number, peptide_sequence, grp) %>% summarize(.groups='keep', mean=mean(value))
# tmt_long %>% filter(gene_symbol=='Hmgn5') %>% group_by(charge, gene_symbol, scan_number, peptide_sequence, grp) %>% summarize(.groups='keep', mean=mean(value))



cmr2535 %>%
  filter(day %in% c(10,24)) %>%
  mutate(x = case_when(tx=='mock' & day==10 ~ 2.5,
                       tx=='mino' & day==10 ~ 1,
                       tx=='mock' & day==24 ~ 3,
                       tx=='mino' & day==24 ~ 1.5)) %>%
  mutate(pch = case_when(day==24 ~ 1,
                         day==10 ~ 20)) -> cmr2535_washout

cmr2535_washout %>%
  distinct(x, tx, day, pch, color) %>%
  mutate(axdisp = case_when(day==10 ~ 'on',
                            day==24 ~ 'off')) %>%
  mutate(legdisp = paste0(gsub('mino','minocycline',tx),ifelse(day==10,'',' washout'))) %>%
  arrange(x) -> xleg


# cmr2535_washout_smry = dviz(cmr2535_washout, xlims=range(xleg$x) + c(-0.5, 0.5), ylims=c(20,20000), log='y',
#                          xvar='jittered_x', yvar='nfl', colorvar='color', xcols=c('cohort','day','tx'),
#                          xbigs=NULL, yats=rep(1:9,4) * rep(10^(1:4),each=9), ybigs=10^(1:4),
#                          xlab='', ylab='NfL (pg/mL)', xaxcex=0.8, ylabcex=0.8, xlwds=0, yaxcex=0.8, 
#                          boxwidth=.17, mar=c(3,3,3,1))
par(mar=c(3,4,3,1))
xlims = c(0.5, 3.5)
ylims = c(20, 20000)
yats=rep(1:9,4) * rep(10^(1:4),each=9)
ybigs=10^(1:4)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, lwd=0, labels=c('10','100','1K','10K'), las=2)
set.seed(2)
cmr2535_washout$jittered_x = jitter(cmr2535_washout$x, amount=.1)
points(cmr2535_washout$jittered_x, cmr2535_washout$nfl, col=cmr2535_washout$color, pch=cmr2535_washout$pch)
cmr2535_washout %>%
  group_by(x, pch, color, tx, day) %>%
  summarize(.groups='keep', 
            n=n(), 
            mean_nfl = mean(nfl), nfl_l95=lower(nfl), nfl_u95=upper(nfl),
            median_nfl = median(nfl), nfl_q25=quantile(nfl,.25), nfl_q75=quantile(nfl,.75)) %>%
  ungroup() -> cmr2535_washout_smry
boxwidth = .17 
rect(xleft=cmr2535_washout_smry$x - boxwidth, xright=cmr2535_washout_smry$x + boxwidth,
     ybottom=cmr2535_washout_smry$nfl_q25, ytop=cmr2535_washout_smry$nfl_q75, 
     border=cmr2535_washout_smry$color, col=NA)
segments(x0=cmr2535_washout_smry$x - boxwidth, x1=cmr2535_washout_smry$x + boxwidth, y0=cmr2535_washout_smry$median_nfl, col=cmr2535_washout_smry$color)
for (anml in unique(cmr2535_washout$animal)) {
  subs = subset(cmr2535_washout, animal==anml)
  points(x=subs$jittered_x, y=subs$nfl, type='l', lwd=.125, col=subs$color)
}
mtext(side=1, line=0.25, at=xleg$x, text=xleg$day, col='#000000', cex=.7)
mtext(side=1, line=1, text='day')
mtext(side=2, line=2.6, text='plasma NfL (pg/mL)', cex=1)
legend(x=2.2, y=15000, xleg$legdisp, pch=xleg$pch, col=xleg$color, cex=.8, bty='n')

mock_p = suppressWarnings(wilcox.test(cmr2535_washout$nfl[cmr2535_washout$tx=='mock' & cmr2535_washout$day==10],
                                          cmr2535_washout$nfl[cmr2535_washout$tx=='mock' & cmr2535_washout$day==24])$p.value)
mino_p = suppressWarnings(wilcox.test(cmr2535_washout$nfl[cmr2535_washout$tx=='mino' & cmr2535_washout$day==10],
                                          cmr2535_washout$nfl[cmr2535_washout$tx=='mino' & cmr2535_washout$day==24])$p.value)

write_stats('Change in NfL among mock treated animals over 14 days "washout" in mouse treatment study 4: P = ',formatC(mock_p,format='g',digits=2),', Wilcoxon test.')
write_stats('Change in NfL among minocycline treated animals over 14 days washout in mouse treatment study 4: P = ',formatC(mino_p,format='g',digits=2),', Wilcoxon test.')

xleg$p = 1
xleg$p[xleg$tx=='mino' & xleg$day==24] = mino_p
xleg$p[xleg$tx=='mock' & xleg$day==24] = mock_p
mtext(side=3, at=xleg$x, line=0, text=p_to_symbol(xleg$p))
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



msc = read_tsv('data/animals/mouse_study_comparison.tsv', col_types=cols())

mp_smry_out %>%
  filter(grouped_day==6) %>%
  mutate(study=1) %>%
  group_by(study) %>%
  summarize(.groups='keep',
            mean_cv = mean(cv),
            fold_change = median[grouped_group=='minocycline']/median[grouped_group=='none']) %>%
  ungroup() -> mp_stats

dr_smry_out %>%
  filter(drug=='saline' | (drug == 'minocycline' & dose_mgkg==50)) %>%
  mutate(study=2) %>%
  group_by(study) %>%
  summarize(.groups='keep',
            mean_cv = mean(cv),
            fold_change = median[drug=='minocycline']/median[drug=='saline']) %>%
  ungroup() -> dr_stats

cmr_1982_smry_out %>%
  filter(day==5) %>%
  mutate(study=3) %>%
  group_by(study) %>%
  summarize(.groups='keep',
            mean_cv = mean(cv),
            fold_change = median[drug=='mino']/median[drug=='none']) %>%
  ungroup() -> cmr1982_stats

cmr2535_day5_smry %>%
  mutate(study=4) %>%
  group_by(study) %>%
  summarize(.groups='keep',
            mean_cv = mean(cv),
            fold_change = median[tx=='mino']/median[tx=='mock']) %>%
  ungroup() -> cmr2535_stats

msc %>%
  inner_join(rbind(mp_stats, dr_stats, cmr1982_stats, cmr2535_stats), by='study') -> msc_with_stats

write_supp_table(msc_with_stats, 'Comparison of experimental parameters and results across 4 mouse studies.')


unnecessary_message_end_figure_4 = dev.off()









######
# Figure S2
######

tell_user('done.\nCreating Figure S2...')

resx=300
png('display_items/figure-s2.png', width=6.5*resx, height=4*resx, res=resx)

layout_matrix = matrix(c(1,2),nrow=1, byrow=T)
layout(layout_matrix)
panel = 1

par(mar=c(3,4,3,1))
mino_color = meta$color[meta$cohort=='mino-5']
weights_rel = weights_rel_study1
weights_rel_avg = weights_rel_avg_study1
xlims = c(1,15)
xats = 1:15
xbigs = c(1,5,10,15)
ylims = c(-.25, .25)
yats = seq(-.25, .25, .05)
ybigs = seq(-.20, .20, .1)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, tck=-0.025, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbigs, cex.axis=0.8)
mtext(side=1, line=1.25, text='study day')
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.3, labels=percent(ybigs), cex.axis=0.8)
mtext(side=2, line=2.25, text='weight change')
abline(h=0, lty=1, col='black')
abline(h=-0.2, lty=3, col='red')
for (anml in unique(weights_rel$animal)) {
  subs = subset(weights_rel, animal==anml & !is.na(delta))
  points(subs$day, subs$delta, col=mino_color, type='l', lwd=0.125)
}        
nominal_threshold = min(deg$padj[deg$pvalue > 0.05], na.rm=T)
# plot(x=clipsym(deg$log2fold_change,4), y=-log10(deg$padj))
clipwidth = 4
par(mar=c(3,3,3,3))
xlims = clipwidth*c(-1.05,1.05)
ylims = c(0, 13.5)
plot(x=NA, y=NA, axes=F, ann=F, xaxs='i', yaxs='i', xlim=xlims, ylim=ylims)
points(x=clipsym(deg$l2fc,clipwidth), 
       y=clipdist(-log10(deg$padj),0,13), 
       pch=20, 
       col=deg$color)
points(x=clipsym(deg$l2fc[deg$mg_enriched & deg$padj < 0.05],clipwidth), 
       y=clipdist(-log10(deg$padj[deg$mg_enriched & deg$padj < 0.05]),0,13), 
       pch=1, 
       col='#000000',lwd=1.5,lty=3)
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=floor(min(xlims)):ceiling(max(xlims)), tck=-0.025, labels=NA)
axis(side=1, line=-1, at=floor(min(xlims)):ceiling(max(xlims)), lwd=0, cex.axis=0.7)
mtext(side=1, line=1, text='log2(fold change)')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
axis(side=2, at=floor(min(ylims)):ceiling(max(ylims)), tck=-0.025, labels=NA)
axis(side=2, at=seq(floor(min(ylims)),ceiling(max(ylims)),5), tck=-0.05, cex.axis=0.7, las=2)
axis(side=2, at=max(ylims), lwd=0, las=2, labels=paste0('≥',max(ylims)), cex.axis=0.7)
mtext(side=2, line=2, text='-log10(P)')
par(xpd=F)
abline(h=-log10(c(0.05, nominal_threshold)), lwd=0.25)
mtext(side=4, at=-log10(c(0.05)), line=0.25, las=2,
      text=c('FDR < 0.05'), cex=0.5)
par(xpd=T)
# select genes to label
callouts = rank(deg$padj) <= 5 | deg$padj < 1e-5 & deg$l2fc < 0 | deg$mg_enriched & deg$padj < 0.05
text(x=clipsym(deg$l2fc[callouts],clipwidth) + ifelse(deg$gene[callouts] %in% c('Ptk2b'), .5, 0), 
     y=clipdist(-log10(deg$padj[callouts]),0,13) + ifelse(deg$gene[callouts] %in% c('Aif1','Jpt1'),.25,0), 
     labels=deg$gene[callouts], 
     pos=ifelse(deg$l2fc[callouts] < 0 & deg$l2fc[callouts] > -clipwidth, 2, 4), 
     font=4, cex=0.6, family='mono')
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

unnecessary_message_end_figure_s2 = dev.off()









######
# Figure 5
######

tell_user('done.\nCreating Figure 5...')



resx=600
png('display_items/figure-5.png', width=6.5*resx, height=2.75*resx, res=resx)

layout_matrix = matrix(1:3,nrow=1, byrow=T)
layout(layout_matrix, widths=c(1.1, 1, 1))
panel = 1

meta = read_tsv('data/cultures/coc_meta.tsv', col_types=cols())
coc = read_tsv('data/cultures/coc.tsv', col_types=cols())

par(mar=c(3,4,2,1))
xlims = c(-0.5,10.5)
ylims = c(0, 5.5)
xats = seq(0, 10, 1)
xbigs = seq(0, 10, 2)
yats = seq(0, 5, 0.1)
ybigs = seq(0, 5, 1)
#yats = rep(1:9,3) * 10^(rep(0:2, each=9))
#ybigs = 10^(0:2)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=xats, tck=-0.025, lwd.ticks=1, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd.ticks=1, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbigs, cex.axis=0.8)
mtext(side=1, line=0.4, adj=-0.15, text='day:', cex=0.7)
axis(side=2, at=yats, tck=-0.025, lwd.ticks=1, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, lwd.ticks=1, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, lwd=0, line=-0.25, labels=ybigs, cex.axis=0.8, las=2)
abline(h=1, lty=3)
mtext(side=2, line=2.0, text='NfL fold change')
coc %>%
  inner_join(meta, by='treatment') %>%
  group_by(plate, treatment, id) %>%
  mutate(fold_change = nfl / nfl[day==0]) -> coc_anno
points(x=coc_anno$day, y=coc_anno$fold_change, col=coc_anno$color, pch=20)
coc_anno %>%
  rename(y=fold_change) %>%
  group_by(day, treatment, color) %>%
  summarize(.groups='keep',
            n = n(),
            mean = mean(y),
            l95 = mean(y) - 1.96 * sd(y) / sqrt(n()),
            u95 = mean(y) + 1.96 * sd(y) / sqrt(n()),
            median = median(y),
            q25 = quantile(y, .25),
            q75 = quantile(y, .75),
            cv = sd(y) / mean(y))  %>%
  ungroup() -> coc_smry
coc_smry$wilcoxon_p = as.numeric(NA)
for (i in which(coc_smry$treatment=='minocycline 25 µM')) {
  treated = subset(coc_anno, treatment=='minocycline 25 µM' & day==coc_smry$day[i])
  control = subset(coc_anno, treatment=='control' & day==coc_smry$day[i])
  if (length(unique(c(treated,control)))==1) {
    coc_smry$wilcoxon_p[i] = 1.00
  } else {
    coc_smry$wilcoxon_p[i] = suppressWarnings(wilcox.test(treated$fold_change, control$fold_change)$p.value)
  }
}
coc_smry$psymb = p_to_symbol(coc_smry$wilcoxon_p)
boxwidth = 0.5
rect(xleft=coc_smry$day-boxwidth, xright=coc_smry$day+boxwidth, ybottom=coc_smry$q25, ytop=coc_smry$q75, border=coc_smry$color, col=NA)
segments(x0=coc_smry$day-boxwidth, x1=coc_smry$day+boxwidth, y0=coc_smry$median, col=coc_smry$color)
for (tx in unique(coc_smry$treatment)) {
  subs = subset(coc_smry, treatment==tx)
  points(subs$day, subs$median, col=subs$color, type='l', lwd=0.5)
}
coc_smry %>%
  filter(!is.na(wilcoxon_p)) -> subs
par(xpd=T)
text(x=subs$day, y=rep(1.02*max(ylims), nrow(subs)), labels=subs$psymb)
par(xpd=F)

washout_color = '#FCE412'
mino_color = meta$color[meta$treatment=='minocycline 25 µM']
ybot = -1.0
ytop = -0.6
par(xpd=T)
rect(xleft=0, xright = 6, ybottom = ybot, ytop=ytop, col=alpha(mino_color, 0.3), lwd=1.5, border=NA)
text(x=(0+6)/2, y=mean(c(ybot,ytop)), labels='treatment', col=mino_color, cex=0.8, lwd=1.5)
rect(xleft=6, xright = 10, ybottom = ybot, ytop=ytop, col=alpha(washout_color, 0.3), lwd=1.5, border=NA)
text(x=(6+10)/2, y=mean(c(ybot,ytop)), labels='washout', col='black', cex=0.8, lwd=1.5)
par(xpd=F)

par(xpd=T)
legend('topleft', meta$treatment, col=meta$color, pch=20, bty='n', cex=0.8)
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = -0.2)
panel = panel + 1


coc_anno %>%
  select(-color) -> coc_anno_out
coc_smry %>%
  select(-color) -> coc_smry_out

write_supp_table(coc_anno_out, 'Individual NfL values in co-culture study 1')
write_supp_table(coc_smry_out, 'Summarized NfL values in co-culture study 1')

par(mar=c(0.5,0,2,0.5))
figure_5b = image_convert(image_read('data/cultures/RGB-MAX_20X_CTRL-3_Iba1_568_TUJ1_488_DAPI002-scalebar.tif'),'png')
plot(as.raster(figure_5b))
mtext(side=3, line = -0.7, text='control', cex=0.9)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.0, line = -0.2)
panel = panel + 1

figure_5c = image_convert(image_read('data/cultures/RGB-MAX_20X_Mino25-1_Iba1_568_TUJ1_488_DAPI001-scalebar.tif'),'png')
plot(as.raster(figure_5c))
mtext(side=3, line = -0.7, text='minocycline 25 µM', cex=0.9)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.0, line = -0.2)
panel = panel + 1

unnecessary_message_end_figure_5 = dev.off()







######
# Figure S3
######

tell_user('done.\nCreating Figure S3...')



resx=600
png('display_items/figure-s3.png', width=6.5*resx, height=6.5*resx, res=resx)

layout_matrix = matrix(1:4,nrow=2, byrow=T)
layout(layout_matrix)
panel = 1

par(mar=c(0.5,1,2,0.1))

panel_paths = tibble(path=c('RGB-MAX_20X_CTRL-3_Tuj1.tif',
                'RGB-MAX_20X_CTRL-3_Iba1.tiff',
                'RGB-MAX_20X_Mino25-1_Tuj1.tif',
                'RGB-MAX_20X_Mino25-1_Iba1.tiff'),
                condition = c('control','control','minocycline 25 µM','minocycline 25 µM'),
                channel = c('Tuj1','Iba1','Tuj1','Iba1'))

for (i in 1:nrow(panel_paths)) {
  figure_content = image_convert(image_read(paste0('data/cultures/',panel_paths$path[i])),'png')
  plot(as.raster(figure_content))
  if (i %% 2 == 1) {
    mtext(side=2, line = -0.1, text=panel_paths$condition[i], cex=1.2)
  }
  mtext(side=3, line = 0.1, text=panel_paths$channel[i], cex=1)
  mtext(LETTERS[panel], side=3, cex=2, adj = 0.1, line = 0)
  panel = panel + 1
  
}


unnecessary_message_end_figure_s3 = dev.off()





#####
# SUPPLEMENT
#####

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))



