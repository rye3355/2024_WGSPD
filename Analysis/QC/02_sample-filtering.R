# Filtering samples based on gnomadv3 HQ sample filter

library(data.table)
library(stringr)
library(ggplot2)
# Read in sample-qc data
d <- fread("files/20240903_subset_sample_qc1.tsv")
hq_d <- fread("files/20240903_subset_hq_sample_qc1.tsv")
names(d) <- gsub("sample_qc1\\.", "", names(d))
names(hq_d) <- gsub("hq_sample_qc1\\.", "", names(hq_d))

# Read in manifest
manifest <- fread("../../subsetting/2024_WGSPD_merged-manifest.tsv")

# Annotate case/con
d <- merge(d, manifest[, c("s", "CASECON")], by = "s",
           all.x = T)
hq_d <- merge(hq_d, manifest[, c("s", "CASECON")], by = "s",
              all.x = T)

# Which were filtered
d$kept <- d$s %in% hq_d$s
d$size <- 0.5
d$size[!(d$s %in% hq_d$s)] <- 1





# Plotting function
create_pretty_boxplots <- function(df, aes, aes_col, threshold=NULL,
                                   threshold_max=NULL, file='file_out', title='', x_label='', y_label='',
                                   key_label='', xlim=NULL, legend=FALSE, save_figure=FALSE, 
                                   width=160, height=90, scaling=1, facet=FALSE, facet_grid=NULL, jitter_size=0.5,
                                   outlier.shape=NA, n_ticks=10, alpha=0.6, title.hjust=0.5, ggplot_theme=theme_classic, violin = T)
{
  p = ggplot(df, aes) +
    #geom_boxplot(outlier.shape=outlier.shape, coef=0, color='grey50', fill='grey95', show.legend=FALSE) + 
    coord_flip(ylim=xlim) +
    labs(title=title, x=y_label, y=x_label, color=key_label) + 
    scale_y_continuous(breaks = scales::pretty_breaks(n=n_ticks)) +
    guides(color = guide_legend(override.aes = list(size=2))) +
    ggplot_theme() +
    theme(axis.title.x = element_text(margin = ggplot2::margin(t=10)),
          plot.title = element_text(hjust=title.hjust))
  if (violin) {
    p <- p + 
      geom_violin(trim=T)
  }
  if (!is.null(threshold)) {
    p <- p + geom_hline(yintercept=threshold, linetype='dashed')
  }
  if (!is.null(threshold_max)) {
    p <- p + geom_hline(yintercept=threshold_max, linetype='dashed')
  }
  if (facet){
    p <- p + facet_grid
  }
  p <- p + geom_jitter(width=0.2, height=0, size=jitter_size, 
                       aes_col, show.legend=legend, alpha=alpha, stroke=0.05) + 
  
  if (save_figure) {
    ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
  }
  return(p)
}


# Stratify by Affected status
y_list <- c("call_rate", "dp_stats.mean", "gq_stats.mean", 
            "n_insertion", "n_deletion", "n_snp", 
            "r_het_hom_var",
            "n_singleton", "r_ti_tv", "r_insertion_deletion")

for (field in y_list) {
  plot <- create_pretty_boxplots(d, aes_string(x = "CASECON", y = field),
                                 aes(color = kept), violin = T, alpha = 0.3) +
    ggtitle(field)
  print(plot)
  ggsave(paste0('figures/02_qc1_', field, '.jpg'),
         plot, width=300, height=250, units='mm')
}


for (field in y_list) {
  plot <- create_pretty_boxplots(hq_d, aes_string(x = "CASECON", y = field), 
                                 aes(color = CASECON), violin = T, alpha = 0.3) +
    ggtitle(field)
  print(plot)
  ggsave(paste0('figures/02_hq_qc1_', field, '.jpg'), 
         plot, width=300, height=250, units='mm')
}
