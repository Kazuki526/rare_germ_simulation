library(tidyverse)
library(pipeR)
library(gridExtra)
loadNamespace('cowplot')
setwd("~/simulate/main")
write_df= function(x, path, delim='\t', na='NA', append=FALSE, col_names=!append, ...) {
  file = if (grepl('gz$', path)) {
    gzfile(path, ...)
  } else if (grepl('bz2$', path)) {
    bzfile(path, ...)
  } else if (grepl('xz$', path)) {
    xzfile(path, ...)
  } else {path}
  utils::write.table(x, file,
                     append=append, quote=FALSE, sep=delim, na=na,
                     row.names=FALSE, col.names=col_names)
}


#tbl = read_tsv("result_test/all_simulation_result.tsv")
tbl = read_tsv("all_simulation_result.tsv")
mut_tbl = read_tsv("all_simulation_result_mutator.tsv")
non_number = 2.8051689860835
syn_number = 1.54870775347913
ns_regression = 0.144
ns_correlation = 0.1988867



#prior distribution
plot = tbl %>>%
  dplyr::select(mutation_rate,mutator_effect,mutator_mutation_rate,
                mutator_damage,non_damage_e) %>>%
  dplyr::rename(`mutation rate (μ)`=mutation_rate,`mutator effect (t)`=mutator_effect,
                `mutator mutation rate (μm)`=mutator_mutation_rate,
                `selection to mutatin rate (smut)`=mutator_damage,
                `selection to nonsynonymous (s)`=non_damage_e)%>>%
  filter(`mutator mutation rate (μm)`<0.002)%>>%
  tidyr::gather() %>>%
  ggplot()+
  geom_histogram(aes(x=value))+
  theme_bw()+
  facet_wrap(~ key,scales = "free")
plot
ggsave("~/Dropbox/work/simulate/one_locus/bef_data/prior_distribution.pdf",plot,width = 12,height = 8 )

#簡単にabc
conf_wid =0.1
plot_para = tbl %>>%
  filter(non_num+syn_num >(non_number+syn_number)*(1-conf_wid),non_num+syn_num<(non_number+syn_number)*(1+conf_wid)) %>>%
  filter(non_num/syn_num >non_number/syn_number*(1-conf_wid),non_num/syn_num<non_number/syn_number*(1+conf_wid)) %>>%(?.)%>>%
  filter(reg_R >ns_regression*(1-conf_wid),reg_R<ns_regression*(1+conf_wid)) %>>%(?.)%>>%
  #filter(correlation > ns_correlation*(1-conf_wid),correlation<ns_correlation*(1+conf_wid)) %>>%(?.)%>>%
  dplyr::select(mutation_rate,mutator_effect,mutator_mutation_rate,mutator_damage,non_damage_e) %>>%
  mutate(mutator_effect_class=ifelse(mutation_rate>2e-8,"high","low"))%>>%
  #mutate(mutator_effect_class=ifelse(mutator_effect*mutator_damage*mutation_rate-mutation_rate>1,"out","ok"))%>>%(?.%>>%count(mutator_effect_class))%>>%
  dplyr::rename(`mutation rate (μ)`=mutation_rate,`mutator effect (t)`=mutator_effect,
                `mutator mutation rate (μm)`=mutator_mutation_rate,
                `selection to mutatin rate (smut)`=mutator_damage,
                `selection to nonsynonymous (s)`=non_damage_e)%>>%
  #filter(`mutator mutation rate (μm)`<0.002)%>>%
  tidyr::gather(key="key",value="value",-mutator_effect_class) %>>%
  ggplot()+
  #geom_histogram(aes(x=value,fill=mutator_effect_class))+
  geom_histogram(aes(x=value))+xlab("Value")+ylab("Count")+
  theme_classic()+
  facet_wrap(~ key,scales = "free")+
  theme(strip.background = element_rect(color="white"),strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),axis.text = element_text(size=12,color="black"))
plot_para
ggsave("~/Dropbox/work/simulate/posterior_distribution.pdf",plot_para,width = 14,height = 6)
plot_mfreq = tbl %>>%
  filter(non_num+syn_num >(non_number+syn_number)*(1-conf_wid),non_num+syn_num<(non_number+syn_number)*(1+conf_wid)) %>>%
  filter(non_num/syn_num >non_number/syn_number*(1-conf_wid),non_num/syn_num<non_number/syn_number*(1+conf_wid)) %>>%
  #filter(reg_R >ns_regression*(1-conf_wid),reg_R<ns_regression*(1+conf_wid)) %>>%
  filter(correlation > ns_correlation*(1-conf_wid),correlation < ns_correlation*(1+conf_wid)) %>>%
  mutate(mutator_effect_class=ifelse(mutation_rate>3e-8,"high","low"))%>>%
  dplyr::select(mutator_freq,mutator_effect_class) %>>%
  tidyr::gather(key="key",value="value",-mutator_effect_class) %>>%
  ggplot()+
  #geom_histogram(aes(x=value,fill=mutator_effect_class),binwidth = 0.001)+
  geom_histogram(aes(x=value),binwidth = 0.001)+
  theme_classic()+
  facet_wrap(~ key,scales = "free")+
  theme(strip.background = element_rect(color="white"),strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),axis.text = element_text(size=12,color="black"))
plot = cowplot::plot_grid(plot_para+theme(legend.position = "none"),
                          plot_mfreq,rel_widths = c(3,1),labels = c("a","b"))
plot
 if(conf_wid==0.05){
  ggsave("~/Dropbox/work/simulate/one_locus/bef_data/conf005_abc_mutatorfreq.pdf",plot)
}else if(conf_wid==0.1){
  ggsave("~/Dropbox/work/simulate/one_locus/conf01_abc_mutatorfreq.pdf",plot)
}else if(conf_wid==0.2){
  ggsave("~/Dropbox/work/simulate/one_locus/bef_data/conf02_abc_mutatorfreq.pdf",plot)
}


tbl %>>%
  filter(non_num+syn_num >(non_number+syn_number)*(1-conf_wid),non_num+syn_num<(non_number+syn_number)*(1+conf_wid)) %>>%
  filter(non_num/syn_num >non_number/syn_number*(1-conf_wid),non_num/syn_num<non_number/syn_number*(1+conf_wid)) %>>%
  filter(reg_R >ns_regression*(1-conf_wid),reg_R<ns_regression*(1+conf_wid)) %>>%
  filter(correlation > ns_correlation*(1-conf_wid),correlation<ns_correlation*(1+conf_wid)) %>>%View
  mutate(cor_rank=ifelse(correlation > ns_correlation*(1-conf_wid)&correlation<ns_correlation*(1+conf_wid),"ok","out"))%>>%
  #filter(mutator_freq<0.02)%>>%
  ggplot()+geom_point(aes(x=mutator_freq*mutator_effect,y=correlation))

do_summary=mut_tbl%>>%group_by(slot,replicate)%>>%summarise(summary="ok")
tbl %>>%
  filter(non_num+syn_num >(non_number+syn_number)*(1-conf_wid),non_num+syn_num<(non_number+syn_number)*(1+conf_wid)) %>>%
  filter(non_num/syn_num >non_number/syn_number*(1-conf_wid),non_num/syn_num<non_number/syn_number*(1+conf_wid)) %>>%
  #filter(reg_R >ns_regression*(1-conf_wid),reg_R<ns_regression*(1+conf_wid)) %>>%
  filter(correlation > ns_correlation*(1-conf_wid),correlation<ns_correlation*(1+conf_wid)) %>>%
  left_join(do_summary)%>>%View

tbl %>>%(?.%>>%count()) %>>% filter(correlation!=0) %>>% (?.%>>%count()) %>>%
  mutate(mutator_selection=ifelse(mutator_effect*mutator_damage*mutation_rate-mutation_rate>1,"out","ok"))%>>%
  count(mutator_selection)
