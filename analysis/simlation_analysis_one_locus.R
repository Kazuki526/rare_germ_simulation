library(tidyverse)
library(pipeR)
library(gridExtra)
loadNamespace('cowplot')
setwd("~/simulate/")
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
if(0){
tbl = read_tsv("one_locus/t1/simulation_result1.tsv")
for(i in 2:100){
  tbl=bind_rows(tbl,
               read_tsv(paste0("one_locus/t1/simulation_result",i,".tsv"),
                col_types = paste0("i",paste0(rep("d",16),collapse = ""))))
}
write_df(tbl,"one_locus/t1/all_result.tsv")
}

tbl =read_tsv("one_locus/t1/all_result.tsv")
tn_num = 0.7642481
ts_num = 0.490959
regression_R = 0.03703
corr = 0.04642742


#tsg_non, synの数
tbl %>>%ggplot()+geom_point(aes(x=tsg_non_num,y=mutation_rate,color=mutation_rate_freq))+
  geom_vline(xintercept = 0.7642481)
tbl %>>%ggplot()+geom_point(aes(x=tsg_syn_num,y=mutation_rate,color=correlation))+
  geom_vline(xintercept = 0.490959)



#regの傾きを評価するのは原点通る場合よりそのままのregressionの方が良さそう
if(0){
ggplot(data=tbl)+geom_point(aes(x=reg_R,y=reg_R_zero))
ggplot(data=tbl)+geom_point(aes(x=reg_R,y=correlation))
ggplot(data=tbl)+geom_point(aes(x=reg_R_zero,y=correlation))
}
#regressionの傾きの分布
tbl %>>%ggplot()+geom_histogram(aes(x=reg_R))+geom_vline(xintercept = 0.03703)
#原点通る場合
#tbl %>>%ggplot()+geom_histogram(aes(x=reg_R_zero))+geom_vline(xintercept = 0.2866)

#mutation_rateの分散(標準偏差)との関係
tbl %>>%ggplot()+geom_point(aes(y=mutation_rate_sd, x=reg_R))+geom_vline(xintercept = 0.03703)
#tbl %>>%ggplot()+geom_point(aes(y=mutation_rate_sd, x=reg_R_zero))+geom_vline(xintercept = 0.2866)

#prior distribution
plot = tbl %>>%
  dplyr::select(mutation_rate,mutater_effect,mutater_mutation_rate,
                mutater_damage,tsg_non_damage_e) %>>%
  tidyr::gather() %>>%
  ggplot()+
  geom_histogram(aes(x=value))+
  theme_bw()+
  facet_wrap(~ key,scales = "free")
ggsave("~/Dropbox/work/simulate/one_locus/bef_data/prior_distribution.pdf",plot,width = 12,height = 8 )

#簡単にabc
conf_wid =0.05
plot_para = tbl %>>%
  filter(tsg_non_num >tn_num*(1-conf_wid),tsg_non_num<tn_num*(1+conf_wid)) %>>%
  filter(tsg_syn_num >ts_num*(1-conf_wid),tsg_syn_num<ts_num*(1+conf_wid)) %>>%
  filter(reg_R >regression_R*(1-conf_wid),reg_R<regression_R*(1+conf_wid)) %>>%(?.)%>>%
  #filter(correlation > corr*(1-conf_wid),correlation<corr*(1+conf_wid)) %>>%
  dplyr::select(mutation_rate,mutater_effect,mutater_mutation_rate,
                mutater_damage,tsg_non_damage_e) %>>%
  tidyr::gather() %>>%
  ggplot()+
  geom_histogram(aes(x=value))+
  theme_bw()+
  facet_wrap(~ key,scales = "free")

plot_mfreq = tbl %>>%
  filter(tsg_non_num >tn_num*(1-conf_wid),tsg_non_num<tn_num*(1+conf_wid)) %>>%
  filter(tsg_syn_num >ts_num*(1-conf_wid),tsg_syn_num<ts_num*(1+conf_wid)) %>>%
  filter(reg_R >regression_R*(1-conf_wid),reg_R<regression_R*(1+conf_wid)) %>>%
  #filter(correlation > corr*(1-conf_wid),correlation<corr*(1+conf_wid)) %>>%
  dplyr::select(mutater_freq) %>>%
  tidyr::gather() %>>%
  ggplot()+
  geom_histogram(aes(x=value))+
  theme_bw()+
  facet_wrap(~ key,scales = "free")
plot = cowplot::plot_grid(plot_para,plot_mfreq,rel_widths = c(3,1),labels = c("a","b"))
if(conf_wid==0.05){
  ggsave("~/Dropbox/work/simulate/one_locus/bef_data/conf005_abc_mutaterfreq.pdf")
}else if(conf_wid==0.1){
  ggsave("~/Dropbox/work/simulate/one_locus/bef_data/conf01_abc_mutaterfreq.pdf")
}else if(conf_wid==0.2){
  ggsave("~/Dropbox/work/simulate/one_locus/bef_data/conf02_abc_mutaterfreq.pdf")
}
