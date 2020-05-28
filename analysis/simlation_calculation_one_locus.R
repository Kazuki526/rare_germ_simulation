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


tbl = read_tsv("result_test/all_simulation_result.tsv")%>>%filter(correlation>0)
mut_tbl = read_tsv("result_test/all_simulation_result_mutator.tsv")
non_number = 2.8051689860835
syn_number = 1.54870775347913
ns_regression = 0.144
ns_correlation = 0.1988867

#tsg_non, synの数
tbl %>>%filter(correlation!=0)%>>%
  ggplot()+geom_point(aes(x=non_num,y=mutation_rate,color=mutation_rate_freq))+
  geom_vline(xintercept = non_number)
tbl %>>%filter(correlation!=0)%>>%
  ggplot()+geom_point(aes(y=syn_num,x=mutation_rate,color=correlation))+
  geom_hline(yintercept = syn_number*1.2)+geom_hline(yintercept = syn_number*0.8)
tbl %>>%filter(correlation!=0)%>>%
  mutate(mutator_s_rank=ifelse(mutation_rate*mutator_effect*mutator_damage>1,"out","ok"))%>>%
  ggplot()+geom_histogram(aes(x=correlation,fill=mutator_s_rank))+
  geom_vline(xintercept = ns_correlation*1.2)+geom_vline(xintercept = ns_correlation*0.8)

#### mutator_freqの予測 ####
tbl %>>%filter(correlation!=0)%>>%

  #mutate(smut=(mutation_rate*mutator_effect*mutator_damage))%>>%ggplot()+geom_histogram(aes(x=smut))
  mutate(e_mutator_freq=mutator_mutation_rate/(mutation_rate*mutator_effect*mutator_damage))%>>%
  filter(e_mutator_freq<0.1,(mutation_rate*mutator_effect*mutator_damage)<1)%>>%
  #filter(mutator_freq>e_mutator_freq)%>>%
  ggplot()+
  geom_point(aes(x=e_mutator_freq,y=mutator_freq,color=mutator_effect))+
  geom_abline(slope=1,intercept = 0,color="red")


########mutator_freqが与えられていると、、、#########
tbl %>>%filter(correlation>0)%>>%
  mutate(e_mutation_rate = mutation_rate*(1-mutator_freq+mutator_freq*mutator_effect))%>>%
  ggplot()+
  geom_point(aes(x=e_mutation_rate,y=syn_num))


#mutation_rateの平均
tbl %>>%filter(mutation_rate_freq<1e-4,mutator_freq<0.01)%>>%
  mutate(p=mutator_freq/2)%>>%
  mutate(e_mutation_rate=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  ggplot()+
  geom_point(aes(x=mutation_rate_freq,y=e_mutation_rate,color=mutator_freq))
#mutation_rateの標準偏差
tbl %>>%
  mutate(p=mutator_freq/2)%>>%
  mutate(e_mutr=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                                    (mutator_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                                    (mutator_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutation_rate_sd*(10**6),y=e_mutation_rate_sd*(10**6),color=correlation))
#予測されるmutation_rateの標準偏差と相関係数の関係
tbl %>>%filter(mutation_rate_sd<0.1,correlation>0)%>>%
  mutate(p=mutator_freq/2)%>>%
  mutate(e_mutr=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                                    (mutator_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                                    (mutator_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  filter(mutation_rate_sd>0.005,e_mutation_rate_sd<0.001)%>>%View
  ggplot()+
  geom_point(aes(x=mutation_rate_sd,y=e_mutation_rate_sd,color=correlation))+
  geom_abline(intercept = 0,slope = 1)

#--------------------------------------------------------------------------------------------------------------------------
#expect mutator_freqとobserveの違い
tbl %>>%filter(mutator_freq<0.1)%>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  ggplot()+
  geom_point(aes(x=mutator_freq,y=e_mutator_freq,color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#mutator の標準偏差
tbl%>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutator_freq, e_mutator_sd=(2*p*(1-p))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutator_sd,y=e_mutator_sd,color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#mutation_rateの平均
tbl %>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutator_freq)%>>%
  mutate(e_mutation_rate=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  ggplot()+
  geom_point(aes(x=mutation_rate_freq,y=e_mutation_rate,color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#mutation_rateの標準偏差
tbl %>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutator_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutator_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutator_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutation_rate_sd*(10**6),y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#予測されるmutation_rateの標準偏差と相関係数の関係
tbl %>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutator_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutator_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutator_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+
  geom_point(aes(x=correlation,y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)+
  geom_vline(xintercept = corr*0.9)+geom_vline(xintercept = corr*1.1)

######## correlationがobserveの10%以内になるe_mutation_rate_sdは？
tbl %>>%
  filter(correlation<corr*1.1,correlation>corr*0.9)%>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutator_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutator_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutator_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutation_rate_sd*(10**6),y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)
tbl %>>%
  filter(correlation<corr*1.1,correlation>corr*0.9)%>>%
  dplyr::rename(s=mutator_damage,m=mutator_mutation_rate) %>>%
  mutate(e_mutator_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutator_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutator_effect*2*p*(1-p)+mutator_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutator_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutator_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+
  geom_point(aes(x=correlation,y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)

