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

tbl =read_tsv("t1/all_result.tsv")%>%bind_rows(read_tsv("t2/all_result.tsv"))
tn_num = 0.7642481
ts_num = 0.490959
regression_R = 0.03703
corr = 0.04642742

########mutater_freqが与えられていると、、、#########
#mutater の標準偏差
tbl%>>%
  mutate(p=mutater_freq/2, e_mutater_sd=(2*p*(1-p))**0.5)%>>%
  ggplot()+geom_point(aes(x=e_mutater_sd,y=mutater_sd,color=correlation))
#mutation_rateの平均
tbl %>>%
  mutate(p=mutater_freq/2)%>>%
  mutate(e_mutation_rate=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  ggplot()+
  geom_point(aes(x=mutation_rate_freq,y=e_mutation_rate,color=correlation))
#mutation_rateの標準偏差
tbl %>>%
  mutate(p=mutater_freq/2)%>>%
  mutate(e_mutr=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                                    (mutater_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                                    (mutater_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutation_rate_sd*(10**6),y=e_mutation_rate_sd*(10**6),color=correlation))
#予測されるmutation_rateの標準偏差と相関係数の関係
tbl %>>%
  mutate(p=mutater_freq/2)%>>%
  mutate(e_mutr=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                                    (mutater_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                                    (mutater_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+
  geom_point(aes(x=correlation,y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)+
  geom_vline(xintercept = corr*0.9)+geom_vline(xintercept = corr*1.1)

#--------------------------------------------------------------------------------------------------------------------------
#expect mutater_freqとobserveの違い
tbl %>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  ggplot()+
  geom_point(aes(x=mutater_freq,y=e_mutater_freq,color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#mutater の標準偏差
tbl%>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutater_freq, e_mutater_sd=(2*p*(1-p))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutater_sd,y=e_mutater_sd,color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#mutation_rateの平均
tbl %>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutater_freq)%>>%
  mutate(e_mutation_rate=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  ggplot()+
  geom_point(aes(x=mutation_rate_freq,y=e_mutation_rate,color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#mutation_rateの標準偏差
tbl %>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutater_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutater_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutater_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutation_rate_sd*(10**6),y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)
#予測されるmutation_rateの標準偏差と相関係数の関係
tbl %>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutater_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutater_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutater_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+
  geom_point(aes(x=correlation,y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)+
  geom_vline(xintercept = corr*0.9)+geom_vline(xintercept = corr*1.1)

######## correlationがobserveの10%以内になるe_mutation_rate_sdは？
tbl %>>%
  filter(correlation<corr*1.1,correlation>corr*0.9)%>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutater_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutater_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutater_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+geom_point(aes(x=mutation_rate_sd*(10**6),y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)
tbl %>>%
  filter(correlation<corr*1.1,correlation>corr*0.9)%>>%
  dplyr::rename(s=mutater_damage,m=mutater_mutation_rate) %>>%
  mutate(e_mutater_freq = m/s+0.5-m/2-((m/s+(1-m)/2)**2-m/s)**0.5) %>>%
  mutate(p=e_mutater_freq)%>>%
  mutate(e_mutr=((1-p)**2+mutater_effect*2*p*(1-p)+mutater_effect**2*(p**2))*mutation_rate)%>>%
  mutate(e_mutation_rate_sd=((mutation_rate-e_mutr)**2*(1-p)**2+
                               (mutater_effect*mutation_rate-e_mutr)**2*2*p*(1-p)+
                               (mutater_effect**2*mutation_rate-e_mutr)**2*(p**2))**0.5)%>>%
  ggplot()+
  geom_point(aes(x=correlation,y=e_mutation_rate_sd*(10**6),color=correlation))+
  geom_abline(intercept = 0,slope = 1)

