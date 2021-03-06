#include "population.hpp"

Population::Population(Constant& nums, const Parameters& param){
  individuals.reserve(nums.N);
  for(std::size_t i=0; i < nums.N; i++){
    Individual ind;
    ind.set_param(nums,param);
    individuals.push_back(ind);
  }
}

void Population::set_common_variant(const Constant& nums){
  double common_num = nums.N *2*0.01;
  for(std::size_t m=0; m < nums.tsg_non_site; m++){
    if(num_tsg_non_mutation[m] > common_num){tsg_non_common.insert(m);}
  }
  for(std::size_t m=0; m < nums.tsg_syn_site; m++){
    if(num_tsg_syn_mutation[m] > common_num){tsg_syn_common.insert(m);}
  }
}

/* make offspring from two Individual */
Individual Population::reproduct(Constant& nums, const Parameters& param, std::size_t i1, std::size_t i2, bool common_check){
  std::vector<std::size_t> m1 = individuals[i1].gamate_mutater(nums,param);
  std::vector<std::size_t> m2 =  individuals[i2].gamate_mutater(nums,param);
  std::vector<std::size_t> mutater(param.mutater_locus,0);
  for(std::size_t i=0; i < param.mutater_locus; i++){
    mutater[i] = m1[i] + m2[i];
  }
  //std::cout <<"m";std::cout.flush();
/* tsg nonsynonymous */
  std::vector<std::size_t> tsg_non=individuals[i1].gamate_tsg_non(nums,param);
  std::vector<std::size_t> tn2=individuals[i2].gamate_tsg_non(nums,param);
  tsg_non.insert(tsg_non.end(), tn2.begin(), tn2.end());
  //std::cout <<"n";std::cout.flush();
/* tsg synonymous */
  std::vector<std::size_t> tsg_syn=individuals[i1].gamate_tsg_syn(nums,param);
  std::vector<std::size_t> ts2=individuals[i2].gamate_tsg_syn(nums,param);
  tsg_syn.insert(tsg_syn.end(), ts2.begin(), ts2.end());
  //std::cout <<"s";std::cout.flush();
  if(common_check){
    Individual ind(mutater, tsg_non, tsg_syn, tsg_non_common, tsg_syn_common);
    ind.set_param(nums, param);
    return(ind);
  }else{
    Individual ind(mutater, tsg_non, tsg_syn);
    ind.set_param(nums, param);
    return(ind);
  }
}

/* get all mutation AC */
void Population::mutation_count(const Constant& nums, const Parameters& param){
  const double rare =nums.N*2.0*0.05*0.01;
/* mutater */
  std::vector<std::size_t> mutater_num; mutater_num.reserve(nums.N);
  for(Individual &ind: individuals){
    mutater_num.push_back(ind.get_mutater());
  }
  mutater_freq = (double)std::accumulate(mutater_num.begin(),mutater_num.end(),0)/nums.N;
/* tsg nonsynonymous */
  num_tsg_non_mutation.clear();
  num_tsg_non_mutation.resize(nums.tsg_non_site);
  for(Individual &ind: individuals){
    for(std::size_t het_mu: ind.get_tsg_non_het()){
      num_tsg_non_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_tsg_non_hom()){
      num_tsg_non_mutation[hom_mu]+=2;
    }
  }
  rare_tsg_non_freq=0;
  for(std::size_t mu=0; mu < nums.tsg_non_site; mu++){
    if(num_tsg_non_mutation[mu] < rare){
      rare_tsg_non_freq += (double) num_tsg_non_mutation[mu] /nums.N;
    }
  }
/* tsg synonymous */
  num_tsg_syn_mutation.clear();
  num_tsg_syn_mutation.resize(nums.tsg_syn_site);
  for(Individual &ind: individuals){
    for(std::size_t het_mu: ind.get_tsg_syn_het()){
      num_tsg_syn_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_tsg_syn_hom()){
      num_tsg_syn_mutation[hom_mu]+=2;
    }
  }
  rare_tsg_syn_freq=0;
  for(std::size_t mu=0; mu < nums.tsg_syn_site; mu++){
    if(num_tsg_syn_mutation[mu] < rare){
      rare_tsg_syn_freq += (double) num_tsg_syn_mutation[mu] /nums.N;
    }
  }
}

void Population::next_generation(Constant& nums, const Parameters& param, bool common_check){
  if(common_check){set_common_variant(nums);}
  /* add_mutations & set fitness */
  fitness.clear(); fitness.reserve(nums.N);
  for(std::size_t i=0; i < nums.N; i++){
    fitness.push_back(individuals[i].get_fitness());
  }
  /* reproduct N individual */
  std::vector<Individual> next_inds;
  next_inds.reserve(nums.N);
  std::discrete_distribution<std::size_t> dist(fitness.begin(), fitness.end());
  for(std::size_t i=0; i < nums.N; i++){
    Individual ind = reproduct(nums, param, dist(nums.mt), dist(nums.mt), common_check);
    if(ind.mutater_mut_r >= 1){stop_or_not=true;break;}
    next_inds.push_back(ind);
  }
  individuals = next_inds;
  /* mutation_count */
  mutation_count(nums, param);
  if(common_check){
    tsg_non_common.clear();
    tsg_syn_common.clear();
  }
}

void Population::correlation_ns(){
  std::cout <<individuals.size() <<" "<<std::flush;
  const int N =individuals.size();
  const int rare = N/1000;
  std::vector<int> non_rare_num(N,0), syn_rare_num(N,0),non_notrare_num(N,0), syn_notrare_num(N,0);
  std::vector<std::size_t> mutater_num; mutater_num.reserve(N);
  std::vector<double> mutation_rate; mutation_rate.reserve(N);
  for(std::size_t i=0; i < N;i++){
    for(std::size_t tn_het: individuals[i].get_tsg_non_het()){
      if(num_tsg_non_mutation[tn_het] < rare){non_rare_num[i]++;}else{non_notrare_num[i]++;}
    }
    for(std::size_t tn_hom: individuals[i].get_tsg_non_hom()){
      if(num_tsg_non_mutation[tn_hom] < rare){non_rare_num[i]+=2;}else{non_notrare_num[i]+=2;}
    }
    for(std::size_t ts_het: individuals[i].get_tsg_syn_het()){
      if(num_tsg_syn_mutation[ts_het] < rare){syn_rare_num[i]++;}else{syn_notrare_num[i]++;}
    }
    for(std::size_t ts_hom: individuals[i].get_tsg_syn_hom()){
      if(num_tsg_syn_mutation[ts_hom] < rare){syn_rare_num[i]+=2;}else{syn_notrare_num[i]+=2;}
    }
    mutater_num.push_back(individuals[i].get_mutater());
    mutation_rate.push_back(individuals[i].get_mut_r());
  }
  double nonav,synav;
  nonav = (double)std::accumulate(non_rare_num.begin(),non_rare_num.end(),0.0)/N;
  synav = (double)std::accumulate(syn_rare_num.begin(),syn_rare_num.end(),0.0)/N;
  // mutater_freq = (double)std::accumulate(mutater_num.begin(),mutater_num.end(),0)/N;
  mutation_rate_ave = (double)std::accumulate(mutation_rate.begin(),mutation_rate.end(),0.0)/N;
  double nonv=0.0, synv=0.0, nonsynv=0.0, mutv=0.0, mutrv=0.0;
  double xy=0.0, xx=0.0, xx_=0.0, xy_=0.0;
  for(std::size_t i=0; i < N; i++){
    nonv+=(double)((non_rare_num[i]-nonav)*(non_rare_num[i]-nonav));
    synv+=(double)((syn_rare_num[i]-synav)*(syn_rare_num[i]-synav));
    nonsynv+=(double)((non_rare_num[i]-nonav)*(syn_rare_num[i]-synav));
    mutv+=(double)((mutater_num[i]-mutater_freq)*(mutater_num[i]-mutater_freq));
    mutrv+=(double)((mutation_rate[i]-mutation_rate_ave)*(mutation_rate[i]-mutation_rate_ave));
    xx+=(non_rare_num[i]-nonav)*(non_rare_num[i]-nonav);
    xy+=(non_rare_num[i]-nonav)*(syn_rare_num[i]-synav);
    xx_+=non_rare_num[i]*non_rare_num[i];
    xy_+=non_rare_num[i]*syn_rare_num[i];
  }
  rare_tsg_non_sd = std::sqrt((double)nonv/N);
  rare_tsg_syn_sd = std::sqrt((double)synv/N);
  rare_non_syn_correlation = (double)(nonsynv/N)/(rare_tsg_non_sd*rare_tsg_syn_sd);
  rare_num_reg = (double) xy/xx;
  rare_num_reg_zero = (double) xy_/xx_;
  mutater_sd = std::sqrt((double)mutv/N);
  mutation_rate_sd = std::sqrt((double)mutrv/N);
}
