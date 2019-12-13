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
Individual Population::reproduct(Constant& nums,const Parameters& param, std::size_t i1, std::size_t i2, bool common_check){
/* mutator */
  std::vector<std::size_t> mutator=individuals[i1].gamate_mutator(nums,param);
  std::vector<std::size_t> mut2=individuals[i2].gamate_mutator(nums,param);
  mutator.insert(mutator.end(), mut2.begin(), mut2.end());
/* tsg nonsynonymous */
  std::vector<std::size_t> tsg_non=individuals[i1].gamate_tsg_non(nums,param);
  std::vector<std::size_t> tn2=individuals[i2].gamate_tsg_non(nums,param);
  tsg_non.insert(tsg_non.end(), tn2.begin(), tn2.end());
/* tsg synonymous */
  std::vector<std::size_t> tsg_syn=individuals[i1].gamate_tsg_syn(nums,param);
  std::vector<std::size_t> ts2=individuals[i2].gamate_tsg_syn(nums,param);
  tsg_syn.insert(tsg_syn.end(), ts2.begin(), ts2.end());
  if(common_check){
    Individual ind(mutator, tsg_non, tsg_syn, tsg_non_common, tsg_syn_common);
    ind.set_param(nums, param);
    return(ind);
  }else{
    Individual ind(mutator, tsg_non, tsg_syn);
    ind.set_param(nums, param);
    return(ind);
  }
}

/* get all mutation AC */
void Population::mutation_count(const Constant& nums, const Parameters& param){
  const double rare =nums.N*2.0*0.05*0.01;
  std::vector<std::size_t> mutator_num;
  mutator_num.reserve(nums.N);
  num_tsg_non_mutation.clear();
  num_tsg_non_mutation.resize(nums.tsg_non_site);
  num_tsg_syn_mutation.clear();
  num_tsg_syn_mutation.resize(nums.tsg_syn_site);
  for(Individual &ind: individuals){
/* mutator */
    mutator_num.push_back(ind.get_mutator_num());
/* tsg nonsynonymous */
    for(std::size_t het_mu: ind.get_tsg_non_het()){
      num_tsg_non_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_tsg_non_hom()){
      num_tsg_non_mutation[hom_mu]+=2;
    }
/* tsg synonymous */
    for(std::size_t het_mu: ind.get_tsg_syn_het()){
      num_tsg_syn_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_tsg_syn_hom()){
      num_tsg_syn_mutation[hom_mu]+=2;
    }
  }
/* count each freq */
  mutator_freq = (double)std::accumulate(mutator_num.begin(),mutator_num.end(),0)/nums.N;
  rare_tsg_non_freq=0;
  for(std::size_t mu=0; mu < nums.tsg_non_site; mu++){
    if(num_tsg_non_mutation[mu] < rare){
      rare_tsg_non_freq += (double) num_tsg_non_mutation[mu] /nums.N;
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
     next_inds.push_back(reproduct(nums, param, dist(nums.mt), dist(nums.mt), common_check));
  }
  individuals = next_inds;
  /* mutation_count */
  mutation_count(nums, param);
  if(common_check){
    tsg_non_common.clear();
    tsg_syn_common.clear();
  }
}

void Population::correlation_ns(const Constant& nums){
  const int N =nums.N;
  const int patient_n =nums.patient_n;
  const double rare = (N-patient_n)*2*0.05*0.01;
/* select cancer patient */
/*0~patient_n-1 => cancer patients & patient_n~N => control case */

/* count mutation */
  std::vector<std::size_t> tsg_non_mutation(nums.tsg_non_site,0);
  std::vector<std::size_t> tsg_syn_mutation(nums.tsg_syn_site,0);
  for(std::size_t i =patient_n; i < N; i++){
    for(std::size_t het_mu: individuals[i].get_tsg_non_het()){
      tsg_non_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: individuals[i].get_tsg_non_hom()){
      tsg_non_mutation[hom_mu]+=2;
    }
    for(std::size_t het_mu: individuals[i].get_tsg_syn_het()){
      tsg_syn_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: individuals[i].get_tsg_syn_hom()){
      tsg_syn_mutation[hom_mu]+=2;
    }
  }
/* count patient rare num */
  std::vector<int> non_rare_num(patient_n,0), syn_rare_num(patient_n,0);
  std::vector<double> mutation_rate; mutation_rate.reserve(patient_n);
  for(std::size_t i=0; i <patient_n; i++){
    for(std::size_t tn_het: individuals[i].get_tsg_non_het()){
      if(tsg_non_mutation[tn_het] <= rare){non_rare_num[i]++;}
    }
    for(std::size_t tn_hom: individuals[i].get_tsg_non_hom()){
      if(tsg_non_mutation[tn_hom] <= rare){non_rare_num[i]+=2;}
    }
    for(std::size_t ts_het: individuals[i].get_tsg_syn_het()){
      if(tsg_syn_mutation[ts_het] <= rare){syn_rare_num[i]++;}
    }
    for(std::size_t ts_hom: individuals[i].get_tsg_syn_hom()){
      if(tsg_syn_mutation[ts_hom] <= rare){syn_rare_num[i]+=2;}
    }
    mutation_rate.push_back(individuals[i].get_mutation_r());
  }
  rare_tsg_non_freq = (double)std::accumulate(non_rare_num.begin(),non_rare_num.end(),0.0)/patient_n;
  rare_tsg_syn_freq = (double)std::accumulate(syn_rare_num.begin(),syn_rare_num.end(),0.0)/patient_n;
  mutation_rate_ave = (double)std::accumulate(mutation_rate.begin(),mutation_rate.end(),0.0)/patient_n;
  double nonv=0.0, synv=0.0, nonsynv=0.0, mutv=0.0, mutrv=0.0;
  double xy=0.0, xx=0.0, xx_=0.0, xy_=0.0;
  for(std::size_t i=0; i < patient_n; i++){
    nonv+=(double)((non_rare_num[i]-rare_tsg_non_freq)*(non_rare_num[i]-rare_tsg_non_freq));
    synv+=(double)((syn_rare_num[i]-rare_tsg_syn_freq)*(syn_rare_num[i]-rare_tsg_syn_freq));
    nonsynv+=(double)((non_rare_num[i]-rare_tsg_non_freq)*(syn_rare_num[i]-rare_tsg_syn_freq));
    mutrv+=(double)((mutation_rate[i]-mutation_rate_ave)*(mutation_rate[i]-mutation_rate_ave));
    xx+=(non_rare_num[i]-rare_tsg_non_freq)*(non_rare_num[i]-rare_tsg_non_freq);
    xy+=(non_rare_num[i]-rare_tsg_non_freq)*(syn_rare_num[i]-rare_tsg_syn_freq);
    xx_+=non_rare_num[i]*non_rare_num[i];
    xy_+=non_rare_num[i]*syn_rare_num[i];
  }
  rare_tsg_non_sd = std::sqrt((double)nonv/patient_n);
  rare_tsg_syn_sd = std::sqrt((double)synv/patient_n);
  rare_non_syn_correlation = (double)(nonsynv/patient_n)/(rare_tsg_non_sd*rare_tsg_syn_sd);
  rare_num_reg = (double) xy/xx;
  rare_num_reg_zero = (double) xy_/xx_;
  mutation_rate_sd = std::sqrt((double)mutrv/patient_n);
}

void Population::out_mutator_state(const Constant& nums, std::ofstream& mutout){
  std::vector<int> mutator_num(nums.new_mutator_id,0);
  double have_multi_mutator=0,have_multi_homo=0, have_single_mutator=0;
  for(Individual &ind: individuals){
    if(ind.get_mutator_num() >0){
      if(ind.get_mutator_num() ==1){have_single_mutator++;}else{have_multi_mutator++;}
      for(std::size_t het: ind.get_mutator_het()){mutator_num[het]++;}
      for(std::size_t hom: ind.get_mutator_hom()){mutator_num[hom]+=2;have_multi_homo++;}
    }
  }

  std::sort(mutator_num.begin(),mutator_num.end(),std::greater<int>());
  mutout << mutator_freq << "\t";
  mutout << have_multi_mutator << "\t";
  mutout << have_multi_homo << "\t"<< have_single_mutator;
  for(std::size_t i=0; i<10; i++){
    mutout << "\t" <<mutator_num[i];
  }
  mutout << "\n";
  mutout.flush();
}
