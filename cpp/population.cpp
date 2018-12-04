#include "population.hpp"

Population::Population(const Constant& nums, const Parameters& param){
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
  for(std::size_t m=0; m < nums.cont_non_site; m++){
    if(num_cont_non_mutation[m] > common_num){cont_non_common.insert(m);}
  }
}

/* make offspring from two Individual */
Individual Population::reproduct(Constant& nums, const Parameters& param, std::size_t i1, std::size_t i2, bool common_check){
  int mutater =
      individuals[i1].gamate_mutater(nums,param) +
      individuals[i2].gamate_mutater(nums,param);
/* tsg nonsynonymous */
  std::vector<std::size_t> tsg_non=individuals[i1].gamate_tsg_non(nums,param);
  std::vector<std::size_t> tn2=individuals[i2].gamate_tsg_non(nums,param);
  tsg_non.insert(tsg_non.end(), tn2.begin(), tn2.end());
/* tsg synonymous */
  std::vector<std::size_t> tsg_syn=individuals[i1].gamate_tsg_syn(nums,param);
  std::vector<std::size_t> ts2=individuals[i2].gamate_tsg_syn(nums,param);
  tsg_syn.insert(tsg_syn.end(), ts2.begin(), ts2.end());
/* control nonsynonymous */
  std::vector<std::size_t> cont_non=individuals[i1].gamate_cont_non(nums,param);
  std::vector<std::size_t> cn2=individuals[i2].gamate_cont_non(nums,param);
  cont_non.insert(cont_non.end(), cn2.begin(), cn2.end());
  if(common_check){
    Individual ind(mutater, tsg_non, tsg_syn, cont_non, tsg_non_common, tsg_syn_common, cont_non_common);
    ind.set_param(nums, param);
    return(ind);
  }else{
    Individual ind(mutater, tsg_non, tsg_syn, cont_non);
    ind.set_param(nums, param);
    return(ind);
  }
}

/* get all mutation AC */
void Population::mutation_count(const Constant& nums, const Parameters& param){
  const double rare =nums.N*2.0*0.05*0.01;
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
/* control nonsynonymous */
  num_cont_non_mutation.clear();
  num_cont_non_mutation.resize(nums.cont_non_site);
  for(Individual &ind: individuals){
    for(std::size_t het_mu: ind.get_cont_non_het()){
      num_cont_non_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_cont_non_hom()){
      num_cont_non_mutation[hom_mu]+=2;
    }
  }
  rare_cont_non_freq=0;
  for(std::size_t mu=0; mu < nums.cont_non_site; mu++){
    if(num_cont_non_mutation[mu] < rare){
      rare_cont_non_freq += (double) num_cont_non_mutation[mu] /nums.N;
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
    cont_non_common.clear();
  }
}

void Population::regression_onset_age(Constant& nums, const Parameters& param){
  const double rare =nums.N*2.0*0.05*0.01;
  std::normal_distribution<> norm_dist(nums.mean_onset_age, nums.onset_age_sd);
/* set patient number(id) for analyzing onset age */
  std::vector<std::size_t> patient_ids(nums.N);
  std::iota(patient_ids.begin(),patient_ids.end(),0);
  std::shuffle(patient_ids.begin(),patient_ids.end(),nums.mt);

/* search rare variant nums of each patient */
  int rare_tsg_non_num=0, rare_tsg_syn_num=0, rare_cont_non_num=0;
  std::vector<int> rare_tsg_non(nums.patient_n, 0);
  std::vector<int> rare_tsg_syn(nums.patient_n, 0);
  std::vector<int> rare_cont_non(nums.patient_n, 0);
  std::vector<double> onset_age_list; onset_age_list.reserve(nums.patient_n);
  for(std::size_t t=0; t < nums.patient_n; t++){
    Individual& ind = individuals[patient_ids[t]];
    /* tsg nonsynonymous rare variant */
    for(std::size_t het: ind.get_tsg_non_het()){
      if(num_tsg_non_mutation[het] < rare){rare_tsg_non[t]++;rare_tsg_non_num++;}
    }
    for(std::size_t hom: ind.get_tsg_non_hom()){
      if(num_tsg_non_mutation[hom] < rare){rare_tsg_non[t]+=2;rare_tsg_non_num+=2;}
    }
    /* tsg synonymous rare variant */
    for(std::size_t het: ind.get_tsg_syn_het()){
      if(num_tsg_syn_mutation[het] < rare){rare_tsg_syn[t]++;rare_tsg_syn_num++;}
    }
    for(std::size_t hom: ind.get_tsg_syn_hom()){
      if(num_tsg_syn_mutation[hom] < rare){rare_tsg_syn[t]+=2;rare_tsg_syn_num+=2;}
    }
    /* control nonsynonymous rare variant */
    for(std::size_t het: ind.get_cont_non_het()){
      if(num_cont_non_mutation[het] < rare){rare_cont_non[t]++;rare_cont_non_num++;}
    }
    for(std::size_t hom: ind.get_cont_non_hom()){
      if(num_cont_non_mutation[hom] < rare){rare_cont_non[t]+=2;rare_cont_non_num+=2;}
    }
    onset_age_list.push_back(norm_dist(nums.mt) - ind.get_damage());
  }
  rare_tsg_non_freq = (double)rare_tsg_non_num/nums.patient_n;
  rare_tsg_syn_freq = (double)rare_tsg_syn_num/nums.patient_n;
  rare_cont_non_freq = (double)rare_cont_non_num/nums.patient_n;
  tsg_non_regression = linear_model(rare_tsg_non, onset_age_list);
  tsg_syn_regression = linear_model(rare_tsg_syn, onset_age_list);
  cont_non_regression = linear_model(rare_cont_non, onset_age_list);
}

bool around_percentn(const double value, const double desired_value, double n=0.2){
  bool result = true;
  if(value < desired_value*(1.0-n)){result=false;}
  if(value > desired_value*(1.0+n)){result=false;}
  return result;
}

bool accept_reject_judge(const Constant& nums, const Population& population){
  bool result=true;
  if(!around_percentn(population.rare_tsg_non_freq, nums.rare_tsg_non_num)){
    result = false;
  }else if(!around_percentn(population.rare_tsg_syn_freq, nums.rare_tsg_syn_num)){
    result = false;
  }else if(!around_percentn(population.rare_cont_non_freq, nums.rare_cont_non_num)){
    result = false;
  }else if(!around_percentn(population.tsg_non_regression, nums.tsg_non_onset_regression)){
    result = false;
  }else if(!around_percentn(population.tsg_syn_regression, nums.tsg_syn_onset_regression)){
    result = false;
  }else if(!around_percentn(population.cont_non_regression, nums.cont_non_onset_regression)){
    result = false;
  }
  return result;
}
