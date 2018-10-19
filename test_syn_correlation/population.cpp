#include "population.hpp"

Population::Population(const Constant& nums){
  individuals.reserve(nums.N);
  for(std::size_t i=0; i < nums.N; i++){
    Individual ind;
    individuals.push_back(ind);
  }
}
void Population::set_params(const Constant& nums, const Parameters& param){
  for(std::size_t i=0; i < nums.N; i++){
    individuals[i].set_param(nums,param);
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
std::vector<std::size_t> mutater=individuals[i1].gamate_mutater(nums);
std::vector<std::size_t> mutater2=individuals[i2].gamate_mutater(nums);
for(std::size_t m=0; m < mutater.size(); m++){
  mutater[m] += mutater2[m];
}
/* tsg nonsynonymous */
  std::vector<std::size_t> tsg_non=individuals[i1].gamate_tsg_non(nums);
  std::vector<std::size_t> tn2=individuals[i2].gamate_tsg_non(nums);
  tsg_non.insert(tsg_non.end(), tn2.begin(), tn2.end());
/* tsg synonymous */
  std::vector<std::size_t> tsg_syn=individuals[i1].gamate_tsg_syn(nums);
  std::vector<std::size_t> ts2=individuals[i2].gamate_tsg_syn(nums);
  tsg_syn.insert(tsg_syn.end(), ts2.begin(), ts2.end());
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
    individuals[i].add_mutations(nums, param);
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
