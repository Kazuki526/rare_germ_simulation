#include"individual.hpp"

Individual::Individual(const std::vector<std::size_t>& mutator,
           const std::vector<std::size_t>& non,
           const std::vector<std::size_t>& syn,
           const std::unordered_set<std::size_t>& non_common,
           const std::unordered_set<std::size_t>& syn_common){
  mutator_num=mutator.size();
  for(std::size_t mu: mutator){
    if(std::count(mutator.begin(), mutator.end(), mu) >2){
      if(std::find(mutator_hom.begin(), mutator_hom.end(), mu) != mutator_hom.end()){
        mutator_hom.push_back(mu);
      }
    }else{
      mutator_het.push_back(mu);
    }
  }
  bool common_focal =false;
  if(! non_common.empty() || ! syn_common.empty()){common_focal=true;}
  if(common_focal){ /* need common variant check */
    for(std::size_t mu: non){
      if( non_common.find(mu) !=  non_common.end()){continue;}
      if(std::count(non.begin(), non.end(), mu) >2){
        if(std::find(non_hom.begin(), non_hom.end(), mu) != non_hom.end()){
          non_hom.push_back(mu);
        }
      }else{
        non_het.push_back(mu);
      }
    }
    for(std::size_t mu: syn){
      if( syn_common.find(mu) !=  syn_common.end()){continue;}
      if(std::count(syn.begin(), syn.end(), mu) >2){
        if(std::find(syn_hom.begin(), syn_hom.end(), mu) != syn_hom.end()){
          syn_hom.push_back(mu);
        }
      }else{
        syn_het.push_back(mu);
      }
    }
  }else{ /* common_focal = false */
    for(std::size_t mu: non){
      if(std::count(non.begin(), non.end(), mu) >2){
        if(std::find(non_hom.begin(), non_hom.end(), mu) != non_hom.end()){
          non_hom.push_back(mu);
        }
      }else{
        non_het.push_back(mu);
      }
    }
    for(std::size_t mu: syn){
      if(std::count(syn.begin(), syn.end(), mu) >2){
        if(std::find(syn_hom.begin(), syn_hom.end(), mu) != syn_hom.end()){
          syn_hom.push_back(mu);
        }
      }else{
        syn_het.push_back(mu);
      }
    }
  }
}

void Individual::set_param(Constant& nums, const Parameters& param){
  /* mutation rate */
  mutation_r = param.mutation_rate;
  /* all mutator affect synergistic */
  for(std::size_t i=1; i <= mutator_num; i++){mutation_r*=param.mutator_effect;}
  /* mutator effect, dif allele additive and hom synergistic */
  //if(mutator_num>0){mutation_r *= param.mutator_effect*mutator_het.size() + param.mutator_effect*param.mutator_effect*mutator_hom.size();}
  /* if mutator is mutation of repair pathway */
  //if(mutator_num >0){mutation_r *= param.repair_power*(1- std::pow(1-t, mutator_num));}
  /* set fitness */
  fitness = 1;
  fitness*=(1 - param.mutator_damage*(mutation_r - param.mutation_rate));
  if(fitness<0){fitness=0;}
  for(std::size_t mu: non_het){fitness*=(1-param.non_damage[mu]);}
  for(std::size_t mu: non_hom){fitness*=((1-param.non_damage[mu])*(1-param.non_damage[mu]));}
}

/* gamate methods */
std::vector<std::size_t> Individual::gamate_mutator(Constant& nums, const Parameters& param){
  std::vector<std::size_t> new_mutator =mutator_hom;
  for(std::size_t mu: mutator_het){
    if(nums.bern(nums.mt)){new_mutator.push_back(mu);}
  }
  /* add mutation */
 std::bernoulli_distribution p_mutator(param.mutator_mutation_rate);
 if(p_mutator(nums.mt)){
   new_mutator.push_back(nums.get_new_mutator());
 }
  return new_mutator;
}

std::vector<std::size_t> Individual::gamate_non(Constant& nums, const Parameters& param){
  std::vector<std::size_t> new_non =non_hom;
  for(std::size_t het_mu: non_het){
    if(nums.bern(nums.mt)){new_non.push_back(het_mu);}
  }
  /* add mutation */
  std::poisson_distribution<> pois_non((nums.non_site/3)*mutation_r);
  std::uniform_int_distribution<> non_mut(0, nums.non_site -1);
  int non_num=pois_non(nums.mt);
  for(std::size_t n=1; n <=non_num; n++){
      new_non.push_back(non_mut(nums.mt));
  }
  return new_non;
}

std::vector<std::size_t> Individual::gamate_syn(Constant& nums, const Parameters& param){
  std::vector<std::size_t> new_syn=syn_hom;
  for(std::size_t het_mu: syn_het){
    if(nums.bern(nums.mt)){new_syn.push_back(het_mu);}
  }
  /* add mutation */
  std::poisson_distribution<> pois_syn((nums.syn_site/3)*mutation_r);
  std::uniform_int_distribution<> syn_mut(0, nums.syn_site -1);
  int syn_num=pois_syn(nums.mt);
  for(std::size_t n=1; n <= syn_num; n++){
      new_syn.push_back(syn_mut(nums.mt));
  }
  return new_syn;
}
