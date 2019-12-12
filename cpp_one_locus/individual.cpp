#include"individual.hpp"

Individual::Individual(const std::vector<std::size_t>& mutator,
           const std::vector<std::size_t>& tsg_non,
           const std::vector<std::size_t>& tsg_syn,
           const std::unordered_set<std::size_t>& tsg_non_common,
           const std::unordered_set<std::size_t>& tsg_syn_common){
  mutator_num=mutator.size();
  for(std::size_t mu: mutator){
    if(std::count(mutator.begin(), mutator.end(), mu) >2){
      if(std::find(mutator_hom.begin(), mutator_hom.end(), mu) != mutater_hom.end()){
        mutator_hom.push_back(mu);
      }
    }else{
      mutator_het.push_back(mu);
    }
  }
  bool common_focal =false;
  if(! tsg_non_common.empty() || ! tsg_syn_common.empty()){common_focal=true;}
  if(common_focal){ /* need common variant check */
    for(std::size_t mu: tsg_non){
      if( tsg_non_common.find(mu) !=  tsg_non_common.end()){continue;}
      if(std::count(tsg_non.begin(), tsg_non.end(), mu) >2){
        if(std::find(tsg_non_hom.begin(), tsg_non_hom.end(), mu) != tsg_non_hom.end()){
          tsg_non_hom.push_back(mu);
        }
      }else{
        tsg_non_het.push_back(mu);
      }
    }
    for(std::size_t mu: tsg_syn){
      if( tsg_syn_common.find(mu) !=  tsg_syn_common.end()){continue;}
      if(std::count(tsg_syn.begin(), tsg_syn.end(), mu) >2){
        if(std::find(tsg_syn_hom.begin(), tsg_syn_hom.end(), mu) != tsg_syn_hom.end()){
          tsg_syn_hom.push_back(mu);
        }
      }else{
        tsg_syn_het.push_back(mu);
      }
    }
  }else{ /* common_focal = false */
    for(std::size_t mu: tsg_non){
      if(std::count(tsg_non.begin(), tsg_non.end(), mu) >2){
        if(std::find(tsg_non_hom.begin(), tsg_non_hom.end(), mu) != tsg_non_hom.end()){
          tsg_non_hom.push_back(mu);
        }
      }else{
        tsg_non_het.push_back(mu);
      }
    }
    for(std::size_t mu: tsg_syn){
      if(std::count(tsg_syn.begin(), tsg_syn.end(), mu) >2){
        if(std::find(tsg_syn_hom.begin(), tsg_syn_hom.end(), mu) != tsg_syn_hom.end()){
          tsg_syn_hom.push_back(mu);
        }
      }else{
        tsg_syn_het.push_back(mu);
      }
    }
  }
}

void Individual::set_param(Constant& nums, const Parameters& param){
  /* set fitness */
  fitness=1;
  for(std::size_t i=1; i <= mutator_num; i++){fitness*=(1-param.mutator_damage);}
  for(std::size_t mu: tsg_non_het){fitness*=(1-param.tsg_non_damage[mu]);}
  for(std::size_t mu: tsg_non_hom){fitness*=((1-param.tsg_non_damage[mu])*(1-param.tsg_non_damage[mu]));}
  //for(std::size_t m=0; m < tsg_non_het.size();m++){fitness*=(1-param.tsg_non_damage_e);}
  //for(std::size_t m=0; m < tsg_non_hom.size();m++){fitness*=((1-param.tsg_non_damage_e)*(1-param.tsg_non_damage_e));}
  /* mutation rate */
  mutation_r = param.mutation_rate;
  /*#######################################################################*/
  /* all mutator affect synergistic */
  for(std::size_t i=1; i <= mutator_num; i++){mutation_r*=param.mutator_effect;}
  /* mutator effect, dif allele additive and hom synergistic */
  //if(mutator_num>0){mutation_r *= param.mutator_effect*mutator_het.size() + param.mutator_effect*param.mutator_effect*mutator_hom.size();}
  /*#######################################################################*/
}

/* gamate methods */
std::vector<std::size_t> Individual::gamate_mutator(Constant& nums, Parameters& param){
  std::vector<std::size_t> new_mutator =mutator_hom;
  for(std::size_t mu: mutator_het){
    if(nums.bern(nums.mt)){new_mutator.push_back(mu);}
  }
  /* add mutation */
 std::bernoulli_distribution p_mutator(param.mutator_mutation_rate);
 if(p_mutator(nums.mt)){
   new_mutator.push_back(param.get_new_mutator());
 }
  return new_mutator;
}

std::vector<std::size_t> Individual::gamate_tsg_non(Constant& nums, const Parameters& param){
  std::vector<std::size_t> new_tsg_non =tsg_non_hom;
  for(std::size_t het_mu: tsg_non_het){
    if(nums.bern(nums.mt)){new_tsg_non.push_back(het_mu);}
  }
  /* add mutation */
  std::poisson_distribution<> pois_tn(nums.tsg_non_site*mutation_r);
  std::uniform_int_distribution<> tn_mut(0, nums.tsg_non_site -1);
  int tn_num=pois_tn(nums.mt);
  for(std::size_t n=1; n <=tn_num; n++){
      new_tsg_non.push_back(tn_mut(nums.mt));
  }
  return new_tsg_non;
}

std::vector<std::size_t> Individual::gamate_tsg_syn(Constant& nums, const Parameters& param){
  std::vector<std::size_t> new_tsg_syn=tsg_syn_hom;
  for(std::size_t het_mu: tsg_syn_het){
    if(nums.bern(nums.mt)){new_tsg_syn.push_back(het_mu);}
  }
  /* add mutation */
  std::poisson_distribution<> pois_ts(nums.tsg_syn_site*mutation_r);
  std::uniform_int_distribution<> ts_mut(0, nums.tsg_syn_site -1);
  int ts_num=pois_ts(nums.mt);
  for(std::size_t n=1; n <= ts_num; n++){
      new_tsg_syn.push_back(ts_mut(nums.mt));
  }
  return new_tsg_syn;
}
