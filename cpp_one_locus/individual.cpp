#include"individual.hpp"

Individual::Individual(const std::vector<std::size_t>& m,
           const std::vector<std::size_t>& tsg_non,
           const std::vector<std::size_t>& tsg_syn,
           const std::unordered_set<std::size_t>& tsg_non_common,
           const std::unordered_set<std::size_t>& tsg_syn_common){
  mutater=m;
  mutater_num=m.size();
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
  for(std::size_t i=1; i <= mutater_num; i++){fitness*=(1-param.mutater_damage);}
  for(std::size_t mu: tsg_non_het){fitness*=(1-param.tsg_non_damage[mu]);}
  for(std::size_t mu: tsg_non_hom){fitness*=((1-param.tsg_non_damage[mu])*(1-param.tsg_non_damage[mu]));}
  //for(std::size_t m=0; m < tsg_non_het.size();m++){fitness*=(1-param.tsg_non_damage_e);}
  //for(std::size_t m=0; m < tsg_non_hom.size();m++){fitness*=((1-param.tsg_non_damage_e)*(1-param.tsg_non_damage_e));}
  /* mutation rate */
  mut_r = param.mutation_rate;
  /*#######################################################################*/
  for(int i=1; i <= mutater_num; i++){mut_r*=param.mutater_effect;}
  /*#######################################################################*/
}

/* gamate methods */
std::vector<std::size_t> Individual::gamate_mutater(Constant& nums, Parameters& param){
  std::vector<std::size_t> new_mutater;
  if(mutater_num > 1){
    for(std::size_t mut_allele: mutater){

    }
  }else{
    for(std::size_t i=1; i<=mutater_num;i++){
      if(nums.bern(nums.mt)){new_mutater.push_back(mutater[i]);}
    }
  }
  /* add mutation */
 std::bernoulli_distribution p_mutater(param.mutater_mutation_rate);
 if(p_mutater(nums.mt)){
   new_mutater.push_back(param.get_new_mutater());
 }
  return new_mutater;
}

std::vector<std::size_t> Individual::gamate_tsg_non(Constant& nums, const Parameters& param){
  std::vector<std::size_t> new_tsg_non =tsg_non_hom;
  for(std::size_t het_mu: tsg_non_het){
    if(nums.bern(nums.mt)){new_tsg_non.push_back(het_mu);}
  }
  /* add mutation */
  std::poisson_distribution<> pois_tn(nums.tsg_non_site*mut_r);
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
  std::poisson_distribution<> pois_ts(nums.tsg_syn_site*mut_r);
  std::uniform_int_distribution<> ts_mut(0, nums.tsg_syn_site -1);
  int ts_num=pois_ts(nums.mt);
  for(std::size_t n=1; n <= ts_num; n++){
      new_tsg_syn.push_back(ts_mut(nums.mt));
  }
  return new_tsg_syn;
}
