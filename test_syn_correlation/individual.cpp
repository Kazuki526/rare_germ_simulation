#include"individual.hpp"

Individual::Individual(const std::vector<std::size_t>& m,
           const std::vector<std::size_t>& tsg_non,
           const std::vector<std::size_t>& tsg_syn,
           const std::unordered_set<std::size_t>& tsg_non_common,
           const std::unordered_set<std::size_t>& tsg_syn_common){
  mutater=m;
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

void Individual::set_param(const Constant& nums, const Parameters& param){
  mutater.resize(param.mutater_locas);
  mut_r = nums.mutation_rate;
  for(const std::size_t mutater_num: mutater){
    for(int i=1; i <= mutater_num; i++){mut_r*=param.mutater_effect;}
  }
}

void Individual::add_mutations(Constant& nums, const Parameters& param){
  /* new mutater mutation */
  std::bernoulli_distribution p_mutater(param.mutater_mutation_rate);
  for(std::size_t mutater_posi=0; mutater_posi < param.mutater_locas; mutater_posi++){
    mutater[mutater_posi] +=p_mutater(nums.mt);
    if(mutater[mutater_posi]>2){mutater[mutater_posi]-=2;}
  }
  /* new TSG nonsynonymous mutation */
  std::poisson_distribution<> pois_tn(nums.tsg_non_site*mut_r);
  std::poisson_distribution<> pois_ts(nums.tsg_syn_site*mut_r);
  std::uniform_int_distribution<> tn_mut(0, nums.tsg_non_site -1);
  std::uniform_int_distribution<> ts_mut(0, nums.tsg_syn_site -1);
  int tn_num=pois_tn(nums.mt), ts_num=pois_ts(nums.mt);
/* tsg nonsynonymous */
  if(tn_num > 0){
    std::vector<std::size_t> new_tsg_non;
    while(tn_num > 0){
      tn_num--;
      new_tsg_non.push_back(tn_mut(nums.mt));
    }
    for(std::size_t &mut: new_tsg_non){
      tsg_non_het.push_back(mut);
    }
  }
/* tsg synonymous */
  if(ts_num > 0){
    std::vector<std::size_t> new_tsg_syn;
    while(ts_num > 0){
      ts_num--;
      new_tsg_syn.push_back(ts_mut(nums.mt));
    }
    for(std::size_t &mut: new_tsg_syn){
      tsg_syn_het.push_back(mut);
    }
  }
  /* set fitness */
  fitness=1;
  for(std::size_t mut: mutater){fitness-=param.mutater_damage*mut;}
  for(std::size_t mu: tsg_non_het){fitness-=param.tsg_non_damage[mu];}
  for(std::size_t mu: tsg_non_hom){fitness-=param.tsg_non_damage[mu]*2;}
  if(fitness < 0){fitness=0;}
}

/* gamate methods */
std::vector<std::size_t> Individual::gamate_mutater(Constant& nums){
  std::vector<std::size_t> new_mutater(mutater.size(),0);
  for(std::size_t mutater_posi=0; mutater_posi < mutater.size();mutater_posi++){
    if(mutater[mutater_posi] ==2){
      new_mutater[mutater_posi]=1;
    }else if(mutater[mutater_posi]==1){
      if(nums.bern(nums.mt)){new_mutater[mutater_posi]++;}
    }
  }
  return new_mutater;
}

std::vector<std::size_t> Individual::gamate_tsg_non(Constant& nums){
  std::vector<std::size_t> new_tsg_non;
  for(std::size_t het_mu: tsg_non_het){
    if(nums.bern(nums.mt)){new_tsg_non.push_back(het_mu);}
  }
  for(std::size_t hom_mu: tsg_non_hom){
    new_tsg_non.push_back(hom_mu);
  }
  return new_tsg_non;
}

std::vector<std::size_t> Individual::gamate_tsg_syn(Constant& nums){
  std::vector<std::size_t> new_tsg_syn;
  for(std::size_t het_mu: tsg_syn_het){
    if(nums.bern(nums.mt)){new_tsg_syn.push_back(het_mu);}
  }
  for(std::size_t hom_mu: tsg_syn_hom){
    new_tsg_syn.push_back(hom_mu);
  }
  return new_tsg_syn;
}
