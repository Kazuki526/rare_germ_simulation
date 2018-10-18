#include"individual.hpp"

Individual::Individual(const int m,
           const std::vector<std::size_t>& tsg_non,
           const std::vector<std::size_t>& tsg_syn,
           const std::vector<std::size_t>& cont_non,
           const std::unordered_set<std::size_t>& tsg_non_common,
           const std::unordered_set<std::size_t>& tsg_syn_common,
           const std::unordered_set<std::size_t>& cont_non_common){
  mutater=m;
  bool common_focal =false;
  if(! tsg_non_common.empty() || ! tsg_syn_common.empty() || ! cont_non_common.empty()){common_focal=true;}
  if(common_focal){ /* need common variant check */
    for(std::size_t mu: tsg_non){
      if( tsg_non_common.find(mu) !=  tsg_non_common.end()){continue;}
      if(std::count(tsg_non.begin(), tsg_non.end(), mu) ==2){
        if(std::find(tsg_non_hom.begin(), tsg_non_hom.end(), mu) != tsg_non_hom.end()){
          tsg_non_hom.push_back(mu);
        }
      }else{
        tsg_non_het.push_back(mu);
      }
    }
    for(std::size_t mu: tsg_syn){
      if( tsg_syn_common.find(mu) !=  tsg_syn_common.end()){continue;}
      if(std::count(tsg_syn.begin(), tsg_syn.end(), mu) ==2){
        if(std::find(tsg_syn_hom.begin(), tsg_syn_hom.end(), mu) != tsg_syn_hom.end()){
          tsg_syn_hom.push_back(mu);
        }
      }else{
        tsg_syn_het.push_back(mu);
      }
    }
    for(std::size_t mu: cont_non){
      if( cont_non_common.find(mu) !=  cont_non_common.end()){continue;}
      if(std::count(cont_non.begin(), cont_non.end(), mu) ==2){
        if(std::find(cont_non_hom.begin(), cont_non_hom.end(), mu) != cont_non_hom.end()){
          cont_non_hom.push_back(mu);
        }
      }else{
        cont_non_het.push_back(mu);
      }
    }
  }else{ /* common_focal = false */
    for(std::size_t mu: tsg_non){
      if(std::count(tsg_non.begin(), tsg_non.end(), mu) ==2){
        if(std::find(tsg_non_hom.begin(), tsg_non_hom.end(), mu) != tsg_non_hom.end()){
          tsg_non_hom.push_back(mu);
        }
      }else{
        tsg_non_het.push_back(mu);
      }
    }
    for(std::size_t mu: tsg_syn){
      if(std::count(tsg_syn.begin(), tsg_syn.end(), mu) ==2){
        if(std::find(tsg_syn_hom.begin(), tsg_syn_hom.end(), mu) != tsg_syn_hom.end()){
          tsg_syn_hom.push_back(mu);
        }
      }else{
        tsg_syn_het.push_back(mu);
      }
    }
    for(std::size_t mu: cont_non){
      if(std::count(cont_non.begin(), cont_non.end(), mu) ==2){
        if(std::find(cont_non_hom.begin(), cont_non_hom.end(), mu) != cont_non_hom.end()){
          cont_non_hom.push_back(mu);
        }
      }else{
        cont_non_het.push_back(mu);
      }
    }
  }
}

void Individual::set_param(const Constant& nums, const Parameters& param){
  mut_r = nums.mutation_rate * param.mutation_rate_coef;
  for(int i=1; i <= mutater; i++){mut_r*=param.mutater_effect;}
  damage=0;
  for(std::size_t mu: tsg_non_het){damage+=param.tsg_non_damage[mu];}
  for(std::size_t mu: tsg_non_hom){damage+=param.tsg_non_damage[mu]*2;}
  double cont_damage=0;
  for(std::size_t mu: cont_non_het){cont_damage+=param.cont_non_damage[mu];}
  for(std::size_t mu: cont_non_hom){cont_damage+=param.cont_non_damage[mu]*2;}
  if(param.cont_non_fitness_only){
    fitness = 1 - (damage+cont_damage)*param.fitness_coef;
  }else{
    damage += cont_damage;
    fitness = 1 - damage*param.fitness_coef;
  }
  if(fitness < 0){fitness=0;}
}

void Individual::add_mutations(Constant& nums, const Parameters& param){
  /* new mutater mutation */
  std::bernoulli_distribution p_mutater(param.mutater_mutation_rate);
  int new_mutater = p_mutater(nums.mt);
  mutater+=new_mutater;
  if(mutater>2){
    mutater-=2;
  }
  /* new TSG nonsynonymous mutation */
  std::poisson_distribution<> pois_tn(nums.tsg_non_site*mut_r);
  std::poisson_distribution<> pois_ts(nums.tsg_syn_site*mut_r);
  std::poisson_distribution<> pois_cn(nums.cont_non_site*mut_r);
  std::uniform_int_distribution<> tn_mut(0, nums.tsg_non_site -1);
  std::uniform_int_distribution<> ts_mut(0, nums.tsg_syn_site -1);
  std::uniform_int_distribution<> cn_mut(0, nums.cont_non_site -1);
  int tn_num=pois_tn(nums.mt), ts_num=pois_ts(nums.mt), cn_num=pois_cn(nums.mt);
/* tsg nonsynonymous */
  if(tn_num > 0){
    std::vector<std::size_t> new_tsg_non;
    while(tn_num > 0){
      tn_num--;
      new_tsg_non.push_back(tn_mut(nums.mt));
    }
    for(std::size_t &mut: new_tsg_non){
      /* mutation occur at homo site */
      int focal=-1;
      for(std::size_t mu=0; mu < tsg_non_hom.size(); mu++){
        if(tsg_non_hom[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        tsg_non_hom.erase(tsg_non_hom.begin()+focal);
        tsg_non_het.push_back(mut);
        continue;
      }
      /* mutation occur at hetero site (focal=-1) */
      for(std::size_t mu=0; mu < tsg_non_het.size(); mu++){
        if(tsg_non_het[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        tsg_non_het.erase(tsg_non_het.begin()+focal);
        tsg_non_hom.push_back(mut);
        continue;
      }
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
      /* mutation occur at homo site */
      int focal=-1;
      for(std::size_t mu=0; mu < tsg_syn_hom.size(); mu++){
        if(tsg_syn_hom[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        tsg_syn_hom.erase(tsg_syn_hom.begin()+focal);
        tsg_syn_het.push_back(mut);
        continue;
      }
      /* mutation occur at hetero site (focal=-1) */
      for(std::size_t mu=0; mu < tsg_syn_het.size(); mu++){
        if(tsg_syn_het[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        tsg_syn_het.erase(tsg_syn_het.begin()+focal);
        tsg_syn_hom.push_back(mut);
        continue;
      }
      tsg_syn_het.push_back(mut);
    }
  }
/* control nonsynonymous */
  if(cn_num > 0){
    std::vector<std::size_t> new_cont_non;
    while(cn_num > 0){
      cn_num--;
      new_cont_non.push_back(cn_mut(nums.mt));
    }
    for(std::size_t &mut: new_cont_non){
      /* mutation occur at homo site */
      int focal=-1;
      for(std::size_t mu=0; mu < cont_non_hom.size(); mu++){
        if(cont_non_hom[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        cont_non_hom.erase(cont_non_hom.begin()+focal);
        cont_non_het.push_back(mut);
        continue;
      }
      /* mutation occur at hetero site (focal=-1) */
      for(std::size_t mu=0; mu < cont_non_het.size(); mu++){
        if(cont_non_het[mu] == mut){focal=mu;}
      }
      if(focal != -1){
        cont_non_het.erase(cont_non_het.begin()+focal);
        cont_non_hom.push_back(mut);
        continue;
      }
      cont_non_het.push_back(mut);
    }
  }
}

/* gamate methods */
int Individual::gamate_mutater(Constant& nums){
  int new_mutater=0;
  if(mutater == 2){
    new_mutater++;
  }else if(mutater == 1){
    if(nums.bern(nums.mt)){new_mutater++;}
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

std::vector<std::size_t> Individual::gamate_cont_non(Constant& nums){
  std::vector<std::size_t> new_cont_non;
  for(std::size_t het_mu: cont_non_het){
    if(nums.bern(nums.mt)){new_cont_non.push_back(het_mu);}
  }
  for(std::size_t hom_mu: cont_non_hom){
    new_cont_non.push_back(hom_mu);
  }
  return new_cont_non;
}
