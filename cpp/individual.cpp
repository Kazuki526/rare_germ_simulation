#include"individual.hpp"

Individual::Individual(const int m, const std::vector<std::size_t>& tsg_non){
  mutater=m;
  if(!tsg_non.empty()){
    for(std::size_t mu: tsg_non){
      if(std::count(tsg_non.begin(), tsg_non.end(), mu) ==2){
        if(std::find(tsg_non_hom.begin(), tsg_non_hom.end(), mu) != tsg_non_hom.end()){
          tsg_non_hom.push_back(mu);
        }
      }else{
        tsg_non_het.push_back(mu);
      }
    }
  }
}

void Individual::set_param(const Constant& nums, const Parameters& param){
  mut_r = nums.mutation_rate * param.mutation_rate_coef;
  for(int i=1; i <= mutater; i++){mut_r*=param.mutater_effect;}
  damage=0;
  for(std::size_t& mu: tsg_non_het){damage+=param.tsg_non_damage[mu];}
  for(std::size_t& mu: tsg_non_hom){damage+=param.tsg_non_damage[mu]*2;}
  //fitness=1 - damage*param.fitness_coef;
}

void Individual::add_mutations(Constant& nums, const Parameters& param){
  /* new mutater mutation */
  std::bernoulli_distribution p_mutater(param.mutater_mutation_rate);
  int new_mutater = p_mutater(nums.mt);
  if(new_mutater>2){new_mutater=2;}
  mutater+=new_mutater;
  if(mutater>2){
    mutater-=2;
    if(mutater==2){mutater=0;}
  }
  /* new TSG nonsynonymous mutation */
  std::poisson_distribution<> pois_tn(nums.tsg_non_site*mut_r);
  std::uniform_int_distribution<> tn_mut(0, nums.tsg_non_site -1);
  int tn_num=pois_tn(nums.mt);
  std::vector<std::size_t> new_tsg_non;
  while(tn_num > 0){
    tn_num--;
    new_tsg_non.push_back(tn_mut(nums.mt));
  }
  if(!new_tsg_non.empty()){
    for(std::size_t &mut: new_tsg_non){
      /* mutation occur at homo site */
      bool focal=false;
      for(std::size_t mu=0; mu < tsg_non_hom.size(); mu++){
        if(tsg_non_hom[mu] == mut){focal=true;}
      }
      if(focal){
        tsg_non_hom.erase(tsg_non_hom.begin()+focal);
        tsg_non_het.push_back(mut);
        continue;
      }
      /* mutation occur at hetero site (focal=false) */
      for(std::size_t mu=0; mu < tsg_non_het.size(); mu++){
        if(tsg_non_het[mu] == mut){focal=true;}
      }
      if(focal){
        tsg_non_het.erase(tsg_non_het.begin()+focal);
        tsg_non_hom.push_back(mut);
        continue;
      }
      tsg_non_het.push_back(mut);
    }
  }
}

/* gamate methods */
int Individual::gamate_mutater(Constant nums){
  int new_mutater=0;
  std::bernoulli_distribution bern(0.5);
  if(mutater == 2){
    new_mutater++;
  }else if(mutater == 1){
    if(bern(nums.mt)){new_mutater++;}
  }
  return new_mutater;
}

std::vector<size_t> Individual::gamate_tsg_non(Constant nums){
  std::vector<std::size_t> new_tsg_non;
  std::bernoulli_distribution bern(0.5);
  for(std::size_t het_mu: tsg_non_het){
    if(bern(nums.mt)){new_tsg_non.push_back(het_mu);}
  }
  for(std::size_t hom_mu: tsg_non_hom){
    new_tsg_non.push_back(hom_mu);
  }
  return new_tsg_non;
}
