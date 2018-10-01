#include"individual.hpp"

Individual::Individual(const int& m, std::vector<int>& tn){
  mutater=m;
  if(!tn.empty()){
    for(int &mu: tn){
      if(std::count(tn.begin(), tn.end(), mu) ==2){
        if(std::find(tsg_non_hom.begin(), tsg_non_hom.end(), mu) != tsg_non_hom.end()){
          tsg_non_hom.push_back(mu);
        }
      }else{
        tsg_non_het.push_back(mu);
      }
    }
  }
}

void Individual::set_param(Parameters& param){
  mut_r = param.mutation_rate * param.mutation_rate_coef;
  for(int i=1; i <= mutater; i++){mut_r*=param.mutater_effect;}
  damage=0;
  for(int& mu: tsg_non_het){damage+=param.tsg_non_damage[mu];}
  for(int& mu: tsg_non_hom){damage+=param.tsg_non_damage[mu]*2;}
  fitness=1 - damage*param.fitness_coef;
}

void Individual::add_mutations(Parameters& param){
  // new mutater mutation
  std::bernoulli_distribution p_mutater(param.mutater_mutation_rate);
  int new_mutater = p_mutater(param.mt);
  if(new_mutater>2){new_mutater=2;}
  mutater+=new_mutater;
  if(mutater>2){
    mutater-=2;
    if(mutater==2){mutater=0;}
  }
  // new TSG nonsynonymous mutation
  std::poisson_distribution<> pois_tn(param.tsg_non_site*mut_r);
  std::uniform_int_distribution<> tn_mut(0, param.tsg_non_site -1);
  int tn_num=pois_tn(param.mt);
  std::vector<int> new_tsg_non;
  while(tn_num > 0){
    tn_num--;
    new_tsg_non.push_back(tn_mut(param.mt));
  }
  if(!new_tsg_non.empty()){
    for(int &mut: new_tsg_non){
      // mutation occur at homo site
      bool focal=false;
      for(int mu=0; mu < tsg_non_hom.size(); mu++){
        if(tsg_non_hom[mu] == mut){focal=true;}
      }
      if(focal){
        tsg_non_hom.erase(tsg_non_hom.begin()+focal);
        tsg_non_het.push_back(mut);
        continue;
      }
      // mutation occur at hetero site (focal=false)
      for(int mu=0; mu < tsg_non_het.size(); mu++){
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
