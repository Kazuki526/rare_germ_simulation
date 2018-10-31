#include "parameters.hpp"

Parameters::Parameters(Constant& nums){
  mutater_effect = get_mutater_effect(nums);
  mutater_mutation_rate = get_mutater_mutation_rate(nums);
  mutater_damage = get_mutater_damage(nums);
  tsg_non_damage_e = get_tsg_non_damage_e(nums);

  expected_mutater_freq =
    std::pow(mutater_effect,
      (mutater_mutation_rate/
        (mutater_damage+(2-mutater_damage)*mutater_mutation_rate)));
}
void Parameters::reset(Constant& nums){
  mutater_effect = get_mutater_effect(nums);
  mutater_mutation_rate = get_mutater_mutation_rate(nums);

  expected_mutater_freq =
    std::pow(mutater_effect,
      (mutater_mutation_rate/
        (mutater_damage+(2-mutater_damage)*mutater_mutation_rate)));
}

void Parameters::set_damage(Constant& nums){
  std::exponential_distribution<> tsg_ex(1.0/tsg_non_damage_e);
  for(std::size_t s=0; s < nums.tsg_non_site; s++){
    tsg_non_damage.push_back(tsg_ex(nums.mt));
  }
}

double Parameters::get_mutater_effect(Constant& nums){
  std::uniform_real_distribution<> dist(2.0, 100.0);
  return dist(nums.mt);
}
double Parameters::get_mutater_mutation_rate(Constant& nums){
  return log_random(0.0001, 0.01, nums.mt);
}
double Parameters::get_mutater_damage(Constant& nums){
  std::uniform_real_distribution<> dist(0.0, 0.01);
  return dist(nums.mt);
}
double Parameters::get_tsg_non_damage_e(Constant& nums){
  std::uniform_real_distribution<> dist(0.015, 0.05);
  return dist(nums.mt);
}


double log_random(double start, double end, std::mt19937& mt){
  start = std::log10(start==0 ? 0.000000001 : start);
  end   = std::log10(end);
  std::uniform_real_distribution<> dist(start, end);
  return std::pow(10, dist(mt));
}
bool equiv_lm(const std::vector<double>& mutation) {
  bool focal=true;
  std::size_t vect_size = mutation.size();
  double xy=0, xx=0 , average=0;
  for(double i: mutation){average+=i;}
  average =(double) average/vect_size;
  double x_ave = (double)(1+vect_size)/2;
  for(std::size_t t=0; t < vect_size; t++){
    xx += (x_ave-t)*(x_ave-t); xy += (x_ave-t)*(average-mutation[t]);
  }
  if((double)average/5000 > xy/xx){focal=false;}
  return focal;
}
