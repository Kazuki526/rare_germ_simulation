#include "parameters.hpp"

Parameters::Parameters(Constant& nums){
  mutater_locas = get_mutater_locas(nums);
  mutater_effect = get_mutater_effect(nums)/mutater_locas;
  mutater_mutation_rate = get_mutater_mutation_rate(nums);
  mutater_damage = get_mutater_damage(nums);
  tsg_non_damage_e = get_tsg_non_damage_e(nums);

  std::exponential_distribution<> tsg_ex(tsg_non_damage_e);
  for(std::size_t s=0; nums.tsg_non_site > s; s++){
    tsg_non_damage.push_back(tsg_ex(nums.mt));
  }
}

  int Parameters::get_mutater_locas(Constant& nums){
    double dbl = log_random(0.1, 100.0, nums.mt);
    return(dbl<1 ? 1 : std::round(dbl));
  }
  double Parameters::get_mutater_effect(Constant& nums){
    std::uniform_real_distribution<> dist(2.0, 100.0);
    return dist(nums.mt);
  }
  double Parameters::get_mutater_mutation_rate(Constant& nums){
    return log_random(0.000001, 0.001, nums.mt);
  }
  double Parameters::get_mutater_damage(Constant& nums){
    std::uniform_real_distribution<> dist(0.0, 0.1);
    return dist(nums.mt);
  }
  double Parameters::get_tsg_non_damage_e(Constant& nums){
    std::uniform_real_distribution<> dist(0.1, 0.5);
    return dist(nums.mt);
  }
  void Parameters::set_tsg_non_damage(Constant& nums){
    std::exponential_distribution<> tsg_ex(tsg_non_damage_e);
    for(std::size_t s=0; nums.tsg_non_site > s; s++){
      tsg_non_damage.push_back(tsg_ex(nums.mt));
    }
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
