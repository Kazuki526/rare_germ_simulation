#include "parameters.hpp"

Constant::Constant(){
  std::random_device rnd;
  std::mt19937 mt_(rnd());
  mt = mt_;
}

void Parameters::set_damage(Constant& nums){
  mutation_rate_coef = get_mutation_rate_coef(nums);
  mutater_effect = get_mutater_effect(nums);
  mutater_mutation_rate = get_mutater_mutation_rate(nums);
  mutater_damage = get_mutater_damage(nums);
  tsg_non_damage_e = get_tsg_non_damage_e(nums);
  cont_non_damage_e = get_cont_non_damage_e(nums);
  fitness_coef = get_fitness_coef(nums);

  std::exponential_distribution<> tsg_ex(tsg_non_damage_e);
  for(std::size_t s=0; nums.tsg_non_site > s; s++){
    tsg_non_damage.push_back(tsg_ex(nums.mt));
  }
  std::exponential_distribution<> cont_ex(cont_non_damage_e);
  for(std::size_t s=0; nums.cont_non_site > s; s++){
    cont_non_damage.push_back(cont_ex(nums.mt));
  }
}

void Parameters::change_param(Constant& nums,const std::size_t time){
  switch (variable_param) {
    case 0: mutation_rate_coef=variable_param_values[time]; break;
    case 1: mutater_effect=variable_param_values[time]; break;
    case 2: mutater_mutation_rate=variable_param_values[time]; break;
    case 3: mutater_damage=variable_param_values[time]; break;
    case 4:
      for(std::size_t m=0; m < nums.tsg_non_site; m++){
        tsg_non_damage[m] *= variable_param_values[time]/tsg_non_damage_e;
      }
      tsg_non_damage_e = variable_param_values[time];
      break;
    case 5:
      for(std::size_t m=0; m < nums.cont_non_site; m++){
        cont_non_damage[m] *= variable_param_values[time]/cont_non_damage_e;
      }
      cont_non_damage_e = variable_param_values[time];
      break;
    case 6:
      if(variable_param_values[time] == 0){
        cont_non_fitness_only = false;
      }else{
        cont_non_fitness_only = true;
      }
      break;
    case 7: fitness_coef=variable_param_values[time]; break;
  }
}

double Parameters::get_mutation_rate_coef(Constant& nums){
  std::uniform_real_distribution<> dist(0.1, 1.0);
  return dist(nums.mt);
}
double Parameters::get_mutater_effect(Constant& nums){
  std::uniform_real_distribution<> dist(5.0, 100.0);
  return dist(nums.mt);
}
double Parameters::get_mutater_mutation_rate(Constant& nums){
  return log_random(0, 0.001, nums.mt);
}
double Parameters::get_mutater_damage(Constant& nums){
  std::uniform_real_distribution<> dist(0.0, 5.0);
  return dist(nums.mt);
}
double Parameters::get_tsg_non_damage_e(Constant& nums){
  std::uniform_real_distribution<> dist(0.1, 5.0);
  return dist(nums.mt);
}
double Parameters::get_cont_non_damage_e(Constant& nums){
  std::uniform_real_distribution<> dist(0.1, 2.0);
  return dist(nums.mt);
}
double Parameters::get_cont_non_fitness_only(Constant& nums){
  std::uniform_int_distribution<> dist(0, 1);
  return dist(nums.mt);
}
double Parameters::get_fitness_coef(Constant& nums){
  return log_random(0.01,0.5,nums.mt);
}


void Parameters::set_tsg_non_damage(Constant& nums){
  std::exponential_distribution<> tsg_ex(tsg_non_damage_e);
  for(std::size_t s=0; nums.tsg_non_site > s; s++){
    tsg_non_damage.push_back(tsg_ex(nums.mt));
  }
}
void Parameters::set_cont_non_damage(Constant& nums){
  std::exponential_distribution<> cont_ex(cont_non_damage_e);
  for(std::size_t s=0; nums.cont_non_site > s; s++){
    cont_non_damage.push_back(cont_ex(nums.mt));
  }
}


double log_random(double start, double end, std::mt19937 mt){
  start = std::log10(start==0 ? 0.000000001 : start);
  end   = std::log10(end);
  std::uniform_real_distribution<> dist(start, end);
  return std::pow(10, dist(mt));
}