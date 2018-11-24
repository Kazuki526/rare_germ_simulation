#include "parameters.hpp"

Parameters::Parameters(Constant& nums){
  mutation_rate = get_mutation_rate(nums);
  mutater_effect = get_mutater_effect(nums);
  mutater_mutation_rate = get_mutater_mutation_rate(nums);
  mutater_damage = get_mutater_damage(nums);
  tsg_non_damage_e = get_tsg_non_damage_e(nums);
  cont_non_damage_e = get_cont_non_damage_e(nums);

  double m=mutater_mutation_rate, s=mutater_damage, e=mutater_effect;
  double p = m/s+0.5-m/2-std::pow((m/s+(1-m)/2)*(m/s+(1-m)/2)-m/s,0.5);
  double mutr = ((1-p)*(1-p)+mutater_effect*2*p*(1-p)+e*e*(p*p))*mutation_rate;
  expected_mutation_sd =
    std::pow((mutation_rate-mutr)*(mutation_rate-mutr)*(1-p)*(1-p)+
             (e*mutation_rate-mutr)*(e*mutation_rate-mutr)*2*p*(1-p)+
             (e*e*mutation_rate-mutr)*(e*e*mutation_rate-mutr)*p*p,0.5);
}
void Parameters::reset(Constant& nums){
  mutater_effect = get_mutater_effect(nums);
  mutater_mutation_rate = get_mutater_mutation_rate(nums);

  double m=mutater_mutation_rate, s=mutater_damage, e=mutater_effect;
  double p = m/s+0.5-m/2-std::pow((m/s+(1-m)/2)*(m/s+(1-m)/2)-m/s,0.5);
  double mutr = ((1-p)*(1-p)+mutater_effect*2*p*(1-p)+e*e*(p*p))*mutation_rate;
  expected_mutation_sd =
    std::pow((mutation_rate-mutr)*(mutation_rate-mutr)*(1-p)*(1-p)+
             (e*mutation_rate-mutr)*(e*mutation_rate-mutr)*2*p*(1-p)+
             (e*e*mutation_rate-mutr)*(e*e*mutation_rate-mutr)*p*p,0.5);
}

void Parameters::set_damage(Constant& nums){
  std::exponential_distribution<> tsg_ex(1.0/tsg_non_damage_e);
  std::exponential_distribution<> cont_ex(1.0/cont_non_damage_e);
  for(std::size_t s=0; s < nums.tsg_non_site; s++){
    tsg_non_damage.push_back(tsg_ex(nums.mt));
  }
  for(std::size_t s=0; s < nums.cont_non_site; s++){
    cont_non_damage.push_back(cont_ex(nums.mt));
  }
}

double Parameters::get_mutation_rate(Constant& nums){
  std::uniform_real_distribution<> dist(0.00000002,0.00000004);
  return dist(nums.mt);
}
double Parameters::get_mutater_effect(Constant& nums){
  std::uniform_real_distribution<> dist(70.0, 120.0);
  return dist(nums.mt);
}
double Parameters::get_mutater_mutation_rate(Constant& nums){
  return log_random(0.0003, 0.0006, nums.mt);
}
double Parameters::get_mutater_damage(Constant& nums){
  std::uniform_real_distribution<> dist(0.001, 0.008);
  return dist(nums.mt);
}
double Parameters::get_tsg_non_damage_e(Constant& nums){
  std::uniform_real_distribution<> dist(0.015, 0.03);
  return dist(nums.mt);
}
double Parameters::get_cont_non_damage_e(Constant& nums){
  std::uniform_real_distribution<> dist(0.01, 0.03);
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
  double xy=0.0, xx=0.0, average=0.0;
  //double yy=0.0;
  for(double i: mutation){average+=i;}
  average =(double) average/vect_size;
  double x_ave = (double)(1+vect_size)/2;
  for(std::size_t t=0; t < vect_size; t++){
    xx += (x_ave-t)*(x_ave-t);
    //yy += (average-mutation[t])*(average-mutation[t]);
    xy += (x_ave-t)*(average-mutation[t]);
  }
  if((double)average/10000.0 > std::fabs(xy/xx)){focal=false;}
  //if(std::fabs(xy/(std::pow(xx,0.5)*std::pow(yy,0.5))) < 0.1){focal=false;}
  return focal;
}
