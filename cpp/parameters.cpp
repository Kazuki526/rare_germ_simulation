#include "parameters.hpp"

Parameters::Parameters(Constant& nums){
  mutation_rate = set_mutation_rate(nums);
  mutator_effect = set_mutator_effect(nums);
  mutator_mutation_rate = set_mutator_mutation_rate(nums);
  mutator_damage = set_mutator_damage(nums);
  non_damage_e = set_non_damage_e(nums);

  /*double m=mutator_mutation_rate, s=mutator_damage, e=mutator_effect;
  double p = m/s+0.5-m/2-std::pow((m/s+(1-m)/2)*(m/s+(1-m)/2)-m/s,0.5);
  double mutr = ((1-p)*(1-p)+mutator_effect*2*p*(1-p)+e*e*(p*p))*mutation_rate;
  expected_mutation_sd =
    std::pow((mutation_rate-mutr)*(mutation_rate-mutr)*(1-p)*(1-p)+
             (e*mutation_rate-mutr)*(e*mutation_rate-mutr)*2*p*(1-p)+
             (e*e*mutation_rate-mutr)*(e*e*mutation_rate-mutr)*p*p,0.5);*/
  mutator_s = mutation_rate*mutator_effect*mutator_damage - mutation_rate;
}
void Parameters::reset(Constant& nums){
  mutation_rate = set_mutation_rate(nums);
  mutator_effect = set_mutator_effect(nums);
  mutator_mutation_rate = set_mutator_mutation_rate(nums);
  mutator_damage = set_mutator_damage(nums);
  non_damage_e = set_non_damage_e(nums);

  double m=mutator_mutation_rate, s=mutator_damage, e=mutator_effect;
  double p = m/s+0.5-m/2-std::pow((m/s+(1-m)/2)*(m/s+(1-m)/2)-m/s,0.5);
  double mutr = ((1-p)*(1-p)+mutator_effect*2*p*(1-p)+e*e*(p*p))*mutation_rate;
  expected_mutation_sd =
    std::pow((mutation_rate-mutr)*(mutation_rate-mutr)*(1-p)*(1-p)+
             (e*mutation_rate-mutr)*(e*mutation_rate-mutr)*2*p*(1-p)+
             (e*e*mutation_rate-mutr)*(e*e*mutation_rate-mutr)*p*p,0.5);
}

/* set each nonsynonymous site damage */
void Parameters::set_damage(Constant& nums){
  std::exponential_distribution<> non_ex(1.0/non_damage_e);
  for(std::size_t s=0; s < nums.non_site; s++){
    non_damage.push_back(non_ex(nums.mt));
  }
}

/*set each parameters */
double Parameters::set_mutation_rate(Constant& nums){
  //std::uniform_real_distribution<> dist(0, 0.0000001); first
  std::uniform_real_distribution<> dist(0, 0.0000001);
  return dist(nums.mt);
}
double Parameters::set_mutator_effect(Constant& nums){
  std::uniform_real_distribution<> dist(2, 10000);
  return dist(nums.mt);
}
double Parameters::set_non_damage_e(Constant& nums){
  //std::uniform_real_distribution<> dist(0, 0.05); first
  std::uniform_real_distribution<> dist(0, 0.06);
  return dist(nums.mt);
}
double Parameters::set_mutator_damage(Constant& nums){
  std::uniform_real_distribution<> dist(1,100000);
  return dist(nums.mt);
}
double Parameters::set_mutator_mutation_rate(Constant& nums){
  std::uniform_real_distribution<> dist(0, 0.1);
  return dist(nums.mt);
}




double log_random(double start, double end, std::mt19937& mt){
  start = std::log10(start==0 ? 0.000000001 : start);
  end   = std::log10(end);
  std::uniform_real_distribution<> dist(start, end);
  return std::pow(10, dist(mt));
}

bool equib_lm(const std::vector<double>& mutation) {
  bool focal=true;
  std::size_t vect_size = mutation.size();
  double xy=0.0, xx=0.0, yy=0.0, average=0.0;
  for(double i: mutation){average+=i;}
  average =(double) average/vect_size;
  double x_ave = (double)(1+vect_size)/2;
  for(std::size_t t=0; t < vect_size; t++){
    xx += (x_ave-t)*(x_ave-t);
    //yy += (average-mutation[t])*(average-mutation[t]);
    xy += (x_ave-t)*(average-mutation[t]);
  }
  if(xy/xx < 0.0){focal=false;}
  //if(std::fabs(xy/(std::pow(xx,0.5)*std::pow(yy,0.5))) < 0.1){focal=false;}
  return focal;
}
