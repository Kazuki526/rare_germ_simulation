#include "parameters.hpp"

void Parameters::set_damage_mt(Constant& nums){
  std::exponential_distribution<> tsg_ex(tsg_non_damage_e);
  std::exponential_distribution<> cont_ex(cont_non_damage_e);
  for(std::size_t s=0; nums.tsg_non_site > s; s++){
    tsg_non_damage.push_back(tsg_ex(nums.mt));
  }
  for(std::size_t s=0; nums.cont_non_site > s; s++){
    cont_non_damage.push_back(cont_ex(nums.mt));
  }
}
void Parameters::change_param(std::size_t p, double doub){

}
