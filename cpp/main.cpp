#include"individual.hpp"
#include"population.hpp"

int main()
{
  using std::cout;
  //const double mean_onset_age = 61.5;
  //const double age_sd = 13.5;
  Constant nums;
  Parameters param;
  param.set_damage_mt(nums);
  Population population(nums, param);
  for(int t=1; t <= 100; t++){
    population.next_generation(nums, param);
    cout << population.rare_tsg_non_freq << std::endl;
  }
  return 0;
}
