#include"individual.hpp"
#include"population.hpp"

int main()
{
  using std::cout;
  //const double mean_onset_age = 61.5;
  //const double age_sd = 13.5;
  Parameters param;
  param.set_damage_mt();
  Population population(param);
  for(int t=1; 100 >= t; t++){
    population.one_generation(param);
    int rare_mut_num = 0;
    for(int &n: population.tsg_non_mutation_count(param)){
      if((n>0) && (n<(param.N *0.05/100))){rare_mut_num+=n;}
    }
    cout << rare_mut_num << std::endl;
  }
  return 0;
}
