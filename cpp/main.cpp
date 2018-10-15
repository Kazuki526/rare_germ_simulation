#include <fstream>
#include "population.hpp"

void one_replicate(Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile){
  population.set_params(nums,param);
  bool tn_equiv=true, ts_equiv=true, cn_equiv=true;
  std::vector<double> tsg_non;
  std::vector<double> tsg_syn;
  std::vector<double> cont_non;
  int t=0;
  while((tn_equiv || ts_equiv || cn_equiv) && t < 500){
    if(t!=0 && t%20==0){
      population.next_generation(nums, param, true);
    }else{
        population.next_generation(nums, param);
    }
    if(t<100){
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_syn.push_back(population.rare_tsg_syn_freq);
      cont_non.push_back(population.rare_cont_non_freq);
    }else{
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_non.erase(tsg_non.begin());
      tsg_syn.push_back(population.rare_tsg_syn_freq);
      tsg_syn.erase(tsg_syn.begin());
      cont_non.push_back(population.rare_cont_non_freq);
      cont_non.erase(cont_non.begin());
      tn_equiv = equiv_lm(tsg_non);
      ts_equiv = equiv_lm(tsg_syn);
      cn_equiv = equiv_lm(cont_non);
    }
    t++;
  }
  outfile << population.rare_tsg_non_freq << "\t";
  outfile << population.rare_tsg_syn_freq << "\t";
  outfile << population.rare_cont_non_freq << "\n";
}

int main()
{
  std::ofstream outfile;
  outfile.open("equiv_out.tsv", std::ios::out);
  Constant nums;
  Parameters param;
  param.set_damage_mt(nums);
  Population population(nums);
  one_replicate(nums, param, population, outfile);

  return 0;
}
