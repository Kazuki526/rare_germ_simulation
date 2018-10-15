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
  population.regression_onset_age(nums, param);
  bool accept_reject = accept_reject_judge(nums, population);
  outfile << param.mutation_rate_coef <<"\t";
  outfile << param.mutater_effect <<"\t";
  outfile << param.mutater_mutation_rate <<"\t";
  outfile << param.mutater_damage <<"\t";
  outfile << param.tsg_non_damage_e <<"\t";
  outfile << param.cont_non_damage_e <<"\t";
  outfile << param.cont_non_fitness_only <<"\t";
  outfile << param.fitness_coef <<"\t";
  outfile << population.rare_tsg_non_freq << "\t";
  outfile << population.tsg_non_regression <<"\t";
  outfile << population.rare_tsg_syn_freq << "\t";
  outfile << population.tsg_syn_regression <<"\t";
  outfile << population.rare_cont_non_freq << "\t";
  outfile << population.cont_non_regression <<"\t";
  outfile << accept_reject << std::endl;
}

void one_set(Constant& nums,std::ofstream& outfile){
  Parameters param(nums);
  Population population(nums);
  for(int time=0; time < 15; time++){
    param.change_param(nums, time);
    one_replicate(nums, param, population, outfile);
    std::cout << "done " << time+1 << " time\n";
  }
}

int main()
{
  std::ofstream outfile;
  outfile.open("simulation_result.tsv", std::ios::out);
  outfile << "mutation_rate_coef\tmutater_effect\tmutater_mutation_rate\t";
  outfile << "mutater_damage\ttsg_non_damage_e\tcont_non_damage_e\t";
  outfile << "cont_non_fitness_only\tfitness_coef\t";
  outfile << "tsg_non_num\ttsg_non_regression\t";
  outfile << "tsg_syn_num\ttsg_syn_regression\t";
  outfile << "cont_non_num\tcont_non_regression\taccept_reject\n";
  Constant nums;
  for(int set_num=0; set_num<6; set_num++){
    one_set(nums, outfile);
  }

  return 0;
}
