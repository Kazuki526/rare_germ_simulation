#include <fstream>
#include "population.hpp"

void print_out(const Parameters& param,
               const Population& population, std::ofstream& outfile){
  outfile << param.mutater_locas << "\t";
  outfile << param.mutater_effect <<"\t";
  outfile << param.mutater_mutation_rate <<"\t";
  outfile << param.mutater_damage <<"\t";
  outfile << param.tsg_non_damage_e <<"\t";
  outfile << population.rare_tsg_non_freq << "\t";
  outfile << population.rare_tsg_non_sd <<"\t";
  outfile << population.rare_tsg_syn_freq << "\t";
  outfile << population.rare_tsg_syn_sd <<"\t";
  outfile << population.rare_non_syn_correlation << "\n";
  outfile.flush();
}

bool one_replicate(Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile){
  population.set_params(nums,param);
  bool tn_equiv=true, ts_equiv=true;
  std::vector<double> tsg_non;
  std::vector<double> tsg_syn;
  int t=0;
  bool over_mutation=false;
  while((tn_equiv || ts_equiv) && t < 1000){
    if(t!=0 && t%20==0){
      population.next_generation(nums, param, true);
    }else{
        population.next_generation(nums, param);
    }
    if(t<100){
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_syn.push_back(population.rare_tsg_syn_freq);
    }else{
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_non.erase(tsg_non.begin());
      tsg_syn.push_back(population.rare_tsg_syn_freq);
      tsg_syn.erase(tsg_syn.begin());
      tn_equiv = equiv_lm(tsg_non);
      ts_equiv = equiv_lm(tsg_syn);
    }
    if((population.rare_tsg_non_freq > nums.rare_tsg_non_num*2)||
       (population.rare_tsg_syn_freq > nums.rare_tsg_syn_num*2)){
         over_mutation=true;break;
       }
    t++;
  }
  if(over_mutation){
    return(!over_mutation); // return false
  }else{
    std::vector<double> tn_num,tn_sd,ts_num,ts_sd,cor;
    for(int time=0; time <= 90; time++){
      if(time %10 ==0){
        population.correlation_ns();
        //print_out(nums,param,population,outfile);
        tn_num.push_back(population.rare_tsg_non_freq);
        tn_sd.push_back(population.rare_tsg_non_sd);
        ts_num.push_back(population.rare_tsg_syn_freq);
        ts_sd.push_back(population.rare_tsg_syn_sd);
        cor.push_back(population.rare_non_syn_correlation);
        population.next_generation(nums,param,true);
      }else{
        population.next_generation(nums,param);
      }
    }
    population.rare_tsg_non_freq = (double)std::accumulate(tn_num.begin(),tn_num.end(),0.0) /tn_num.size();
    population.rare_tsg_non_sd = (double)std::accumulate(tn_sd.begin(),tn_sd.end(),0.0) /tn_sd.size();
    population.rare_tsg_syn_freq = (double)std::accumulate(ts_num.begin(),ts_num.end(),0.0) /ts_num.size();
    population.rare_tsg_syn_sd = (double)std::accumulate(ts_sd.begin(),ts_sd.end(),0.0) /ts_sd.size();
    population.rare_non_syn_correlation = (double)std::accumulate(cor.begin(),cor.end(),0.0) /cor.size();
    print_out(param,population,outfile);
    return(!over_mutation); // return true
  }
}
int main()
{
  std::ofstream outfile;
  outfile.open("simulation_result.tsv", std::ios::out);
  outfile << "mutation_locas\tmutater_effect\tmutater_mutation_rate\t";
  outfile << "mutater_damage\ttsg_non_damage_e\t";
  outfile << "tsg_non_num\ttsg_non_sd\t";
  outfile << "tsg_syn_num\ttsg_syn_sd\t";
  outfile << "correlation\n";
  Constant nums;
  int time=0;
  while(time <4000){
    Parameters param(nums);
    Population population(nums);
    bool replicate_result = one_replicate(nums, param, population, outfile);
    if(replicate_result){
      std::cout << "done" << time+1 << "time\n";
      time++;
    }
  }

  return 0;
}
