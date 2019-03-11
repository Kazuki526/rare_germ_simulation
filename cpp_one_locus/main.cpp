#include <fstream>
#include <stdio.h>
#include <string>
#include "population.hpp"

void print_out(const Parameters& param,
               const Population& population, std::ofstream& outfile){
  outfile << param.mutation_rate <<"\t";
  outfile << param.mutater_effect <<"\t";
  outfile << param.mutater_mutation_rate <<"\t";
  outfile << param.mutater_damage <<"\t";
  outfile << param.tsg_non_damage_e <<"\t";
  outfile << population.rare_tsg_non_freq << "\t";
  outfile << population.rare_tsg_non_sd <<"\t";
  outfile << population.rare_tsg_syn_freq << "\t";
  outfile << population.rare_tsg_syn_sd <<"\t";
  outfile << population.mutater_freq << "\t";
  outfile << population.mutater_sd << "\t";
  outfile << population.mutation_rate_ave << "\t";
  outfile << population.mutation_rate_sd << "\t";
  outfile << population.rare_non_syn_correlation << "\t";
  outfile << population.rare_num_reg << "\t";
  outfile << population.rare_num_reg_zero << "\n";
  outfile.flush();
}

bool one_replicate(Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile){
  bool mut_equiv=true, tn_equiv=true, ts_equiv=true;
  std::vector<double> mutater;
  std::vector<double> tsg_non;
  std::vector<double> tsg_syn;
  int t=0;
  bool over_mutation=false;
  while((mut_equiv || tn_equiv || ts_equiv) && t < 2000){
    t++;
    if(t%20==0){
      population.next_generation(nums, param, true);
    }else{
        population.next_generation(nums, param);
    }

    if(t<200){
      mutater.push_back(population.mutater_freq);
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_syn.push_back(population.rare_tsg_syn_freq);
    }else{
      mutater.push_back(population.mutater_freq);
      mutater.erase(mutater.begin());
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_non.erase(tsg_non.begin());
      tsg_syn.push_back(population.rare_tsg_syn_freq);
      tsg_syn.erase(tsg_syn.begin());
      mut_equiv = equib_lm(mutater);
      tn_equiv = equib_lm(tsg_non);
      ts_equiv = equib_lm(tsg_syn);
    }
    if((population.rare_tsg_non_freq > nums.rare_tsg_non_num*2)||
       (population.rare_tsg_syn_freq > nums.rare_tsg_syn_num*2)){
         over_mutation=true;break;
    }
    //std::cout <<t<<" "<<population.rare_tsg_syn_freq<<" "<<population.mutater_freq<<"\n";std::cout.flush();
  }
  if(over_mutation){
    return(!over_mutation); // return false
  }else{
    std::vector<double> tn_num,tn_sd,ts_num,ts_sd,cor,mut,mut_sd,mutr,mutr_sd,rare_num_reg, rare_num_reg_zero;
    std::vector<double> tn_0r,tn_1r,tn_2r,tn_0nr,tn_1nr,tn_2nr,ts_0r,ts_1r,ts_2r,ts_0nr,ts_1nr,ts_2nr;
    for(int time=0; time <= 500; time++){
      if(time %10 ==0){
        population.correlation_ns();
        //print_out(nums,param,population,outfile);
        tn_num.push_back(population.rare_tsg_non_freq);
        tn_sd.push_back(population.rare_tsg_non_sd);
        ts_num.push_back(population.rare_tsg_syn_freq);
        ts_sd.push_back(population.rare_tsg_syn_sd);
        cor.push_back(population.rare_non_syn_correlation);
        rare_num_reg.push_back(population.rare_num_reg);
        rare_num_reg_zero.push_back(population.rare_num_reg_zero);
        mut.push_back(population.mutater_freq);
        mut_sd.push_back(population.mutater_sd);
        mutr.push_back(population.mutation_rate_ave);
        mutr_sd.push_back(population.mutation_rate_sd);
        population.next_generation(nums,param,true);
      }else{
        population.next_generation(nums,param);
      }
    }
    outfile << t <<"\t";
    if(tn_2r.empty()){tn_2r.push_back(0);tn_2nr.push_back(0);ts_2r.push_back(0);ts_2nr.push_back(0);}
    population.rare_tsg_non_freq = (double)std::accumulate(tn_num.begin(),tn_num.end(),0.0) /tn_num.size();
    population.rare_tsg_non_sd = (double)std::accumulate(tn_sd.begin(),tn_sd.end(),0.0) /tn_sd.size();
    population.rare_tsg_syn_freq = (double)std::accumulate(ts_num.begin(),ts_num.end(),0.0) /ts_num.size();
    population.rare_tsg_syn_sd = (double)std::accumulate(ts_sd.begin(),ts_sd.end(),0.0) /ts_sd.size();
    population.rare_non_syn_correlation = (double)std::accumulate(cor.begin(),cor.end(),0.0) /cor.size();
    population.rare_num_reg =(double)std::accumulate(rare_num_reg.begin(),rare_num_reg.end(),0.0) /rare_num_reg.size();
    population.rare_num_reg_zero = (double)std::accumulate(rare_num_reg_zero.begin(),rare_num_reg_zero.end(),0.0) /rare_num_reg_zero.size();
    population.mutater_freq = (double)std::accumulate(mut.begin(),mut.end(),0.0) /mut.size();
    population.mutater_sd = (double)std::accumulate(mut_sd.begin(),mut_sd.end(),0.0) /mut_sd.size();
    population.mutation_rate_ave = (double)std::accumulate(mutr.begin(),mutr.end(),0.0) /mutr.size();
    population.mutation_rate_sd = (double)std::accumulate(mutr_sd.begin(),mutr_sd.end(),0.0) /mutr_sd.size();
    print_out(param,population,outfile);
    return(!over_mutation); // return true
  }
}
int main(int argc,char *argv[])
{
  std::ofstream outfile;
  outfile.open("simulation_result"+std::string(argv[1])+".tsv", std::ios::out);
  outfile << "generation\tmutation_rate\tmutater_effect\t";
  outfile << "mutater_mutation_rate\tmutater_damage\ttsg_non_damage_e\t";
  outfile << "tsg_non_num\ttsg_non_sd\t";
  outfile << "tsg_syn_num\ttsg_syn_sd\t";
  outfile << "mutater_freq\tmutater_sd\tmutation_rate_freq\tmutation_rate_sd\t";
  outfile << "correlation\treg_R\treg_R_zero\n";
  Constant nums;
  int time=1;
  while(time <=20000){
    Parameters param(nums);
    int t=0;
    while(param.expected_mutation_sd<0.000002){param.reset(nums);t++;}
    std::cout << t << " time reparam\t";
    param.set_damage(nums);
    Population population(nums,param);
    bool replicate_result = one_replicate(nums, param, population, outfile);
    if(replicate_result){
      std::cout << "done" << time << "time\n";
      time++;
    }else{std::cout<<"restart: "<<std::flush;}
  }

  return 0;
}
