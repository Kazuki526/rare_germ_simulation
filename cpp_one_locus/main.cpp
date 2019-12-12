#include <fstream>
#include <stdio.h>
#include <string>
#include "population.hpp"

void print_out(const Parameters& param,
               const Population& population, std::ofstream& outfile){
  outfile << param.mutation_rate <<"\t";
  outfile << param.mutator_effect <<"\t";
  outfile << param.mutator_mutation_rate <<"\t";
  outfile << param.mutator_damage <<"\t";
  outfile << param.tsg_non_damage_e <<"\t";
  outfile << population.rare_tsg_non_freq << "\t";
  outfile << population.rare_tsg_non_sd <<"\t";
  outfile << population.rare_tsg_syn_freq << "\t";
  outfile << population.rare_tsg_syn_sd <<"\t";
  //outfile << population.mutator_freq << "\t";
  //outfile << population.mutator_sd << "\t";
  outfile << population.mutation_rate_ave << "\t";
  outfile << population.mutation_rate_sd << "\t";
  outfile << population.rare_non_syn_correlation << "\t";
  outfile << population.rare_num_reg << "\t";
  outfile << population.rare_num_reg_zero << "\n";
  outfile.flush();
}

void one_replicate(Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile){
  bool mut_equiv=true, tn_equiv=true, ts_equiv=true;
  std::vector<double> mutator;
  std::vector<double> tsg_non;
  std::vector<double> tsg_syn;
  int t=0;
  while(mut_equiv || tn_equiv || ts_equiv){
    t++;
    if(t%20==0){
      population.next_generation(nums, param, true);
    }else{
        population.next_generation(nums, param);
    }

    if(t<=1000){
      mutator.push_back(population.mutator_freq);
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_syn.push_back(population.rare_tsg_syn_freq);
    }else{
      mutator.push_back(population.mutator_freq);
      mutator.erase(mutator.begin());
      tsg_non.push_back(population.rare_tsg_non_freq);
      tsg_non.erase(tsg_non.begin());
      tsg_syn.push_back(population.rare_tsg_syn_freq);
      tsg_syn.erase(tsg_syn.begin());
      if(mut_equiv){mut_equiv = equib_lm(mutator);}
      if(tn_equiv) {tn_equiv = equib_lm(tsg_non);}
      if(ts_equiv) {ts_equiv = equib_lm(tsg_syn);}
    }
    if((population.rare_tsg_non_freq > nums.rare_tsg_non_num*2)||
       (population.rare_tsg_syn_freq > nums.rare_tsg_syn_num*2)){break;}
    //std::cout <<t<<" "<<population.rare_tsg_syn_freq<<" "<<population.mutator_freq<<"\n";std::cout.flush();
    //if(t % 100 == 0){std::cout << "now " << t << " generation " << mut_equiv << tn_equiv << ts_equiv << "\n";std::cout.flush();}
  }
  //std::cout << t << " generation => finish and go sumary\n";std::cout.flush();
  //std::vector<double> tn_num,tn_sd,ts_num,ts_sd,cor,mut,mut_sd,mutr,mutr_sd,rare_num_reg, rare_num_reg_zero;
  std::vector<double> tn_num,tn_sd,ts_num,ts_sd,cor,mutr,mutr_sd,rare_num_reg, rare_num_reg_zero;
  std::vector<double> tn_0r,tn_1r,tn_2r,tn_0nr,tn_1nr,tn_2nr,ts_0r,ts_1r,ts_2r,ts_0nr,ts_1nr,ts_2nr;
  for(int time=0; time <= 500; time++){
    if(time %10 ==0){
      population.correlation_ns(nums);
      //print_out(nums,param,population,outfile);
      tn_num.push_back(population.rare_tsg_non_freq);
      tn_sd.push_back(population.rare_tsg_non_sd);
      ts_num.push_back(population.rare_tsg_syn_freq);
      ts_sd.push_back(population.rare_tsg_syn_sd);
      cor.push_back(population.rare_non_syn_correlation);
      rare_num_reg.push_back(population.rare_num_reg);
      rare_num_reg_zero.push_back(population.rare_num_reg_zero);
      //mut.push_back(population.mutator_freq);
      //mut_sd.push_back(population.mutator_sd);
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
  //population.mutator_freq = (double)std::accumulate(mut.begin(),mut.end(),0.0) /mut.size();
  //population.mutator_sd = (double)std::accumulate(mut_sd.begin(),mut_sd.end(),0.0) /mut_sd.size();
  population.mutation_rate_ave = (double)std::accumulate(mutr.begin(),mutr.end(),0.0) /mutr.size();
  population.mutation_rate_sd = (double)std::accumulate(mutr_sd.begin(),mutr_sd.end(),0.0) /mutr_sd.size();
  print_out(param,population,outfile);
}
int main(int argc,char *argv[])
{
  std::ofstream outfile;
  outfile.open("simulation_result"+std::string(argv[1])+".tsv", std::ios::out);
  outfile << "generation\tmutation_rate\tmutator_effect\t";
  outfile << "mutator_mutation_rate\tmutator_damage\ttsg_non_damage_e\t";
  outfile << "tsg_non_num\ttsg_non_sd\t";
  outfile << "tsg_syn_num\ttsg_syn_sd\t";
  //outfile << "mutator_freq\tmutator_sd\tmutation_rate_freq\tmutation_rate_sd\t";
  outfile << "mutation_rate_freq\tmutation_rate_sd\t";
  outfile << "correlation\treg_R\treg_R_zero\n";
  Constant nums;
  int time=1;
  while(time <=100){
    Parameters param(nums);
    //int t=0;
    //while(param.expected_mutation_sd<0.0000005){param.reset(nums);t++;}
    //std::cout << t << " time reparam\t";
    param.set_damage(nums);
    Population population(nums,param);
    one_replicate(nums, param, population, outfile);
    std::cout << "done " << time << " time\n";std::cout.flush();
    time++;
  }

  return 0;
}
