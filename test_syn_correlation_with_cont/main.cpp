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
  outfile << param.cont_non_damage_e <<"\t";
  outfile << population.rare_tsg_non_freq << "\t";
  outfile << population.rare_tsg_non_sd <<"\t";
  outfile << population.rare_tsg_syn_freq << "\t";
  outfile << population.rare_tsg_syn_sd <<"\t";
  outfile << population.rare_cont_non_freq << "\t";
  outfile << population.rare_cont_non_sd <<"\t";
  outfile << population.mut0_rare_non_num<<"\t";
  outfile << population.mut1_rare_non_num<<"\t";
  outfile << population.mut2_rare_non_num<<"\t";
  outfile << population.mut0_notrare_non_num<<"\t";
  outfile << population.mut1_notrare_non_num<<"\t";
  outfile << population.mut2_notrare_non_num<<"\t";
  outfile << population.mut0_rare_syn_num<<"\t";
  outfile << population.mut1_rare_syn_num<<"\t";
  outfile << population.mut2_rare_syn_num<<"\t";
  outfile << population.mut0_notrare_syn_num<<"\t";
  outfile << population.mut1_notrare_syn_num<<"\t";
  outfile << population.mut2_notrare_syn_num<<"\t";
  outfile << population.mut0_rare_cont_num<<"\t";
  outfile << population.mut1_rare_cont_num<<"\t";
  outfile << population.mut2_rare_cont_num<<"\t";
  outfile << population.mut0_notrare_cont_num<<"\t";
  outfile << population.mut1_notrare_cont_num<<"\t";
  outfile << population.mut2_notrare_cont_num<<"\t";
  outfile << population.mutater_freq << "\t";
  outfile << population.mutater_sd << "\t";
  outfile << population.mutation_rate_ave << "\t";
  outfile << population.mutation_rate_sd << "\t";
  outfile << population.rare_non_syn_correlation << "\t";
  outfile << population.rare_non_cont_correlation << "\n";
  outfile.flush();
}

bool one_replicate(Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile){
  bool tn_equiv=true, ts_equiv=true, cn_equiv=true;
  std::vector<double> tsg_non;
  std::vector<double> tsg_syn;
  std::vector<double> cont_non;
  int t=0;
  bool over_mutation=false;
  while((tn_equiv || ts_equiv || cn_equiv) && t < 2000){
    t++;
    if(t%20==0){
      population.next_generation(nums, param, true);
    }else{
        population.next_generation(nums, param);
    }

    if(t<500){
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
    if((population.rare_tsg_non_freq > nums.rare_tsg_non_num*2)||
       (population.rare_tsg_syn_freq > nums.rare_tsg_syn_num*2)||
       (population.rare_cont_non_freq > nums.rare_cont_non_num*2)){
         over_mutation=true;break;
       }
  }
  if(over_mutation){
    return(!over_mutation); // return false
  }else{
    std::vector<double> tn_num,tn_sd,ts_num,ts_sd,cn_num,cn_sd;
    std::vector<double> nscor,nccor,mut,mut_sd,mutr,mutr_sd;
    std::vector<double> tn_0r,tn_1r,tn_2r,tn_0nr,tn_1nr,tn_2nr;
    std::vector<double> ts_0r,ts_1r,ts_2r,ts_0nr,ts_1nr,ts_2nr;
    std::vector<double> cn_0r,cn_1r,cn_2r,cn_0nr,cn_1nr,cn_2nr;
    for(int time=0; time <= 90; time++){
      if(time %10 ==0){
        population.correlation_ns();
        //print_out(nums,param,population,outfile);
        tn_num.push_back(population.rare_tsg_non_freq);
        tn_sd.push_back(population.rare_tsg_non_sd);
        ts_num.push_back(population.rare_tsg_syn_freq);
        ts_sd.push_back(population.rare_tsg_syn_sd);
        cn_num.push_back(population.rare_cont_non_freq);
        cn_sd.push_back(population.rare_cont_non_sd);
        nscor.push_back(population.rare_non_syn_correlation);
        nccor.push_back(population.rare_non_cont_correlation);
        mut.push_back(population.mutater_freq);
        mut_sd.push_back(population.mutater_sd);
        mutr.push_back(population.mutation_rate_ave);
        mutr_sd.push_back(population.mutation_rate_sd);
        tn_0r.push_back(population.mut0_rare_non_num);
        tn_1r.push_back(population.mut1_rare_non_num);
        tn_2r.push_back(population.mut2_rare_non_num);
        tn_0nr.push_back(population.mut0_notrare_non_num);
        tn_1nr.push_back(population.mut1_notrare_non_num);
        tn_2nr.push_back(population.mut2_notrare_non_num);
        ts_0r.push_back(population.mut0_rare_syn_num);
        ts_1r.push_back(population.mut1_rare_syn_num);
        ts_2r.push_back(population.mut2_rare_syn_num);
        ts_0nr.push_back(population.mut0_notrare_syn_num);
        ts_1nr.push_back(population.mut1_notrare_syn_num);
        cn_0r.push_back(population.mut0_rare_cont_num);
        cn_1r.push_back(population.mut1_rare_cont_num);
        cn_2r.push_back(population.mut2_rare_cont_num);
        cn_0nr.push_back(population.mut0_notrare_cont_num);
        cn_1nr.push_back(population.mut1_notrare_cont_num);
        cn_2nr.push_back(population.mut2_notrare_cont_num);
        if(population.mut2_rare_non_num != 0.0 &&
           population.mut2_notrare_non_num != 0.0 &&
           population.mut2_rare_cont_num != 0.0 &&
           population.mut2_notrare_cont_num != 0.0 &&
           population.mut2_rare_syn_num != 0.0 &&
           population.mut2_notrare_syn_num != 0.0){
             tn_2r.push_back(population.mut2_rare_non_num);
             tn_2nr.push_back(population.mut2_notrare_non_num);
             ts_2r.push_back(population.mut2_rare_syn_num);
             ts_2nr.push_back(population.mut2_notrare_syn_num);
             cn_2r.push_back(population.mut2_rare_cont_num);
             cn_2nr.push_back(population.mut2_notrare_cont_num);
           }
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
    population.rare_cont_non_freq = (double)std::accumulate(cn_num.begin(),cn_num.end(),0.0) /cn_num.size();
    population.rare_cont_non_sd = (double)std::accumulate(cn_sd.begin(),cn_sd.end(),0.0) /cn_sd.size();
    population.rare_non_syn_correlation = (double)std::accumulate(nscor.begin(),nscor.end(),0.0) /nscor.size();
    population.rare_non_cont_correlation = (double)std::accumulate(nccor.begin(),nccor.end(),0.0) /nccor.size();
    population.mutater_freq = (double)std::accumulate(mut.begin(),mut.end(),0.0) /mut.size();
    population.mutater_sd = (double)std::accumulate(mut_sd.begin(),mut_sd.end(),0.0) /mut_sd.size();
    population.mutation_rate_ave = (double)std::accumulate(mutr.begin(),mutr.end(),0.0) /mutr.size();
    population.mutation_rate_sd = (double)std::accumulate(mutr_sd.begin(),mutr_sd.end(),0.0) /mutr_sd.size();
    population.mut0_rare_non_num = (double)std::accumulate(tn_0r.begin(),tn_0r.end(),0.0) /tn_0r.size();
    population.mut1_rare_non_num = (double)std::accumulate(tn_1r.begin(),tn_1r.end(),0.0) /tn_1r.size();
    population.mut2_rare_non_num = (double)std::accumulate(tn_2r.begin(),tn_2r.end(),0.0) /tn_2r.size();
    population.mut0_notrare_non_num = (double)std::accumulate(tn_0nr.begin(),tn_0nr.end(),0.0) /tn_0nr.size();
    population.mut1_notrare_non_num = (double)std::accumulate(tn_1nr.begin(),tn_1nr.end(),0.0) /tn_1nr.size();
    population.mut2_notrare_non_num = (double)std::accumulate(tn_2nr.begin(),tn_2nr.end(),0.0) /tn_2nr.size();
    population.mut0_rare_syn_num = (double)std::accumulate(ts_0r.begin(),ts_0r.end(),0.0) /ts_0r.size();
    population.mut1_rare_syn_num = (double)std::accumulate(ts_1r.begin(),ts_1r.end(),0.0) /ts_1r.size();
    population.mut2_rare_syn_num = (double)std::accumulate(ts_2r.begin(),ts_2r.end(),0.0) /ts_2r.size();
    population.mut0_notrare_syn_num = (double)std::accumulate(ts_0nr.begin(),ts_0nr.end(),0.0) /ts_0nr.size();
    population.mut1_notrare_syn_num = (double)std::accumulate(ts_1nr.begin(),ts_1nr.end(),0.0) /ts_1nr.size();
    population.mut2_notrare_syn_num = (double)std::accumulate(ts_2nr.begin(),ts_2nr.end(),0.0) /ts_2nr.size();
    population.mut0_rare_cont_num = (double)std::accumulate(cn_0r.begin(),cn_0r.end(),0.0) /cn_0r.size();
    population.mut1_rare_cont_num = (double)std::accumulate(cn_1r.begin(),cn_1r.end(),0.0) /cn_1r.size();
    population.mut2_rare_cont_num = (double)std::accumulate(cn_2r.begin(),cn_2r.end(),0.0) /cn_2r.size();
    population.mut0_notrare_cont_num = (double)std::accumulate(cn_0nr.begin(),cn_0nr.end(),0.0) /cn_0nr.size();
    population.mut1_notrare_cont_num = (double)std::accumulate(cn_1nr.begin(),cn_1nr.end(),0.0) /cn_1nr.size();
    population.mut2_notrare_cont_num = (double)std::accumulate(cn_2nr.begin(),cn_2nr.end(),0.0) /cn_2nr.size();
    print_out(param,population,outfile);
    return(!over_mutation); // return true
  }
}
int main(int argc,char *argv[]){
  std::ofstream outfile;
  outfile.open("simulation_result"+std::string(argv[1])+".tsv", std::ios::out);
  outfile << "generation\tmutation_rate\tmutater_effect\t";
  outfile << "mutater_mutation_rate\tmutater_damage\ttsg_non_damage_e\tcont_non_damage_e\t";
  outfile << "tsg_non_num\ttsg_non_sd\ttsg_syn_num\ttsg_syn_sd\tcont_non_num\tcont_non_sd\t";
  outfile << "mut0_rare_non_num\tmut1_rare_non_num\tmut2_rare_non_num\t";
  outfile << "mut0_notrare_non_num\tmut1_notrare_non_num\tmut2_notrare_non_num\t";
  outfile << "mut0_rare_syn_num\tmut1_rare_syn_num\tmut2_rare_syn_num\t";
  outfile << "mut0_notrare_syn_num\tmut1_notrare_syn_num\tmut2_notrare_syn_num\t";
  outfile << "mut0_rare_cont_num\tmut1_rare_cont_num\tmut2_rare_cont_num\t";
  outfile << "mut0_notrare_cont_num\tmut1_notrare_cont_num\tmut2_notrare_cont_num\t";
  outfile << "mutater_freq\tmutater_sd\tmutation_rate_freq\tmutation_rate_sd\t";
  outfile << "non_syn_correlation\ttsg_cont_correlation\n";
  Constant nums;
  int time=1;
  while(time <= 2500){
    Parameters param(nums);
    while(param.expected_mutation_sd<0.0000015){param.reset(nums);}
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
