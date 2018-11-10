#include <fstream>
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
  outfile << population.mutater_freq << "\t";
  outfile << population.mutater_sd << "\t";
  outfile << population.mutation_rate_ave << "\t";
  outfile << population.mutation_rate_sd << "\t";
  outfile << population.rare_non_syn_correlation << "\n";
  outfile.flush();
}

bool one_replicate(Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile){
  bool tn_equiv=true, ts_equiv=true;
  std::vector<double> tsg_non;
  std::vector<double> tsg_syn;
  int t=0;
  bool over_mutation=false;
  while((tn_equiv || ts_equiv) && t < 2000){
    t++;
    if(t%20==0){
      population.next_generation(nums, param, true);
    }else{
        population.next_generation(nums, param);
    }

    if(t<200){
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
  }
  if(over_mutation){
    return(!over_mutation); // return false
  }else{
    std::vector<double> tn_num,tn_sd,ts_num,ts_sd,cor,mut,mut_sd,mutr,mutr_sd;
    std::vector<double> tn_0r,tn_1r,tn_2r,tn_0nr,tn_1nr,tn_2nr,ts_0r,ts_1r,ts_2r,ts_0nr,ts_1nr,ts_2nr;
    for(int time=0; time <= 90; time++){
      if(time %10 ==0){
        population.correlation_ns();
        //print_out(nums,param,population,outfile);
        tn_num.push_back(population.rare_tsg_non_freq);
        tn_sd.push_back(population.rare_tsg_non_sd);
        ts_num.push_back(population.rare_tsg_syn_freq);
        ts_sd.push_back(population.rare_tsg_syn_sd);
        cor.push_back(population.rare_non_syn_correlation);
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
        if(population.mut2_rare_non_num != 0.0 &&
           population.mut2_notrare_non_num != 0.0 &&
           population.mut2_rare_syn_num != 0.0 &&
           population.mut2_notrare_syn_num != 0.0){
             tn_2r.push_back(population.mut2_rare_non_num);
             tn_2nr.push_back(population.mut2_notrare_non_num);
             ts_2r.push_back(population.mut2_rare_syn_num);
             ts_2nr.push_back(population.mut2_notrare_syn_num);
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
    population.rare_non_syn_correlation = (double)std::accumulate(cor.begin(),cor.end(),0.0) /cor.size();
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
    print_out(param,population,outfile);
    return(!over_mutation); // return true
  }
}
int main()
{
  std::ofstream outfile;
  outfile.open("simulation_result.tsv", std::ios::out);
  outfile << "generation\tmutation_rate\tmutater_effect\t";
  outfile << "mutater_mutation_rate\tmutater_damage\ttsg_non_damage_e\t";
  outfile << "tsg_non_num\ttsg_non_sd\t";
  outfile << "tsg_syn_num\ttsg_syn_sd\t";
  outfile << "mut0_rare_non_num\tmut1_rare_non_num\tmut2_rare_non_num\t";
  outfile << "mut0_notrare_non_num\tmut1_notrare_non_num\tmut2_notrare_non_num\t";
  outfile << "mut0_rare_syn_num\tmut1_rare_syn_num\tmut2_rare_syn_num\t";
  outfile << "mut0_notrare_syn_num\tmut1_notrare_syn_num\tmut2_notrare_syn_num\t";
  outfile << "mutater_freq\tmutater_sd\tmutation_rate_freq\tmutation_rate_sd\t";
  outfile << "correlation\n";
  Constant nums;
  int time=1;
  while(time <=2000){
    Parameters param(nums);
    while(param.expected_mutation_sd<0.0000012){param.reset(nums);}
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
