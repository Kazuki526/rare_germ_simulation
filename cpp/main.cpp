#include "population.hpp"

void print_out(const Parameters& param,
               const Population& population, std::ofstream& outfile){
  outfile << param.mutation_rate <<"\t";
  outfile << param.mutator_effect <<"\t";
  outfile << param.mutator_mutation_rate <<"\t";
  outfile << param.mutator_damage <<"\t";
  outfile << param.non_damage_e <<"\t";
  outfile << population.rare_non_freq << "\t";
  outfile << population.rare_non_sd <<"\t";
  outfile << population.rare_syn_freq << "\t";
  outfile << population.rare_syn_sd <<"\t";
  outfile << population.mutator_freq << "\t";
  outfile << population.mutation_rate_ave << "\t";
  outfile << population.mutation_rate_sd << "\t";
  outfile << population.rare_non_syn_correlation << "\t";
  outfile << population.rare_num_reg << "\n";
  outfile.flush();
}

void one_replicate(int replicate, Constant& nums, const Parameters& param,
                   Population& population, std::ofstream& outfile,
                   std::ofstream& mutout, std::ofstream& logout){
  logout <<param.mutation_rate<<":"<<param.mutator_effect<<":"<<param.mutator_mutation_rate;
  logout <<":"<<param.mutator_damage<<":"<<param.non_damage_e<<" => ";logout.flush();
  bool mut_equiv=true, non_equiv=true, syn_equiv=true;
  std::vector<double> mutator_freq_list;
  std::vector<double> non_freq_list;
  std::vector<double> syn_freq_list;
  int generation=0;
  while(mut_equiv || non_equiv || syn_equiv){
    generation++;
    if(generation%50==0){
      population.next_generation(nums, param, true);
      logout << generation << "> ";logout.flush();
    }else{
        population.next_generation(nums, param);
    }

    if(generation<=500){
      mutator_freq_list.push_back(population.mutator_freq);
      non_freq_list.push_back(population.rare_non_freq);
      syn_freq_list.push_back(population.rare_syn_freq);
    }else{
      mutator_freq_list.push_back(population.mutator_freq);
      mutator_freq_list.erase(mutator_freq_list.begin());
      non_freq_list.push_back(population.rare_non_freq);
      non_freq_list.erase(non_freq_list.begin());
      syn_freq_list.push_back(population.rare_syn_freq);
      syn_freq_list.erase(syn_freq_list.begin());
      if(mut_equiv){mut_equiv = equib_lm(mutator_freq_list);}
      if(non_equiv) {non_equiv = equib_lm(non_freq_list);}
      if(syn_equiv) {syn_equiv = equib_lm(syn_freq_list);}
    }
    if((population.rare_non_freq > nums.rare_non_num*1.5)||
       (population.rare_syn_freq > nums.rare_syn_num*1.5)){logout<<"stop";logout.flush();break;}
    if((generation>=500)&&
       ((population.rare_non_freq < nums.rare_non_num*0.5)||
        (population.rare_syn_freq < nums.rare_syn_num*0.5))){logout<<"stop";logout.flush();break;}
  }
  //logout << " finish and go sumary\n";logout.flush();
  if(!mut_equiv && !non_equiv && !syn_equiv){ /* not stoped by too much rare variant */
    std::vector<double> n_num,n_sd,s_num,s_sd,cor,mut,mutr,mutr_sd,rare_num_reg;
    for(int aft_gene=1; aft_gene <= 500; aft_gene++){
      generation++;
      if(aft_gene %50 ==0){
        population.next_generation(nums, param, true);
        mutout << replicate << "\t" << generation << "\t";
        population.out_mutator_state(nums, mutout);
      }else{
        population.next_generation(nums, param);
      }
      if(aft_gene %10 ==0){
        logout << generation << "=> ";logout.flush();
        population.correlation_ns(nums);
        //print_out(nums,param,population,outfile);
        n_num.push_back(population.rare_non_freq);
        n_sd.push_back(population.rare_non_sd);
        s_num.push_back(population.rare_syn_freq);
        s_sd.push_back(population.rare_syn_sd);
        cor.push_back(population.rare_non_syn_correlation);
        rare_num_reg.push_back(population.rare_num_reg);
        mut.push_back(population.mutator_freq);
        mutr.push_back(population.mutation_rate_ave);
        mutr_sd.push_back(population.mutation_rate_sd);
      }
    }
    outfile << replicate << "\t";
    outfile << generation <<"\t";
    population.rare_non_freq = (double)std::accumulate(n_num.begin(),n_num.end(),0.0) /n_num.size();
    population.rare_non_sd = (double)std::accumulate(n_sd.begin(),n_sd.end(),0.0) /n_sd.size();
    population.rare_syn_freq = (double)std::accumulate(s_num.begin(),s_num.end(),0.0) /s_num.size();
    population.rare_syn_sd = (double)std::accumulate(s_sd.begin(),s_sd.end(),0.0) /s_sd.size();
    population.rare_non_syn_correlation = (double)std::accumulate(cor.begin(),cor.end(),0.0) /cor.size();
    population.rare_num_reg =(double)std::accumulate(rare_num_reg.begin(),rare_num_reg.end(),0.0) /rare_num_reg.size();
    population.mutator_freq = (double)std::accumulate(mut.begin(),mut.end(),0.0) /mut.size();
    population.mutation_rate_ave = (double)std::accumulate(mutr.begin(),mutr.end(),0.0) /mutr.size();
    population.mutation_rate_sd = (double)std::accumulate(mutr_sd.begin(),mutr_sd.end(),0.0) /mutr_sd.size();
    print_out(param,population,outfile);
  }else{ /* when stoped by too much rare variant */
    population.correlation_ns(nums);
    outfile << replicate << "\t";
    outfile << generation <<"\t";
    population.rare_non_syn_correlation=0;
    population.rare_num_reg=0;
    print_out(param,population,outfile);
  }
}
int main(int argc,char *argv[])
{
  std::ofstream outfile;
  outfile.open("simulation_result"+std::string(argv[1])+".tsv", std::ios::out);
  outfile << "replicate\tgeneration\tmutation_rate\tmutator_effect\t";
  outfile << "mutator_mutation_rate\tmutator_damage\tnon_damage_e\t";
  outfile << "non_num\tnon_sd\t";
  outfile << "syn_num\tsyn_sd\t";
  outfile << "mutator_freq\tmutation_rate_freq\tmutation_rate_sd\t";
  outfile << "correlation\treg_R\n";

  std::ofstream mutout;
  mutout.open("simulation_result_mutator"+std::string(argv[1])+".tsv",std::ios::out);
  mutout << "replicate\tgeneration\t";
  mutout << "mutator_freq\thave_multi_mutator\thave_multi_homo\thave_single_mutator\t";
  mutout << "AF1\tAF2\tAF3\tAF4\tAF5\tAF6\tAF7\tAF8\tAF9\tAF10\n";

  std::ofstream logout;
  logout.open("log/outlog"+std::string(argv[1])+".txt",std::ios::out);

  Constant nums;
  int replicate=1;
  while(replicate <=5000){
    Parameters param(nums);
    //int t=0;
    //while(param.expected_mutation_sd<0.0000005){param.reset(nums);t++;}
    //logout << t << " time reparam\t";
    param.set_damage(nums);
    Population population(nums,param);
    if(param.mutator_s>1 || param.mutation_rate>0.00000008 || param.non_damage_e>0.04){
      outfile << replicate << "\t0\t";
      print_out(param,population,outfile);
      logout << "done " << replicate << " time\n";logout.flush();
    }else{
      one_replicate(replicate, nums, param, population, outfile, mutout,logout);
      logout << "done " << replicate << " time\n";logout.flush();
    }
    replicate++;
  }

  return 0;
}
