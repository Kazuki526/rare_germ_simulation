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
      population.mutation_count(nums, param);
      logout << generation << "> ";logout.flush();
    }else{
        population.next_generation(nums, param);
        population.mutation_count(nums, param);
    }

    if(generation<=100){
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
    if((population.rare_non_freq > nums.rare_non_num*2)||
       (population.rare_syn_freq > nums.rare_syn_num*2)){logout<<"stop";logout.flush();break;}
//    if((generation>=500)&&
//       ((population.rare_non_freq < nums.rare_non_num*0.5)||
//        (population.rare_syn_freq < nums.rare_syn_num*0.5))){logout<<"stop";logout.flush();break;}
  }

  logout << generation << "=> summary >";logout.flush();
  if(!mut_equiv && !non_equiv && !syn_equiv){ /* not stoped by too much rare variant */
    if(generation<500){
      while(generation==500){
        generation++;
        if(generation%50==0){
          population.next_generation(nums, param, true);
        }else{
          population.next_generation(nums, param);
        }
      }
    }else{
      for(int gene=1; gene<=50; gene++){
        generation++;
        population.next_generation(nums, param);
      }
    }
    logout <<" >";logout.flush();
    population.mutation_count(nums, param);
    logout <<" >";logout.flush();
    population.correlation_ns(nums);
    outfile << replicate << "\t";
    outfile << generation <<"\t";outfile.flush();
    print_out(param,population,outfile);
    mutout << replicate << "\t" << generation << "\t";mutout.flush();
    population.out_mutator_state(nums, mutout);
  }else{ /* when stoped by too much rare variant */
    logout <<" >";logout.flush();
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
  mutout << "replicate\tgeneration\tnum_id\t";
  mutout << "mutator_freq\thave_multi_mutator\thave_multi_homo\thave_single_mutator\t";
  mutout << "AF1\tAF2\tAF3\tAF4\tAF5\tAF6\tAF7\tAF8\tAF9\tAF10\n";

  std::ofstream logout;
  logout.open("log/outlog"+std::string(argv[1])+".txt",std::ios::out);

  Constant nums;
  int replicate=1;
  while(replicate <=1000000){
    nums.new_mutator_id=0;
    Parameters param(nums);
    param.set_damage(nums);
    Population population(nums,param);
    //if(param.mutator_s>1 || param.mutation_rate>0.00000008 || param.non_damage_e>0.04){
    if(param.mutator_s>1 || param.mutation_rate>0.00000008){
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
