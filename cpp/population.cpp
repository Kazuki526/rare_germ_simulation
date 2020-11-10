#include "population.hpp"

Population::Population(Constant& nums, const Parameters& param){
  individuals.reserve(nums.N);
  for(std::size_t i=0; i < nums.N; i++){
    Individual ind;
    ind.set_param(nums,param);
    individuals.push_back(ind);
  }
}

void Population::set_common_variant(const Constant& nums){
  double common_num = nums.N *2*0.01;
  for(std::size_t m=0; m < nums.non_site; m++){
    if(num_non_mutation[m] > common_num){non_common.insert(m);}
  }
  for(std::size_t m=0; m < nums.syn_site; m++){
    if(num_syn_mutation[m] > common_num){syn_common.insert(m);}
  }
}

/* make offspring from two Individual */
Individual Population::reproduct(Constant& nums,const Parameters& param, std::size_t i1, std::size_t i2, bool common_check){
/* mutator */
  std::vector<std::size_t> mutator=individuals[i1].gamate_mutator(nums,param);
  std::vector<std::size_t> mut2=individuals[i2].gamate_mutator(nums,param);
  mutator.insert(mutator.end(), mut2.begin(), mut2.end());
/* nonsynonymous */
  std::vector<std::size_t> non=individuals[i1].gamate_non(nums,param);
  std::vector<std::size_t> non2=individuals[i2].gamate_non(nums,param);
  non.insert(non.end(), non2.begin(), non2.end());
/* synonymous */
  std::vector<std::size_t> syn=individuals[i1].gamate_syn(nums,param);
  std::vector<std::size_t> syn2=individuals[i2].gamate_syn(nums,param);
  syn.insert(syn.end(), syn2.begin(), syn2.end());
  if(common_check){
    Individual ind(mutator, non, syn, non_common, syn_common);
    ind.set_param(nums, param);
    return(ind);
  }else{
    Individual ind(mutator, non, syn);
    ind.set_param(nums, param);
    return(ind);
  }
}

/* get all mutation AC */
void Population::mutation_count(const Constant& nums, const Parameters& param){
  const double rare =nums.N*2.0*0.05*0.01;
  std::vector<std::size_t> mutator_num;
  mutator_num.reserve(nums.N);
  num_non_mutation.clear();
  num_non_mutation.resize(nums.non_site);
  num_syn_mutation.clear();
  num_syn_mutation.resize(nums.syn_site);
  for(Individual &ind: individuals){
/* mutator */
    mutator_num.push_back(ind.get_mutator_num());
/* nonsynonymous */
    for(std::size_t het_mu: ind.get_non_het()){
      num_non_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_non_hom()){
      num_non_mutation[hom_mu]+=2;
    }
/* synonymous */
    for(std::size_t het_mu: ind.get_syn_het()){
      num_syn_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_syn_hom()){
      num_syn_mutation[hom_mu]+=2;
    }
  }
/* count each freq */
  mutator_freq = (double)std::accumulate(mutator_num.begin(),mutator_num.end(),0)/nums.N;
  rare_non_freq=0;
  for(std::size_t mu=0; mu < nums.non_site; mu++){
    if(num_non_mutation[mu] < rare){
      rare_non_freq += (double) num_non_mutation[mu] /nums.N;
    }
  }
  rare_syn_freq=0;
  for(std::size_t mu=0; mu < nums.syn_site; mu++){
    if(num_syn_mutation[mu] < rare){
      rare_syn_freq += (double) num_syn_mutation[mu] /nums.N;
    }
  }
}

void Population::next_generation(Constant& nums, const Parameters& param, bool common_check){
  if(common_check){set_common_variant(nums);}
  /* add_mutations & set fitness */
  fitness.clear(); fitness.reserve(nums.N);
  for(std::size_t i=0; i < nums.N; i++){
    fitness.push_back(individuals[i].get_fitness());
  }
  /* reproduct N individual */
  std::vector<Individual> next_inds;
  next_inds.reserve(nums.N);
  std::discrete_distribution<std::size_t> dist(fitness.begin(), fitness.end());
  for(std::size_t i=0; i < nums.N; i++){
     next_inds.push_back(reproduct(nums, param, dist(nums.mt), dist(nums.mt), common_check));
  }
  individuals = next_inds;
  if(common_check){
    non_common.clear();
    syn_common.clear();
  }
}

void Population::correlation_ns(const Constant& nums){
  const int N =nums.N;
  const int sample_n =nums.sample_n;
  const double rare = (N-sample_n)*2*0.05*0.01;
/* select sample */
/* individual is randomly arranged */
/*0 ~ sample_n-1 => sample */

/* count sample rare num */
  std::vector<int> non_rare_num(sample_n,0), syn_rare_num(sample_n,0);
  std::vector<double> mutation_rate; mutation_rate.reserve(sample_n);
  for(std::size_t i=0; i <sample_n; i++){
    for(std::size_t non_het: individuals[i].get_non_het()){
      if(num_non_mutation[non_het] <= rare){non_rare_num[i]++;}
    }
    for(std::size_t non_hom: individuals[i].get_non_hom()){
      if(num_non_mutation[non_hom] <= rare){non_rare_num[i]+=2;}
    }
    for(std::size_t syn_het: individuals[i].get_syn_het()){
      if(num_syn_mutation[syn_het] <= rare){syn_rare_num[i]++;}
    }
    for(std::size_t syn_hom: individuals[i].get_syn_hom()){
      if(num_syn_mutation[syn_hom] <= rare){syn_rare_num[i]+=2;}
    }
    mutation_rate.push_back(individuals[i].get_mutation_r());
  }

  rare_non_freq = (double)std::accumulate(non_rare_num.begin(),non_rare_num.end(),0.0)/sample_n;
  rare_syn_freq = (double)std::accumulate(syn_rare_num.begin(),syn_rare_num.end(),0.0)/sample_n;
  mutation_rate_ave = (double)std::accumulate(mutation_rate.begin(),mutation_rate.end(),0.0)/sample_n;
  double nonv=0.0, synv=0.0, nonsynv=0.0, mutv=0.0, mutrv=0.0;
  for(std::size_t i=0; i < sample_n; i++){
    nonv+=(double)((non_rare_num[i]-rare_non_freq)*(non_rare_num[i]-rare_non_freq));
    synv+=(double)((syn_rare_num[i]-rare_syn_freq)*(syn_rare_num[i]-rare_syn_freq));
    nonsynv+=(double)((non_rare_num[i]-rare_non_freq)*(syn_rare_num[i]-rare_syn_freq));
    mutrv+=(double)((mutation_rate[i]-mutation_rate_ave)*(mutation_rate[i]-mutation_rate_ave));

  }
  rare_non_sd = std::sqrt((double)nonv/sample_n);
  rare_syn_sd = std::sqrt((double)synv/sample_n);
  rare_non_syn_correlation = (double)(nonsynv/sample_n)/(rare_non_sd*rare_syn_sd);
  rare_num_reg = (double) nonsynv/nonv;
  mutation_rate_sd = std::sqrt((double)mutrv/sample_n);
}

void Population::out_mutator_state(const Constant& nums, std::ofstream& mutout){
  mutout << nums.new_mutator_id; mutout.flush();
  std::vector<int> mutator_ac(nums.new_mutator_id,0);
  mutout <<"\t"; mutout.flush();
  double have_multi_mutator=0,have_multi_homo=0, have_single_mutator=0;
  for(Individual &ind: individuals){
    if(ind.get_mutator_num() >0){
      if(ind.get_mutator_num() ==1){have_single_mutator++;}else{have_multi_mutator++;}
      for(std::size_t het: ind.get_mutator_het()){
        mutator_ac[het-1]++;
      }
      for(std::size_t hom: ind.get_mutator_hom()){
        mutator_ac[hom-1]+=2;
        have_multi_homo++;
      }
    }
  }
  std::sort(mutator_ac.begin(),mutator_ac.end(),std::greater<int>());
  mutout << mutator_freq << "\t";
  mutout << have_multi_mutator/nums.N << "\t";
  mutout << have_multi_homo/nums.N << "\t"<< have_single_mutator/nums.N;
  for(std::size_t i=0; i<10; i++){
    mutout << "\t" <<(double)mutator_ac[i]/nums.N;
  }
  mutout << "\n";
  mutout.flush();
}
