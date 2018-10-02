#include"individual.hpp"
#include"population.hpp"

Population::Population(const Constant& nums, const Parameters& param){
  for(std::size_t i=0; i < nums.N; i++){
    Individual ind;
    ind.set_param(nums, param);
    individuals.push_back(ind);
  }
}

/* make offspring from two Individual */
Individual Population::reproduct(Constant& nums, const Parameters& param, std::size_t i1, std::size_t i2){
  int mutater =
      individuals[i1].gamate_mutater(nums) +
      individuals[i2].gamate_mutater(nums);
  std::vector<std::size_t> tsg_non=individuals[i1].gamate_tsg_non(nums);
  std::vector<std::size_t> tn2=individuals[i2].gamate_tsg_non(nums);
  tsg_non.insert(tsg_non.end(), tn2.begin(), tn2.end());
  Individual ind(mutater, tsg_non);
  ind.set_param(nums, param);
  return(ind);
}

/* get all mutation AC */
void Population::mutation_count(const Constant& nums, const Parameters& param){
  num_tsg_non_mutation.clear();
  num_tsg_non_mutation.resize(nums.tsg_non_site);
  for(Individual &ind: individuals){
    for(std::size_t het_mu: ind.get_tsg_non_het()){
      num_tsg_non_mutation[het_mu]++;
    }
    for(std::size_t hom_mu: ind.get_tsg_non_hom()){
      num_tsg_non_mutation[hom_mu]+=2;
    }
  }
  rare_tsg_non_freq=0;
  for(std::size_t mu=0; mu < nums.tsg_non_site; mu++){
    if(num_tsg_non_mutation[mu] < nums.N*0.05/100){
      rare_tsg_non_freq += (double) num_tsg_non_mutation[mu] /nums.N;
    }
  }
}

void Population::next_generation(Constant& nums, const Parameters& param){
  /* add_mutations & set fitness */
  fitness.clear(); fitness.reserve(nums.N);
  for(std::size_t i=0; i < nums.N; i++){
    individuals[i].add_mutations(nums, param);
    //fitness.push_back(individuals[i].get_fitness());
    fitness.push_back((1- individuals[i].get_damage() * param.fitness_coef));
  }
  /* reproduct N individual */
  std::vector<Individual> next_inds;
  next_inds.reserve(nums.N);
  std::discrete_distribution<std::size_t> dist(fitness.begin(), fitness.end());
  for(std::size_t i=0; i < nums.N; i++){
     next_inds.push_back(reproduct(nums, param, dist(nums.mt), dist(nums.mt)));
  }
  individuals = next_inds;
  /* mutation_count */
  mutation_count(nums, param);
}
