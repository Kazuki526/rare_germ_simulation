#include"individual.hpp"
#include"population.hpp"

Population::Population(Parameters& param){
  for(size_t i=0; i < param.N; i++){
    Individual ind;
    ind.set_param(param);
    individuals.push_back(ind);
  }
}

std::vector<double> Population::add_new_mutations(Parameters& param){
  std::vector<double> fitness;
  for(size_t i=0; i < param.N; i++){
    individuals[i].add_mutations(param);
    fitness.push_back(individuals[i].get_fitness());
  }
  return fitness;
}

// make offspring from two Individual
Individual Population::reproduct(Parameters& param, size_t i1, size_t i2){
  int mutater=0;
  std::vector<int> tsg_non;
  std::bernoulli_distribution bern(0.5);
  // individuals[i1] 's gamate
  int muta = individuals[i1].get_mutater();
  if(muta == 2){
    mutater++;
  }else if(muta == 1){
    if(bern(param.mt)){mutater++;}
  }
  for(int &het_mu: individuals[i1].get_tsg_non_het()){
    if(bern(param.mt)){tsg_non.push_back(het_mu);}
  }
  for(int &hom_mu: individuals[i1].get_tsg_non_hom()){
    tsg_non.push_back(hom_mu);
  }
  // individuals[i2] 's gamate
  muta = individuals[i2].get_mutater();
  if(muta == 2){
    mutater++;
  }else if(muta == 1){
    if(bern(param.mt)){mutater++;}
  }
  for(int &het_mu: individuals[i2].get_tsg_non_het()){
    if(bern(param.mt)){tsg_non.push_back(het_mu);}
  }
  for(int &hom_mu: individuals[i2].get_tsg_non_hom()){
    tsg_non.push_back(hom_mu);
  }
  Individual ind(mutater, tsg_non);
  ind.set_param(param);
  return(ind);
}

// reproduct N time
void Population::next_generation(Parameters& param, const std::vector<double>& fitness){
  std::vector<Individual> next_inds;
  std::discrete_distribution<std::size_t> dist(fitness.begin(), fitness.end());
  for(size_t i=0; i < param.N; i++){
     next_inds.push_back(reproduct(param, dist(param.mt), dist(param.mt)));
  }
  individuals = next_inds;
}

// add_new_mutations + next_generation
void Population::one_generation(Parameters& param){
  std::vector<double> fitness = add_new_mutations(param);
  next_generation(param, fitness);
}

// get all mutation AC
std::vector<int> Population::tsg_non_mutation_count(Parameters& param){
  std::vector<int> tsg_non_mutation(param.tsg_non_site);
  for(Individual &ind: individuals){
    for(int &het_mu: ind.get_tsg_non_het()){
      tsg_non_mutation[het_mu]++;
    }
    for(int &hom_mu: ind.get_tsg_non_hom()){
      tsg_non_mutation[hom_mu]+=2;
    }
  }
  return(tsg_non_mutation);
}
