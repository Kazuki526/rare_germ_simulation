#pragma once
#ifndef POP_HPP
#define POP_HPP

#include"individual.hpp"

class Population : public Individual
{
private:
  std::vector<Individual> individuals;
public:
  Population(Parameters& param);
  std::vector<int> get_mutater_list();
  std::vector<double> add_new_mutations(Parameters& param); //return fitness_vector
  Individual reproduct(Parameters& param, size_t i1, size_t i2);
  void next_generation(Parameters& param, const std::vector<double>& fitness);
  void one_generation(Parameters& param);
  std::vector<int> tsg_non_mutation_count(Parameters& param);
};


#endif
