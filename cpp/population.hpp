#pragma once
#ifndef POP_HPP
#define POP_HPP

#include"individual.hpp"

class Population
{
private:
  std::vector<Individual> individuals;
  std::vector<double> fitness;
  std::vector<std::size_t> num_tsg_non_mutation;

  Individual reproduct(Constant& nums, const Parameters& param, std::size_t i1, std::size_t i2);
  void mutation_count(const Constant& nums, const Parameters& param);
public:
  double rare_tsg_non_freq;

  Population(const Constant& nums, const Parameters& param);
  std::vector<int> get_mutater_list();
  void next_generation(Constant& nums, const Parameters& param);
};


#endif
