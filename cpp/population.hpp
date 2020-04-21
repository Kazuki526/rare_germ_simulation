#pragma once
#ifndef POP_HPP
#define POP_HPP

#include <numeric>
#include <string>
#include <functional>
#include <fstream>
#include <stdio.h>
#include "individual.hpp"

class Population
{
private:
  std::vector<Individual> individuals;
  std::vector<double> fitness;
  std::vector<std::size_t> num_non_mutation;
  std::vector<std::size_t> num_syn_mutation;
  std::unordered_set<std::size_t> non_common;
  std::unordered_set<std::size_t> syn_common;

  Individual reproduct(Constant& nums, const Parameters& param, std::size_t i1, std::size_t i2, bool common_check);
  void mutation_count(const Constant& nums, const Parameters& param);
public:
  double mutation_rate_ave;
  double mutation_rate_sd;
  double mutator_freq;
  //double mutator_sd;
  double rare_non_freq;
  double rare_syn_freq;
  double rare_non_sd;
  double rare_syn_sd;
  double rare_non_syn_correlation;
  double rare_num_reg;

  Population(Constant& nums, const Parameters& param);
  void set_common_variant(const Constant& nums);
  void next_generation(Constant& nums, const Parameters& param, bool common_variant=false);
  void correlation_ns(const Constant& nums);
  void out_mutator_state(const Constant& nums, std::ofstream& mutout);
};


#endif