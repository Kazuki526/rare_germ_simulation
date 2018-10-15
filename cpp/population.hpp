#pragma once
#ifndef POP_HPP
#define POP_HPP

#include <numeric>
#include "individual.hpp"
#include "linear_model.hpp"

class Population
{
private:
  std::vector<Individual> individuals;
  std::vector<double> fitness;
  std::vector<std::size_t> num_tsg_non_mutation;
  std::vector<std::size_t> num_tsg_syn_mutation;
  std::vector<std::size_t> num_cont_non_mutation;
  std::unordered_set<std::size_t> tsg_non_common;
  std::unordered_set<std::size_t> tsg_syn_common;
  std::unordered_set<std::size_t> cont_non_common;

  Individual reproduct(Constant& nums, const Parameters& param, std::size_t i1, std::size_t i2, bool common_check);
  void mutation_count(const Constant& nums, const Parameters& param);
public:
  double rare_tsg_non_freq;
  double rare_tsg_syn_freq;
  double rare_cont_non_freq;
  double tsg_non_regression;
  double tsg_syn_regression;
  double cont_non_regression;

  Population(const Constant& nums);
  void set_params(const Constant& nums, const Parameters& param);
  void set_common_variant(const Constant& nums);
  void next_generation(Constant& nums, const Parameters& param, bool common_variant=false);
  void regression_onset_age(Constant& nums, const Parameters& param);
};

bool accept_reject_judge(const Constant& nums, const Population population);

#endif
