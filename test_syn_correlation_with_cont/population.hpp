#pragma once
#ifndef POP_HPP
#define POP_HPP

#include <numeric>
#include "individual.hpp"

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
  double mutation_rate_ave;
  double mutation_rate_sd;
  double mutater_freq;
  double mutater_sd;
  double rare_tsg_non_freq;
  double rare_tsg_syn_freq;
  double rare_cont_non_freq;
  double rare_tsg_non_sd;
  double rare_tsg_syn_sd;
  double rare_cont_non_sd;
  double mut0_rare_non_num;
  double mut1_rare_non_num;
  double mut2_rare_non_num;
  double mut0_notrare_non_num;
  double mut1_notrare_non_num;
  double mut2_notrare_non_num;
  double mut0_rare_syn_num;
  double mut1_rare_syn_num;
  double mut2_rare_syn_num;
  double mut0_notrare_syn_num;
  double mut1_notrare_syn_num;
  double mut2_notrare_syn_num;
  double mut0_rare_cont_num;
  double mut1_rare_cont_num;
  double mut2_rare_cont_num;
  double mut0_notrare_cont_num;
  double mut1_notrare_cont_num;
  double mut2_notrare_cont_num;
  double rare_non_syn_correlation;
  double rare_non_cont_correlation;

  Population(const Constant& nums, const Parameters& param);
  void set_common_variant(const Constant& nums);
  void next_generation(Constant& nums, const Parameters& param, bool common_variant=false);
  void correlation_ns();
};


#endif