#pragma once
#ifndef IND_HPP
#define IND_HPP

#include "parameters.hpp"

class Individual
{
private:
  std::size_t mutator_num;
  std::vector<std::size_t> mutator_het;
  std::vector<std::size_t> mutator_hom;
  std::vector<std::size_t> tsg_non_het;
  std::vector<std::size_t> tsg_non_hom;
  std::vector<std::size_t> tsg_syn_het;
  std::vector<std::size_t> tsg_syn_hom;
  double mutation_r;
  double fitness;
public:
  Individual(){mutator_num=0;}
  Individual(const std::vector<std::size_t>& mutator,
             const std::vector<std::size_t>& tsg_non,
             const std::vector<std::size_t>& tsg_syn,
             const std::unordered_set<std::size_t>& tsg_non_common={},
             const std::unordered_set<std::size_t>& tsg_syn_common={});
  const std::size_t get_mutator_num(){return mutator_num;}
  const std::vector<std::size_t>& get_mutator_het(){return mutator_het;}
  const std::vector<std::size_t>& get_mutator_hom(){return mutator_hom;}
  const std::vector<std::size_t>& get_tsg_non_het(){return tsg_non_het;}
  const std::vector<std::size_t>& get_tsg_non_hom(){return tsg_non_hom;}
  const std::vector<std::size_t>& get_tsg_syn_het(){return tsg_syn_het;}
  const std::vector<std::size_t>& get_tsg_syn_hom(){return tsg_syn_hom;}
  double get_mutation_r(){return mutation_r;}
  double get_fitness(){return fitness;}
  void set_param(Constant& nums, const Parameters& param);
  /* gamate */
  std::vector<std::size_t> gamate_mutator(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_tsg_non(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_tsg_syn(Constant& nums, const Parameters& param);
};

#endif
