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
  std::vector<std::size_t> non_het;
  std::vector<std::size_t> non_hom;
  std::vector<std::size_t> syn_het;
  std::vector<std::size_t> syn_hom;
  double mutation_r;
  double fitness;
public:
  Individual(){mutator_num=0;}
  Individual(const std::vector<std::size_t>& mutator,
             const std::vector<std::size_t>& non,
             const std::vector<std::size_t>& syn,
             const std::unordered_set<std::size_t>& non_common={},
             const std::unordered_set<std::size_t>& syn_common={});
  const std::size_t get_mutator_num(){return mutator_num;}
  const std::vector<std::size_t>& get_mutator_het(){return mutator_het;}
  const std::vector<std::size_t>& get_mutator_hom(){return mutator_hom;}
  const std::vector<std::size_t>& get_non_het(){return non_het;}
  const std::vector<std::size_t>& get_non_hom(){return non_hom;}
  const std::vector<std::size_t>& get_syn_het(){return syn_het;}
  const std::vector<std::size_t>& get_syn_hom(){return syn_hom;}
  double get_mutation_r(){return mutation_r;}
  double get_fitness(){return fitness;}
  void set_param(Constant& nums, const Parameters& param);
  /* gamate */
  std::vector<std::size_t> gamate_mutator(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_non(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_syn(Constant& nums, const Parameters& param);
};

#endif
