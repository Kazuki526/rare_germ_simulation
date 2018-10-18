#pragma once
#ifndef IND_HPP
#define IND_HPP

#include "parameters.hpp"

class Individual
{
private:
  int mutater;
  std::vector<std::size_t> tsg_non_het;
  std::vector<std::size_t> tsg_non_hom;
  std::vector<std::size_t> tsg_syn_het;
  std::vector<std::size_t> tsg_syn_hom;
  std::vector<std::size_t> cont_non_het;
  std::vector<std::size_t> cont_non_hom;

  double damage;
  double fitness;
public:
  double mut_r;
  Individual(): mutater(0){};
  Individual(const int m,
             const std::vector<std::size_t>& tsg_non,
             const std::vector<std::size_t>& tsg_syn,
             const std::vector<std::size_t>& cont_non,
             const std::unordered_set<std::size_t>& tsg_non_common={},
             const std::unordered_set<std::size_t>& tsg_syn_common={},
             const std::unordered_set<std::size_t>& cont_non_common={});
  int get_mutater(){return mutater;}
  const std::vector<std::size_t>& get_tsg_non_het(){return tsg_non_het;}
  const std::vector<std::size_t>& get_tsg_non_hom(){return tsg_non_hom;}
  const std::vector<std::size_t>& get_tsg_syn_het(){return tsg_syn_het;}
  const std::vector<std::size_t>& get_tsg_syn_hom(){return tsg_syn_hom;}
  const std::vector<std::size_t>& get_cont_non_het(){return cont_non_het;}
  const std::vector<std::size_t>& get_cont_non_hom(){return cont_non_hom;}
  double get_damage(){return damage;}
  double get_fitness(){return fitness;}
  void set_param(const Constant& nums, const Parameters& param);
  void add_mutations(Constant& nums, const Parameters& param);
  /* gamate */
  int gamate_mutater(Constant& nums);
  std::vector<std::size_t> gamate_tsg_non(Constant& nums);
  std::vector<std::size_t> gamate_tsg_syn(Constant& nums);
  std::vector<std::size_t> gamate_cont_non(Constant& nums);
};

#endif
