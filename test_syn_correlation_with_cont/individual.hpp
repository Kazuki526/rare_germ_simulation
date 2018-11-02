#pragma once
#ifndef IND_HPP
#define IND_HPP

#include "parameters.hpp"

class Individual
{
private:
  std::size_t mutater;
  std::vector<std::size_t> tsg_non_het;
  std::vector<std::size_t> tsg_non_hom;
  std::vector<std::size_t> tsg_syn_het;
  std::vector<std::size_t> tsg_syn_hom;
  std::vector<std::size_t> cont_non_het;
  std::vector<std::size_t> cont_non_hom;

  double fitness;
public:
  double mut_r;
  Individual(): mutater(0){};
  Individual(const std::size_t& m,
             const std::vector<std::size_t>& tsg_non,
             const std::vector<std::size_t>& tsg_syn,
             const std::vector<std::size_t>& cont_non,
             const std::unordered_set<std::size_t>& tsg_non_common={},
             const std::unordered_set<std::size_t>& tsg_syn_common={},
             const std::unordered_set<std::size_t>& cont_non_common={});
  const std::size_t& get_mutater(){return mutater;}
  const std::vector<std::size_t>& get_tsg_non_het(){return tsg_non_het;}
  const std::vector<std::size_t>& get_tsg_non_hom(){return tsg_non_hom;}
  const std::vector<std::size_t>& get_tsg_syn_het(){return tsg_syn_het;}
  const std::vector<std::size_t>& get_tsg_syn_hom(){return tsg_syn_hom;}
  const std::vector<std::size_t>& get_cont_non_het(){return cont_non_het;}
  const std::vector<std::size_t>& get_cont_non_hom(){return cont_non_hom;}
  double get_mut_r(){return mut_r;}
  double get_fitness(){return fitness;}
  void set_param(const Constant& nums, const Parameters& param);
  /* gamate */
  std::size_t gamate_mutater(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_tsg_non(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_tsg_syn(Constant& nums, const Parameters& param);
  std::vector<std::size_t> gamate_cont_non(Constant& nums, const Parameters& param);
};

#endif
