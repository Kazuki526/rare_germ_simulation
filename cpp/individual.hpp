#pragma once
#ifndef IND_HPP
#define IND_HPP

#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <unordered_set>

struct Constant{
  const std::size_t N = 50000;
  const std::size_t patient_n =6418;
  const double mutation_rate = 0.00000001;
  const std::size_t tsg_non_site = 191624;
  const std::size_t tsg_syn_site = 61470;
  const std::size_t cont_non_site = 923307;
  const double mean_onset_age = 61.5;
  const double onset_age_sd = 13.5;
  std::mt19937 mt;
  Constant(){
    std::random_device rnd;
    std::mt19937 mt_(rnd());
    mt = mt_;
  }
};

struct Parameters
{
  double mutation_rate_coef=10;
  double mutater_effect=100;
  double mutater_mutation_rate=0.001;
  double mutater_damage=0.3;
  double tsg_non_damage_e=0.5;
  double cont_non_damage_e=0.1;
  std::vector<double> tsg_non_damage;
  std::vector<double> cont_non_damage;
  bool cont_non_fitness_only=true;
  double fitness_coef=0.1;
  double cancer_prob_coef=1;
  void set_damage_mt(Constant& nums);
  void change_param(std::size_t p, double doub);
};

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
  double mut_r;
  double damage;
  double fitness;
public:
  Individual(): mutater(0){};
  Individual(const int m,
             const std::vector<std::size_t>& tsg_non,
             const std::vector<std::size_t>& tsg_syn,
             const std::vector<std::size_t>& cont_non,
             const std::unordered_set<std::size_t> tsg_non_common={},
             const std::unordered_set<std::size_t> tsg_syn_common={},
             const std::unordered_set<std::size_t> cont_non_common={});
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
  int gamate_mutater(Constant nums);
  std::vector<size_t> gamate_tsg_non(Constant nums);
  std::vector<size_t> gamate_tsg_syn(Constant nums);
  std::vector<size_t> gamate_cont_non(Constant nums);
};

#endif
