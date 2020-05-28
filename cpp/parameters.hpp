#pragma once
#ifndef PARAM_HPP
#define PARAM_HPP

#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <unordered_set>
#include <cmath>

// non = nonsynonymous, syn = synonymous; real site num = site_num/3
struct Constant{
  const std::size_t N = 50000;
  const std::size_t sample_n = 503;
  // tsg_non_site = , control_non_site = 1270003
  const double non_site = 1270003;
  // tsg_syn_site = , control_syn_site = 414929
  const double syn_site = 414929;
  // rare_tsg_non_num = , rare_cg_non_num = 2.8051689860835
  const double rare_non_num = 2.8051689860835;
  //rare_tsg_non_sd = , rare_cg_non_sd = 1.90663893622778
  const double rare_non_sd = 1.90663893622778;
  //rare_tsg_syn_num = , rare_cg_syn_num = 1.54870775347913
  const double rare_syn_num = 1.54870775347913;
  //rare_tsg_syn_sd = , rare_cg_syn_sd = 1.38039536172101
  const double rare_syn_sd = 1.38039536172101;
  std::size_t new_mutator_id;
  std::size_t get_new_mutator(){new_mutator_id++; return new_mutator_id;}

  std::mt19937 mt;
  std::bernoulli_distribution bern;
  Constant(){
    new_mutator_id = 0;
    std::random_device rnd;
    std::mt19937 mt_(rnd());
    mt = mt_;
    std::bernoulli_distribution bern_(0.5);
    bern = bern_;
  };
};

struct Parameters
{
  double mutation_rate;
  //double repair_power=10000;
  double mutator_effect;
  double mutator_mutation_rate;
  double mutator_damage;
  double non_damage_e;
  double mutator_s;
  std::vector<double> non_damage;
  double expected_mutation_sd;
  Parameters(Constant& nums);
  void reset(Constant& nums);
  void set_damage(Constant& nums);

private:
  double set_mutation_rate(Constant& nums);
  double set_mutator_effect(Constant& nums);
  double set_mutator_mutation_rate(Constant& nums);
  double set_mutator_damage(Constant& nums);
  double set_non_damage_e(Constant& nums);
  void set_non_damage(Constant& nums);
};

double log_random(double start, double end, std::mt19937& mt);
bool equib_lm(const std::vector<double>& mutation);

#endif
