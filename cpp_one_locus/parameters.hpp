#pragma once
#ifndef PARAM_HPP
#define PARAM_HPP

#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <unordered_set>
#include <cmath>

struct Constant{
  const std::size_t N = 56418;
  const std::size_t patient_n =6418;
  const std::size_t tsg_non_site = 191624;
  const std::size_t tsg_syn_site = 61470;
  const double rare_tsg_non_num = 0.7642481;
  const double rare_tsg_non_sd = 0.9125167;
  const double rare_tsg_syn_num = 0.490959;
  const double rare_tsg_syn_sd = 0.7278306;
  std::mt19937 mt;
  std::bernoulli_distribution bern;
  Constant(){
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
  double mutater_effect;
  double mutater_mutation_rate;
  double mutater_damage;
  double tsg_non_damage_e;
  std::vector<double> tsg_non_damage;
  double expected_mutation_sd;
  Parameters(Constant& nums);
  void reset(Constant& nums);
  void set_damage(Constant& nums);
  int get_new_mutater(){new_mutater_id++;return(new_mutater_id);};

private:]
  int new_mutater_id = 0;
  double get_mutation_rate(Constant& nums);
  double get_mutater_effect(Constant& nums);
  double get_mutater_mutation_rate(Constant& nums);
  double get_mutater_damage(Constant& nums);
  double get_tsg_non_damage_e(Constant& nums);
  void set_tsg_non_damage(Constant& nums);
};

double log_random(double start, double end, std::mt19937& mt);
bool equib_lm(const std::vector<double>& mutation);

#endif
