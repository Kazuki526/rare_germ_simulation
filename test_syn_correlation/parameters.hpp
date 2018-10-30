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
  const std::size_t N = 50000;
  const std::size_t patient_n =6418;
  const double mutation_rate = 0.00000003;
  const std::size_t tsg_non_site = 191624;
  const std::size_t tsg_syn_site = 61470;
  const double rare_tsg_non_num = 0.92474;
  const double rare_tsg_syn_num = 0.6611;
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
  double mutater_effect;
  double mutater_mutation_rate;
  double mutater_damage;
  double tsg_non_damage_e;
  std::vector<double> tsg_non_damage;
  double expected_mutater_freq;
  Parameters(Constant& nums);
  void reset(Constant& nums);
  void set_damage(Constant& nums);

private:
  double get_mutater_effect(Constant& nums);
  double get_mutater_mutation_rate(Constant& nums);
  double get_mutater_damage(Constant& nums);
  double get_tsg_non_damage_e(Constant& nums);
  void set_tsg_non_damage(Constant& nums);
};

double log_random(double start, double end, std::mt19937& mt);
bool equiv_lm(const std::vector<double>& mutation);

#endif
