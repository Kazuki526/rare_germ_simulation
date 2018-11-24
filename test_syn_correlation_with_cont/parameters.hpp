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
  const std::size_t tsg_non_site = 191624;
  const std::size_t tsg_syn_site = 61470;
  const std::size_t syn_site = 319944+61470;
  const std::size_t cont_non_site = 923307;
  const double rare_tsg_non_num = 0.92474;
  const double rare_tsg_non_sd = 1.004913;
  const double rare_tsg_syn_num = 0.6611;
  const double rare_tsg_syn_sd = 0.8668836;
  const double rare_syn_num = 1.585385+0.6611;
  const double rare_cont_non_num = 4.24805;
  const double tsg_non_tsg_syn_cor = 0.119441;
  const double tsg_non_syn_cor = 0.195513;
  const double tsg_non_cont_non_cor = 0.18858;
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
  double cont_non_damage_e;
  std::vector<double> tsg_non_damage;
  std::vector<double> cont_non_damage;
  double expected_mutation_sd;
  Parameters(Constant& nums);
  void reset(Constant& nums);
  void set_damage(Constant& nums);

private:
  double get_mutation_rate(Constant& nums);
  double get_mutater_effect(Constant& nums);
  double get_mutater_mutation_rate(Constant& nums);
  double get_mutater_damage(Constant& nums);
  double get_tsg_non_damage_e(Constant& nums);
  double get_cont_non_damage_e(Constant& nums);
};

double log_random(double start, double end, std::mt19937& mt);
bool equiv_lm(const std::vector<double>& mutation);

#endif
