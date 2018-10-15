#pragma once
#ifndef PARAM_HPP
#define PARAM_HPP

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

#endif
