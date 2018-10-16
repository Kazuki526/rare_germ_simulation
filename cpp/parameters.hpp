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
  const double rare_tsg_non_num = 0.92474;
  const double rare_tsg_syn_num = 0.6611;
  const double rare_cont_non_num = 4.24805;
  const double tsg_non_onset_regression = -0.70911;
  const double tsg_syn_onset_regression = -0.69074;
  const double cont_non_onset_regression = -0.29612;
  std::mt19937 mt;
  Constant();
};

struct Parameters
{                                           // variable_param num
  double mutation_rate_coef;                // 0
  double mutater_effect;                    // 1
  double mutater_mutation_rate;             // 2
  double mutater_damage;                    // 3
  double tsg_non_damage_e;                  // 4
  double cont_non_damage_e;                 // 5
  std::vector<double> tsg_non_damage;
  std::vector<double> cont_non_damage;
  bool cont_non_fitness_only;               // 6
  double fitness_coef;                      // 7
  Parameters(Constant& nums);
  //void change_param(Constant& nums,const std::size_t time);

private:
  //std::size_t variable_param;
  //std::vector<double> variable_param_values;

  double get_mutation_rate_coef(Constant& nums);
  double get_mutater_effect(Constant& nums);
  double get_mutater_mutation_rate(Constant& nums);
  double get_mutater_damage(Constant& nums);
  double get_tsg_non_damage_e(Constant& nums);
  double get_cont_non_damage_e(Constant& nums);
  double get_cont_non_fitness_only(Constant& nums);
  double get_fitness_coef(Constant& nums);
  void set_tsg_non_damage(Constant& nums);
  void set_cont_non_damage(Constant& nums);
};

double log_random(double start, double end, std::mt19937 mt);


#endif
