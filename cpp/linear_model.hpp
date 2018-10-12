#pragma once
#ifndef LM_HPP
#define LM_HPP

#include <iostream>
#include <vector>


bool equiv_lm(const std::vector<double>& mutation);

double linear_model(const std::vector<int>& rare_variant, const std::vector<double>& onset_age);

#endif
