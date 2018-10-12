#include "linear_model.hpp"

bool equiv_lm(const std::vector<double>& mutation) {
  bool focal=true;
  std::size_t vect_size = mutation.size();
  double xy=0, xx=0 , average=0;
  for(double i: mutation){average+=i;}
  average =(double) average/vect_size;
  double x_ave = (double)(1+vect_size)/2;
  for(std::size_t t=0; t < vect_size; t++){
    xx += (x_ave-t)*(x_ave-t); xy += (x_ave-t)*(average-mutation[t]);
  }
  if((double)average/5000 > xy/xx){focal=false;}
  return focal;
}

double linear_model(const std::vector<int>& rare_variant, const std::vector<double>& onset_age){
  double xx=0, xy=0;
  double x_=0, y_=0;
  std::size_t p_num =rare_variant.size();
  for(std::size_t t=0; t < p_num; t++){
    x_ += rare_variant[t]; y_ += onset_age[t];
  }
  x_ /= p_num; y_ /=p_num;
  for(std::size_t t=0; t < p_num; t++){
    xy += rare_variant[t] * (onset_age[t] - y_);
    xx += rare_variant[t] * (rare_variant[t] - x_);
  }
  return xy/xx;
}
