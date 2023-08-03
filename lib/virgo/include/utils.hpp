//
// Created by 69029 on 6/23/2021.
//

#ifndef UNLAYERED_LIBRA_UTILS_HPP
#define UNLAYERED_LIBRA_UTILS_HPP

#include "config_pc.hpp"

#define MIN(x, y)	((x) > (y) ? (y) : (x))

#define MAX(x, y)	((x) > (y) ? (x) : (y))

int mylog(long long x);

void initBetaTable(vector<F> &beta_g, u8 gLength, const vector<F>::const_iterator &r, const F &init);

void
initBetaTable(vector<F> &beta_g, int gLength, const vector<F>::const_iterator &r_0,
              const vector<F>::const_iterator &r_1, const F &alpha, const F &beta);

template <class T>
void myResize(vector<T> &vec, u64 sz) {
    if (vec.size() < sz) vec.resize(sz);
}
#endif //UNLAYERED_LIBRA_UTILS_HPP
