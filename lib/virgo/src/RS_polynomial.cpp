#include "RS_polynomial.h"
#include <iostream>

using namespace std;

namespace virgo {
    fieldElement *__dst[3];
    fieldElementPacked *__avx_dst[2];
    fieldElement *twiddle_factor;

    void init_scratch_pad(int order) {
        __dst[0] = new fieldElement[order];//(fieldElement*)malloc(order * sizeof(fieldElement));
        __dst[1] = new fieldElement[order];//(fieldElement*)malloc(order * sizeof(fieldElement));
        __dst[2] = new fieldElement[order];
        __avx_dst[0] = new fieldElementPacked[order / packed_size];
        __avx_dst[1] = new fieldElementPacked[order / packed_size];
        twiddle_factor = new fieldElement[order];
    }

    void delete_scratch_pad() {
        delete[] __dst[0];
        delete[] __dst[1];
        delete[] __avx_dst[0];
        delete[] __avx_dst[1];
        delete[] twiddle_factor;
    }

    void fast_fourier_transform(const fieldElement *coefficients, int coef_len, int order, fieldElement root_of_unity,
                                fieldElement *result) {
        fieldElement rot_mul[fieldElement::maxOrder()];
        //note: malloc and free will not call the constructor and destructor, not recommended unless for efficiency
        //assert(sizeof(fieldElement) * 2 == sizeof(__hhash_digest));
        //In sake of both memory and time efficiency, use the non-recursive version
        int lg_order = -1;
        rot_mul[0] = root_of_unity;
        for (int i = 0; i < fieldElement::maxOrder(); ++i) {
            if (i > 0)
                rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
            if ((1LL << i) == order) {
                lg_order = i;
            }
        }
        int lg_coef = -1;
        for (int i = 0; i < fieldElement::maxOrder(); ++i) {
            if ((1LL << i) == coef_len) {
                lg_coef = i;
            }
        }
        assert(lg_order != -1 && lg_coef != -1);
        assert(rot_mul[lg_order] == 1);
        //we can merge both cases, but I just don't want to do so since it's easy to make mistake
        if (lg_coef > lg_order) {
            assert(false);
        } else {
            //initialize leaves
            int blk_sz = (order / coef_len);
            for (int j = 0; j < blk_sz; ++j) {
                for (int i = 0; i < coef_len; ++i) {
                    __dst[lg_coef & 1][(j << lg_coef) | i] = coefficients[i];
                }
            }
            fieldElement *x_arr = new fieldElement[1 << lg_order];
            {
				for (int dep = lg_coef - 1; dep >= 0; --dep) {
					int blk_size = 1 << (lg_order - dep);
					int half_blk_size = blk_size >> 1;
					int cur = dep & 1;
					int pre = cur ^ 1;

					x_arr[0] = fieldElement(1);
					for (int j = 1; j < blk_size; ++j)
						x_arr[j] = x_arr[j - 1] * rot_mul[dep];
					for (int k = 0; k < blk_size / 2; ++k) {
						int double_k = (k) & (half_blk_size - 1);
						for (int j = 0; j < (1 << dep); ++j) {
							auto l_value = __dst[pre][double_k << (dep + 1) | j], r_value =
									x_arr[k] * __dst[pre][double_k << (dep + 1) | (1 << dep) | j];
							__dst[cur][k << dep | j] = l_value + r_value;
							__dst[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
						}
					}
				}
            }
            delete[] x_arr;
        }

        for (int i = 0; i < order; ++i)
            result[i] = __dst[0][i];
    }

    void inverse_fast_fourier_transform(fieldElement *evaluations, int coef_len, int order, fieldElement root_of_unity,
                                        fieldElement *dst) {
        if (coef_len > order) {
            //more coefficient than evaluation
            fprintf(stderr, "Warning, Request do inverse fft with inefficient number of evaluations.");
            fprintf(stderr, "Will construct a polynomial with less order than required.");
            coef_len = order;
        }
        //assume coef_len <= order

        //subsample evalutions

        fieldElement *sub_eval;
        bool need_free = false;
        if (coef_len != order) {
            need_free = true;
            sub_eval = (fieldElement *) malloc(coef_len * sizeof(fieldElement));
            for (int i = 0; i < coef_len; ++i) {
                sub_eval[i] = evaluations[i * (order / coef_len)];
            }
        } else
            sub_eval = evaluations;

        fieldElement new_rou = fieldElement(1);
        for (int i = 0; i < order / coef_len; ++i)
            new_rou = new_rou * root_of_unity;
        order = coef_len;

        fieldElement inv_rou = fieldElement(1), tmp = new_rou;
        int lg_order = -1;
        for (int i = 0; i < fieldElement::maxOrder(); ++i) {
            if ((1LL << i) == order) {
                lg_order = i;
            }
        }
        int lg_coef = -1;
        for (int i = 0; i < fieldElement::maxOrder(); ++i) {
            if ((1LL << i) == coef_len) {
                lg_coef = i;
            }
        }
        assert(lg_order != -1 && lg_coef != -1);

        for (int i = 0; i < lg_order; ++i) {
            inv_rou = inv_rou * tmp;
            tmp = tmp * tmp;
        }
        assert(inv_rou * new_rou == fieldElement(1));

        fast_fourier_transform(sub_eval, order, coef_len, inv_rou, dst);

        if (need_free)
            free(sub_eval);

        fieldElement inv_n = fieldElement(order).inv();
        assert(inv_n * fieldElement(order) == fieldElement(1));

        for (int i = 0; i < coef_len; ++i) {
            dst[i] = dst[i] * inv_n;
        }
    }
}