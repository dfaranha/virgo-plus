#pragma once

#include <algorithm>
#include "RS_polynomial.h"
#include <vector>
#include "fri.h"
#include "vpd_prover.h"
#include <chrono>

using namespace std;

namespace virgo_ext {
		extern uint64_t __glb_c1, __glb_c2, __glb_c3, __glb_c4, __glb_c5;
    namespace poly_commit {
        extern fieldElement *all_pri_msk_arr;
        extern fieldElement *all_pub_msk_arr;
        extern fieldElement *inner_prod_evals;
        extern fieldElement *l_coef, *l_eval, *q_coef, *q_eval; //l for private input, q for public randomness
        extern fieldElement *lq_eval, *h_coef, *lq_coef, *h_eval;
        extern fieldElement *h_eval_arr;
        extern int l_coef_len, l_eval_len, q_coef_len, q_eval_len;
        extern int slice_size, slice_count, slice_real_ele_cnt;
        extern int mask_position_gap; //masks are positioned in specific way for efficiency
        extern bool pre_prepare_executed;

        class ldt_commitment {
			public:
				__hhash_digest *commitment_hhash;
				fieldElement *randomness;
				fieldElement *final_rs_code;
				int mx_depth;
				int repeat_no;
			};

        class poly_commit_prover {
			private:

			public:
				double total_time;

				std::vector<fieldElement> all_pri_mask;

				fieldElement inner_prod(const fieldElement *a, const fieldElement *b, u64 l);
				__hhash_digest commit_private_array(fieldElement *private_array, int log_array_length, std::vector<fieldElement> private_mask_array);
				__hhash_digest commit_public_array(std::vector<fieldElement> &all_pub_msk, fieldElement *public_array, int r_0_len, fieldElement target_sum, fieldElement *all_sum);
				ldt_commitment commit_phase(int log_length);
			};

			class poly_commit_verifier {
				public:
					poly_commit_prover *p;

					bool verify_poly_commitment(fieldElement *all_sum, int log_length, fieldElement *public_array,
												std::vector<fieldElement> &all_pub_mask, double &v_time, int &proof_size,
												double &p_time, __hhash_digest merkle_tree_l, __hhash_digest merkle_tree_h);
			};
    }

}