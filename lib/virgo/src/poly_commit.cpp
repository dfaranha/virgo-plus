#include "poly_commit.h"
#include "utils.hpp"

namespace virgo_ext {
	fieldElement *poly_commit::all_pri_msk_arr;
	fieldElement *poly_commit::all_pub_msk_arr;
	fieldElement *poly_commit::inner_prod_evals;
	fieldElement *poly_commit::l_coef, *poly_commit::l_eval, *poly_commit::q_coef, *poly_commit::q_eval; //l for private input, q for public randomness
	fieldElement *poly_commit::lq_eval, *poly_commit::h_coef, *poly_commit::lq_coef, *poly_commit::h_eval;
	fieldElement *poly_commit::h_eval_arr;
	int poly_commit::l_coef_len, poly_commit::l_eval_len, poly_commit::q_coef_len, poly_commit::q_eval_len;
	int poly_commit::slice_size, poly_commit::slice_count, poly_commit::slice_real_ele_cnt;
	int poly_commit::mask_position_gap; //masks are positioned in specific way for efficiency
	bool poly_commit::pre_prepare_executed = false;

	namespace poly_commit {
		fieldElement poly_commit_prover::inner_prod(const fieldElement *a, const fieldElement *b, u64 l) {
			fieldElement ret = 0;
			for (int i = 0; i < l; ++i)
				ret = ret + a[i] * b[i];
			return ret;
		}

		__hhash_digest poly_commit_prover::commit_private_array(fieldElement *private_array, int log_array_length, std::vector<fieldElement> private_mask_array) {
			total_time = 0;
			std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
			pre_prepare_executed = true;
			slice_count = 1 << log_slice_number;
			slice_count++; //for masks
			slice_size = (1 << (log_array_length + rs_code_rate - log_slice_number));
			slice_real_ele_cnt = slice_size >> rs_code_rate;
			//prepare polynomial division/decomposition
			//things are going to be sliced
			//get evaluations
			l_eval_len = slice_count * slice_size;
			l_eval = new fieldElement[l_eval_len];
			mask_position_gap = slice_size / private_mask_array.size();
			for (int j = 0; j < 64; ++j) {
				if ((1LL << j) <= mask_position_gap && (1LL << (j + 1)) > mask_position_gap) {
					mask_position_gap = 1LL << j;
					break;
				}
			}
			//padding masks
			int mask_size_after_padding = slice_size / mask_position_gap;

			while (private_mask_array.size() < mask_size_after_padding)
				private_mask_array.push_back(fieldElement(0));
			all_pri_mask = private_mask_array;

			fieldElement *tmp = new fieldElement[std::max(mask_size_after_padding, slice_real_ele_cnt)];

			auto fft_t0 = std::chrono::high_resolution_clock::now();

			init_scratch_pad(slice_size);
			for (int i = 0; i < slice_count; ++i) {
				if (i == slice_count - 1) {
					assert(mask_size_after_padding <= slice_size);
					all_pri_msk_arr = new fieldElement[mask_size_after_padding];
					for (int j = 0; j < mask_size_after_padding; ++j)
						all_pri_msk_arr[j] = private_mask_array[j];
					inverse_fast_fourier_transform(all_pri_msk_arr, mask_size_after_padding,
													mask_size_after_padding,
													fieldElement::getRootOfUnity(mylog(mask_size_after_padding)),
													tmp);
					fast_fourier_transform(tmp, mask_size_after_padding, slice_size,
											fieldElement::getRootOfUnity(mylog(slice_size)),
											&l_eval[(i * slice_size)]);
				} else {
					// bool all_zero = true;
					// auto zero = fieldElement::zero();
					// for (int j = 0; j < slice_real_ele_cnt; ++j) {
					// 	if (private_array[i * slice_real_ele_cnt + j] == zero)
					// 		continue;
					// 	all_zero = false;
					// 	break;
					// }
					// if (all_zero) {
					// 	for (int j = 0; j < slice_size; ++j)
					// 		l_eval[i * slice_size + j] = zero;
					// } else {
						inverse_fast_fourier_transform(&private_array[i * slice_real_ele_cnt], slice_real_ele_cnt,
														slice_real_ele_cnt,
														fieldElement::getRootOfUnity(mylog(slice_real_ele_cnt)),
														tmp);
						fast_fourier_transform(tmp, slice_real_ele_cnt, slice_size,
												fieldElement::getRootOfUnity(mylog(slice_size)),
												&l_eval[i * slice_size]);
					// }
				}
			}

			auto fft_t1 = std::chrono::high_resolution_clock::now();
			//printf("FFT Prepare time %lf\n", std::chrono::duration_cast<std::chrono::duration<double>>(fft_t1 - fft_t0).count());
			delete[] tmp;

			auto ret = vpd_prover_init(l_eval, l_coef, log_array_length, slice_size, slice_count);

			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
					t1 - t0);
			total_time += time_span.count();
			//printf("VPD prepare time %lf\n", time_span.count());
			return ret;
		}

		__hhash_digest poly_commit_prover::commit_public_array(std::vector<fieldElement> &all_pub_msk, fieldElement *public_array, int r_0_len, fieldElement target_sum, fieldElement *all_sum) {
			std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
			assert(pre_prepare_executed);
			fri::virtual_oracle_witness = new fieldElement[slice_size * slice_count];
			fri::virtual_oracle_witness_msk = new fieldElement[slice_size];
			fri::virtual_oracle_witness_msk_mapping = new int[slice_size];
			fri::virtual_oracle_witness_mapping = new int[slice_size * slice_count];
			q_eval_len = l_eval_len;
			q_eval = new fieldElement[q_eval_len];

			int mask_size_after_padding = slice_size / mask_position_gap;
			assert(mask_size_after_padding >= all_pri_mask.size());
			while (all_pub_msk.size() < mask_size_after_padding)
				all_pub_msk.push_back(fieldElement(0));
			all_pub_msk_arr = new fieldElement[all_pub_msk.size()];
			for (int i = 0; i < all_pub_msk.size(); ++i)
				all_pub_msk_arr[i] = all_pub_msk[i];
			fieldElement *tmp;
			tmp = new fieldElement[std::max((int) all_pub_msk.size(), slice_size)];
			fieldElement *all_pub_msk_eval = new fieldElement[slice_size];

			double fft_time = 0;
			double re_mapping_time = 0;

			auto fft_t0 = std::chrono::high_resolution_clock::now();

			inverse_fast_fourier_transform(all_pub_msk_arr, all_pub_msk.size(), all_pub_msk.size(),
											fieldElement::getRootOfUnity(mylog(all_pub_msk.size())), tmp);
			fast_fourier_transform(tmp, all_pub_msk.size(), slice_size,
									fieldElement::getRootOfUnity(mylog(slice_size)), all_pub_msk_eval);


			for (int i = 0; i < all_pub_msk.size(); ++i) {
				assert(all_pub_msk[i] == all_pub_msk_eval[i * mask_position_gap]);
			}
			for (int i = 0; i < slice_count; ++i) {
				if (i == slice_count - 1) {
					for (int j = 0; j < slice_size; ++j) {
						q_eval[j + (i * slice_size)] = all_pub_msk_eval[j];
					}
				} else {
					inverse_fast_fourier_transform(&public_array[i * slice_real_ele_cnt], slice_real_ele_cnt,
													slice_real_ele_cnt,
													fieldElement::getRootOfUnity(mylog(slice_real_ele_cnt)), tmp);
					fast_fourier_transform(tmp, slice_real_ele_cnt, slice_size,
											fieldElement::getRootOfUnity(mylog(slice_size)),
											&q_eval[i * slice_size]);
				}
			}

			auto fft_t1 = std::chrono::high_resolution_clock::now();
			fft_time += std::chrono::duration_cast<std::chrono::duration<double>>(fft_t1 - fft_t0).count();

			fieldElement sum = fieldElement(0);
			assert((slice_count - 1) * slice_real_ele_cnt == (1 << r_0_len));
			for (int i = 0; i < (slice_count - 1) * slice_real_ele_cnt; ++i) {
				assert((i << rs_code_rate) < q_eval_len);
				sum = sum + q_eval[i << rs_code_rate] * l_eval[i << rs_code_rate];
			}
			assert(mylog(slice_size / all_pub_msk.size()) == mylog(mask_position_gap));
			for (int i = 0; i < slice_size; i += mask_position_gap) {
				sum = sum + all_pri_mask[i / mask_position_gap] * all_pub_msk[i / mask_position_gap];
				assert(all_pub_msk[i / mask_position_gap] == q_eval[i + (slice_count - 1) * slice_size]);
				assert(all_pri_mask[i / mask_position_gap] == l_eval[i + (slice_count - 1) * slice_size]);
			}
			//assert(sum == target_sum);
			//do fft for q_eval
			//experiment
			//poly mul first
			lq_eval = new fieldElement[std::max(2 * slice_size / mask_position_gap, 2 * slice_real_ele_cnt)];
			h_coef = new fieldElement[std::max(slice_real_ele_cnt, slice_size / mask_position_gap)];
			lq_coef = new fieldElement[std::max(2 * slice_size / mask_position_gap, 2 * slice_real_ele_cnt)];
			h_eval = new fieldElement[std::max(slice_size,
												std::max(slice_size / mask_position_gap, slice_real_ele_cnt))];
			h_eval_arr = new fieldElement[slice_count * slice_size];
			const int log_leaf_size = log_slice_number + 1;
			for (int i = 0; i < slice_count; ++i) {
				if (i == slice_count - 1) {
					assert(mask_position_gap != 1);
					int id = 0;
					bool all_zero = true;
					auto zero = fieldElement(0);
					for (int j = 0; j < slice_size; j += mask_position_gap / 2, id++) {
						lq_eval[id] = l_eval[i * slice_size + j] * q_eval[i * slice_size + j];
						if (!(lq_eval[id] == zero)) {
							all_zero = false;
						}
					}
					if (all_zero) {
						for (int j = 0; j < 2 * slice_size / mask_position_gap; ++j)
							lq_coef[j] = zero;
						for (int j = 0; j < slice_size / mask_position_gap; ++j)
							h_coef[j] = zero;
						for (int j = 0; j < slice_size; ++j)
							h_eval[j] = zero;
					} else {
						fft_t0 = std::chrono::high_resolution_clock::now();
						inverse_fast_fourier_transform(lq_eval, 2 * slice_size / mask_position_gap,
														2 * slice_size / mask_position_gap,
														fieldElement::getRootOfUnity(
																mylog(2 * slice_size / mask_position_gap)), lq_coef);
						for (int j = 0; j < slice_size / mask_position_gap; ++j) {
							h_coef[j] = lq_coef[j + slice_size / mask_position_gap];
						}
						fast_fourier_transform(h_coef, slice_size / mask_position_gap, slice_size,
												fieldElement::getRootOfUnity(mylog(slice_size)), h_eval);
						fft_t1 = std::chrono::high_resolution_clock::now();
						fft_time += std::chrono::duration_cast<std::chrono::duration<double>>(
								fft_t1 - fft_t0).count();
					}
					//test
					auto rou = fieldElement::getRootOfUnity(mylog(slice_size));
					auto inv_x = fieldElement(1);
					auto rou_n = fieldElement::fastPow(rou, slice_size / mask_position_gap);
					auto x_n = fieldElement(1);
					auto inv_rou = rou.inv();
					for (int j = 0; j < slice_size; ++j) {
						//	assert(l_eval[i * slice_size + j] * q_eval[i * slice_size + j] == g_eval[i][j] + (fieldElement::fastPow(x, slice_size / mask_position_gap) - 1) * h_eval[i][j]);
						auto g = l_eval[i * slice_size + j] * q_eval[i * slice_size + j] - (x_n - 1) * h_eval[j];
						if (j < slice_size / 2) {
							fri::virtual_oracle_witness_msk[j << 1] =
									(g - (lq_coef[0] + h_coef[0])) * fieldElement(slice_size / mask_position_gap) *
									inv_x;
							fri::virtual_oracle_witness_msk_mapping[j] = j << 1;
						} else {
							int jj = j - slice_size / 2;
							fri::virtual_oracle_witness_msk[jj << 1 | 1] =
									(g - (lq_coef[0] + h_coef[0])) * fieldElement(slice_size / mask_position_gap) *
									inv_x;
							fri::virtual_oracle_witness_msk_mapping[j] = jj << 1;
						}
						inv_x = inv_x * inv_rou;
						x_n = x_n * rou_n;
					}
					all_sum[i] = (lq_coef[0] + h_coef[0]) * fieldElement(slice_size / mask_position_gap);
				} else {
					assert(2 * slice_real_ele_cnt <= slice_size);
					// bool all_zero = true;
					// auto zero = fieldElement(0);
					for (int j = 0; j < 2 * slice_real_ele_cnt; ++j) {
						lq_eval[j] = l_eval[i * slice_size + j * (slice_size / (2 * slice_real_ele_cnt))] *
										q_eval[i * slice_size + j * (slice_size / (2 * slice_real_ele_cnt))];
						// if (lq_eval[j] != zero) {
						// 	all_zero = false;
						// }
					}
					// if (all_zero) {
					// 	for (int j = 0; j < 2 * slice_real_ele_cnt; ++j)
					// 		lq_coef[j] = zero;
					// 	for (int j = 0; j < slice_real_ele_cnt; ++j)
					// 		h_coef[j] = zero;
					// 	for (int j = 0; j < slice_size; ++j)
					// 		h_eval[j] = zero;
					// } else {
						fft_t0 = std::chrono::high_resolution_clock::now();
						inverse_fast_fourier_transform(lq_eval, 2 * slice_real_ele_cnt, 2 * slice_real_ele_cnt,
														fieldElement::getRootOfUnity(mylog(2 * slice_real_ele_cnt)),
														lq_coef);
						for (int j = 0; j < slice_real_ele_cnt; ++j)
							h_coef[j] = lq_coef[j + slice_real_ele_cnt];
						fast_fourier_transform(h_coef, slice_real_ele_cnt, slice_size,
												fieldElement::getRootOfUnity(mylog(slice_size)), h_eval);
						fft_t1 = std::chrono::high_resolution_clock::now();
						fft_time += std::chrono::duration_cast<std::chrono::duration<double>>(
								fft_t1 - fft_t0).count();
					// }
					auto rou = fieldElement::getRootOfUnity(mylog(slice_size));
					auto inv_rou = rou.inv();
					auto rou_n = fieldElement::fastPow(rou, slice_real_ele_cnt);
					auto x_n = fieldElement(1);
					auto inv_x = fieldElement(slice_real_ele_cnt); //absorb the constant factor
					auto remap_t0 = std::chrono::high_resolution_clock::now();
					auto const_sum = fieldElement(0) - (lq_coef[0] + h_coef[0]);
					for (int j = 0; j < slice_size; ++j) {
						//	assert(l_eval[i * slice_size + j] * q_eval[i * slice_size + j] == g_eval[i][j] + (fieldElement::fastPow(x, slice_real_ele_cnt) - 1) * h_eval[i][j]);
						auto g = l_eval[i * slice_size + j] * q_eval[i * slice_size + j] - (x_n - 1) * h_eval[j];
						if (j < slice_size / 2) {
							fri::virtual_oracle_witness[(j) << log_leaf_size | (i << 1) | 0] =
									(g + const_sum) * inv_x;
							fri::virtual_oracle_witness_mapping[j << log_slice_number | i] =
									j << log_leaf_size | (i << 1) | 0;
						} else {
							int jj = j - slice_size / 2;
							fri::virtual_oracle_witness[(jj) << log_leaf_size | (i << 1) | 1] =
									(g + const_sum) * inv_x;
							fri::virtual_oracle_witness_mapping[jj << log_slice_number | i] =
									jj << log_leaf_size | (i << 1) | 0;
						}
						inv_x = inv_x * inv_rou;
						x_n = x_n * rou_n;
					}
					auto remap_t1 = std::chrono::high_resolution_clock::now();
					re_mapping_time += std::chrono::duration_cast<std::chrono::duration<double>>(
							remap_t1 - remap_t0).count();
					assert(i < slice_number + 1);
					all_sum[i] = (lq_coef[0] + h_coef[0]) * fieldElement(slice_real_ele_cnt);
				}
				for (int j = 0; j < slice_size; ++j)
					h_eval_arr[i * slice_size + j] = h_eval[j];
			}
			delete[] tmp;
			delete[] lq_eval;
			delete[] h_coef;
			delete[] lq_coef;
			delete[] h_eval;
			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
					t1 - t0);
			total_time += time_span.count();
			//printf("PostGKR FFT time %lf\n", fft_time);
			//printf("PostGKR remap time %lf\n", re_mapping_time);
			//printf("PostGKR prepare time 0 %lf\n", time_span.count());

			t0 = std::chrono::high_resolution_clock::now();
			auto ret = fri::request_init_commit(r_0_len, 1);
			t1 = std::chrono::high_resolution_clock::now();
			time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
			total_time += time_span.count();

			//printf("PostGKR prepare time 1 %lf\n", time_span.count());
			return ret;
		}
	}
}
