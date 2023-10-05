#ifndef __vpd_verifierfpd
#define __vpd_verifierfpd
#include <vector>
#include "fieldElement.hpp"
namespace virgo_ext {
    bool vpd_verify(prime_field::field_element all_mask_sum, double &v_time);
}
#endif