#pragma once

#include <vector>
#include "fieldElement.hpp"
namespace virgo_ext {
    bool vpd_verify(prime_field::field_element all_mask_sum, double &v_time);
}