//
// Created by 69029 on 6/25/2021.
//

#ifndef UNLAYERED_LIBRA_CONFIG_PC_HPPfpd
#define UNLAYERED_LIBRA_CONFIG_PC_HPPfpd
#define USE_VIRGO

#ifdef USE_VIRGO
#include "poly_commit.h"
#include "timer.hpp"
#include "fieldElement.hpp"
#include "polynomial.h"
#include "RS_polynomial.h"
#include "constants.h"
#include "my_hhash.h"
#define F   virgo_ext::fieldElement
#define F_ONE   virgo_ext::fieldElement::one()
#define F_ZERO  virgo_ext::fieldElement::zero()
#endif

#ifdef USE_HYRAX_P224
#include <hyrax-p224/src/polyVerifier.hpp>
#define F   hyrax_p224::fieldElement
#define G   hyrax_p224::groupElement
#define F_ONE   hyrax_p224::fieldElement::one()
#define F_ZERO  hyrax_p224::fieldElement::zero()
#endif

#endif //UNLAYERED_LIBRA_CONFIG_PC_HPP
