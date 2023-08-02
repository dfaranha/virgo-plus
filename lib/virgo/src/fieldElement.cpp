#undef NDEBUG
#include "assert.h"
#include "string.h"

#include "fieldElement.hpp"
#include <cstdio>
#include <chrono>
#include <immintrin.h>
#include <cstdint>

using namespace std;

namespace virgo {
    mp_limb_t fieldElement::mod[MOD_LIMBS] = {0};
	mp_limb_t fieldElement::rou[MOD_LIMBS] = {0};
	unsigned int fieldElement::len = 0;
	unsigned int fieldElement::__max_order = 0;
    bool fieldElement::initialized = false;
    int fieldElement::multCounter, fieldElement::addCounter;
    bool fieldElement::isCounting;
    bool fieldElement::isSumchecking;

    fieldElement::fieldElement() {
		memset(elem, 0, MOD_LIMBS * sizeof(mp_limb_t));
    }

    fieldElement::fieldElement(const fieldElement &b) {
		memcpy(&elem, &b.elem, MOD_LIMBS * sizeof(mp_limb_t));
    }

    fieldElement::fieldElement(long long x) {
		memset(elem, 0, MOD_LIMBS * sizeof(mp_limb_t));
		if (x >= 0) {
			mpn_add_1(elem, elem, MOD_LIMBS, x);
		} else {
			mpn_sub_1(elem, mod, MOD_LIMBS, -x);
		}
    }

    fieldElement fieldElement::operator+(const fieldElement &other) const {
        if (isCounting)
        {
            ++addCounter;
        }
        fieldElement ret;
		mp_limb_t carry = mpn_add_n(ret.elem, elem, other.elem, MOD_LIMBS);
		if (carry || (mpn_cmp(ret.elem, mod, MOD_LIMBS) >= 0)) {
			carry = mpn_sub_n(ret.elem, ret.elem, mod, MOD_LIMBS);
		}
        return ret;
    }

    fieldElement fieldElement::operator*(const fieldElement &other) const {
        if (isCounting)
        {
            ++multCounter;
        }
        fieldElement ret;
		mp_limb_t t[2 * MOD_LIMBS] = {0};
		mp_limb_t q[2 * MOD_LIMBS] = {0};
		mpn_mul_n(t, elem, other.elem, MOD_LIMBS);
		mpn_tdiv_qr(q, t, 0, t, 2 * MOD_LIMBS, mod, MOD_LIMBS);
		memcpy(ret.elem, t, MOD_LIMBS * sizeof(mp_limb_t));
        return ret;
    }

    fieldElement fieldElement::operator-(const fieldElement &other) const {
        if (isCounting)
        {
            ++addCounter;
        }
		mp_limb_t t[MOD_LIMBS + 1] = {0};
        fieldElement ret;
		if (mpn_sub_n(ret.elem, elem, other.elem, MOD_LIMBS)) {
			mpn_add_n(ret.elem, ret.elem, mod, MOD_LIMBS);
		}
        return ret;
    }

    fieldElement fieldElement::operator-() const {
        if (isCounting)
        {
            ++addCounter;
        }
        return zero() - *this;
    }

    void fieldElement::init(unsigned long long prime, unsigned long long root) {
        initialized = true;
        srand(3396);
        isCounting = false;
		mod[0] = prime;
		rou[0] = root;
		len = 64 * MOD_LIMBS - __builtin_clzll(mod[MOD_LIMBS - 1]);

		__max_order = 0;
		prime = prime - 1;
		while (prime % 2 == 0) {
			__max_order++;
			prime = prime >> 1;
		}
    }

	unsigned int fieldElement::maxOrder() {
		return __max_order;
	}

    fieldElement fieldElement::random() {
        fieldElement ret;
		mpn_random(ret.elem, MOD_LIMBS);
		while (mpn_cmp(ret.elem, mod, MOD_LIMBS) > 0) {
			mpn_sub_n(ret.elem, ret.elem, mod, MOD_LIMBS);
		}
        return ret;
    }

    bool fieldElement::operator!=(const fieldElement &other) const {
        return mpn_cmp(elem, other.elem, MOD_LIMBS) != 0;
    }

    bool fieldElement::operator==(const fieldElement &other) const {
        return !(*this != other);
    }

    fieldElement &fieldElement::operator=(const fieldElement &other) {
		memcpy(elem, other.elem, MOD_LIMBS * sizeof(mp_limb_t));
        return *this;
    }

    fieldElement &fieldElement::operator+=(const fieldElement &other) {
        *this = *this + other;
        return *this;
    }

    fieldElement &fieldElement::operator-=(const fieldElement &other) {
        *this = *this - other;
        return *this;
    }

    fieldElement &fieldElement::operator*=(const fieldElement &other) {
        *this = *this * other;
        return *this;
    }

    fieldElement::operator bool() const {
        return elem;
    }

    bool fieldElement::isNegative() const {
		mp_limb_t t[MOD_LIMBS] = {0};
		mpn_rshift(t, mod, MOD_LIMBS, 1);
        return (mpn_cmp(elem, t, MOD_LIMBS) > 0);
    }

    unsigned char fieldElement::getBitWidth() const {
		return 64 * MOD_LIMBS - __builtin_clzll(elem[MOD_LIMBS - 1]);
    }

    unsigned char fieldElement::getBit(unsigned int i) const {
        return (elem[i / 64] >> (i % 64)) & 1;
    }

    bool fieldElement::operator<(const fieldElement &other) const {
		return (mpn_cmp(elem, other.elem, MOD_LIMBS) < 0);
    }

    bool fieldElement::isZero() {
		for (int i = MOD_LIMBS - 1; i >= 0; i++) {
			if (elem[i]) {
				return false;
			}
		}
        return true;
    }

    fieldElement fieldElement::abs() const {
		fieldElement res = -*this;
        return (mpn_cmp(res.elem, this->elem, MOD_LIMBS) < 0) ? res : *this;
    }

    fieldElement fieldElement::sqr() const {
        return (*this) * (*this);
    }

    fieldElement fieldElement::inv() const {
		fieldElement ret;
		mp_size_t cn;
		mp_limb_t s[1], t[2];

		memcpy(s, &mod, sizeof(mp_limb_t));

		mpn_gcdext(t, ret.elem, &cn, (mp_ptr)elem, MOD_LIMBS, s, MOD_LIMBS);
		if (cn < 0) {
			memset(ret.elem - cn, 0, (1 + cn)*sizeof(mp_limb_t));
			mpn_sub_n(ret.elem, mod, ret.elem, MOD_LIMBS);
		} else {
			memset(ret.elem + cn, 0, (1 - cn)*sizeof(mp_limb_t));
		}

		return ret;
    }

    void fieldElement::setAbs() {
        *this = this -> abs();
    }

    void fieldElement::setSqr() {
        *this = this -> sqr();
    }

    void fieldElement::setInv() {
        *this = this -> inv();
    }

    void fieldElement::print(FILE *fileno) const {
		int j = 0;
		for (int i = MOD_LIMBS - 1; i >= 0; i--) {
			j += fprintf(fileno, "%llu", (uint64_t) elem[i]);
		}
		fprintf(fileno, "\n");
    }

    fieldElement fieldElement::getRootOfUnity(int log_order) {
        fieldElement root;
        //general root of unity, have log_order 2^30
		memcpy(root.elem, rou, MOD_LIMBS * sizeof(mp_limb_t));

        assert(log_order <= __max_order);

        for (int i = 0; i < __max_order - log_order; ++i) {
            root = root * root;
		}


        return root;
    }

    fieldElement fieldElement::zero() {
        return fieldElement(0);
    }

    fieldElement fieldElement::one() {
        return fieldElement(1);
    }

    fieldElement fieldElement::innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n) {
        fieldElement res = fieldElement::zero();
        for (int i = 0; i < n; ++i)
            res += a[i] * b[i];
        return res;
    }

    char *fieldElement::toString() const {
		int j = 0;
        char *s = new char[100];
		for (int i = MOD_LIMBS - 1; i >= 0; i--) {
			j += sprintf(s + j, "%llu", (uint64_t) elem[i]);
		}
        return s;
    }

    ostream &operator << (ostream &out, const fieldElement &c) {
        out << '(' << c.elem[0] << ')';
        return out;
    }

    fieldElement fieldElement::fastPow(fieldElement x, unsigned long long p) {
        fieldElement ret = fieldElement(1), tmp = x;
        while (p)
        {
            if (p & 1)
            {
                ret = ret * tmp;
            }
            tmp = tmp * tmp;
            p >>= 1;
        }
        return ret;
    }
}
