#undef NDEBUG
#include "assert.h"

#include "fieldElement.hpp"
#include <cstdio>
#include <chrono>
#include <immintrin.h>
#include <cstdint>

//#define PRIME31

using namespace std;

namespace virgo {
    mp_limb_t fieldElement::mod = 0;
	mp_limb_t fieldElement::rou = 0;
	mp_limb_t fieldElement::rcp = 0;
	unsigned int fieldElement::len = 0;
	unsigned int fieldElement::__max_order = 0;
    bool fieldElement::initialized = false;
    int fieldElement::multCounter, fieldElement::addCounter;
    bool fieldElement::isCounting;
    bool fieldElement::isSumchecking;

    fieldElement::fieldElement() {
        elem = 0;
    }

    fieldElement::fieldElement(const fieldElement &b) {
		elem = b.elem;
    }

    fieldElement::fieldElement(long long x) {
        elem = x >= 0 ? x : mod + x;
    }

    fieldElement fieldElement::operator+(const fieldElement &other) const {
        if (isCounting)
        {
            ++addCounter;
        }
        fieldElement ret;
        ret.elem = other.elem + elem;
        if (mod <= ret.elem)
            ret.elem = ret.elem - mod;
        return ret;
    }

    fieldElement fieldElement::operator*(const fieldElement &other) const {
        if (isCounting)
        {
            ++multCounter;
        }
        fieldElement ret;
		auto prod = mymult(elem, other.elem);
        while (prod >= mod)
        	prod-= mod;
		ret.elem = prod;
        return ret;
    }

    fieldElement fieldElement::operator-(const fieldElement &other) const {
        if (isCounting)
        {
            ++addCounter;
        }
        fieldElement ret;
        auto tmp_r = mod - other.elem; //tmp_r == -b.real is true in this prime field
        ret.elem = elem + tmp_r;
        if (ret.elem >= mod)
            ret.elem -= mod;

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
		mod = prime;
		rou = root;
		len = 64 - __builtin_clzll(mod);
		rcp = 1 + ((__int128_t)1 << (2 * len)) / mod;

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
        ret.elem = randomNumber() % mod;
        return ret;
    }

    bool fieldElement::operator!=(const fieldElement &other) const {
        return elem != other.elem;
    }

    bool fieldElement::operator==(const fieldElement &other) const {
        return !(*this != other);
    }

    fieldElement &fieldElement::operator=(const fieldElement &other) {
        elem = other.elem;
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
        return (elem > (mod >> 1));
    }

    unsigned char fieldElement::getBitWidth() const {
        auto dat = elem;
        unsigned char res = 0;
        for (int i = 32; i && dat; i >>= 1) {
            if (dat >> i) {
                res += i;
                dat >>= i;
            }
        }
        return res + 1;
    }

    unsigned char fieldElement::getBit(unsigned int i) const {
        return (elem >> i) & 1;
    }

    bool fieldElement::operator<(const fieldElement &other) const {
        return elem < other.elem;
    }

    bool fieldElement::isZero() {
        return !elem;
    }

    fieldElement fieldElement::abs() const {
		fieldElement res = -*this;
        return res.elem < this -> elem ? res : elem;
    }

    fieldElement fieldElement::sqr() const {
        return (*this) * (*this);
    }

    fieldElement fieldElement::inv() const {
        mp_limb_t p = mod;
        return fastPow(*this, p - 2);
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
        fprintf(fileno, "(%llu)\n", elem);
    }

    fieldElement fieldElement::maxWithZero(const fieldElement &a, const fieldElement &b) {
        if (a.isNegative() && b.isNegative()) return fieldElement::zero();
        return fieldElement(a.isNegative() ? b : b.isNegative() ? a : std::max(a.elem, b.elem));
    }

    fieldElement fieldElement::maxUnsigned(const fieldElement &a, const fieldElement &b) {
        return a < b ? b : a;
    }

    fieldElement fieldElement::getRootOfUnity(int log_order) {
        fieldElement root;
        //general root of unity, have log_order 2^30
        root.elem = rou;
		
        assert(log_order <= __max_order);

        for (int i = 0; i < __max_order - log_order; ++i)
            root = root * root;

        return root;
    }

    fieldElement fieldElement::zero() {
        return fieldElement(0);
    }

    fieldElement fieldElement::one() {
        return fieldElement(1);
    }

    vector<fieldElement> fieldElement::generateRandomness(u64 size) {
        int k = size;
        vector<fieldElement> ret(k);

        for (int i = 0; i < k; ++i)
            ret[i] = fieldElement::random();
        return ret;
    }

    fieldElement fieldElement::innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n) {
        fieldElement res = fieldElement::zero();
        for (int i = 0; i < n; ++i)
            res += a[i] * b[i];
        return res;
    }

    double fieldElement::self_speed_test_add(int repeat) {
        fieldElement a, b;
        a = random();
        b = random();
        fieldElement c;
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < repeat; ++i)
        {
            c = a + b;
            b = c + a;
            a = c + b;
        }
        printf("%llu\n", c.elem);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    }

    double fieldElement::self_speed_test_mult(int repeat) {
        fieldElement a, b;
        a = random();
        b = random();
        fieldElement c;
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < repeat; ++i)
        {
            c = a * b;
            b = c * a;
            a = c * b;
        }
        printf("%llu\n", c.elem);
        auto t1 = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();
    }

    char *fieldElement::toString() const {
        char *s = new char[50];
        sprintf(s, "(%llu)", this -> elem);
        return s;
    }

    ostream &operator << (ostream &out, const fieldElement &c) {
        out << '(' << c.elem << ')';
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

    unsigned long long fieldElement::mymult(const unsigned long long int x, const unsigned long long int y) {
        //return a value between [0, 2PRIME) = x * y mod PRIME
        unsigned long long lo, hi;
		lo = _mulx_u64(x, y, &hi);
        hi = (hi << (64 - len)) | (lo >> len);
		lo = lo & ((1LL << len) - 1);

		uint64_t q = ((__uint128_t)hi * rcp) >> len;
		uint64_t q0 = ((uint64_t)q*mod) & ((1L << len) - 1);
		uint64_t q1 = ((__uint128_t)q*mod)>>len;
		q0 = (lo - q0);
		q1 = (hi - q1) - (lo < q0);
		return (q0 - q1*mod) & ((1L << len) - 1);
    }

    unsigned long long fieldElement::randomNumber() {
        unsigned long long ret = ::random() % 10;
        for (int i = 1; i < 20; ++i)
            ret = (ret * 10ull + (unsigned long long)(::random() % 10)) % mod;
        return ret;
    }
}
