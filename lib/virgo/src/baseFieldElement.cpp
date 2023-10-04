#undef NDEBUG
#include "assert.h"

#include "baseFieldElement.hpp"
#include <cstdio>
#include <chrono>
#include <immintrin.h>
#include <cstdint>

using namespace std;

namespace virgo_ext {
	unsigned long long baseFieldElement::mod = 0;
	unsigned long long baseFieldElement::rou = 0;
	unsigned long long baseFieldElement::rcp = 0;
	unsigned int baseFieldElement::len = 0;
	unsigned int baseFieldElement::__max_order = 0;
	bool baseFieldElement::initialized = false;
	int baseFieldElement::multCounter, baseFieldElement::addCounter;
	bool baseFieldElement::isCounting;
	bool baseFieldElement::isSumchecking;

	baseFieldElement::baseFieldElement() {
		elem[0] = elem[1] = 0;
	}
	
	baseFieldElement::baseFieldElement(const baseFieldElement & b) {
		elem[0] = b.elem[0];
		elem[1] = b.elem[1];
	}

	baseFieldElement::baseFieldElement(long long x) {
		elem[0] = x >= 0 ? x : mod + x;
		elem[1] = 0;
	}

	baseFieldElement::baseFieldElement(long long int x, long long int y) {
		elem[0] = x >= 0 ? x : mod + x;
		elem[1] = y >= 0 ? y : mod + y;
	}

	baseFieldElement baseFieldElement::operator+(const baseFieldElement & other) const {
		if (isCounting) {
			++addCounter;
		}
		baseFieldElement ret;
		for (size_t i = 0; i < 2; i++) {
			ret.elem[i] = other.elem[i] + elem[i];
			if (mod <= ret.elem[i])
				ret.elem[i] = ret.elem[i] - mod;
		}
		return ret;
	}

	baseFieldElement baseFieldElement::operator-(const baseFieldElement & other) const {
		if (isCounting) {
			++addCounter;
		}
		baseFieldElement ret;
		for (size_t i = 0; i < 2; i++) {
			auto t = mod - other.elem[i];	//tmp_r == -b.real is true in this prime field
			ret.elem[i] = elem[i] + t;
			if (ret.elem[i] >= mod)
				ret.elem[i] -= mod;
		}
		return ret;
	}

	baseFieldElement baseFieldElement::operator*(const baseFieldElement & other) const {
		if (isCounting) {
			++multCounter;
		}
		baseFieldElement ret;
		auto t0 = elem[0] + elem[1];
		if (t0 >= mod)
			t0 -= mod;
		auto t1 = other.elem[0] + other.elem[1];
		if (t1 >= mod)
			t1 -= mod;
		auto all_prod = mymult(t0, t1);
		auto ac = mymult(elem[0], other.elem[0]);
		auto bd = mymult(elem[1], other.elem[1]);
		auto nac = mod - ac;
		auto nbd = mod - bd;

		t0 = ac + nbd + nbd + nbd;
		while (t0 >= mod)
			t0 -= mod;
		ret.elem[0] = t0;

		t1 = all_prod + nac + nbd;
		while (t1 >= mod)
			t1 -= mod;
		ret.elem[1] = t1;
		return ret;
	}

	baseFieldElement baseFieldElement::operator-() const {
		if (isCounting) {
			++addCounter;
		}
		return zero() - *this;
	}

	void baseFieldElement::init(unsigned long long prime, unsigned long long root) {
		initialized = true;
		srand(3396);
		isCounting = false;
		mod = prime;
		rou = root;
		len = 64 - __builtin_clzll(mod);
		rcp = 1 + ((__int128_t) 1 << (2 * len)) / mod;

		__max_order = 0;
		unsigned long long order = prime - 1;
		while (order % 2 == 0) {
			__max_order++;
			order = order >> 1;
		}
	}

	unsigned int baseFieldElement::maxOrder() {
		return __max_order;
	}

	baseFieldElement baseFieldElement::random() {
		baseFieldElement ret;
		ret.elem[0] = randomNumber() % mod;
		ret.elem[1] = randomNumber() % mod;
		return ret;
	}

	bool baseFieldElement::operator!=(const baseFieldElement & other) const {
		return elem[0] != other.elem[0] || elem[1] != other.elem[1];
	}
	
	bool baseFieldElement::operator==(const baseFieldElement & other) const {
		return !(*this != other);
	}
	
	baseFieldElement & baseFieldElement::operator=(const baseFieldElement & other) {
		elem[0] = other.elem[0];
		elem[1] = other.elem[1];
		return *this;
	}

	baseFieldElement & baseFieldElement::operator+=(const baseFieldElement & other) {
		*this = *this + other;
		return *this;
	}

	baseFieldElement & baseFieldElement::operator-=(const baseFieldElement & other) {
		*this = *this - other;
		return *this;
	}

	baseFieldElement & baseFieldElement::operator*=(const baseFieldElement & other) {
		*this = *this * other;
		return *this;
	}

	baseFieldElement::operator  bool () const {
		return elem[0] || elem[1];
	}
	
	bool baseFieldElement::isNegative() const {
		return (elem[0] > (mod >> 1)) && (elem[1] == 0);
	}
	
	unsigned char baseFieldElement::getBit(unsigned int i) const {
		assert(elem[1] == 0);
		return (elem[0] >> i) & 1;
	}
	
	bool baseFieldElement::operator<(const baseFieldElement & other) const {
		assert(elem[1] == 0);
		return elem[0] < other.elem[0];
	}
	
	bool baseFieldElement::isZero() {
		return !elem[0] && !elem[1];
	}

	baseFieldElement baseFieldElement::abs() const {
		assert(elem[1] == 0);
		baseFieldElement res = -*this;
		 return res.elem[0] < this->elem[0] ? res : elem[0];
	}
	
	baseFieldElement baseFieldElement::sqr() const {
		return (*this) * (*this);
	}
	
	baseFieldElement baseFieldElement::inv() const {
		baseFieldElement ret;
		auto t0 = mymult(elem[0], elem[0]);
		auto t1 = mymult(elem[1], elem[1]);

		 t0 = t0 + t1 + t1 + t1;
		while (t0 > mod)
			 t0 -= mod;

		 t1 = 1;
		unsigned long long p = mod - 2;
		while (p) {
			if (p & 1) {
				t1 = mymult(t1, t0);
			}
			t0 = mymult(t0, t0);
			p >>= 1;
		}
		ret.elem[0] = mymult(elem[0], t1);
		ret.elem[1] = mod - mymult(elem[1], t1);
		return ret;
	}

	void baseFieldElement::setAbs() {
		*this = this->abs();
	}

	void baseFieldElement::setSqr() {
		*this = this->sqr();
	}

	void baseFieldElement::setInv() {
		*this = this->inv();
	}

	void baseFieldElement::print(FILE *fileno) const {
		fprintf(fileno, "(%llu, %llu)\n", elem[0], elem[1]);
	}
	
	baseFieldElement baseFieldElement::maxWithZero(const baseFieldElement & a,
			const baseFieldElement & b) {
		if (a.isNegative() && b.isNegative())
			return baseFieldElement::zero();
		return baseFieldElement(a.isNegative()? b : b.isNegative()? a : std::max(a.
						elem[0], b.elem[0]));
	}

	baseFieldElement baseFieldElement::maxUnsigned(const baseFieldElement & a,
			const baseFieldElement & b) {
		return a < b ? b : a;
	}

	baseFieldElement baseFieldElement::getRootOfUnity(int log_order) {
		baseFieldElement root;
		root.elem[0] = rou;
		root.elem[1] = 0;

		assert(log_order <= __max_order);

		for (int i = 0; i < __max_order - log_order; ++i)
			root = root * root;

		return root;
	}

	baseFieldElement baseFieldElement::zero() {
		return baseFieldElement(0);
	}

	baseFieldElement baseFieldElement::one() {
		return baseFieldElement(1);
	}

	vector < baseFieldElement > baseFieldElement::generateRandomness(u64 size) {
		int k = size;
		vector < baseFieldElement > ret(k);

		for (int i = 0; i < k; ++i)
			ret[i] = baseFieldElement::random();
		return ret;
	}

	baseFieldElement baseFieldElement::innerProd(vector < baseFieldElement >::iterator a,
			vector < baseFieldElement >::iterator b, u64 n) {
		baseFieldElement res = baseFieldElement::zero();
		for (int i = 0; i < n; ++i)
			res += a[i] * b[i];
		return res;
	}

	double baseFieldElement::self_speed_test_add(int repeat) {
		baseFieldElement a, b;
		a = random();
		b = random();
		baseFieldElement c;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < repeat; ++i) {
			c = a + b;
			b = c + a;
			a = c + b;
		}
		printf("%llu %llu\n", c.elem[0], c.elem[1]);
		auto t1 = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast < std::chrono::duration <
				double >>(t1 - t0).count();
	}

	double baseFieldElement::self_speed_test_mult(int repeat) {
		baseFieldElement a, b;
		a = random();
		b = random();
		baseFieldElement c;
		auto t0 = std::chrono::high_resolution_clock::now();
		for (int i = 0; i < repeat; ++i) {
			c = a * b;
			b = c * a;
			a = c * b;
		}
		printf("%llu %llu\n", c.elem[0], c.elem[1]);
		auto t1 = std::chrono::high_resolution_clock::now();
		return std::chrono::duration_cast < std::chrono::duration <
				double >>(t1 - t0).count();
	}

	char *baseFieldElement::toString() const {
		char *s = new char[50];
		 sprintf(s, "(%llu, %llu)", this->elem[0], this->elem[1]);
		 return s;
	}
	
	ostream & operator<<(ostream & out, const baseFieldElement & c) {
		out << '(' << c.elem[0] << ',' << c.elem[1] << ')';
		return out;
	}

	baseFieldElement baseFieldElement::fastPow(baseFieldElement x, unsigned long long p) {
		baseFieldElement ret = baseFieldElement(1), tmp = x;
		while (p) {
			if (p & 1) {
				ret = ret * tmp;
			}
			tmp = tmp * tmp;
			p >>= 1;
		}
		return ret;
	}

	unsigned long long baseFieldElement::mymult(const unsigned long long int x,
			const unsigned long long int y) {
		//return a value between [0, PRIME) = x * y mod PRIME
		unsigned long long lo, hi;
		lo = _mulx_u64(x, y, &hi);
		hi = (hi << (64 - len)) | (lo >> len);
		lo = lo & ((1LL << len) - 1);

		uint64_t q = ((__uint128_t) hi * rcp) >> len;
		uint64_t q0 = ((uint64_t) q * mod) & ((1L << len) - 1);
		uint64_t q1 = ((__uint128_t) q * mod) >> len;
		q0 = (lo - q0);
		q1 = (hi - q1) - (lo < q0);
		q0 = (q0 - q1*mod) & ((1L << len) - 1);
		while (q0 >= mod) q0 -= mod;
		return q0;
    }

	unsigned long long baseFieldElement::randomNumber() {
		unsigned long long ret =::random() % 10;
		for (int i = 1; i < 20; ++i)
			ret = (ret * 10ull + (unsigned long long)(::random() % 10)) % mod;
		return ret;
	}

	baseFieldElement baseFieldElement::mulNor() {
		unsigned long long lo, hi, neg;

		neg = (mod - elem[0]) + (mod - elem[1]);
		lo = neg + neg + neg;
		while (lo > mod)
			lo -= mod;
		hi = mod - elem[1];
		hi = (elem[0] + hi + hi + hi);
		while (hi > mod)
			hi -= mod;
		baseFieldElement ret;
		ret.elem[0] = lo;
		ret.elem[1] = hi;
		return ret;
	}
}