#undef NDEBUG
#include "assert.h"

#include "fieldElement.hpp"
#include <cstdio>
#include <chrono>
#include <immintrin.h>
#include <cstdint>

using namespace std;

namespace virgo_ext {
	unsigned long long fieldElement::mod = 0;
	unsigned long long fieldElement::rou = 0;
	unsigned long long fieldElement::rcp = 0;
	unsigned int fieldElement::len = 0;
	unsigned int fieldElement::__max_order = 0;
	bool fieldElement::initialized = false;
	int fieldElement::multCounter, fieldElement::addCounter;
	bool fieldElement::isCounting;
	bool fieldElement::isSumchecking;

	fieldElement::fieldElement() {
		elem[0] = elem[1] = 0;
	}
	
	fieldElement::fieldElement(const fieldElement & b) {
		elem[0] = b.elem[0];
		elem[1] = b.elem[1];
	}

	fieldElement::fieldElement(long long x) {
		elem[0] = x >= 0 ? x : mod + x;
		elem[1] = 0;
	}

	fieldElement::fieldElement(long long int x, long long int y) {
		elem[0] = x >= 0 ? x : mod + x;
		elem[1] = y >= 0 ? y : mod + y;
	}

	fieldElement fieldElement::operator+(const fieldElement & other) const {
		if (isCounting) {
			++addCounter;
		}
		fieldElement ret;
		for (size_t i = 0; i < 2; i++) {
			ret.elem[i] = other.elem[i] + elem[i];
			if (mod <= ret.elem[i])
				ret.elem[i] = ret.elem[i] - mod;
		}
		return ret;
	}

	fieldElement fieldElement::operator-(const fieldElement & other) const {
		if (isCounting) {
			++addCounter;
		}
		fieldElement ret;
		for (size_t i = 0; i < 2; i++) {
			auto t = mod - other.elem[i];	//tmp_r == -b.real is true in this prime field
			ret.elem[i] = elem[i] + t;
			if (ret.elem[i] >= mod)
				ret.elem[i] -= mod;
		}
		return ret;
	}

	fieldElement fieldElement::operator*(const fieldElement & other) const {
		if (isCounting) {
			++multCounter;
		}
		fieldElement ret;
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

	fieldElement fieldElement::operator-() const {
		if (isCounting) {
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
		rcp = 1 + ((__int128_t) 1 << (2 * len)) / mod;

		__max_order = 0;
		unsigned long long order = prime * prime - 1;
		while (order % 2 == 0) {
			__max_order++;
			order = order >> 1;
		}
	}

	unsigned int fieldElement::maxOrder() {
		return __max_order;
	}

	fieldElement fieldElement::random() {
		fieldElement ret;
		ret.elem[0] = randomNumber() % mod;
		ret.elem[1] = randomNumber() % mod;
		return ret;
	}

	bool fieldElement::operator!=(const fieldElement & other) const {
		return elem[0] != other.elem[0] || elem[1] != other.elem[1];
	}
	
	bool fieldElement::operator==(const fieldElement & other) const {
		return !(*this != other);
	}
	
	fieldElement & fieldElement::operator=(const fieldElement & other) {
		elem[0] = other.elem[0];
		elem[1] = other.elem[1];
		return *this;
	}

	fieldElement & fieldElement::operator+=(const fieldElement & other) {
		*this = *this + other;
		return *this;
	}

	fieldElement & fieldElement::operator-=(const fieldElement & other) {
		*this = *this - other;
		return *this;
	}

	fieldElement & fieldElement::operator*=(const fieldElement & other) {
		*this = *this * other;
		return *this;
	}

	fieldElement::operator  bool () const {
		return elem[0] || elem[1];
	}
	
	bool fieldElement::isNegative() const {
		return (elem[0] > (mod >> 1)) && (elem[1] == 0);
	}
	
	unsigned char fieldElement::getBit(unsigned int i) const {
		assert(elem[1] == 0);
		return (elem[0] >> i) & 1;
	}
	
	bool fieldElement::operator<(const fieldElement & other) const {
		assert(elem[1] == 0);
		return elem[0] < other.elem[0];
	}
	
	bool fieldElement::isZero() {
		return !elem[0] && !elem[1];
	}

	fieldElement fieldElement::abs() const {
		assert(elem[1] == 0);
		fieldElement res = -*this;
		 return res.elem[0] < this->elem[0] ? res : elem[0];
	}
	
	fieldElement fieldElement::sqr() const {
		return (*this) * (*this);
	}
	
	fieldElement fieldElement::inv() const {
		fieldElement ret;
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

	void fieldElement::setAbs() {
		*this = this->abs();
	}

	void fieldElement::setSqr() {
		*this = this->sqr();
	}

	void fieldElement::setInv() {
		*this = this->inv();
	}

	void fieldElement::print(FILE *fileno) const {
		fprintf(fileno, "(%llu, %llu)\n", elem[0], elem[1]);
	}
	
	fieldElement fieldElement::maxWithZero(const fieldElement & a,
			const fieldElement & b) {
		if (a.isNegative() && b.isNegative())
			return fieldElement::zero();
		return fieldElement(a.isNegative()? b : b.isNegative()? a : std::max(a.
						elem[0], b.elem[0]));
	}

	fieldElement fieldElement::maxUnsigned(const fieldElement & a,
			const fieldElement & b) {
		return a < b ? b : a;
	}

	fieldElement fieldElement::getRootOfUnity(int log_order) {
		fieldElement root;
		root.elem[0] = 0;
		root.elem[1] = rou;

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

	vector < fieldElement > fieldElement::generateRandomness(u64 size) {
		int k = size;
		vector < fieldElement > ret(k);

		for (int i = 0; i < k; ++i)
			ret[i] = fieldElement::random();
		return ret;
	}

	fieldElement fieldElement::innerProd(vector < fieldElement >::iterator a,
			vector < fieldElement >::iterator b, u64 n) {
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

	double fieldElement::self_speed_test_mult(int repeat) {
		fieldElement a, b;
		a = random();
		b = random();
		fieldElement c;
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

	char *fieldElement::toString() const {
		char *s = new char[50];
		 sprintf(s, "(%llu, %llu)", this->elem[0], this->elem[1]);
		 return s;
	}
	
	ostream & operator<<(ostream & out, const fieldElement & c) {
		out << '(' << c.elem[0] << ',' << c.elem[1] << ')';
		return out;
	}

	fieldElement fieldElement::fastPow(fieldElement x, unsigned long long p) {
		fieldElement ret = fieldElement(1), tmp = x;
		while (p) {
			if (p & 1) {
				ret = ret * tmp;
			}
			tmp = tmp * tmp;
			p >>= 1;
		}
		return ret;
	}

	unsigned long long fieldElement::mymult(const unsigned long long int x,
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

	unsigned long long fieldElement::randomNumber() {
		unsigned long long ret =::random() % 10;
		for (int i = 1; i < 20; ++i)
			ret = (ret * 10ull + (unsigned long long)(::random() % 10)) % mod;
		return ret;
	}
}
