#undef NDEBUG
#include "assert.h"

#include "fieldElement.hpp"
#include <cstdio>
#include <chrono>
#include <immintrin.h>
#include <cstdint>

using namespace std;

namespace virgo_ext {
	bool fieldElement::initialized = false;
	int fieldElement::multCounter, fieldElement::addCounter;
	bool fieldElement::isCounting;
	bool fieldElement::isSumchecking;

	fieldElement::fieldElement() {
		elem[0] = elem[1] = elem[2] = baseFieldElement(0);
		assert(elem[0].cleartext);
		assert(elem[1].cleartext);
		assert(elem[2].cleartext);
	}
	
	fieldElement::fieldElement(const fieldElement & b) {
		elem[0] = b.elem[0];
		elem[1] = b.elem[1];
		elem[2] = b.elem[2];
	}

	fieldElement::fieldElement(long long x) {
		elem[0] = baseFieldElement(x);
		elem[1] = elem[2] = baseFieldElement(0);
	}

	fieldElement::fieldElement(baseFieldElement x) {
		elem[0] = x;
		elem[1] = elem[2] = baseFieldElement(0);
	}

	fieldElement::fieldElement(helib::Ctxt &x) {
		elem[0] = baseFieldElement(x);
		elem[1] = elem[2] = baseFieldElement(0);
	}

	fieldElement::fieldElement(long long int x, long long int y) {
		elem[0] = baseFieldElement(x, y);
		elem[1] = elem[2] = baseFieldElement(0);
	}

	fieldElement fieldElement::operator+(const fieldElement & other) const {
		if (isCounting) {
			++addCounter;
		}
		fieldElement ret;
		for (size_t i = 0; i < 3; i++) {
			ret.elem[i] = elem[i] + other.elem[i];
		}
		return ret;
	}

	fieldElement fieldElement::operator-(const fieldElement & other) const {
		if (isCounting) {
			++addCounter;
		}
		fieldElement ret;
		for (size_t i = 0; i < 3; i++) {
			ret.elem[i] = elem[i] - other.elem[i];
		}
		return ret;
	}

	fieldElement fieldElement::operator*(const fieldElement & other) const {
		if (isCounting) {
			++multCounter;
		}
		fieldElement ret;

		baseFieldElement v0, v1, v2, t0, t1;
		v0 = elem[0] * other.elem[0];
		v1 = elem[1] * other.elem[1];
		v2 = elem[2] * other.elem[2];

		t0 = elem[1] + elem[2];
		t1 = other.elem[1] + other.elem[2];
		t0 = t0 * t1 - v1 - v2;
		t0 = t0.mulNor();
		ret.elem[0] = t0 + v0;

		t0 = elem[0] + elem[1];
		t1 = other.elem[0] + other.elem[1];
		ret.elem[1] = t0 * t1 - v0 - v1;
		t0 = v2.mulNor();
		ret.elem[1] = ret.elem[1] + t0;

		t0 = elem[0] + elem[2];
		t1 = other.elem[0] + other.elem[2];
		ret. elem[2] = t0 * t1 - v0 + v1 - v2;

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
		baseFieldElement::init(prime, root);
	}

	unsigned int fieldElement::maxOrder() {
		return baseFieldElement::maxOrder();
	}

	fieldElement fieldElement::random() {
		fieldElement ret;
		for (size_t i = 0; i < 3; i++) {
			ret.elem[i] = baseFieldElement::random();
		}
		return ret;
	}

	bool fieldElement::operator!=(const fieldElement & other) const {
		return elem[0] != other.elem[0] || elem[1] != other.elem[1] || elem[2] != other.elem[2];
	}
	
	bool fieldElement::operator==(const fieldElement & other) const {
		return !(*this != other);
	}
	
	fieldElement & fieldElement::operator=(const fieldElement & other) {
		elem[0] = other.elem[0];
		elem[1] = other.elem[1];
		elem[2] = other.elem[2];
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
		return elem[0] || elem[1] || elem[2];
	}
	
	bool fieldElement::isNegative() const {
		assert(elem[1] == 0);
		assert(elem[2] == 0);
		return elem[0].isNegative();
	}
	
	unsigned char fieldElement::getBit(unsigned int i) const {
		assert(elem[1] == 0);
		assert(elem[2] == 0);
		return (elem[0].getBit(i));
	}
	
	bool fieldElement::operator<(const fieldElement & other) const {
		assert(elem[1] == 0);
		assert(elem[2] == 0);
		return elem[0] < other.elem[0];
	}
	
	bool fieldElement::isZero() {
		return elem[0].isZero() && elem[1].isZero() && elem[2].isZero();
	}

	fieldElement fieldElement::abs() const {
		fieldElement res = *this;
		return res.elem[0].abs();
	}
	
	fieldElement fieldElement::sqr() const {
		return (*this) * (*this);
	}
	
	fieldElement fieldElement::inv() const {
		baseFieldElement v0, v1, v2, t0;
		fieldElement ret;

		t0 = elem[0].sqr();
		v0 = elem[1] * elem[2];
		v2 = v0.mulNor();
		v0 = t0 - v2;

		t0 = elem[2].sqr();
		v2 = t0.mulNor();
		v1 = elem[0] * elem[1];
		v1 = v2 - v1;

		t0 = elem[1].sqr();
		v2 = elem[0] * elem[2];
		v2 = t0 - v2;

		t0 = elem[1] * v2;
		ret.elem[1] = t0.mulNor();

		ret.elem[0] = elem[0] * v0;

		t0 = elem[2] * v1;
		ret.elem[2] = t0.mulNor();

		t0 = ret.elem[0] + ret.elem[1] + ret.elem[2];
		t0 = t0.inv();

		ret.elem[0] = t0 * v0;
		ret.elem[1] = t0 * v1;
		ret.elem[2] = t0 * v2;
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
		for (size_t i = 0; i < 3; i++) {
			elem[i].print(fileno);
		}
	}
	
	fieldElement fieldElement::maxWithZero(const fieldElement & a,
			const fieldElement & b) {
		if (a.isNegative() && b.isNegative())
			return fieldElement::zero();
		return fieldElement(a.isNegative()? b : b.isNegative()? a : std::max(a.
						elem[0].elem[0], b.elem[0].elem[0]));
	}

	fieldElement fieldElement::maxUnsigned(const fieldElement & a,
			const fieldElement & b) {
		return a < b ? b : a;
	}

	fieldElement fieldElement::getRootOfUnity(int log_order) {
		return fieldElement(baseFieldElement::getRootOfUnity(log_order));
	}

	fieldElement fieldElement::zero() {
		return fieldElement(0);
	}

	fieldElement fieldElement::one() {
		return fieldElement(1);
	}

	vector <fieldElement> fieldElement::generateRandomness(u64 size) {
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
		char *s = new char[150];
		 sprintf(s, "(%llu, %llu, %llu, %llu, %llu, %llu)", this->elem[0].elem[0], this->elem[0].elem[1], this->elem[1].elem[0], this->elem[1].elem[1], this->elem[2].elem[0], this->elem[2].elem[1]);
		 return s;
	}
	
	ostream & operator<<(ostream & out, const fieldElement & c) {
		out << '(' << c.elem[0] << ',' << c.elem[1] << ',' << c.elem[2] << ')';
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
}
