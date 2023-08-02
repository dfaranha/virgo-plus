#pragma once

#include <gmp.h>
#include "constants.h"
#include "typedef.hpp"
#include <cassert>
#include <immintrin.h>
#include <vector>
#include <cstdio>
#include <iostream>

using std::vector;
using std::ostream;

#define MOD_LIMBS	2

namespace virgo {
    class fieldElementPacked;

    class fieldElement {
    public:
        fieldElement();

        fieldElement(const fieldElement &b);

        fieldElement(long long x);

        fieldElement(long long x, long long y);

        fieldElement operator+(const fieldElement &other) const;

        fieldElement operator-(const fieldElement &other) const;

        fieldElement operator-() const;

        fieldElement operator*(const fieldElement &other) const;

        bool operator==(const fieldElement &other) const;

        bool operator!=(const fieldElement &other) const;

        fieldElement &operator=(const fieldElement &other);

        fieldElement &operator+=(const fieldElement &other);

        fieldElement &operator-=(const fieldElement &other);

        fieldElement &operator*=(const fieldElement &other);

        friend ostream &operator << (ostream &out, const fieldElement &c);

        bool operator < (const fieldElement &other) const;

        explicit operator bool () const;

        [[nodiscard]] bool isNegative() const;

        [[nodiscard]] unsigned char getBitWidth() const;

        [[nodiscard]] unsigned char getBit(unsigned int i) const;

        [[nodiscard]] __int128_t toint128() const;

        [[nodiscard]] bool isZero();

        [[nodiscard]] fieldElement abs() const;
        [[nodiscard]] fieldElement sqr() const;
        [[nodiscard]] fieldElement inv() const;
        void setAbs();
        void setSqr();
        void setInv();

        void print(::FILE *fileno) const;
        char *toString() const;

        static void init(unsigned long long prime, unsigned long long root);
		static void init(const unsigned char *prime, const unsigned char *root, int size);
        static fieldElement maxWithZero(const fieldElement &a, const fieldElement &b);
        static fieldElement maxUnsigned(const fieldElement &a, const fieldElement &b);
        static fieldElement getRootOfUnity(int log_order); //return a root of unity with log_order 2^[log_order]
        static fieldElement random();

        static fieldElement zero();
        static fieldElement one();
        static vector<fieldElement> generateRandomness(u64 size);
        static fieldElement innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n);

        static fieldElement fastPow(fieldElement x, unsigned long long p);

		static unsigned int maxOrder();

        static bool initialized;
        static int multCounter, addCounter;
        static bool isCounting;
        static bool isSumchecking;
        static mp_limb_t mod[MOD_LIMBS];
		static mp_limb_t rou[MOD_LIMBS];
		static unsigned int len;
        static unsigned int __max_order;

        mp_limb_t elem[MOD_LIMBS];

        static double self_speed_test_mult(int repeat);
        static double self_speed_test_add(int repeat);

    protected:
        static unsigned long long myMod(unsigned long long x);
        static unsigned long long mymult(const unsigned long long x, const unsigned long long y);
        static unsigned long long randomNumber();

        friend fieldElementPacked;
    };
}