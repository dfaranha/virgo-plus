#pragma once

#include "constants.h"
#include "typedef.hpp"
#include <cassert>
#include <immintrin.h>
#include <vector>
#include <cstdio>
#include <iostream>
#include <helib/helib.h>
#include <helib/zeroValue.h>
#include "my_hhash.h"

using std::vector;
using std::ostream;

namespace virgo_ext {
    class baseFieldElementPacked;

    class baseFieldElement {
    public:
        baseFieldElement();

        baseFieldElement(const baseFieldElement &b);

        baseFieldElement(long long x);

        baseFieldElement(helib::Ctxt &x);

        baseFieldElement(long long x, long long y);

        baseFieldElement operator+(const baseFieldElement &other) const;

        baseFieldElement operator-(const baseFieldElement &other) const;

        baseFieldElement operator-() const;

        baseFieldElement operator*(const baseFieldElement &other) const;

        bool operator==(const baseFieldElement &other) const;

        bool operator!=(const baseFieldElement &other) const;

        baseFieldElement &operator=(const baseFieldElement &other);

        baseFieldElement &operator+=(const baseFieldElement &other);

        baseFieldElement &operator-=(const baseFieldElement &other);

        baseFieldElement &operator*=(const baseFieldElement &other);

        friend ostream &operator << (ostream &out, const baseFieldElement &c);

        bool operator < (const baseFieldElement &other) const;

        explicit operator bool () const;

        [[nodiscard]] bool isNegative() const;

        [[nodiscard]] unsigned char getBit(unsigned int i) const;

        [[nodiscard]] __int128_t toint128() const;

        [[nodiscard]] bool isZero();

        [[nodiscard]] baseFieldElement abs() const;
        [[nodiscard]] baseFieldElement sqr() const;
        [[nodiscard]] baseFieldElement inv() const;
        void setAbs();
        void setSqr();
        void setInv();

        void hash(void * buffer);

        void print(::FILE *fileno) const;
        char *toString() const;

        static void init(unsigned long long prime, unsigned long long root);
        static baseFieldElement maxWithZero(const baseFieldElement &a, const baseFieldElement &b);
        static baseFieldElement maxUnsigned(const baseFieldElement &a, const baseFieldElement &b);
        static baseFieldElement getRootOfUnity(int log_order); //return a root of unity with log_order 2^[log_order]
        static baseFieldElement random();

        static baseFieldElement zero();
        static baseFieldElement one();
        static vector<baseFieldElement> generateRandomness(u64 size);
        static baseFieldElement innerProd(vector<baseFieldElement>::iterator a, vector<baseFieldElement>::iterator b, u64 n);

        static baseFieldElement fastPow(baseFieldElement x, unsigned long long p);
		baseFieldElement mulNor();

		static unsigned int maxOrder();

        static bool initialized;
        static int multCounter, addCounter;
        static bool isCounting;
        static bool isSumchecking;

        unsigned long long elem[2];
        helib::Ctxt * elemHE[2];
        helib::Ptxt<helib::BGV> * elemCT[2];

        bool cleartext = true;
        void check_verifier_mode(const baseFieldElement & other) const;

        static double self_speed_test_mult(int repeat);
        static double self_speed_test_add(int repeat);
        static bool verifier_mode;
        static helib::SecKey * sk;
        static unsigned long long randomNumber();

    protected:
        static unsigned long long myMod(unsigned long long x);
        static unsigned long long mymult(const unsigned long long x, const unsigned long long y);
        static unsigned long long mod;
        static unsigned long long rou;
        static unsigned long long rcp;
        static unsigned int len;
        static unsigned int __max_order;

        friend baseFieldElementPacked;
    };
}














