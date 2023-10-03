#include "polynomial.h"

namespace virgo_ext {
    quintuple_poly::quintuple_poly() {}

    quintuple_poly::quintuple_poly(const virgo_ext::fieldElement &aa, const virgo_ext::fieldElement &bb,
                                   const virgo_ext::fieldElement &cc, const virgo_ext::fieldElement &dd,
                                   const virgo_ext::fieldElement &ee, const virgo_ext::fieldElement &ff) {
        a = aa;
        b = bb;
        c = cc;
        d = dd;
        e = ee;
        f = ff;
    }

    quintuple_poly quintuple_poly::operator+(const quintuple_poly &x) const {
        return quintuple_poly(a + x.a, b + x.b, c + x.c, d + x.d, e + x.e, f + x.f);
    }

    virgo_ext::fieldElement quintuple_poly::eval(const virgo_ext::fieldElement &x) const {
        return (((((a * x) + b) * x + c) * x + d) * x + e) * x + f;
    }

    quadruple_poly::quadruple_poly() {}

    quadruple_poly::quadruple_poly(const virgo_ext::fieldElement &aa, const virgo_ext::fieldElement &bb,
                                   const virgo_ext::fieldElement &cc, const virgo_ext::fieldElement &dd,
                                   const virgo_ext::fieldElement &ee) {
        a = aa;
        b = bb;
        c = cc;
        d = dd;
        e = ee;
    }

    quadruple_poly quadruple_poly::operator+(const quadruple_poly &x) const {
        return quadruple_poly(a + x.a, b + x.b, c + x.c, d + x.d, e + x.e);
    }

    virgo_ext::fieldElement quadruple_poly::eval(const virgo_ext::fieldElement &x) const {
        return ((((a * x) + b) * x + c) * x + d) * x + e;
    }

    cubic_poly::cubic_poly() {}

    cubic_poly::cubic_poly(const virgo_ext::fieldElement &aa, const virgo_ext::fieldElement &bb, const virgo_ext::fieldElement &cc,
                           const virgo_ext::fieldElement &dd) {
        a = aa;
        b = bb;
        c = cc;
        d = dd;
    }

    cubic_poly cubic_poly::operator+(const cubic_poly &x) const {
        return cubic_poly(a + x.a, b + x.b, c + x.c, d + x.d);
    }

    virgo_ext::fieldElement cubic_poly::eval(const virgo_ext::fieldElement &x) const {
        return (((a * x) + b) * x + c) * x + d;
    }


    quadratic_poly::quadratic_poly() {}

    quadratic_poly::quadratic_poly(const virgo_ext::fieldElement &aa, const virgo_ext::fieldElement &bb,
                                   const virgo_ext::fieldElement &cc) {
        a = aa;
        b = bb;
        c = cc;
    }

    quadratic_poly quadratic_poly::operator+(const quadratic_poly &x) const {
        return quadratic_poly(a + x.a, b + x.b, c + x.c);
    }

    quadratic_poly quadratic_poly::operator+(const linear_poly &x) const {
        return quadratic_poly(a, b + x.a, c + x.b);
    }

    cubic_poly quadratic_poly::operator*(const linear_poly &x) const {
        return cubic_poly(a * x.a, a * x.b + b * x.a, b * x.b + c * x.a, c * x.b);
    }

    quadratic_poly quadratic_poly::operator*(const virgo_ext::fieldElement &x) const {
        return quadratic_poly(a * x, b * x, c * x);
    }

    virgo_ext::fieldElement quadratic_poly::eval(const virgo_ext::fieldElement &x) const {
        return ((a * x) + b) * x + c;
    }


    linear_poly::linear_poly() {
        a = fieldElement::zero();
        b = fieldElement::zero();
    }

    linear_poly::linear_poly(const virgo_ext::fieldElement &aa, const virgo_ext::fieldElement &bb) {
        a = aa;
        b = bb;
    }

    linear_poly::linear_poly(const virgo_ext::fieldElement &x) {
        a = fieldElement::zero();
        b = x;
    }

    linear_poly linear_poly::operator+(const linear_poly &x) const {
        return linear_poly(a + x.a, b + x.b);
    }

    quadratic_poly linear_poly::operator*(const linear_poly &x) const {
        return quadratic_poly(a * x.a, a * x.b + b * x.a, b * x.b);
    }

    linear_poly linear_poly::operator*(const virgo_ext::fieldElement &x) const {
        return linear_poly(a * x, b * x);
    }

    virgo_ext::fieldElement linear_poly::eval(const virgo_ext::fieldElement &x) const {
        return a * x + b;
    }
}