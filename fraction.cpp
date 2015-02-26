#include <iostream>
#include <string>
#include <sstream>
#include "fraction.h"
using namespace std;

/*
Reference: Creating object arrays properly in C++.
http://www.devx.com/cplus/10MinuteSolution/30508/0/page/2
*/

Fraction::Fraction() {
    // Invalid
    this->n = 0;
    this->d = 0;
}
    
Fraction::Fraction(int n, int d) {
    init(n, d);
}

Fraction::Fraction(int n) {
    init(n, 1);
}

// Copy constructor
Fraction::Fraction (const Fraction &o) {
    if (o.d == 0) {
        // Invalid fraction copy.
        this->n = 0;
        this->d = 0;
        return;
    }
    init(o.n, o.d);
}

void Fraction::init(int n, int d) {
    if (d == 0) {
        cout << "\nERROR: Zero denominator: " << n << + "/" << d << "\n";
        return;
    }
    int g = gcd(n,d);
    n /= g;
    d /= g;
    if (d < 0) { // Denominators are strictly positive.
        d *= -1;
        n *= -1;
    }
    
    this->n = n;
    this->d = d;
}

float Fraction::toFloat() {
    return float(n)/d;
}

string Fraction::toString() {
    stringstream a;
    a << n;
    if (d != 1)
        a << "/" << d;
    return a.str();
}

int Fraction::gcd(int a, int b) {
    return a == 0 ? b : gcd(b%a, a);
}

    
Fraction operator*(const Fraction& o, const Fraction& o2) {
    return Fraction(o.n*o2.n, o.d*o2.d);
}

Fraction operator/(const Fraction& o, const Fraction& o2) {
    return Fraction(o.n*o2.d, o.d*o2.n);
}

Fraction operator+(const Fraction& o, const Fraction& o2) {
    return Fraction(o.n*o2.d + o2.n*o.d, o.d*o2.d);
}

Fraction operator-(const Fraction& o, const Fraction& o2) {
    return Fraction(o.n*o2.d - o2.n*o.d, o.d*o2.d);
}

ostream& operator<< (ostream& stream, Fraction& obj) {
    return stream << obj.toString();
}

Fraction& Fraction::operator*=(const Fraction& o) {
    init(n*o.n, d*o.d);
    return *this;
}

Fraction& Fraction::operator/=(const Fraction& o) {
    init(n*o.d, d*o.n);
    return *this;
}

Fraction& Fraction::operator+=(const Fraction& o) {
    init(n*o.d + o.n*d, d*o.d);
    return *this;
}

Fraction& Fraction::operator-=(const Fraction& o) {
    init(n*o.d - o.n*d, d*o.d);
    return *this;
}

bool Fraction::operator==(const Fraction &o) const {
    // simplest form is guaranteed unique.
    return (n == o.n) && (d == o.d);
}

bool Fraction::operator!=(const Fraction &o) const {
    return !(*this == o);
}

bool Fraction::operator<(const Fraction &o) const {
    // Note: denominators are POSITIVE.
    return o.d*n < o.n*d;
}

bool Fraction::operator>(const Fraction &o) const {
    return o.d*n > o.n*d;
}

bool Fraction::operator<=(const Fraction &o) const {
    return !(*this > o);
}

bool Fraction::operator>=(const Fraction &o) const {
    return !(*this < o);
}

bool Fraction::isPositive() {
    return n > 0;
}

bool Fraction::isNegative() {
    return n < 0;
}

bool Fraction::isZero() {
    return n == 0;
}

bool Fraction::isInvalid() {
    return d == 0;
}

Fraction parseFraction(string s) {
    int i=s.find_first_of('/');
    if (i == string::npos) {
        int n;
        stringstream(s) >> n;

        return Fraction(n);
    } else {
        int n, d;
        stringstream(s.substr(0,i)) >> n;
        stringstream(s.substr(i+1)) >> d;

        return Fraction(n,d);
    }
}