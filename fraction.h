#ifndef FRACTION_H
#define FRACTION_H

#include <string>
using namespace std;

struct Fraction {
    int n, d;

    Fraction();
    Fraction(int n, int d);
    Fraction(int n);
    Fraction (const Fraction &o);

    void init(int n, int d);

    int gcd(int a, int b);
    float toFloat();
    string toString();

    Fraction& operator*=(const Fraction &o);
    Fraction& operator/=(const Fraction &o);
    Fraction& operator+=(const Fraction &o);
    Fraction& operator-=(const Fraction &o);
    bool operator==(const Fraction &o) const;
    bool operator!=(const Fraction &o) const;
    bool operator<(const Fraction &o) const;
    bool operator>(const Fraction &o) const;
    bool operator<=(const Fraction &o) const;
    bool operator>=(const Fraction &o) const;

    bool isPositive();
    bool isNegative();
    bool isZero();
    bool isInvalid();

};

Fraction operator*(const Fraction& o, const Fraction& o2);
Fraction operator/(const Fraction& o, const Fraction& o2);
Fraction operator+(const Fraction& o, const Fraction& o2);
Fraction operator-(const Fraction& o, const Fraction& o2);
ostream& operator<< (ostream& stream, Fraction& obj);

Fraction parseFraction(string s);

#endif