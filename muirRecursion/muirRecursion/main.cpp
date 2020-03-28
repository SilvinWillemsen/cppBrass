//
//  main.cpp
//  muirRecursion
//
//  Created by Silvin Willemsen on 26/03/2020.
//  Copyright © 2020 Silvin Willemsen. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>

template <class T>
static inline int sign (T val)
{
     return (0 < val) - (val < 0);
}

class Z
{
public:
    Z (double coeff = 0, int power = 0) : coeff (coeff), power (power) {};
    const int& getPower() { return power; };
    const double& getCoeff() { return coeff; };
    void setCoeff (double val) { coeff = val; };
    void flipSign() { power *= -1; };
    
    
    void operator= (Z z)
    {
        coeff = z.getCoeff();
        power = z.getPower();
    }
    
    friend Z operator+ (Z z1, Z z2)
    {
        if (z1.getPower() != z2.getPower())
        {
            std::cout << "not the same power" << std::endl;
            return Z(1000, 1000);
        }
        return Z (z1.getCoeff() + z2.getCoeff(), z1.getPower());
    }
    
    friend Z operator* (Z z1, Z z2)
    {
        return Z (z1.getCoeff() * z2.getCoeff(), z1.getPower() + z2.getPower());
    }
    
    friend Z operator* (const double val, Z z)
    {
        return Z (val * z.getCoeff(), z.getPower());
    }
        
    friend Z operator* (Z& z, const double val)
    {
        return val * z;
    }
    
    void operator^ (const int val)
    {
        power += val;
    }
private:
    double coeff;
    int power;
};

template <int M>
class Equation
{
public:
    Equation (int zSign) : zs (M+1, Z()), zSign (zSign)
    {
        for (int i = 0; i <= M; ++i)
        {
            zs[i] = Z (0, zSign * i);
        }
    };
    ~Equation() {};
    
//    Equation A (int n, double r)
//    {
//        if (n == 0)
//            zs[0] = 1;
//            return;
//        if (n % 2 == 0)
//            return A(n-1, r);
//        else
//        {
//            Z z (-r/n, n);
//            return A(n-1, r) + z * A(n-1, r);
//        }
//    };
    
    std::vector<Z>& getZs() { return zs; };
    void shiftCoeffs (int amount)
    {
        std::vector<double> curCoeffs (M+1, 0.0);
        for (int i = 0; i <= M; ++i)
        {
            curCoeffs[i] = zs[i].getCoeff();
        }
        // flip coefficients (this only works if the multiplication is always higher than the amount of non-zero coeffs
        for (int i = 0; i <= amount; ++i)
            zs[i].setCoeff (curCoeffs[amount-i]);
    }
    
    void printCoeffs()
    {
        for (int i = 0; i <= M; ++i)
            std::cout << zs[i].getCoeff() << " z^" << -i << std::endl;
        
    }
    
    friend Equation operator+ (Equation eq1, Equation eq2)
    {
//        std::vector<Z> zs1 = eq1.getZs();
//        std::vector<Z> zs2 = eq2.getZs();
//        Equation<M> tmp (eq1.getZSign()); //not sure if this is right
        for (int i = 0; i <= M; ++i)
        {
            eq1.getZs()[i] = eq1.getZs()[i] + eq2.getZs()[i];
//                if (zs1[i].getPower() == zs2[j].getPower())
//                {
//                    eq1.getZs()[i] * zs2[j].getCoeff();
//                }
        }
        return eq1;
    }
    
    friend Equation operator+ (Z z, Equation eq)
    {
        Equation<M> tmp = eq;
        tmp.getZs()[abs(z.getPower())] = eq.getZs()[abs(z.getPower())] + z; // assuming that the coeffs are sorted by power
        return tmp;
    }
    
    friend Equation operator+ (Equation eq, Z z)
    {
        return z + eq;
    }
    
    friend Equation operator* (const double val, Equation eq)
    {
        Equation<M> tmp = eq;
        for (int i = 0; i <= M; ++i)
        {
            tmp.getZs()[i] = tmp.getZs()[i] * val;
        }
        return tmp;
    }
    
    friend Equation operator* (Equation eq, const double val)
    {
        return val * eq;
    }
    
    friend Equation operator* (Z z, Equation eq)
    {
        if (sign (z.getPower()) != eq.getZSign()) // assuming any multiplication with powers of different signs is going to flip all signs of the equation
        {
            eq.flipZSign();
            for (int i = 0; i <= M; ++i)
                eq.getZs()[i] = eq.getZs()[i] * z.getCoeff();
            eq.shiftCoeffs (abs (z.getPower()));
        }
//
//        for (int i = 0; i <= M; ++i)
//        {
//            eq.shiftCoeffs (z.getPower());
////            eq.getZs()[i]^(z.getPower());
//        }
        return eq;
    };
    
    friend Equation operator* (Equation eq, Z z)
    {
        return z * eq;
    };

    int getZSign() { return zSign; };
    void flipZSign() {
        zSign = -1 * zSign;
        for (int i = 0; i <= M; ++i)
            zs[i].flipSign();
    };
    
private:
    std::vector<Z> zs;
    int zSign;
};

//// RECURSION ////

template <int M>
Equation<M> A (int n, double r, int zPow, long& numRecursion)
{
    ++numRecursion;
    if (n == 0)
    {
        return Equation<M> (-1 * zPow) + Z(1, 0);
    }
    if (n % 2 == 0)
    {
        return A<M>(n-1, r, -1 * zPow, numRecursion);
    }
    else
    {
        return A<M>(n-1, r, -1 * zPow, numRecursion) + Z (-r/n, n * zPow) * A<M>(n-1, r, zPow, numRecursion);
    }
};

int main(int argc, const char * argv[]) {

    
    const int M = 21;
    double r = 0.5;
    int zPow = M % 2 == 1 ? -1 : 1;
    long numRecursion = 0;
    Equation<M> numerator = pow (2.0 / 0.001, r) * A<M> (M, r, zPow, numRecursion);
//    Equation<M> denominator = pow (2.0 / 0.001, r) * A<M> (M, -r, zPow, numRecursion);
   
    
    numerator.printCoeffs();
//    denominator.printCoeffs();
    
    std::cout << "Number of recursions: " << numRecursion << std::endl;
    return 0;
}
