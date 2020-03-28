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

/*
    Class for symbolic variable z, or c * z^p. Z has member variables:
    - coeff (c)
    - power (p)
 */

class Z
{
public:
    // Default constructor sets z to 0 * z^0
    Z (double coeff = 0, int power = 0) : coeff (coeff), power (power) {};
    ~Z() {};
    
    // Getters for coefficients
    const int& getPower() { return power; };
    const double& getCoeff() { return coeff; };
    
    // Set coefficient used when swapping (see Equation::swapCoeffs())
    void setCoeff (double val) { coeff = val; };
    
    // flip sign when multiplied onto with a z with a power with an opposite sign
    void flipSign() { power *= -1; };
    
    
    //// OPERATOR OVERLOADS ////
   
    void operator= (Z z) // z1 = z2;
    {
        coeff = z.getCoeff();
        power = z.getPower();
    }
    
    void operator+= (Z z) // z1 += z2;
    {
        // Assuming that 'this' and z always have the same power
        coeff += z.getCoeff();
    }
    
    void operator*= (const double& val) // z1 *= value;
    {
        coeff *= val;
    }
    
private:
    // in c * z^p, the following variables are
    double coeff;   // c, and
    int power;      // p
};




/*
    The Equation class holds a vector of Z's (see class above). This vector is structured according to the (absolute value of the) power of the z's it stores, i.e. zs[0] contains c * z^0, zs[1] contains c * z^(-1) or c * z^1 etc.
 */

template <int M> // M is used for allocating the size of the z-vector
class Equation
{
public:
    // The constructor fills the zs vector with 0 * z^(i*zSign). For example, if zSign = 1, zs = {0 * z^0, 0 * z^1, 0 * z^2, ...., 0 * z^M}
    
    Equation (int zSign) : zs (M+1, Z()), zSign (zSign)
    {
        for (int i = 0; i <= M; ++i)
        {
            zs[i] = Z (0, zSign * i);
        }
    };
    
    ~Equation() {};
    
    // Getter for the z-vector
    std::vector<Z>& getZs() { return zs; };
    
    /*
     Swaps coefficients in the z-vector around an amount determined by the power of the z multiplied onto the equation, specifically abs(p) / 2. In the recursion, a multiplication with z is always with a z of a higher (absolute) power than whatever equation it is multiplied onto. For example, the following could occur in the recursion:
            z^3 * (A_1(z^(-1), r)) = z^3 * (1 - rz^(-1)) = z^3 - rz^2
     So whatever z^3 was multiplied onto (in this case z^0 - rz^(-1)) had lower (absolute) powers than z^3 itself AND of opposite sign. The coefficient-vector of the equation {1, -r, 0, 0} can now be flipped around index abs(3) / 2 = 1.5, so between -r and 0. If this is done, the coefficient-vector of the solution looks like {0, 0, -r, 1}, which is shown in the equation above (i.e., 0z^0 + 0z^1 - rz^2 + 1z^3).
     */
    
    void swapCoeffs (int amount)
    {
        std::vector<double> curCoeffs = getCoeffs();
        for (int i = 0; i <= amount; ++i)
            zs[i].setCoeff (curCoeffs[amount-i]);

    }
    
    //// Functions used after recursion ////
    
    // Returns the coefficients in vector form
    std::vector<double> getCoeffs() {
        std::vector<double> coeffs (M+1, 0);
        
        for (int i = 0; i <= M; ++ i)
            coeffs[i] = zs[i].getCoeff();
        
        return coeffs;
    }
    
    // Debug: print the coefficients
    void printCoeffs()
    {
        for (int i = 0; i <= M; ++i)
            std::cout << zs[i].getCoeff() << " z^" << -i << std::endl;
    }

    //// OPERATOR OVERLOADS ////
    
    // Addition of two equations. Adds the coefficients of the z's in eq2 to the coefficients in eq1 according to the vector structure described above
    friend Equation operator+ (Equation eq1, Equation eq2) // eq1 + eq2;
    {
        for (int i = 0; i <= M; ++i)
        {
            eq1.getZs()[i] += eq2.getZs()[i];
        }
        return eq1;
    }
    
    // Addition of z and equation.
    friend Equation operator+ (Equation eq, Z z) // eq + z;
    {
        // The indexing here is possible because the z-vector is sorted according to the z-power.
        eq.getZs()[abs(z.getPower())] += z;
        return eq;
    }
    
    // Multiplication of double value and equation. Multiplies all coefficients in the z-vector with the value.
    friend Equation operator* (const double& val, Equation eq) // val * eq
    {
        for (int i = 0; i <= M; ++i)
        {
            eq.getZs()[i] *= val;
        }
        return eq;
    }
    
    // Multiplication of z and equation. Triggers swapCoeffs function
    friend Equation operator* (Z z, Equation eq) // z * eq
    {
        // i'm assuming that any multiplication (between Z and an equation) is with powers of different signs, where the power of Z is bigger than the biggest power in eq. I.e. it is going to flip all signs of the equation it is multiplied on
        eq.flipZSign();
        for (int i = 0; i <= M; ++i)
            eq.getZs()[i] *= z.getCoeff();
        eq.swapCoeffs (abs (z.getPower()));
        return eq;
    };
    
    int getZSign() { return zSign; };
    
    // Flip the signs of the z's stored
    void flipZSign() {
        zSign = -zSign;
        
        for (auto it = zs.begin(); it != zs.end(); ++it)
            it->flipSign();
    };
    
private:
    // Vector storing z's
    std::vector<Z> zs;
    
    // Stores what the current sign of the powers of the z's is
    int zSign;
};


//// RECURSION ////

/* The Muir-recursion:
    A_n(z^{-1}, r) = A_{n-1}(z^{-1},r) - c_nz^nA_{n-1}(z,r)
            { r/n,  if n is odd
    c_n  =  {
            { 0,    if n is even
    A_0 = 1
    The function takes in as arguments
    - n,
    - r,
    - the power of z (-1 in the first term after the equals-sign and 1 in the last term), and
    - the number of the recursion to keep track of how many there were
 */

template <int M>
Equation<M> A (int n, double r, int zPow, long& numRecursion)
{
    ++numRecursion; // increment recursion number
    
    // if n = 0 return 1 (or 1 * z^0)
    if (n == 0)
    {
        return Equation<M> (-zPow) + Z (1, 0);
    }
    
    // if n is even return A_{n-1}(z^{-1}, r)
    if (n % 2 == 0)
    {
        return A<M>(n-1, r, -zPow, numRecursion);
    }
    
    // if n is odd return A_{n-1}(z^{-1}, r) - r/n * z^n * A_{n-1}(z, r)
    else
    {
        return A<M>(n-1, r, -zPow, numRecursion) + Z (-r/n, n * zPow) * A<M>(n-1, r, zPow, numRecursion);
    }
};

int main (int argc, const char * argv[]) {

    const int M = 31;               // The starting number for n in the recursion (A_M(z^{-1},r))
    double r = 0.5;                 // The power of the fractional differentiation we want to approxiate
    int zPow = M % 2 == 1 ? -1 : 1; // Either start with z^1 (if M is even) or z^{-1} (if M is odd)
    long numRecursion = 0;          // The number of recursions it takes to solve
    
    //// SOLVE ////
    Equation<M> numerator = A<M> (M, r, zPow, numRecursion);
   
    // Print the coefficients
    numerator.printCoeffs();
    
    // Save the coefficients in a vector
    std::vector<double> coeffs = numerator.getCoeffs();
    
    // Print how many recursions it took
    std::cout << "Number of recursions: " << numRecursion << std::endl;
    
    return 0;
}
