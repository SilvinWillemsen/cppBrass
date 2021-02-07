/*
  ==============================================================================

    Tube.h
    Created: 5 Sep 2020 1:11:57pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Global.h"

//==============================================================================
/*
*/
class Tube  : public juce::Component
{
public:
    Tube (NamedValueSet& parameters, double k);
    ~Tube() override;

    void drawGeometry (Graphics& g);
    Path visualiseState (Graphics& g, double visualScaling);
    void paint (juce::Graphics&) override;
    void resized() override;

    void calculateThermodynamicConstants();
    void calculateGeometry (NamedValueSet& parameters);
//    void calculateRadii();
    void calculateVelocity();
    void calculatePressure();
    void calculateRadiation();

    void setFlowVelocities (double UbIn, double UrIn) { Ub = UbIn; Ur = UrIn; };
    float getOutput() { return getP (1, N-1); };
    void updateStates();
    
    double getP (int n, int l) { return p[n][l]; };
    double getV (int n, int l) { return v[n][l]; };

    double getH() { return h; };
    double getRho() { return rho; };
    double getC() { return c; };

    double getS (int idx) { return S[idx]; };
    double getSHalf (int idx) { return SHalf[idx]; };
    double getSBar (int idx) { return SBar[idx]; };
    
    double getKinEnergy();
    double getPotEnergy();
    double getRadEnergy();
    double getRadDampEnergy();
    
    double getKinEnergy1() { return kinEnergy1; };
    double getPotEnergy1() { return potEnergy1; };
    double getRadEnergy1() { return radEnergy1; };

private:
    double k, h, c, lambda, rho, L, T;
    int N, NnonExtended;
    
    // Radiation vars
    double R1, rL, Lr, R2, Cr, z1, z2, z3, z4;
    double p1Next, p1, v1Next, v1;
    double oORadTerm;
    
    double Ub, Ur;
    
    double lambdaOverRhoC;
    std::vector<std::vector<double>> vVecs;
    std::vector<std::vector<double>> pVecs;

    // pointers to states
    std::vector<double*> v;
    std::vector<double*> p;

    // states
    std::vector<std::vector<double>> uVecs;
    
    // tube geometry
    std::vector<double> S, SHalf, SBar, oOSBar, rLVec;
    
    double* vTmp = nullptr;
    double* pTmp = nullptr;
    
    double kinEnergy1, potEnergy1, radEnergy1 = -1;

    bool raisedCos = false;
    bool init = true;
    
    double qHRadPrev = 0;

    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Tube)
};
