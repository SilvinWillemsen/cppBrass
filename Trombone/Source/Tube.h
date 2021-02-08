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
    Tube (NamedValueSet& parameters, double k, std::vector<std::vector<double>>& geometry);
    ~Tube() override;

    Path drawGeometry (Graphics& g, int topOrBottom);
    Path visualiseState (Graphics& g, double visualScaling);
    void paint (juce::Graphics&) override;
    void resized() override;

    void calculateThermodynamicConstants();
    int calculateGeometry (std::vector<std::vector<double>>& geometry, NamedValueSet& parameters);
    void calculateRadii();
    void calculateVelocity();
    void calculatePressure();
    void calculateRadiation();

    void setFlowVelocities (double UbIn, double UrIn)
    {
        Ub = UbIn;
        Ur = UrIn;
    };
    float getOutput() { return getP (1, N-1); };
    
    void updateStates();
    
    double getP (int n, int l) {
        if (l <= M)
            return up[n][l];
        else
            return wp[n][l-M-1];
    };
    double getV (int n, int l) {
        if (l <= M-1)
            return uv[n][l];
        else
            return wv[n][l-M-1];
    };
    
    int getNint() { return Nint; };
    float getN() { return N; };

    int getM() { return M; };
    int getMw() { return Mw; };

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
    int Nint, M, Mw;
    int NnonExtended;
    float N;
    
    // Radiation vars
    double R1, rL, Lr, R2, Cr, z1, z2, z3, z4;
    double p1Next, p1, v1Next, v1;
    double oORadTerm;
    
    double Ub, Ur;
    
    double lambdaOverRhoC;
    std::vector<std::vector<double>> uvVecs;
    std::vector<std::vector<double>> upVecs;
    
    std::vector<std::vector<double>> wvVecs;
    std::vector<std::vector<double>> wpVecs;

    
    double upMP1, wpm1, uvNextMPh, uvMPh, wvNextmh, wvmh;

    // pointers to states
    std::vector<double*> uv;
    std::vector<double*> up;
    
    std::vector<double*> wv;
    std::vector<double*> wp;

    // states
    std::vector<std::vector<double>> uVecs;
    std::vector<std::vector<double>> wVecs;

    // tube geometry
    std::vector<double> S, SHalf, SBar, oOSBar, radii;
    
    double* uvTmp = nullptr;
    double* upTmp = nullptr;
    double* wvTmp = nullptr;
    double* wpTmp = nullptr;

    double kinEnergy1 = -1;
    double potEnergy1 = -1;
    double radEnergy1 = -1;

    bool raisedCos = false;
    bool init = true;
    
    double qHRadPrev = 0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Tube)
};
