/*
  ==============================================================================

    Tube.cpp
    Created: 5 Sep 2020 1:11:57pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Tube.h"

//==============================================================================
Tube::Tube (NamedValueSet& parameters, double k) : k (k),
L (*parameters.getVarPointer("L")),
T (*parameters.getVarPointer("T"))
{
    calculateThermodynamicConstants();
    
    h = c * k;
    N = floor(L / h);
    h = L / static_cast<double> (N);
    NnonExtended = N; // assuming that the initialisation happens with L = 2.658
    lambda = c * k / h;
    lambdaOverRhoC = lambda / (rho * c);
    
    // initialise state vectors
    vVecs.resize (2);
    pVecs.resize (2);

    for (int i = 0; i < 2; ++i)
    {
        vVecs[i] = std::vector<double> (N, 0);
        pVecs[i] = std::vector<double> (N, 0);
    }
    
    if (raisedCos)
    {
        int start = N * 0.5 - 5;
        int end = N * 0.5 + 5;
        double scaling = 100000.0;
        for (int n = 0; n < 2; ++n)
        {
            for (int l = start; l < end; ++l)
            {
                pVecs[n][l] = scaling * (1.0 - cos (2.0 * double_Pi * (l-start) / static_cast<float>(end - start))) * 0.5;
            }
        }
    }
    
    v.resize (2);
    p.resize (2);

    for (int i = 0; i < 2; ++i)
    {
        v[i] = &vVecs[i][0];
        p[i] = &pVecs[i][0];
    }
    
    calculateGeometry (parameters);
//    calculateRadii();
    
    // Radiation
    R1 = rho * c;
    rL = rLVec[N-1];
    Lr = 0.613 * rho * rL;
    R2 = 0.505 * rho * c;
    Cr = 1.111 * rL / (rho * c * c);
    
    double zDiv = 2.0 * R1 * R2 * Cr + k * (R1 + R2);
    if (zDiv == 0)
    {
        z1 = 0;
        z2 = 0;
    } else {
        z1 = 2 * R2 * k / zDiv;
        z2 = (2 * R1 * R2 * Cr - k * (R1 + R2)) / zDiv;
    }
    
    z3 = k / (2.0 *Lr) + z1 / (2.0 * R2) + Cr * z1 / k;
    z4 = (z2 + 1.0) / (2.0 * R2) + (Cr * z2 - Cr) / k;
    
    oORadTerm = 1.0 / (1.0 + rho * c * lambda * z3);
    
    p1 = 0;
    v1 = 0;
}

Tube::~Tube()
{
}

void Tube::paint (juce::Graphics& g)
{
    /* This demo code just fills the component's background and
       draws some placeholder text to get you started.

       You should replace everything in this method with your own
       drawing code..
    */

    g.setColour (Colours::gold);
    drawGeometry (g);
//        Path stringPathBottom = drawGeometry (g, 1);
//        g.fillPath (stringPathTop);
//        g.strokePath (stringPathBottom, PathStrokeType(2.0f));
//    }
    g.setColour (Colours::cyan);
    Path state = visualiseState (g, 0.01 * Global::oOPressureMultiplier);
    g.strokePath (state, PathStrokeType (2.0f));

    
//    std::cout << "repainted" << std::endl;

}

void Tube::drawGeometry (Graphics& g)
{
    double visualScaling = 1000.0;
    double colourScaling = 1.0 / 10.0;

    Path stringPath;
    stringPath.startNewSubPath (0, -rLVec[0] * visualScaling + getHeight() * 0.5);
    auto spacing = getWidth() / static_cast<double>(N - 1);
    auto x = spacing;
    auto oONm1 = 1.0 / static_cast<double>(N);
    auto xRat = oONm1;
    int flag = 0;
    int topOrBottom = 1;
    ColourGradient cg (Colour (Global::clamp(255 + p[1][0] * colourScaling, 0, 255),
                               Global::clamp(255 - abs(p[1][0] * colourScaling), 0, 255),
                               Global::clamp(255 - p[1][0] * colourScaling, 0, 255)), getX(), getY(),
                       Colours::white, getWidth(), getY(),
                       false);

    
    for (int y = 1; y < 2 * N; y++)
    {
        if (y < N)
        {
            topOrBottom = 1;
        } else {
            flag = 1;
            topOrBottom = -1;
        }
        
        stringPath.lineTo (x, -topOrBottom * rLVec[2 * (N - 0.5) * flag + topOrBottom * y] * visualScaling + getHeight() * 0.5);
        x += topOrBottom * spacing;

        if (y < N)
        {
            cg.addColour (xRat, Colour (Global::clamp (255 + p[1][y] * colourScaling, 0, 255),
                                        Global::clamp (255 - abs(p[1][y] * colourScaling), 0, 255),
                                        Global::clamp (255 - p[1][y] * colourScaling, 0, 255)));
            xRat += oONm1;
        }
            
    }
    stringPath.lineTo (0, rLVec[0] * visualScaling + getHeight() * 0.5);
    stringPath.closeSubPath();

    g.setColour (Colours::gold);
    g.strokePath (stringPath, PathStrokeType(2.0f));
    g.setGradientFill (cg);
    g.fillPath (stringPath);
    
}

Path Tube::visualiseState (Graphics& g, double visualScaling)
{
    auto stringBounds = getHeight() / 2.0;
    Path stringPath;
    stringPath.startNewSubPath (0, -p[1][0] * visualScaling + stringBounds);

    auto spacing = getWidth() / static_cast<double>(N - 1);
    auto x = spacing;
    
    for (int y = 1; y < N; y++)
    {
        float newY = -p[1][y] * visualScaling + stringBounds; // Needs to be -p, because a positive p would visually go down
        
        if (isnan(x) || isinf(abs(x) || isnan(p[1][y]) || isinf(abs(p[1][y]))))
        {
            std::cout << "Wait" << std::endl;
        };
        
        if (isnan(newY))
            newY = 0;
        stringPath.lineTo(x, newY);
        //        g.drawEllipse(x, newY, 2, 2, 5);
        x += spacing;
    }
    return stringPath;
}

void Tube::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..

}

void Tube::calculateThermodynamicConstants()
{
    double deltaT = T - 26.85;
    c = 3.4723e2 * (1 + 0.00166 * deltaT);      // Speed of sound in air [m/s]
    rho = 1.1769 * (1 - 0.00335 * deltaT);      // Density of air [kg·m^{-3}]
//    eta = 1.846 * (1 + 0.0025 * deltaT);        // Shear viscosity [kg·s^{-1}·m^{-1}]
//    nu = 0.8410 * (1 - 0.0002 * deltaT);        // Root of Prandtl number [-]
//    gamma = 1.4017 * (1 - 0.00002 * deltaT);    // Ratio of specific heats [-]
    
}

void Tube::calculateVelocity()
{
    for (int l = 0; l < N - 1; ++l)
        v[0][l] = v[1][l] - lambdaOverRhoC * (p[1][l+1] - p[1][l]);
}

void Tube::calculatePressure()
{
    for (int l = 1; l < N - 1; ++l)
        p[0][l] = p[1][l] - rho * c * lambda * oOSBar[l] * (SHalf[l] * v[0][l] - SHalf[l-1] * v[0][l-1]);
    
    p[0][0] = p[1][0] - rho * c * lambda * oOSBar[0] * (-2.0 * (Ub + Ur) + 2.0 * SHalf[0] * v[0][0]);
}

void Tube::calculateRadiation()
{
    p[0][N-1] = ((1.0 - rho * c * lambda * z3) * p[1][N-1] - 2.0 * rho * c * lambda * (v1 + z4 * p1 - (SHalf[N-2] * v[0][N-2]) * oOSBar[N-1])) * oORadTerm;

    v1Next = v1 + k / (2.0 * Lr) * (p[0][N-1] + p[1][N-1]);
    p1Next = z1 * 0.5 * (p[0][N-1] + p[1][N-1]) + z2 * p1;
    
}

void Tube::updateStates()
{
    vTmp = v[1];
    v[1] = v[0];
    v[0] = vTmp;
    
    pTmp = p[1];
    p[1] = p[0];
    p[0] = pTmp;
    
    p1 = p1Next;
    v1 = v1Next;
}

void Tube::calculateGeometry (NamedValueSet& parameters)
{
    S.resize (N, 0);
    SHalf.resize (N-1, 0);
    SBar.resize (N, 0);
    oOSBar.resize (N, 0);
    rLVec.resize (N, 0);
    std::vector<double> lengths;
    std::vector<double> lengthN;
    std::vector<double> radii;
    
    radii.push_back (*parameters.getVarPointer ("inner1Rad")); // Rad first inner tube
    radii.push_back (*parameters.getVarPointer ("slideRad"));  // Rad slide
    radii.push_back (*parameters.getVarPointer ("inner2Rad")); // Rad second inner tube
    radii.push_back (*parameters.getVarPointer ("gooseNeckRad")); // Rad gooseneck
    radii.push_back (*parameters.getVarPointer ("tuningRad1")); // First rad tuning slide
    radii.push_back (*parameters.getVarPointer ("tuningRad2")); // Second rad tuning slide
    
    lengths.push_back (*parameters.getVarPointer ("inner1L")); // Length first inner tube
    lengths.push_back (*parameters.getVarPointer ("slideL"));  // Slide length (non-extended)
    lengths.push_back (*parameters.getVarPointer ("inner2L")); // Length first inner tube
    lengths.push_back (*parameters.getVarPointer ("gooseNeckL")); // Rad gooseneck
    lengths.push_back (*parameters.getVarPointer ("tuningL")); // First rad tuning slide
    lengths.push_back (*parameters.getVarPointer ("bellL")); // Second rad tuning slide
    
    double sumLengths = 0;
    for (int i = 0; i < lengths.size(); ++i)
    {
        sumLengths += lengths[i];
    }
    
//    for (int i = 0; i < radii.size(); ++i)
//    {
//        radii[i] = radii[i] * radii[i] * double_Pi; // convert to S
//    }
    
    for (int i = 0; i < lengths.size(); ++i)
    {
        lengthN.push_back (round ((NnonExtended * lengths[i]) / sumLengths));
    }
    
    lengthN[1] += (N - NnonExtended);
    lengthN[5] += 1;
    
    double flare = *parameters.getVarPointer ("flare");
    double x0 = *parameters.getVarPointer ("x0");
    double b = *parameters.getVarPointer ("b");
    int idx = 0;
    int j = 0;
    while (j < lengths.size())
    {
        for (int i = 0; i < lengthN[j]; ++i)
        {
            if (j < 4)
                rLVec[idx] = radii[j];
            else if (j == 4)
                rLVec[idx] = Global::linspace (radii[j], radii[j+1], lengthN[j], i);
            else
                rLVec[N - 1 - i] = b * pow((i / (lengthN[j]) * lengths[j] + x0), -flare);
            
            ++idx;
        }
        ++j;
    }

    for (int i = 0; i < N; ++i)
        S[i] = rLVec[i] * rLVec[i] * double_Pi;
    
    for (int i = 0; i < N - 1; ++i)
        SHalf[i] = (S[i] + S[i+1]) * 0.5;

    SBar[0] = S[0];
    
    for (int i = 0; i < N - 2; ++i)
        SBar[i+1] = (SHalf[i] + SHalf[i+1]) * 0.5;

    SBar[N-1] = S[N-1];
    for (int i = 0; i < N; ++i)
        oOSBar[i] = 1.0 / SBar[i];
}

//void Tube::calculateRadii()
//{
//    radii.resize (N, 0);
//    for (int i = 0; i < N; ++i)
//        radii[i] = sqrt (S[i]) / double_Pi;
//  
//}

double Tube::getKinEnergy()
{
    double kinEnergy = 0;
    for (int i = 0; i < N-1; ++i)
        kinEnergy += rho * 0.5 * h * (SHalf[i] * v[0][i] * v[1][i]);
    
    if (kinEnergy1 < 0)
        kinEnergy1 = kinEnergy;
    
    return kinEnergy;
}


double Tube::getPotEnergy()
{
    double potEnergy = 0;
    for (int i = 0; i < N; ++i)
    {
        potEnergy += 1.0 / (2.0 * rho * c * c) * h * (SBar[i] * p[1][i] * p[1][i] * (i == 0 || i == N-1 ? 0.5 : 1));
    }
    if (potEnergy1 < 0)
        potEnergy1 = potEnergy;
    return potEnergy;

}

double Tube::getRadEnergy()
{
    double radEnergy = SBar[N-1] / 2.0 * (Lr * v1 * v1 + Cr * p1 * p1);
    
    if (radEnergy1 < 0)
        radEnergy1 = radEnergy;
    
    return radEnergy;
}

double Tube::getRadDampEnergy()
{
    double pBar = 0.5 * (p[0][N-1] + p[1][N-1]);
    double muTPv2 = (pBar - 0.5 * (p1Next + p1)) / R1;
    
    double qRad = SBar[N-1] * (R1 * muTPv2 * muTPv2 + R2 * (0.5 * ((p1Next + p1) / R2) * (0.5 * ((p1Next + p1) / R2))));
    double qHRad = k * qRad + qHRadPrev;

    double qHRadPrevTmp = qHRadPrev;
    qHRadPrev = qHRad;
    return qHRadPrevTmp;
                                                            

}
