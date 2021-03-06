/*
  ==============================================================================

    Trombone.h
    Created: 5 Sep 2020 1:12:46pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "Global.h"
#include "Tube.h"
#include "LipModel.h"
//==============================================================================
/*
*/
class Trombone  : public juce::Component
{
public:
    Trombone (NamedValueSet& parameters, double k, std::vector<std::vector<double>>& geometry);
    ~Trombone() override;

    void paint (juce::Graphics&) override;
    void resized() override;

    void calculate();
    void calculateEnergy();
    
    float getOutput() { return tube->getOutput(); };
    float getLipOutput() { return lipModel->getY(); };
    
    void saveToFiles();
    void closeFiles();
    void updateStates();

    void refreshLipModelInputParams() { lipModel->refreshInputParams(); };
    
    void setExtVals (double pVal, double lFVal, double LVal) {
        lipModel->setExtVals (pVal, lFVal);
        tube->setExtVals (LVal);
    };
private:
    std::unique_ptr<Tube> tube;
    std::unique_ptr<LipModel> lipModel;
    
    double k, Pm;
    
    double scaledTotEnergy = 0;
    
    std::ofstream massState, pState, vState, MSave, MwSave, alfSave, energySave, scaledTotEnergySave,
        maxMSave, maxMwSave, Ssave, output;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Trombone)
};
