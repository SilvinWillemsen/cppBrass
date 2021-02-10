/*
  ==============================================================================

    Trombone.cpp
    Created: 5 Sep 2020 1:12:46pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#include <JuceHeader.h>
#include "Trombone.h"

//==============================================================================
Trombone::Trombone (NamedValueSet& parameters, double k, std::vector<std::vector<double>>& geometry) : k (k),
Pm (*parameters.getVarPointer ("Pm"))
{
    // In your constructor, you should add any child components, and
    // initialise any special settings that your component needs.
    
    tube = std::make_unique<Tube> (parameters, k, geometry);
    addAndMakeVisible (tube.get());
    lipModel = std::make_unique<LipModel> (parameters, k);
    lipModel->setTubeParameters (tube->getH(),
                                 tube->getRho(),
                                 tube->getC(),
                                 tube->getSBar(0),
                                 tube->getSHalf(0));
    addAndMakeVisible (lipModel.get());

    massState.open ("massState.csv");
    pState.open ("pState.csv");
    vState.open ("vState.csv");
    alfSave.open ("alfSave.csv");
    MSave.open ("MSave.csv");
    MwSave.open ("MwSave.csv");
    energySave.open ("energySave.csv");
    scaledTotEnergySave.open ("scaledTotEnergySave.csv");
    
    maxMSave.open ("maxM.csv");
    maxMwSave.open ("maxMw.csv");

    Ssave.open ("SSave.csv");
    output.open ("output.csv");

    maxMSave << tube->getMaxM();
    maxMwSave << tube->getMaxMw();

    maxMSave.close();
    maxMwSave.close();

    
}

Trombone::~Trombone()
{
    closeFiles();
}

void Trombone::paint (juce::Graphics& g)
{
    /* This demo code just fills the component's background and
       draws some placeholder text to get you started.

       You should replace everything in this method with your own
       drawing code..
    */
}

void Trombone::resized()
{
    // This method is where you should set the bounds of any child
    // components that your component contains..
//    Rectangle<int> totArea = getLocalBounds();
//    tube->setBounds (totArea.removeFromTop(getHeight));
//    lipModel->setBounds (totArea);
    tube->setBounds (getLocalBounds());

}

void Trombone::calculate()
{
    if (!Global::fixedNonInterpolatedL)
        tube->updateL();
    tube->calculateVelocity();
    lipModel->setTubeStates (tube->getP (1, 0), tube->getV (0, 0));
    lipModel->calculateCollision();
    lipModel->calculateDeltaP();
    lipModel->calculate();
    tube->setFlowVelocities (lipModel->getUb(), lipModel->getUr());
    tube->calculatePressure();
    tube->calculateRadiation();
}

void Trombone::calculateEnergy()
{
    bool excludeLip = !Global::connectedToLip;
//    bool excludeLip = false;

    double kinEnergy = tube->getKinEnergy();
    double potEnergy = tube->getPotEnergy();
    double radEnergy = tube->getRadEnergy();
    double radDamp = tube->getRadDampEnergy();
    double lipEnergy = lipModel->getLipEnergy();
    double lipCollisionEnergy = lipModel->getCollisionEnergy();
    double lipPower = lipModel->getPower();
    double lipDamp = lipModel->getDampEnergy();
    
    double totEnergy = kinEnergy + potEnergy + radEnergy + (excludeLip ? 0 : (lipEnergy + lipCollisionEnergy));
    double energy1 = tube->getKinEnergy1() + tube->getPotEnergy1() + tube->getRadEnergy1() + (excludeLip ? 0 : (lipModel->getLipEnergy1() + lipModel->getCollisionEnergy1()));
    
    energySave << kinEnergy << ", ";
    energySave << potEnergy << ", ";
    energySave << radEnergy << ", ";
    energySave << radDamp << ", ";
    energySave << lipEnergy << ", ";
    energySave << lipCollisionEnergy << ", ";
    energySave << lipPower << ", ";
    energySave << lipDamp << ", ";
    energySave << energy1 << ";\n";

//     scaledTotEnergy = (totEnergy + lipModel->getPower() + lipModel->getDampEnergy() + tube->getRadDampEnergy() - energy1) / energy1;
    scaledTotEnergy = (totEnergy + lipPower + lipDamp + radDamp - energy1) / pow(2, floor (log2 (energy1)));
    scaledTotEnergySave << scaledTotEnergy << ";\n";
}

void Trombone::updateStates()
{
    tube->updateStates();
    lipModel->updateStates();
}

void Trombone::saveToFiles()
{
    if (!Global::onlyWriteOutput)
    {
        massState << getLipOutput() << ";\n";
        
        for (int l = 0; l <= tube->getMaxN() + 1; ++l)
        {
            pState << tube->getP (1, l) << ", ";
            if (l < tube->getMaxN())
                vState << tube->getV (1, l) << ", ";
        }
        pState << ";\n";
        vState << ";\n";
        alfSave << tube->getAlf() << ";\n";
        MSave << tube->getM() << ";\n";
        MwSave << tube->getMw() << ";\n";
        
        for (int l = 0; l <= tube->getMaxN() + 1; ++l)
        {
            Ssave << tube->getS (l);
            if (l == tube->getMaxN() + 1)
                Ssave << ";\n";
            else
                Ssave << ",";
        }
    }
    output << tube->getOutput() << ";\n";
    
}

void Trombone::closeFiles()
{
    massState.close();
    pState.close();
    vState.close();
    alfSave.close();
    MSave.close();
    MwSave.close();
    energySave.close();
    scaledTotEnergySave.close();
    Ssave.close();
    output.close();
    tube->closeFiles();

}
