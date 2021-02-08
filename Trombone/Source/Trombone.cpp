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
    MSave.open ("MSave.csv");
    MwSave.open ("MwSave.csv");
    energySave.open ("energySave.csv");

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
    Rectangle<int> totArea = getLocalBounds();
    tube->setBounds (totArea.removeFromTop (getHeight() * 0.5));
    lipModel->setBounds (totArea);

}

void Trombone::calculate()
{
    tube->calculateVelocity();
    lipModel->setTubeStates (tube->getP (1, 0), tube->getV (0, 0));
    lipModel->calculateCollision();
    lipModel->calculateDeltaP();
    lipModel->calculate();
    tube->setFlowVelocities (lipModel->getUb(), lipModel->getUr());
    tube->calculatePressure();
    tube->calculateRadiation();
    
    calculateEnergy();
}

void Trombone::calculateEnergy()
{
    bool excludeLip = true;
    double totEnergy = tube->getKinEnergy() + tube->getPotEnergy() + tube->getRadEnergy() + (excludeLip ? 0 : (lipModel->getLipEnergy() + lipModel->getCollisionEnergy()));
    double energy1 = tube->getKinEnergy1() + tube->getPotEnergy1() + tube->getRadEnergy1() + (excludeLip ? 0 : (lipModel->getLipEnergy1() + lipModel->getCollisionEnergy1()));
    
    scaledTotEnergy = (totEnergy + lipModel->getPower() + lipModel->getDampEnergy() + tube->getRadDampEnergy() - energy1) / energy1;
//    std::cout << scaledTotEnergy << std::endl;
}

void Trombone::updateStates()
{
    tube->updateStates();
    lipModel->updateStates();
}

void Trombone::saveToFiles()
{
    massState << getLipOutput() << ";\n";
    
    for (int l = 0; l <= tube->getNint() + 1; ++l)
    {
        pState << tube->getP (1, l) << ", ";
        if (l < tube->getNint())
            vState << tube->getV (1, l) << ", ";
    }
    pState << ";\n";
    vState << ";\n";
    
    MSave << tube->getM() << ";\n";
    MwSave << tube->getMw() << ";\n";
    
    energySave << scaledTotEnergy << ";\n";
}

void Trombone::closeFiles()
{
    massState.close();
    pState.close();
    vState.close();
    MSave.close();
    MwSave.close();
    energySave.close();

}
