/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent()
{
    // Make sure you set the size of the component after
    // you add any child components.
    
    // specify the number of input and output channels that we want to open
    
    setAudioChannels (0, 2);
    startTimerHz (15);
    
}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    stopTimer();
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
//    auto test = deviceManager.getAudioDeviceSetup();
//    std::cout << test.sampleRate << std::endl;
//    test.sampleRate = 44100;
//    std::cout << deviceManager.setAudioDeviceSetup (test, false) << std::endl;
//    auto test2 = deviceManager.getAudioDeviceSetup();
//    std::cout << test2.sampleRate << std::endl;

    fs = sampleRate;
    NamedValueSet parameters;
    
    //// Tube ////
    parameters.set ("T", 26.85);
    parameters.set ("LnonExtended", 2.658);
    parameters.set ("Lextended", 3.718);
    parameters.set ("L", 3.718);


    // Geometry
//    parameters.set ("mp", 0.015 * 0.015 * double_Pi);                   // mouthpiece cross-sectional area
//    parameters.set ("mpL", 0.01);                   // mouthpiece length (length ratio)
//    parameters.set ("m2tL", 0.01);                  // mouthpiece to tube (length ratio)
//    parameters.set ("tubeS", 0.01 * 0.01 * double_Pi);                 // tube cross-sectional area
    
    // Geometric information including formula from bell taken from T. Smyth "Trombone synthesis by model and measurement"
    geometry = {
        {0.708, 0.177, 0.711, 0.306, 0.254, 0.502},         // lengths (changed fourth entry to account for bell length "error" in paper)
        {0.0069, 0.0072, 0.0069, 0.0071, 0.0075, 0.0107}    // radii
    };
    
    parameters.set ("flare", 0.7);                 // flare (exponent coeff)
    parameters.set ("x0", 0.0174);                    // position of bell mouth (exponent coeff)
    parameters.set ("b", 0.0063);                   // fitting parameter
    parameters.set ("bellL", 0.21);                  // bell (length ratio)
    
    //// Lip ////
    double f0 = 300.0;
    double H0 = 2.9e-4;
    parameters.set("f0", f0);                       // fundamental freq lips
    parameters.set("Mr", 5.37e-5);                  // mass lips
    parameters.set("omega0", 2.0 * double_Pi * f0); // angular freq
    
    parameters.set("sigmaR", 5);                    // damping
    parameters.set("H0", H0);                       // equilibrium
    parameters.set("barrier", -H0);                 // equilibrium

    parameters.set("w", 1e-2);                      // lip width
    parameters.set("Sr", 1.46e-5);                  // lip area
    
    parameters.set ("Kcol", 10000);
    parameters.set ("alphaCol", 3);
    
    //// Input ////
    parameters.set ("Pm", 300 * Global::pressureMultiplier);
    pressureVal = (*parameters.getVarPointer ("Pm"));
    lipFreqVal = (*parameters.getVarPointer ("f0"));
    LVal = (*parameters.getVarPointer ("LnonExtended")); // start by contracting
    
    trombone = std::make_unique<Trombone> (parameters, 1.0 / fs, geometry);
    addAndMakeVisible (trombone.get());
    
    trombone->setExtVals (pressureVal, lipFreqVal, LVal);
    
    setSize (800, 600);

}

void MainComponent::getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    
    float* const channelData1 = bufferToFill.buffer->getWritePointer (0, bufferToFill.startSample);
    float* const channelData2 = bufferToFill.buffer->getWritePointer (1, bufferToFill.startSample);
    
    float output = 0.0;
    float output2 = 0.0;
    
    for (int i = 0; i < bufferToFill.numSamples; ++i)
    {
        trombone->calculate();
        output = trombone->getOutput() * 0.001 * Global::oOPressureMultiplier;
        if (!done && Global::saveToFiles)
            trombone->saveToFiles();
        
        trombone->updateStates();
        channelData1[i] = Global::outputClamp (output);
        channelData2[i] = Global::outputClamp (output);
        ++t;
    }
    if (Global::saveToFiles && t > 10000 && !done)
    {
        done = true;
        trombone->closeFiles();
        std::cout << "done" << std::endl;
    }
    trombone->refreshLipModelInputParams();
}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (Colours::white);
    g.drawText("Pressure: " + String (pressureVal) + "(Pa)   LipFrequency: " +
               String (lipFreqVal) + "(Hz)   Length:" +
               String (LVal) + "(m)",
               getWidth() - 400, getHeight() - 50, 400, 50, Justification::centredRight);
    // You can add your drawing code here!
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
    controlHeight = 0.2 * getHeight();
    controlY = getHeight() - controlHeight;
    trombone->setBounds (getLocalBounds().withHeight (controlY));
}

void MainComponent::timerCallback()
{
    repaint();
}


void MainComponent::mouseDown (const MouseEvent& e)
{
//    pressureVal = e.y * Global::pressureMultiplier;
//    lipFreqVal = e.x;
//    lipFreqVal = e.y;
    pressureVal = 300 * Global::pressureMultiplier;
    double xRatio = e.x / static_cast<double> (getWidth());
    double fineTuneRange = 0.1;
    double fineTune = fineTuneRange * 2 * (e.y - controlY - controlHeight * 0.5) / controlHeight;
    lipFreqVal = ((1-xRatio) + xRatio * 2.658/3.718) * Global::nonExtendedLipFreq * (1 + fineTune);

    LVal = 2.658 + 1.06 * e.x / static_cast<double> (getWidth());
    
    trombone->setExtVals (pressureVal, lipFreqVal, LVal);
}

void MainComponent::mouseDrag (const MouseEvent& e)
{
//    pressureVal = e.y * Global::pressureMultiplier;
    pressureVal = 300 * Global::pressureMultiplier;
//    lipFreqVal = e.x;
    double xRatio = e.x / static_cast<double> (getWidth());
    
    double fineTuneRange = 0.1;
    double fineTune = fineTuneRange * 2 * (e.y - controlY - controlHeight * 0.5) / controlHeight;
    lipFreqVal = ((1-xRatio) + xRatio * 2.658/3.718) * Global::nonExtendedLipFreq * (1 + fineTune);
    
    LVal = 2.658 + 1.06 * e.x / static_cast<double> (getWidth());

    trombone->setExtVals (pressureVal, lipFreqVal, LVal);
}

void MainComponent::mouseUp (const MouseEvent& e)
{
    pressureVal = 0;
//    lipFreqVal = e.x;
//    lipFreqVal = e.y;
//    double xRatio = e.x / static_cast<double> (getWidth());
//    lipFreqVal = ((1-xRatio) + xRatio * 2.658/3.718) * 520;
//    LVal = 2.658 + 1.06 * e.x / static_cast<double> (getWidth());

    trombone->setExtVals (pressureVal, lipFreqVal, LVal);
}
