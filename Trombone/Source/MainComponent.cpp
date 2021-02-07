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
    Timer::stopTimer();
    HighResolutionTimer::stopTimer();
    // start the hi-res timer
    shutdownAudio();
    
    for (auto sensel : sensels)
    {
        if (sensel->senselDetected)
        {
            sensel->shutDown();
        }
    }
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    sensels.add (new Sensel (0)); // chooses the device in the sensel device list
    
    fs = sampleRate;
    NamedValueSet parameters;

    //// Tube ////
    parameters.set ("T", 26.85);
    parameters.set ("L", 2.658);

    // Geometry
    parameters.set ("inner1Rad", 0.0069); // Rad first inner tube
    parameters.set ("slideRad", 0.0072);  // Rad slide
    parameters.set ("inner2Rad", 0.0069); // Rad second inner tube
    parameters.set ("gooseNeckRad", 0.0071); // Rad gooseneck
    parameters.set ("tuningRad1", 0.0075); // First rad tuning slide
    parameters.set ("tuningRad2", 0.0107); // Second rad tuning slide
    
    parameters.set ("inner1L", 0.708); // Length first inner tube
    parameters.set ("slideL", 0.177);  // Slide length (non-extended)
    parameters.set ("inner2L", 0.711); // Length first inner tube
    parameters.set ("gooseNeckL", 0.241); // Rad gooseneck
    parameters.set ("tuningL", 0.254); // First rad tuning slide
    parameters.set ("bellL", 0.502); // Second rad tuning slide
    
    // formula from bell taken from T. Smyth "Trombone synthesis by model and measurement"
    
    parameters.set ("flare", 0.7);                 // flare (exponent coeff)
    parameters.set ("x0", 0.0174);                    // position of bell mouth (exponent coeff)
    parameters.set ("b", 0.0063);                   // fitting parameter
    
    //// Lip ////
    double f0 = 300.0;
    double H0 = 2.9e-4;
    parameters.set("f0", f0);                       // fundamental freq lips
    parameters.set("Mr", 5.37e-5);                  // mass lips
    parameters.set("omega0", 2.0 * double_Pi * f0); // angular freq
    
    parameters.set("sigmaR", 5);                    // damping
    parameters.set("H0", H0);                   // equilibrium
    parameters.set("barrier", -H0);                   // equilibrium

    parameters.set("w", 1e-2);                      // lip width
    parameters.set("Sr", 1.46e-5);                  // lip area
    
    parameters.set ("Kcol", 100);
    parameters.set ("alphaCol", 5);
    
    //// Input ////
    parameters.set ("Pm", 300 * Global::pressureMultiplier);
    
    trombone = std::make_unique<Trombone> (parameters, 1.0 / fs);
    addAndMakeVisible (trombone.get());
    
    if (sensels.size() != 0)
        if (sensels[0]->senselDetected)
            HighResolutionTimer::startTimer (1000.0 / 150.0); // 150 Hz
    
    setSize (1000, 800);

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
        trombone->updateStates();
        channelData1[i] = Global::outputClamp (output);
        channelData2[i] = Global::outputClamp (output);

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
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    // You can add your drawing code here!
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
    trombone->setBounds (getLocalBounds());
}

void MainComponent::timerCallback()
{
    repaint();
}

void MainComponent::hiResTimerCallback()
{
    double maxPm = 3000.0 * Global::pressureMultiplier;
    double maxf0 = 1000.0;
    for (auto sensel : sensels)
    {
        double finger0X = 0;
        double finger0Y = 0;
        if (sensel->senselDetected)
        {
            sensel->check();
            unsigned int fingerCount = sensel->contactAmount;

            for (int f = 0; f < fingerCount; f++)
            {
                bool state = sensel->fingers[f].state;
                float x = sensel->fingers[f].x;
                float y = sensel->fingers[f].y;
                float Pm = Global::clamp (sensel->fingers[f].force * 10000.0 * Global::pressureMultiplier, 0, maxPm);
                
                int fingerID = sensel->fingers[f].fingerID;
                
                if (fingerID == 0 && state)
                {
                    finger0X = x;
                    finger0Y = y;
                    double f0 = Global::clamp (x * maxf0, 0, maxf0);
                    trombone->setParams (Pm, f0);
                }
            }
            
            if (fingerCount == 0)
            {
                trombone->setParams (0, trombone->getLipFreq());
            }
        }
    }
}
