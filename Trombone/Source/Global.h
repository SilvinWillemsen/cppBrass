/*
  ==============================================================================

    Global.h
    Created: 5 Sep 2020 1:13:49pm
    Author:  Silvin Willemsen

  ==============================================================================
*/

#pragma once

namespace Global {
    
    static double pressureMultiplier = 10.0;
    static double oOPressureMultiplier = 1.0 / pressureMultiplier;

    static bool setTubeTo1 = false;
    static bool connectedToLip = false;
    static bool dontInterpolateAtStart = true;
    
    static std::vector<double> linspace (double start, double finish, int N)
    {
        std::vector<double> res (N, 0);
        for (int i = 0; i < N; ++i)
        {
            res[i] = start + i * (finish - start) / static_cast<double> (N - 1);
        }
        return res;
    }
    static double linspace (double start, double finish, int N, int idx)
    {
        if (idx >= N)
        {
            std::cout << "Idx is outside of range" << std::endl;
            return -1;
            
        }
        return start + idx * (finish - start) / static_cast<double> (N - 1);
    }
    
    static inline double subplus (double val) { return (val + abs(val)) * 0.5; };
    
    static inline int sgn (double val) { return (0 < val) - (val < 0); };

    static double outputClamp (double val)
    {
        if (val < -1.0)
        {
            val = -1.0;
            std::cout << "outputClamped!!" << std::endl;
            return val;
        }
        else if (val > 1.0)
        {
            val = 1.0;
            std::cout << "outputClamped!!" << std::endl;
            return val;
        }
        return val;
    }
}
