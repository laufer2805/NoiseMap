#include <iostream>
#include <iostream>
#include <iostream>
#include </media/data/programming/GitHub/NoiseMap/Scripts/OpenSimplexNoise.cpp>
#include </media/data/programming/GitHub/NoiseMap/Scripts/MathLibrary.cpp>

class PseudoRandomNumberGenerator
{
    private:
        OpenSimplexNoise noise;
    public:
        double* GenerateSequence2D (double seedX, double seedY, int length)
        {
            double* seq = new double[length];
            for (int i = 0; i < length; i++)
            {
                double k = (double)i / length;
                seq[i] = noise.noise2D (seedX + k, seedY + k);
            }
            return seq;
        }
        double Range2D (double seedX, double seedY, double min, double max)
        {
            double random = noise.noise2D (seedX, seedY);
            random = Clamp (random, min, max);
            return random;
        }
        double Range3D (double seedX, double seedY, double seedZ, double min, double max)
        {
            double random = noise.noise3D (seedX, seedY, seedZ);
            random = Clamp (random, min, max);
            return random;
        }
        // Generates random int int range [min, max>
        int IntRange2D (int seedX, int seedY, int min, int max)
        {
            return IntClamp((int)noise.noise2D ((double)seedX, (double)seedY), min, max);
        }
};
