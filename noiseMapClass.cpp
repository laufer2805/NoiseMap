#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;


#define STRETCH_CONSTANT (-0.211324865405187)       /* (1 / sqrt(2 + 1) - 1 ) / 2 */
#define SQUISH_CONSTANT (0.366025403784439)         /* (sqrt(2 + 1) - 1) / 2 */
#define NORM_CONSTANT (47)
#define INT_MAX_VALUE (2147483647)                  /* max value of int type */
#define e (2.718281828)
#define DEFAULT_SEED (0)

typedef struct {
    double x, y, z, w;
} Point;

typedef struct
{
    int r, g, b;
} Color;

struct complex
{
    double x;
    double y;
};

static int grad3[12][3] = {
    {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
    {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
    {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}
};

static int p[256];
static int perm[512];
static int permMod12[512];

class OpenSimplexNoise {
public:
    OpenSimplexNoise(int seed = DEFAULT_SEED) {
        init(seed);
    }

    void init(int seed) {
        for (int i = 0; i < 256; i++) {
            p[i] = i;
        }

        for (int i = 0; i < 256; i++) {
            int j = rand() % (i + 1);
            std::swap(p[i], p[j]);
        }

        for (int i = 0; i < 512; i++) {
            perm[i] = p[i & 255];
            permMod12[i] = perm[i] % 12;
        }
    }

    double noise2D(double x, double y) {
        double stretchOffset = (x + y) * STRETCH_CONSTANT;
        double xs = x + stretchOffset;
        double ys = y + stretchOffset;

        int xsb = fastFloor(xs);
        int ysb = fastFloor(ys);

        double squishOffset = (xsb + ysb) * SQUISH_CONSTANT;
        double xb = xsb + squishOffset;
        double yb = ysb + squishOffset;

        double xins = xs - xsb;
        double yins = ys - ysb;

        double inSum = xins + yins;

        int hash =
            permMod12[xsb & 0xFF] +
            permMod12[ysb & 0xFF];

        double dx0 = x - xb;
        double dy0 = y - yb;

        double value = 0.0;

        if (inSum <= 1) {
            double zins = 1 - inSum;
            if (zins > xins || zins > yins) {
                if (xins > yins) {
                    xb += 1;
                    yb += 1;
                    double dx3 = dx0 - 1 - SQUISH_CONSTANT;
                    double dy3 = dy0 - 1 - SQUISH_CONSTANT;
                    value = dot2D(grad3[hash & 0x0B], dx3, dy3);
                } else {
                    xb += 1;
                    double dx2 = dx0 - 1 - SQUISH_CONSTANT;
                    double dy2 = dy0 - SQUISH_CONSTANT;
                    value = dot2D(grad3[hash & 0x08], dx2, dy2);
                }
            } else {
                yb += 1;
                double dx1 = dx0 - SQUISH_CONSTANT;
                double dy1 = dy0 - 1 - SQUISH_CONSTANT;
                value = dot2D(grad3[hash & 0x04], dx1, dy1);
            }
        } else {
            double zins = 2 - inSum;
            if (zins < xins || zins < yins) {
                if (xins > yins) {
                    xb += 2;
                    yb += 2;
                    double dx3 = dx0 - 2 - 2 * SQUISH_CONSTANT;
                    double dy3 = dy0 - 2 - 2 * SQUISH_CONSTANT;
                    value = dot2D(grad3[hash & 0x03], dx3, dy3);
                } else {
                    xb += 2;
                    double dx2 = dx0 - 2 - 2 * SQUISH_CONSTANT;
                    double dy2 = dy0 - 2 * SQUISH_CONSTANT;
                    value = dot2D(grad3[hash & 0x02], dx2, dy2);
                }
            } else {
                yb += 2;
                double dx1 = dx0 - 2 * SQUISH_CONSTANT;
                double dy1 = dy0 - 2 - 2 * SQUISH_CONSTANT;
                value = dot2D(grad3[hash & 0x01], dx1, dy1);
            }
        }

        return value * NORM_CONSTANT;
    }

    double noise3D(double x, double y, double z) {
        double stretchOffset = (x + y + z) * STRETCH_CONSTANT;
        double xs = x + stretchOffset;
        double ys = y + stretchOffset;
        double zs = z + stretchOffset;

        int xsb = fastFloor(xs);
        int ysb = fastFloor(ys);
        int zsb = fastFloor(zs);

        double squishOffset = (xsb + ysb + zsb) * SQUISH_CONSTANT;
        double xb = xsb + squishOffset;
        double yb = ysb + squishOffset;
        double zb = zsb + squishOffset;

        double xins = xs - xsb;
        double yins = ys - ysb;
        double zins = zs - zsb;

        double inSum = xins + yins + zins;

        int hash =
            permMod12[xsb & 0xFF] +
            permMod12[ysb & 0xFF] +
            permMod12[zsb & 0xFF];

        double dx0 = x - xb;
        double dy0 = y - yb;
        double dz0 = z - zb;

        double value = 0.0;

        if (inSum <= 1) {
            double squishOffsetXY = 0.0;
            double squishOffsetXZ = 0.0;
            double squishOffsetYZ = 0.0;
            if (xins >= yins && xins >= zins) {
                squishOffsetXY = 1;
                squishOffsetXZ = 1;
            } else if (yins >= xins && yins >= zins) {
                squishOffsetXY = 1;
                squishOffsetYZ = 1;
            } else {
                squishOffsetXZ = 1;
                squishOffsetYZ = 1;
            }
            value = singleSimplex(grad3[hash & 0x0F], dx0, dy0, dz0, xins, yins, zins, squishOffsetXY, squishOffsetXZ, squishOffsetYZ);
        } else if (inSum >= 2) {
            double squishOffsetXY = 1;
            double squishOffsetXZ = 1;
            double squishOffsetYZ = 1;
            double dx3 = dx0 - 1 - 3 * SQUISH_CONSTANT;
            double dy3 = dy0 - 1 - 3 * SQUISH_CONSTANT;
            double dz3 = dz0 - 1 - 3 * SQUISH_CONSTANT;
            value = singleSimplex(grad3[hash & 0x0F], dx3, dy3, dz3, xins, yins, zins, squishOffsetXY, squishOffsetXZ, squishOffsetYZ);
        } else {
            double squishOffsetXY = 1;
            double squishOffsetXZ = 0;
            double squishOffsetYZ = 0;
            double dx2 = dx0 - 2 * SQUISH_CONSTANT;
            double dy2 = dy0 - 2 * SQUISH_CONSTANT;
            double dz2 = dz0 - 0 * SQUISH_CONSTANT;
            value = singleSimplex(grad3[hash & 0x0F], dx2, dy2, dz2, xins, yins, zins, squishOffsetXY, squishOffsetXZ, squishOffsetYZ);
        }

        return value * NORM_CONSTANT;
    }

private:
    int fastFloor(double x) {
        int xi = (int)x;
        return x < xi ? xi - 1 : xi;
    }

    double dot2D(int* g, double x, double y) {
        return g[0] * x + g[1] * y;
    }

    double dot3D(int* g, double x, double y, double z) {
        return g[0] * x + g[1] * y + g[2] * z;
    }

    double dot4D(int* g, double x, double y, double z, double w) {
        return g[0] * x + g[1] * y + g[2] * z + g[3] * w;
    }

    double singleSimplex(int* g, double dx, double dy, double dz, double xins, double yins, double zins, double squishOffsetXY, double squishOffsetXZ, double squishOffsetYZ) {
        double value = 0.0;
        double attn = 2 - dx * dx - dy * dy - dz * dz;
        if (attn > 0) {
	        double vec[3] = {xins, yins, zins};
	        value = attn * attn * attn * attn * dot3D(g, vec[0], vec[1], vec[2]);
        }
        return value;
    }
};

// Generates and returns palette as array of Colors
Color* GeneratePalette (int length)
{
    Color* array = new Color[length];
    return array;
}

// fills array with colors from Lospec palette SLSO8
Color* FillPalette(Color* palette, int length)
{
    palette[0] = { 13, 43, 69 };
    palette[1] = { 32, 60, 86 };
    palette[2] = { 84, 78, 104 };
    palette[3] = { 141, 105, 122 };
    palette[4] = { 208, 129, 89 };
    palette[5] = { 255, 170, 94 };
    palette[6] = { 255, 212, 163 };
    palette[7] = { 255, 236, 214 };
    return palette;
}

class NoiseMap
{
    /*
        DrawCircle -> 4 random numbers -> 8 seed numbers
    */
private:
    // helper functions
    bool IsPrime(int number)
    {
        if (number == 1 || number == 0) return false;
        for (int i = 0; i < sqrt(number); i++)
        {
            if (number % i == 0) return false;
        }
        return true;
    }
    int Factorial(int n)
    {
        if (n == 0) return 1;
        return n * Factorial(n - 1);
    }
    double Power(double x, int n)
    {
        if (n == 0) return 1;
        return x * Power(x, n - 1);
    }
    double Exp(double x)
    {
        int approx = 15;
        double sum = 0;
        for (int i = 0; i < approx; i++)
        {
            sum += Power(x, i) / (double)Factorial(i);
        }
        return sum;
    }
    
    double Range (double number,double min, double max)
    {
		double decimals = number - (int)number;
		int numberInt = (int)number;
		int range = (int)max - (int)min;
		numberInt = numberInt % range;
		return (double)numberInt + min + decimals;
	}
    
    double* GenerateSequence (double seed, int size)
    {
		// for size < 10
		double* sequence = new double[size];
		int iterator = 0;
		for (double i = 0; i < size; i++)
		{
			double temp1 = sin(Exp(Range(i, ), temp2 = Exp(sin(seed + (double)i/(double)size));
			sequence[iterator] = rng.noise2D(temp1, temp2);
			iterator++;
		}
		return sequence;
	}

public:
    string seedWord;
    int sizeX, sizeY, startX, startY, paletteSize;
    OpenSimplexNoise rng;
    int** image;                                                                                    // double array of ints contains index of color in palette
    Color* palette;
    
    // Generates and returns palette as array of Colors
    void GeneratePalette ()
    {
        palette = new Color[paletteSize];
    }
    // fills array with colors from Lospec palette SLSO8
    void FillPalette ()
    {
        palette[0] = { 13, 43, 69 };
        palette[1] = { 32, 60, 86 };
        palette[2] = { 84, 78, 104 };
        palette[3] = { 141, 105, 122 };
        palette[4] = { 208, 129, 89 };
        palette[5] = { 255, 170, 94 };
        palette[6] = { 255, 212, 163 };
        palette[7] = { 255, 236, 214 };
    }


    void GenerateImageArray ()
    {
        image = new int* [sizeX];
        for (int i = 0; i < sizeX; i++)
        {
            image[i] = new int[sizeY];
        }
    }

    void DrawCircle (int centerX, int centerY, int radius, int colorIndex)
    {
        for (int i = centerX - radius; i < centerX + radius; i++)
        {
            for (int j = centerY - radius; i < centerY + radius; j++)
            {
                if (pow(i - centerX, 2) + pow(j - centerY, 2) < pow(radius, 2))
                {
                    image[i][j] = colorIndex;
                }
            }
        }
    }

    void FillArray ()
    {
        int seedWordLength = seedWord.length();
        int sizeMin;
        if (sizeX < sizeY) sizeMin = sizeX;
        sizeMin = sizeY;
        int* seedWordArray = new int[seedWordLength];
        for (int i = 0; i < seedWordLength; i++)
        {
            seedWordArray[i] = (int)seedWord[i];
        }
        
        for (int i = 0; i < sizeX; i++)
        {
            for (int j = 0; j < sizeX; j++)
            {
                int temp = (int)rng.noise2D(Exp((double)(i / sizeX)), Exp((double)(j / sizeY)));
                if (temp < 0) temp *= (-1);
                if (IsPrime(temp) == true)																		// ovaj section treba prepraviti - tu se poziva noise2D
                {
                    DrawCircle(i, j, temp % sizeMin, seedWordArray[temp % seedWordLength] % paletteSize);
                }
                else 
                {
                    image[i][j] = temp % paletteSize;
                }
            }
        }
    }

    void PrintArrayToFile ()
    {
        ofstream outFile;
        outFile.open(seedWord + ".ppm");
        outFile << "P3" << endl << sizeX << " " << sizeY << endl << 255 << endl;
        for (int i = 0; i < sizeX; i++)
        {
            for (int j = 0;j < sizeY; j++)
            {
                outFile << palette[image[i][j]].r << " " << palette[image[i][j]].g << " " << palette[image[i][j]].b << endl;
            }
        }
    }

    // Run this function!
    void GenerateNoiseMap()
    {
        // Input values from user
        cout << "Input seed word: ";
        std::getline(std::cin, seedWord);
        cout << "Input sizeX: ";
        cin >> sizeX;
        cout << "Input sizeY: ";
        cin >> sizeY;
        cout << "Input startX: ";
        cin >> startX;
        cout << "Input startY: ";
        cin >> startY;

        // Generating palette
        paletteSize = 8;
        GeneratePalette();
        FillPalette();

        cout << "Generating initial image!" << endl;
        GenerateImageArray();
        cout << "Filling array!" << endl;
        FillArray();
        cout << "Printing to file!" << endl;
        PrintArrayToFile();
        cout << "File done!" << endl;
    }
    void PrintSequence (double seed, int length)
    {
		double* seq = new double[length];
		seq = GenerateSequence(seed, length);
		for (int i = 0; i < length; i++)
		{
			cout << seq[i] << "\t";
		}
	}
};

int main (void)
{
    NoiseMap nm;
    nm.PrintSequence(5, 2);
	return 0;
}
