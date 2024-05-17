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

void GenerateNoiseMap (string seedWord, int sizeX, int sizeY, int startX, int startY)
{
	OpenSimplexNoise RNG;                                                                                                       // Instantiates OpenSimplexNoise class
	int i, j;
    Color* palette = GeneratePalette(8);                                                                                        // Generates palette as dynamic array
    palette = FillPalette(palette, 8);                                                                                          // Fills palette and saves to same array

    int seedWordLength = seedWord.length();                                                                                     // Saves seed word to int array
    int* seedWordArray = new int[seedWordLength];
    for (i = 0; i < seedWordLength; i++)
    {
        seedWordArray[i] = (int)seedWord[i];                                                                                    // Saves each letter as int to array
        if (seedWordArray[i] < 0) seedWordArray[i] *= (-1);                                                                     // if nonstandard letter used int will be negative
    }

	ofstream outFile;
	outFile.open (seedWord + ".ppm");
	
	outFile << "P3" << endl << sizeX << " " << sizeY << endl << 255 << endl;
	
	

	for (i = 0; i < sizeX; i++)
	{
		for (j = 0; j < sizeY; j++)
		{ 
            int temp = (int)RNG.noise2D((double)i, (double)j);
            if (temp < 0) temp *= (-1);
            int randomLetter = seedWordArray[temp % seedWordLength];															// random letter from seedWord as int
            int randomValue = (int)RNG.noise3D((double)(startX + i), (double)(startY + j), (double)randomLetter);				// get random int
            if (randomValue < 0) randomValue *= (-1);																			// check if negative number
            randomValue = randomValue % 8;                                                                                      // limit to [0, 7] 
            outFile << palette[randomValue].r << " " << palette[randomValue].g << " " << palette[randomValue].b << endl;        // save random member of palette to pixel at coodinates i, j
		}
	}

	cout << "File Written!" << endl;                                                                                            // Close file and indicate that algo is done
	outFile.close();
}

int main (void)
{
	// Instagram native resolution: 1080 x 1080
    
    string seedWord;
    int sizeX, sizeY, startX, startY;
	cout << "Enter seed word: ";
	getline(cin, seedWord);
    cout << "Enter resolution: \nWidth: ";
    cin >> sizeX;
    cout << "Height: ";
    cin >> sizeY;
    cout << "Enter starting point: \nx = ";
    cin >> startX;
    cout << "y = ";
    cin >> startY;
	
	cout << "Program started!" << endl;
	
	GenerateNoiseMap(seedWord, sizeX, sizeY, startX, startY);
	
	
	return 0;
}
