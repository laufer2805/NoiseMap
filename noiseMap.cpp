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

double f (double x)
{
	return (exp(x) + exp(-x)) / 2;
}

double g (double x)
{
	return (exp(x) - exp(-x)) / 2;
}

complex MandelbrotHelperfunction (complex z, complex c)
{
    complex z0;
    z0.x = z.x*z.x - z.y*z.y + c.x;
    z0.y = 2*z.x*z.y + 2*c.y;
    return z0;
}

double abs(complex z)
{
    return sqrt(z.x*z.x + z.y*z.y);
}

bool isMandelbrot (complex c, int maxIterations = 100, double threshold = 4)
{
    int i;
    complex z;
    for (i = 0; i < 100; i++)
    {
        z = MandelbrotHelperfunction(z, c);
        if (abs(z) > threshold)
        {
            return false;
        }
    }
    return true;
}


void GenerateNoiseMap (string seedWord, int sizeX, int sizeY, double startX, double startY)
{
	
	OpenSimplexNoise RNG;
	int i, j, seedWordLength = seedWord.length();
	double rValue, gValue, bValue, seedDouble = 0;
	
	for (i = 0; i < seedWordLength; i++)
	{
		seedDouble += (double)seedWord[i];
	}
	seedDouble /= (double)seedWordLength;
	
	ofstream outFile;
	outFile.open (seedWord + ".ppm");
	outFile << "P3" << endl << sizeX << " " << sizeY << endl << 255 << endl;

	for (i = 0; i < sizeX; i++)
	{
		for (j = 0; j < sizeY; j++)
		{
			// grayscale
			/*
			double temp = RNG.noise2D (startX + (double)i, startY + (double)j);
			if (temp < 0) temp *= (-1);
			matrix[i][j] = (int)temp % 256;
			*/
			//rgb
			rValue = RNG.noise2D (20 * (double)i / sizeX, 20 * (double)j / sizeY);
			
			// conversion to pixel values
			if (rValue < 0) rValue *= (-1);
			if (gValue < 0) gValue *= (-1);
			if (bValue < 0) bValue *= (-1);
			outFile << (int)rValue % 256 << " " << (int)gValue % 256 << " " << (int)bValue % 256 << endl;
		}
	}

	cout << "File Written!" << endl;
	outFile.close();
}


void GenerateMandelbrotMap (int sizeX, int sizeY)
{
	ofstream outFile;
	cout << "File Opened!" << endl;
	outFile.open ("MandelBrotMap.pbm");
	complex c;
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			cout << "i = " << i << endl << "j = " << j << endl;
			c.x = (double)(i - sizeX/2);
			c.y = (double)(j - sizeY/2);
			if (isMandelbrot(c) == true)
			{
				outFile << "1" << " ";
			}
			else
			{
				outFile << "0" << " ";
			}
		}
		outFile << endl;
	}
	cout << "File Generated" << endl;
	outFile.close();
}

int main (void)
{
	/*
	 * define extents of matrix
	 * generate matrix and fill it with pseudorandom numbers
	 * limit generated numbers to range [0, 255]
	 * save matrix to file
	 * */
	
	// generate and fill matrix
	cout << "Program Started!" << endl;
	// seedWordLength 6
	GenerateNoiseMap ("123456", 1920, 1080, 0, 0);
	
	
	return 0;
}
