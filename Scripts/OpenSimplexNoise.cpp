#include <iostream>
#include <cmath>

#define STRETCH_CONSTANT (-0.211324865405187)    /* (1 / sqrt(2 + 1) - 1 ) / 2 */
#define SQUISH_CONSTANT (0.366025403784439)      /* (sqrt(2 + 1) - 1) / 2 */
#define NORM_CONSTANT (47)

#define DEFAULT_SEED (0)

typedef struct {
    double x, y, z, w;
} Point;

static int grad3[12][3] = {
    {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0},
    {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1},
    {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}
};

static int grad4[32][4] = {
    {0, 1, 1, 1}, {0, 1, 1, -1}, {0, 1, -1, 1}, {0, 1, -1, -1},
    {0, -1, 1, 1}, {0, -1, 1, -1}, {0, -1, -1, 1}, {0, -1, -1, -1},
    {1, 0, 1, 1}, {1, 0, 1, -1}, {1, 0, -1, 1}, {1, 0, -1, -1},
    {-1, 0, 1, 1}, {-1, 0, 1, -1}, {-1, 0, -1, 1}, {-1, 0, -1, -1},
    {1, 1, 0, 1}, {1, 1, 0, -1}, {1, -1, 0, 1}, {1, -1, 0, -1},
    {-1, 1, 0, 1}, {-1, 1, 0, -1}, {-1, -1, 0, 1}, {-1, -1, 0, -1},
    {1, 1, 1, 0}, {1, 1, -1, 0}, {1, -1, 1, 0}, {1, -1, -1, 0},
    {-1, 1, 1, 0}, {-1, 1, -1, 0}, {-1, -1, 1, 0}, {-1, -1, -1, 0}
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
	        int xsb = fastFloor(dx);
	        int ysb = fastFloor(dy);
	        int zsb = fastFloor(dz);
	        int index = perm[(xsb & 0xFF) + perm[(ysb & 0xFF) + perm[(zsb & 0xFF)]]]; // Corrected line
	        double vec[3] = {xins, yins, zins};
	        value = attn * attn * attn * attn * dot3D(g, vec[0], vec[1], vec[2]);
        }
        return value;
    }
};
