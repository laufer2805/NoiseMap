#include <iostream>
#include <iostream>
#include <iostream>
#include <cmath>

#define Pi	3.14159265358979323846

int Factorial (int n)
{
	if (n == 0) return 1;
	return n * Factorial (n - 1);
}

double Power (double x, int n)
{
	if (n == 0) return 1;
	return x * Power (x, n - 1);
}

double Exp (double x)
{
	int approx = 15;
	double sum = 0;
	for (int i = 0; i < approx; i++)
	{
		sum += Power (x, i) / (double)Factorial (i);
	}
	return sum;
}

double Clamp (double value, double min, double max)
{
	double range = (max - min) / 2;
	double c = min + range;
	return range * sin (value) + c;
}

int IntClamp (int value, int min, int max)
{
	int range = max - min;
	if (value < 0) value *= -1;
	return value % range + min;
}
