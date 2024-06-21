#include <iostream>
using namespace std;

class Vector
{
    public:
        int dimension;
        float* v;
        void GenerateVector ()
        {
            v = new float[dimension];
        }
        void InputVector()
        {
            for (int i = 0; i < dimension; i++)
            {
                cout << "v[" << i << "] = ";
                cin >> v[i];
            }
        }
        static void PrintVector (Vector v)
        {
            cout << "(";
            for (int i = 0; i < v.dimension; i++)
            {
                cout << v.v[i];
                if (i != v.dimension - 1)
                {
                    cout << ", ";
                }
                else
                {
                    cout << ")";
                }
            }
        }
};
float Power (float x, int n)
{
    if (n == 0) return 1;
    return x * Power(x, n - 1);
}
class Matrix
{
    public:
        int size;
        float** m;
        
        void GenerateMatrix ()
        {
            m = new float*[size];
            for (int i = 0; i < size; i++)
            {
                m[i] = new float[size];
            }
        }
        
        void InputMatrix ()
        {
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < size; j++)
                {
                    cout << "m[" << i << "][" << j << "] = ";
                    cin >> m[i][j];
                }
            }
        }
        
        static void PrintMatrix (Matrix A)
        {
            for (int i = 0; i < A.size; i++)
            {
                for (int j = 0; j < A.size; j++)
                {
                    cout << A.m[i][j] << "\t";
                }
                cout << endl;
            }
        }
        
        static float det(Matrix A)
        {
            if (A.size == 2)
            {
                return A.m[0][0] * A.m[1][1] - A.m[0][1] * A.m[1][0];
            }
            int l = 0;
            float sum = 0;
            for (int k = 0; k < A.size; k++)
            {
                float sign = Power(-1, l + k);
                Matrix temp;
                temp.size = A.size - 1;
                temp.GenerateMatrix();
                int tempI = 0, tempJ = 0;
                for (int i = 0; i < A.size; i++)
                {
                    for (int j = 0; j < A.size; j++)
                    {
                        if (i != l && j != k)
                        {
                            temp.m[tempI][tempJ] = A.m[i][j];
                            tempJ++;
                            if (tempJ == temp.size)
                            {
                                tempJ = 0;
                                tempI++;
                            }
                        }
                    }
                }
                sum += sign * A.m[l][k] * det(temp);
            }
            return sum;
        }
};

int main (void)
{
    Matrix a;
    a.size = 3;
    a.GenerateMatrix ();
    a.InputMatrix ();
    cout << Matrix::det(a) << endl;
    return 0;
}