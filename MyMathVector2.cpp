// Online C++ compiler to run C++ program online
#include <iostream>
#include <cmath>

class Vector2
{
    public:
        float x, y;
        void define(float x, float y)
        {
            Vector2::x = x;
            Vector2::y = y;
        }
        
        void Add (Vector2 a, Vector2 b)
        {
            Vector2::x = a.x + b.x;
            Vector2::y = a.y + b.y;
        }
        void ScalarMultiply (float scalar)
        {
            Vector2::x = scalar * x;
            Vector2::y = scalar * y;
        }
        float ScalarProduct (Vector2 v1, Vector2 v2)
        {
            return v1.x * v2.x + v1.y * v2.y;
        }
        float Magnitude(void)
        {
            return sqrt(Vector2::x * Vector2::x + Vector2::y * Vector2::y);
        }
        void Normalized(void)
        {
            float x, y, magnitude;
            magnitude = Vector2::Magnitude();
            Vector2::ScalarMultiply (1/magnitude);
        }
        /*float Angle (Vector2 from, Vector2 to)
        {
            return acos(ScalarProduct(from, to)/(from.Magnitude()*to.Magnitude()));
        }*/
        
    static float ScalarProd (Vector2 v1, Vector2 v2)
    {
            return v1.x * v2.x + v1.y * v2.y;
    }
    static float Magn(Vector2 v)
    {
        return sqrt(v.x * v.x + v.y * v.y);
    }
    static float Angle (Vector2 from, Vector2 to)
        {
            return acos(ScalarProd(from, to) / (Magn (from) * Magn (to)));
        }
};

int main() {
    Vector2 a, b;
    a.define (1, 2);
    b.define (3, 4);
    Vector2 c;
    c.Add(a, b);
    c.ScalarMultiply (5);
    c.Normalized();
    std::cout << Vector2::Angle(a, b) << std::endl;

    return 0;
}