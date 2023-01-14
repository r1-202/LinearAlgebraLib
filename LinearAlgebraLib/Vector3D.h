#include<math.h>
#include<iostream>

//the 3d vector (homogenous coordinate = 0)
class Vector3D
{
    public:
    float x,y,z;
    Vector3D(): x(0), y(0), z(0){}
    Vector3D(float x, float y, float z): x(x), y(y), z(z) {}
    //addition and subtraction
    Vector3D operator+(const Vector3D &vec)
    {
        return Vector3D(x+vec.x,y+vec.y,z+vec.z);
    }
    friend Vector3D operator-(const Vector3D &vec);
    Vector3D operator-(const Vector3D &vec)
    {
        return (*this)+(-vec);
    }
    //the dot product
    float operator*(const Vector3D &vec)
    {
        return x*vec.x+y*vec.y+z*vec.z;
    }
    //the scalar product
    Vector3D operator^(const Vector3D &vec)
    {
        return Vector3D(y*vec.z-z*vec.y,
                        z*vec.x-x*vec.z,
                        x*vec.y-y*vec.x);
    }
    //multiplication and division by a scalar
    friend Vector3D operator*(float s, const Vector3D &vec);
    Vector3D operator/(float s)
    {
        return(Vector3D(x/s,y/s,z/s));
    }
    //normalizing
    float normSquared()
    {
        return x*x+y*y+z*z;
    }
    float norm()
    {
        return sqrt((*this).normSquared());
    }
    void normalize()
    {
        (*this)=Vector3D((*this)/this->norm());
    }
};
//scalar and vector multiplication
Vector3D operator*(float s,const Vector3D &vec)
{
    return Vector3D(s*vec.x,s*vec.y,s*vec.z);
}
Vector3D operator-(const Vector3D &vec)
{
    return Vector3D(-vec.x,-vec.y,-vec.z);
}