//the 3d point (homogenous coordinate = 1)
class Point3D
{
    public:
    float x,y,z,h;
    Point3D(): x(0), y(0), z(0), h(1){}
    Point3D(float x, float y, float z): x(x), y(y), z(z),h(1) {}
    //addition and subtraction   
};