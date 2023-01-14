#include "Vector3D.h"
//submatrix class with some utility functions
class Matrix3
{
    public:
    float mat[3][3];
    Matrix3(float const data[3][3])
    {
        for(int i =0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                mat[i][j]=data[i][j];
            }
        }
    }
    Matrix3()
    {
        for(unsigned int i = 0;i<3;i++)
        {
            mat[i][i]=1;
        }
    }
    Vector3D operator*(const Vector3D &vec)
    {
        return Vector3D(mat[0][0]*vec.x+mat[0][1]*vec.y+mat[0][2]*vec.z,
                        mat[1][0]*vec.x+mat[1][1]*vec.y+mat[1][2]*vec.z,
                        mat[2][0]*vec.x+mat[2][1]*vec.y+mat[2][2]*vec.z);
    }
    Matrix3 inverse()
    {
        float invdet = mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1] * mat[1][2]) -
             mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
             mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
        if(invdet==0)
        {
            std::cout
            <<"Matrix Not Invertible";
        }
        invdet = 1/invdet;
        float data[3][3];
        data[0][0] = (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) * invdet;
        data[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * invdet;
        data[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * invdet;
        data[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * invdet;
        data[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * invdet;
        data[1][2] = (mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2]) * invdet;
        data[2][0] = (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]) * invdet;
        data[2][1] = (mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1]) * invdet;
        data[2][2] = (mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]) * invdet;

        return Matrix3(data);
    }
        void insertColumnVector(const Vector3D &vec,int column)
    {
        mat[0][column]=vec.x;
        mat[1][column]=vec.y;
        mat[2][column]=vec.z;
    }
    void standardBasisToBasis(const Vector3D &u,const Vector3D &v,const Vector3D &w)
    {
        this->insertColumnVector(u, 0);
        this->insertColumnVector(v, 1);
        this->insertColumnVector(w, 2);
    }
    Matrix3 operator*(const Matrix3 &mat3)
    {
        float data[3][3];
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                float partialResult=0;
                for(int k=0;k<3;k++)
                {
                    partialResult += mat[i][k]*mat3.mat[k][j];
                }
                data[i][j]=partialResult;
            }
        }
        return Matrix3(data);
    }
    float determinant()
    {
        return mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1] * mat[1][2]) -
             mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
             mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    }
    float trace()
    {
        return mat[0][0]+mat[1][1]+mat[2][2];
    }
};