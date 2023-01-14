#include "Matrix3.h"
#include "Point3D.h"
//the 3d matrix (4 by 4 to allow translations and projective transformations)
//this class can be used to perform standard linear, affine, and general projective transformations of space
class Matrix4
{
    float mat[4][4];
    //build a matrix from data
    Matrix4(float const data[4][4])
    {
        for(int i = 0;i<4;i++)
        {
            for(int j = 0;j<4;j++)
            {
                mat[i][j]=data[i][j];
            }
        }
    }
    //identity
    Matrix4()
    {
        for(unsigned int i = 0;i<4;i++)
        {
            mat[i][i]=1;
        }
    }
    //build a matrix from linear and translational parts
    Matrix4(const Matrix3 &mat3, const Vector3D &trans)
    {
        for(int i=0;i<3;i++)
        {
            for(int j =0;j<3;j++)
            {
                mat[i][j]=mat3.mat[i][j];
            }
            mat[0][3]=trans.x;
            mat[1][3]=trans.y;
            mat[2][3]=trans.z;
            mat[3][3]=1;
        }
    }
    //invert the matrix as an affine transformation
    Matrix4 invert()
    {
        float subdata[3][3];
        for(int i =0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                subdata[i][j]=mat[i][j];
            }
        }
        Matrix3 linearPart = Matrix3(subdata);
        linearPart = linearPart.inverse();
        Vector3D trans = Vector3D(mat[0][3],mat[1][3],mat[2][3]);
        trans = -(linearPart*trans);
        return Matrix4(linearPart,trans);
    }
    //the transpose
    Matrix4 transpose()
    {
        float data[4][4];
        for(int i = 0; i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                data[i][j]=mat[j][i];
            }
        }
        return Matrix4(data);
    }
    //the matrix product
    Matrix4 operator*(const Matrix4 &mat4)
    {
        float data[4][4];
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                float partialResult=0;
                for(int k=0;k<4;k++)
                {
                    partialResult += mat[i][k]*mat4.mat[k][j];
                }
                data[i][j]=partialResult;
            }
        }
        return Matrix4(data);
    }
    //linear transformation of basis
    static Matrix4 basisToBasis(const Vector3D &b11, const Vector3D &b12, const Vector3D &b13,
    const Vector3D &b21, const Vector3D &b22, const Vector3D &b23)
    {
        Matrix3 B1 = Matrix3();
        B1.standardBasisToBasis(b11,b12,b13);
        Matrix3 B2 = Matrix3();
        B2.standardBasisToBasis(b21,b22,b23);
        Matrix3 B12 = B1*(B2.inverse());
        return Matrix4(B12, Vector3D(0,0,0));
    }
    //Matrix4 acting on a point
    Point3D operator*(const Point3D &pt)
    {
        float data[4];
        float point[4];//copy the point to an array to be able to use indexes in a for loop
        point[0]=pt.x;point[1]=pt.y;point[2]=pt.z;point[4]=1;//last coord is always 1
        for(int i = 0;i<4;i++)
        {
            float partialResult = 0;
            for(int k =0;k<4;k++)
            {
                partialResult+=mat[i][k]*point[k];
            }
            data[i]=partialResult;
        }
        return Point3D(data[0]/data[3],data[1]/data[3],data[2]/data[3]);
    }
    //matrix to transform the standard projective frame (e1,e2,e3,e4,u(=e1+e2+e3+e4)) to any four points
    static Matrix4 standardProjectiveFrameToPoints(Point3D pt0, Point3D pt1, Point3D pt2, Point3D pt3, Point3D pt4)
    {
        //"enable" projective transformations
        Matrix4 M;
        for(int i =0;i<4;i++)
        {
            M.mat[3][i]=1;
        }
        //send pt0,pt1,pt2,pt3 to e1,e2,e3,e4 (transformation of 4-space)
        M.mat[0][0]=pt0.x;M.mat[0][1]=pt1.x;M.mat[0][2]=pt2.x;M.mat[0][3]=pt3.x;
        M.mat[1][0]=pt0.y;M.mat[1][1]=pt1.y;M.mat[1][2]=pt2.y;M.mat[1][3]=pt3.y;
        M.mat[2][0]=pt0.z;M.mat[2][1]=pt1.z;M.mat[2][2]=pt2.z;M.mat[2][3]=pt3.z;
        Matrix4 N;
        
        //find q such that Nq = pt4
        N=M.invert();
        float pt4h[4];
        pt4h[0]=pt4.x;pt4h[1]=pt4.y;pt4h[2]=pt4.z;pt4h[3]=1;
        float q[4];
        for (int i = 0; i < 4; i++)
        {
            double partialResult = 0;
            for (int j = 0; j < 4; j++)                {
               partialResult += N.mat[i][j] * pt4h[j];
            }
            q[i] = partialResult;
        }

        //matrix that sends e1,e2,e3,e4,u=e1+e2+e3+e4 to pt0,pt1,pt3,pt3,pt4
        Matrix4 L;
        L.mat[0][0]=q[0]*pt0.x;L.mat[0][1]=q[1]*pt1.x;L.mat[0][2]=q[2]*pt2.x;L.mat[0][3]=q[3]*pt3.x;
        L.mat[1][0]=q[0]*pt0.y;L.mat[1][1]=q[1]*pt1.y;L.mat[1][2]=q[2]*pt2.y;L.mat[1][3]=q[3]*pt3.y;
        L.mat[2][0]=q[0]*pt0.z;L.mat[2][1]=q[1]*pt1.z;L.mat[2][2]=q[2]*pt2.z;L.mat[2][3]=q[3]*pt3.z;
        L.mat[3][0]=q[0]*pt0.h;L.mat[3][1]=q[1]*pt1.h;L.mat[3][2]=q[2]*pt2.h;L.mat[3][3]=q[3]*pt3.h;
        return L;
    }
    //matrix to transform a projective frame to any four points
    static Matrix4 projectiveFrameToPoints
    (Point3D p0, Point3D p1, Point3D p2, Point3D p3, Point3D p4,
    Point3D q0, Point3D q1, Point3D q2, Point3D q3, Point3D q4)
    {
        Matrix4 A = (standardProjectiveFrameToPoints(p0,p1,p2,p3,p4)).invert();
        Matrix4 B = standardProjectiveFrameToPoints(q0,q1,q2,q3,q4);
        return B*A;
    }
};