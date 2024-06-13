/*==========================================================================
Copyright by:
    Macrograph, Research

    For further information, contact:
        Dr. Seungtaik Oh

        E-mail : stoh@macrograph.co.kr
        Phone  : +82 10 2923 6290
        Address : Seoul, Korea

This copyright notice must be included with all copies of the source codes.
============================================================================*/
#pragma once

#include <iostream>
#include <cmath>

namespace mg
{
template <typename T>
class Vector3
{
public:
    //! Default constructor
    Vector3() : x(0), y(0), z(0) {}
    //! Copy constructor
    Vector3(const Vector3<T>& v) : x(v.x), y(v.y), z(v.z) {}
    Vector3(T tx, T ty, T tz) : x(tx), y(ty), z(tz) {}
    Vector3(T* v) : x(v[0]), y(v[1]), z(v[2]) {}

    // Get data
    void get(T* array) { array[0] = x, array[1] = y, array[2] = z; }

    //! Data pointer
    T* data() { return (T*)&x; }

    //! Dot product
    T dotProduct(const Vector3<T>& V) const
    {
        return x*V.x + y*V.y + z*V.z;
    }

    //! Squared norm
    T norm2() const
    {
        return x*x + y*y + z*z;
    }

    //! Norm
    T norm() const
    {
        return sqrt(x*x + y*y + z*z);
    }

    //! Normalize
    void normalize()
    {
        const T length = this->norm();
        x /= length;
        y /= length;
        z /= length;
    }

    //! Squared distance to a given vector
    T dist2To(const Vector3<T>& V) const
    {
        return (x - V.x)*(x - V.x) + (y - V.y)*(y - V.y) + (z - V.z)*(z - V.z);
    }

    //! Distance to a given vector
    T distTo(const Vector3<T>& V) const
    {
        return sqrt((x - V.x)*(x - V.x) + (y - V.y)*(y - V.y) + (z - V.z)*(z - V.z));
    }

    //! Cross product
    const Vector3<T> crossProduct(const Vector3<T>& V) const
    {
        return Vector3<T>(y*V.z - z*V.y,
            z*V.x - x*V.z,
            x*V.y - y*V.x);
    }

    // operator
    const   Vector3<T>&    operator += (const Vector3<T>& V)
    {
        x += V.x;
        y += V.y;
        z += V.z;

        return *this;
    }
    const   Vector3<T>&    operator -= (const Vector3<T>& V)
    {
        x -= V.x;
        y -= V.y;
        z -= V.z;

        return *this;
    }
    const   Vector3<T>&    operator *= (T k)
    {
        x *= k;
        y *= k;
        z *= k;

        return *this;
    }
    const   Vector3<T>&    operator /= (T k)
    {
        x /= k;
        y /= k;
        z /= k;

        return *this;
    }
    friend const Vector3<T> operator+ (const Vector3<T>& L, const Vector3<T>& R) { return Vector3<T>(L) += R; }
    friend const Vector3<T> operator- (const Vector3<T>& L, const Vector3<T>& R) { return Vector3<T>(L) -= R; }
    friend const Vector3<T> operator* (const Vector3<T>& L, T k)         { return Vector3<T>(L) *= k; }
    friend const Vector3<T> operator/ (const Vector3<T>& L, T k)         { return Vector3<T>(L) /= k; }
    friend const Vector3<T> operator* (T k, const Vector3<T>& R)         { return Vector3<T>(R) *= k; }

    friend const Vector3<T> operator- (const Vector3<T>& U)              { return Vector3<T>(-U.x, -U.y, -U.z); }

    friend std::ostream& operator << (std::ostream& os, const Vector3<T>& V)
    {
        os << "(" << V.x << ", " << V.y << ", " << V.z << ")";
        return os;
    }

public:
    T x;
    T y;
    T z;
};

typedef Vector3<float>      Vector3f;
typedef Vector3<double>     Vector3d;
typedef Vector3<int>        Vector3i;
typedef Vector3<unsigned int> Vector3ui;
typedef Vector3<unsigned char> Vector3uc;

template <typename T>
class Vector2
{
public:
    //! Default constructor
    Vector2() : x(0), y(0) {}
    //! Copy constructor
    Vector2(const Vector2<T>& v) : x(v.x), y(v.y) {}
    Vector2(T tx, T ty) : x(tx), y(ty) {}
    Vector2(T* v) : x(v[0]), y(v[1]) {}

    //! Data pointer
    T* data() { return (T*)&x; }

    //! Dot product
    T dotProduct(const Vector2<T>& V) const
    {
        return x*V.x + y*V.y;
    }

    //! Squared norm
    T norm2() const
    {
        return x*x + y*y;
    }

    //! Norm
    T norm() const
    {
        return sqrt(x*x + y*y);
    }

    //! Normalize
    void normalize()
    {
        const T length = this->norm();
        x /= length;
        y /= length;
    }

    //! Squared distance to a given vector
    T dist2To(const Vector2<T>& V) const
    {
        return (x - V.x)*(x - V.x) + (y - V.y)*(y - V.y);
    }

    //! Distance to a given vector
    T distTo(const Vector2<T>& V) const
    {
        return sqrt((x - V.x)*(x - V.x) + (y - V.y)*(y - V.y));
    }

    // operator
    const   Vector2<T>&    operator += (const Vector2<T>& V)
    {
        x += V.x;
        y += V.y;

        return *this;
    }
    const   Vector2<T>&    operator -= (const Vector2<T>& V)
    {
        x -= V.x;
        y -= V.y;

        return *this;
    }
    const   Vector2<T>&    operator *= (T k)
    {
        x *= k;
        y *= k;

        return *this;
    }
    const   Vector2<T>&    operator /= (T k)
    {
        x /= k;
        y /= k;

        return *this;
    }
    friend const Vector2<T> operator+ (const Vector2<T>& L, const Vector2<T>& R) { return Vector2<T>(L) += R; }
    friend const Vector2<T> operator- (const Vector2<T>& L, const Vector2<T>& R) { return Vector2<T>(L) -= R; }
    friend const Vector2<T> operator* (const Vector2<T>& L, T k)         { return Vector2<T>(L) *= k; }
    friend const Vector2<T> operator/ (const Vector2<T>& L, T k)         { return Vector2<T>(L) /= k; }
    friend const Vector2<T> operator* (T k, const Vector2<T>& R)         { return Vector2<T>(R) *= k; }

    friend const Vector2<T> operator- (const Vector2<T>& U)              { return Vector2<T>(-U.x, -U.y); }

    friend std::ostream& operator << (std::ostream& os, const Vector2<T>& V)
    {
        os << "(" << V.x << ", " << V.y  << ")";
        return os;
    }
public:
    T x;
    T y;
};

typedef Vector2<float>      Vector2f;
typedef Vector2<double>     Vector2d;
typedef Vector2<int>        Vector2i;

template <typename T>
void computeTriangleNormal(const Vector3<T>* v, Vector3<T>& n)
{
    Vector3<T> p, q;
    p.x = v[1].x - v[0].x;
    p.y = v[1].y - v[0].y;
    p.z = v[1].z - v[0].z;

    q.x = v[2].x - v[0].x;
    q.y = v[2].y - v[0].y;
    q.z = v[2].z - v[0].z;

    n = p.crossProduct(q);

    n.normalize();
}

template <typename T>
T computeTriangleArea(const Vector3<T>* v)
{
    Vector3<T> p, q;
    p.x = v[1].x - v[0].x;
    p.y = v[1].y - v[0].y;
    p.z = v[1].z - v[0].z;

    q.x = v[2].x - v[0].x;
    q.y = v[2].y - v[0].y;
    q.z = v[2].z - v[0].z;

    return (T)0.5*(p.crossProduct(q).norm());
}

template <typename T>
Vector3<T> transform3x4(T* M, Vector3<T> v)
{
    T tx, ty, tz;

    tx = M[0] * v.x + M[1] * v.y + M[2] * v.z + M[3];
    ty = M[4] * v.x + M[5] * v.y + M[6] * v.z + M[7];
    tz = M[8] * v.x + M[9] * v.y + M[10] * v.z + M[11];

    return Vector3<T>(tx, ty, tz);
}

template <typename T>
Vector3<T> transform3x3(T* M, Vector3<T> v)
{
    T tx, ty, tz;

    tx = M[0] * v.x + M[1] * v.y + M[2] * v.z;
    ty = M[3] * v.x + M[4] * v.y + M[5] * v.z;
    tz = M[6] * v.x + M[7] * v.y + M[8] * v.z;

    return Vector3<T>(tx, ty, tz);
}

}// end of namespace