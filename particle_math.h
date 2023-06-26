#include <math.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#ifndef PHYSICS_MATH_H

#define FLT_MAX 3.402823466e+38F
#define FLT_MIN 1.175494351e-38F

struct V3 {
    float x, y, z;

    inline V3 &operator+=(V3 a);
    inline V3 &operator-=(V3 a);
};

struct Vertex4 {
    __m128 x;
    __m128 y;
    __m128 z;
};

struct Quaternion {
    float scalar; 
    V3 vector;
          
};

struct Quaternion4 {
    __m128 w;
    Vertex4 vector;
};


 
static inline V3 v3(float x, float y, float z) {
    V3 result;
    result.x = x;
    result.y = y;
    result.z = z;

    return result;
}

static inline V3 v3_zero() {
    V3 result = {};
    return result;
}

static inline V3 operator+(V3 a, V3 b) {
    V3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

static inline V3 operator-(V3 a, V3 b) {
    V3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

static inline V3 operator-(V3 a) {
    V3 result;

    result.x = -a.x;
    result.y = -a.y;
    result.z = -a.z;

    return result;
}

inline V3 &V3:: operator+=(V3 a) {
    *this = *this + a;

    return *this;
}

inline V3 &V3:: operator-=(V3 a) {
    *this = *this - a;

    return *this;
}
 
inline V3 operator*(float a, V3 b) {
    V3 result;
    result.x = a*b.x;
    result.y = a*b.y;
    result.z = a*b.z;

    return result;
}

inline V3 operator/(V3 a, float b) {
    V3 result;
    result.x = a.x/b;
    result.y = a.y/b;
    result.z = a.z/b;

    return result;
}

static inline float v3_mag(V3 a) {
    return sqrtf((a.x * a.x) + (a.y * a.y) + (a.z * a.z));
}

static inline V3 v3_norm(V3 a) {
    V3 result;
    if (a.x == 0 && a.y == 0 && a.z == 0) return v3_zero();
    result = a/v3_mag(a);
    return result;
}

static inline float v3_dot(V3 a, V3 b) {
    return ((a.x*b.x) + (a.y*b.y) + (a.z*b.z));
}

static inline V3 v3_cross(V3 a, V3 b) {
    V3 result;
    result.x = a.y*b.z - a.z*b.y;
    result.y = a.z*b.x - a.x*b.z;
    result.z = a.x*b.y - a.y*b.x;

    return result;
}

static inline V3 v3_pairwise_mul(V3 a, V3 b) {
    V3 result;
    result.x = a.x*b.x;
    result.y = a.y*b.y;
    result.z = a.z*b.z;
    return result;
}

inline V3 v3_pariwise_div(V3 a, V3 b) {
    V3 result;
    result.x = a.x/b.x;
    result.y = a.y/b.y;
    result.z = a.z/b.z;

    return result;
}

static void v3_swap(V3 *a, V3 *b) {
    V3 temp = *a;
    *a = *b;
    *b = temp;
}

static inline V3 v3_lerp(float delta_x, float lerp_factor, V3 y0, V3 y1) {
    return y0 + (lerp_factor)*((y1-y0)/(delta_x));
}



static Quaternion quat4(float w, float x, float y, float z) {
    Quaternion result;
    result.scalar = w;
    result.vector.x = x;
    result.vector.y = y;
    result.vector.z = z;

    return result;
}

static Quaternion q4(float w, V3 vec) {
    Quaternion result;
    result.scalar = w;
    result.vector = vec;

    return result;
}

static Quaternion q4_identity() {
    Quaternion result;
    result.scalar = 1;
    result.vector = {};
    return result;
}

static Quaternion rotation_q4(float theta, V3 axis) {
    Quaternion result;
    result.scalar = cosf(theta/2.0f);
    result.vector = sinf(theta/2.0f)*axis;

    return result;
}

static Quaternion operator*(Quaternion a, Quaternion b) {
    Quaternion result; 
    result.vector = 
        a.scalar*b.vector + b.scalar*a.vector + (v3_cross(a.vector, b.vector));
    result.scalar = (a.scalar*b.scalar) - v3_dot(a.vector, b.vector);

    return result;
}

static Quaternion operator*(Quaternion a, V3 b) {
    Quaternion result; 
    result.vector = 
        a.scalar*b + (v3_cross(a.vector, b));
    result.scalar = -v3_dot(a.vector, b);

    return result;
}

static Quaternion q4_conj(Quaternion a) {
    Quaternion result; 

    result.scalar = a.scalar;
    result.vector = -a.vector;
    return result;
}

static Quaternion q4_norm(Quaternion a) {
    Quaternion result;
    float mag = sqrtf(a.scalar * a.scalar +
          a.vector.x*a.vector.x + 
          a.vector.y*a.vector.y +
          a.vector.z*a.vector.z);

    result.scalar = a.scalar / mag;
    result.vector = a.vector / mag;
    return result;
}

static V3 v3_rotate_on_axis(V3 axis, float theta, V3 vec) {
    Quaternion result;
    Quaternion rotation = q4_norm(q4(cosf(theta/2), sinf(theta/2)*axis));
    result = (rotation*vec*q4_conj(rotation));
    return result.vector;
}

static inline V3 v3_rotate_q4(V3 vec, Quaternion q) {
    Quaternion result;
    result = q*vec;
    V3 conjugate_vec = -q.vector;
    //result.vector = result.scalar*conjugate_vec + (v3_cross(result.vector, conjugate_vec));
    result.vector = result.scalar*conjugate_vec + q.scalar*result.vector + (v3_cross(result.vector, conjugate_vec));
    return result.vector;
}

static inline float f_lerp(float delta_x, float lerp_factor, float y0, float y1) {
    return y0 + (lerp_factor)*((y1-y0)/(delta_x));
}


static inline Vertex4 vertex4_scale(Vertex4 _a, __m128 _scalar) {
    Vertex4 result;
    result.x = _mm_mul_ps(_a.x, _scalar);
    result.y = _mm_mul_ps(_a.y, _scalar);
    result.z = _mm_mul_ps(_a.z, _scalar);

    return result;
}

static inline Vertex4 vertex4_pairwise_mul(Vertex4 _a, Vertex4 _b) {
    Vertex4 result;
    result.x = _mm_mul_ps(_a.x, _b.x);
    result.y = _mm_mul_ps(_a.y, _b.y);
    result.z = _mm_mul_ps(_a.z, _b.z);

    return result;
}

static inline Vertex4 vertex4_cross(Vertex4 _a, Vertex4 _b) {
    Vertex4 result;
    result.x = _mm_sub_ps(_mm_mul_ps(_a.y, _b.z), _mm_mul_ps(_a.z, _b.y));
    result.y = _mm_sub_ps(_mm_mul_ps(_a.z, _b.x), _mm_mul_ps(_a.x, _b.z));
    result.z = _mm_sub_ps(_mm_mul_ps(_a.x, _b.y), _mm_mul_ps(_a.y, _b.x));
    return result;

}

static inline __m128 vertex4_dot(Vertex4 _a, Vertex4 _b) {
    __m128 result;
    result = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_a.x, _b.x), _mm_mul_ps(_a.y, _b.y)),
                        _mm_mul_ps(_a.z, _b.z));
    return result;
}

static inline Vertex4 vertex4_norm(Vertex4 _a) {
    Vertex4 result;
    __m128 _inv_magnitude = _mm_rsqrt_ps(vertex4_dot(_a, _a));
    result = vertex4_scale(_a, _inv_magnitude);
    return result;
}


static inline Vertex4 vertex4_add(Vertex4 _a, Vertex4 _b) {
    Vertex4 result;

    result.x = _mm_add_ps(_a.x, _b.x);
    result.y = _mm_add_ps(_a.y, _b.y);
    result.z = _mm_add_ps(_a.z, _b.z);

    return result;
}

static inline Vertex4 vertex4_sub(Vertex4 _a, Vertex4 _b) {
    Vertex4 result;

    result.x = _mm_sub_ps(_a.x, _b.x);
    result.y = _mm_sub_ps(_a.y, _b.y);
    result.z = _mm_sub_ps(_a.z, _b.z);

    return result;
}

#define PHYSICS_MATH_H
#endif
