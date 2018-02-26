#pragma once

////////////////////////////////////////////////////////////////////
bool rg_equals(float a, float b, float epsilon = 1e-5);

float rg_abs(float a);
float rg_min(float a, float b);
float rg_max(float a, float b);
float rg_sin(float a);

float rg_sin(float a);
float rg_asin(float a);
float rg_cos(float a);
float rg_acos(float a);
float rg_tan(float a);
float rg_atan(float a);

////////////////////////////////////////////////////////////////////
union rg_vec2 {
    struct { float x, y; };
    struct { float u, v; };

    rg_vec2();
    rg_vec2(float x);
    rg_vec2(float x, float y);

    float& operator[](int i);
    float  operator[](int i) const;

    rg_vec2& operator+=(rg_vec2 v);
    rg_vec2& operator-=(rg_vec2 v);
    rg_vec2& operator*=(rg_vec2 v);
    rg_vec2& operator/=(rg_vec2 v);

    rg_vec2& operator+=(float s);
    rg_vec2& operator-=(float s);
    rg_vec2& operator*=(float s);
    rg_vec2& operator/=(float s);

    float* data();
    const float* data() const;

    static rg_vec2 axis_x();
    static rg_vec2 axis_y();

private:
    float _data[2];
};

bool operator==(rg_vec2 a, rg_vec2 b);
bool operator!=(rg_vec2 a, rg_vec2 b);

rg_vec2 operator+(rg_vec2 a, rg_vec2 b);
rg_vec2 operator-(rg_vec2 a, rg_vec2 b);
rg_vec2 operator*(rg_vec2 a, rg_vec2 b);
rg_vec2 operator/(rg_vec2 a, rg_vec2 b);

rg_vec2 operator+(float s, rg_vec2 v);
rg_vec2 operator-(float s, rg_vec2 v);
rg_vec2 operator*(float s, rg_vec2 v);
rg_vec2 operator/(float s, rg_vec2 v);

rg_vec2 operator+(rg_vec2 v, float s);
rg_vec2 operator-(rg_vec2 v, float s);
rg_vec2 operator*(rg_vec2 v, float s);
rg_vec2 operator/(rg_vec2 a, float s);

float rg_dot(rg_vec2 a, rg_vec2 b);
float rg_length_squared(rg_vec2 v);
float rg_length(rg_vec2 v);
float rg_inv_length(rg_vec2 v);
rg_vec2 rg_normalize(rg_vec2 v);

union rg_vec3 {
    struct { float x, y, z; };
    struct { float u, v, s; };
    struct { float r, g, b; };
    struct { rg_vec2 xy; float _z; };
    struct { rg_vec2 uv; float _s; };
    struct { float _x; rg_vec2 yz; };
    struct { float _u; rg_vec2 vs; };

    rg_vec3();
    rg_vec3(float x);
    rg_vec3(float x, float y, float z);

    float& operator[](int i);
    float  operator[](int i) const;

    rg_vec3& operator+=(rg_vec3 v);
    rg_vec3& operator-=(rg_vec3 v);
    rg_vec3& operator*=(rg_vec3 v);
    rg_vec3& operator/=(rg_vec3 v);

    rg_vec3& operator+=(float s);
    rg_vec3& operator-=(float s);
    rg_vec3& operator*=(float s);
    rg_vec3& operator/=(float s);

    float* data();
    const float* data() const;

    static rg_vec3 axis_x();
    static rg_vec3 axis_y();
    static rg_vec3 axis_z();

private:
    float _data[3];
};

rg_vec3 operator+(rg_vec3 a, rg_vec3 b);
rg_vec3 operator-(rg_vec3 a, rg_vec3 b);
rg_vec3 operator*(rg_vec3 a, rg_vec3 b);
rg_vec3 operator/(rg_vec3 a, rg_vec3 b);

rg_vec3 operator+(rg_vec3 v, float s);
rg_vec3 operator-(rg_vec3 v, float s);
rg_vec3 operator*(rg_vec3 v, float s);
rg_vec3 operator/(rg_vec3 v, float s);

rg_vec3 operator+(float s, rg_vec3 v);
rg_vec3 operator-(float s, rg_vec3 v);
rg_vec3 operator*(float s, rg_vec3 v);
rg_vec3 operator/(float s, rg_vec3 v);

float rg_dot(rg_vec3 a, rg_vec3 b);
float rg_length_squared(rg_vec3 v);
float rg_length(rg_vec3 v);
float rg_inv_length(rg_vec3 v);
rg_vec3 rg_normalize(rg_vec3 v);
rg_vec3 rg_cross(rg_vec3 a, rg_vec3 b);

union rg_vec4 {
    struct { float x, y, z, w; };
    struct { float u, v, s, t; };
    struct { float r, g, b, a; };
    struct { rg_vec2 xy; rg_vec2 zw; };
    struct { rg_vec2 uv; rg_vec2 st; };
    struct { float _x; rg_vec2 yz; float _w; };
    struct { float _u; rg_vec2 vs; float _t; };
    struct { rg_vec3 xyz; float _w; };
    struct { rg_vec3 uvs; float _t; };
    struct { rg_vec3 rgb; float _a; };
    struct { float _x; rg_vec3 yzw; };
    struct { float _u; rg_vec3 vst; };

    rg_vec4();
    rg_vec4(float x);
    rg_vec4(float x, float y, float z, float w);

    float& operator[](int i);
    float  operator[](int i) const;

    bool operator==(rg_vec4 v);
    bool operator!=(rg_vec4 v);

    rg_vec4& operator+=(rg_vec4 v);
    rg_vec4& operator-=(rg_vec4 v);
    rg_vec4& operator*=(rg_vec4 v);
    rg_vec4& operator/=(rg_vec4 v);

    rg_vec4& operator+=(float s);
    rg_vec4& operator-=(float s);
    rg_vec4& operator*=(float s);
    rg_vec4& operator/=(float s);

    float* data();
    const float* data() const;

    static rg_vec4 axis_x();
    static rg_vec4 axis_y();
    static rg_vec4 axis_z();
    static rg_vec4 axis_w();

private:
    float _data[4];
};

rg_vec4 operator+(rg_vec4 a, rg_vec4 b);
rg_vec4 operator-(rg_vec4 a, rg_vec4 b);
rg_vec4 operator*(rg_vec4 a, rg_vec4 b);
rg_vec4 operator/(rg_vec4 a, rg_vec4 b);

rg_vec4 operator+(rg_vec4 v, float s);
rg_vec4 operator-(rg_vec4 v, float s);
rg_vec4 operator*(rg_vec4 v, float s);
rg_vec4 operator/(rg_vec4 v, float s);

rg_vec4 operator+(float s, rg_vec4 v);
rg_vec4 operator-(float s, rg_vec4 v);
rg_vec4 operator*(float s, rg_vec4 v);
rg_vec4 operator/(float s, rg_vec4 v);

float rg_dot(rg_vec4 a, rg_vec4 b);
float rg_length_squared(rg_vec4 v);
float rg_length(rg_vec4 v);
float rg_inv_length(rg_vec4 v);
rg_vec4 rg_normalize(rg_vec4 v);

struct rg_quat {
    float x, y, z, w;

    rg_quat();
    rg_quat(float x, float y, float z, float w);
};

struct rg_mat2 {
    rg_mat2();
    rg_mat2(float x);
    rg_mat2(rg_vec2 a, rg_vec2 b);
    rg_mat2(float m11, float m21, float m12, float m22);

    rg_vec2& operator[](int i);
    rg_vec2  operator[](int i) const;

    float& operator()(int i, int j);
    float  operator()(int i, int j) const;

    float* data();
    const float* data() const;

    static rg_mat2 identity();

private:
    rg_vec2 _data[2];
};

struct rg_mat3 {
    rg_mat3();
    rg_mat3(float x);
    rg_mat3(rg_vec3 a, rg_vec3 b, rg_vec3 c);
    rg_mat3(float m11, float m21, float m31,
            float m12, float m22, float m32,
            float m13, float m23, float m33);

    rg_vec3& operator[](int i);
    rg_vec3  operator[](int i) const;

    float& operator()(int i, int j);
    float  operator()(int i, int j) const;

    float* data();
    const float* data() const;

    static rg_mat3 identity();

private:
    rg_vec3 _data[3];
};

struct rg_mat4 {
    rg_mat4();
    rg_mat4(float x);
    rg_mat4(rg_vec4 a, rg_vec4 b, rg_vec4 c, rg_vec4 d);
    rg_mat4(float m11, float m21, float m31, float m41,
            float m12, float m22, float m32, float m42,
            float m13, float m23, float m33, float m43,
            float m14, float m24, float m34, float m44);

    rg_vec4& operator[](int i);
    rg_vec4  operator[](int i) const;

    float& operator()(int i, int j);
    float  operator()(int i, int j) const;

    float* data();
    const float* data() const;

    static rg_mat4 identity();
    static rg_mat4 from_translation(rg_vec3 t);
    static rg_mat4 from_rotation(rg_vec4 r);
    static rg_mat4 from_scale(rg_vec3 s);
    static rg_mat4 from_trs(rg_vec3 t, rg_vec3 r, rg_vec3 s);

private:
    rg_vec4 _data[4];
};

rg_mat4 rg_look_at(rg_vec3 eye, rg_vec3 target, rg_vec3 up);
rg_mat4 rg_perspective(float fovy, float aspect, float znear, float zfar);

rg_mat4 operator*(const rg_mat4& a, const rg_mat4& b);

#define RG_MATH_IMPLEMENTATION
#if defined(RG_MATH_IMPLEMENTATION)

#include <cassert>
#include <cmath>
#include <limits>

////////////////////////////////////////////////////////////////////
bool rg_equals(float a, float b, float epsilon) {
    bool result = (rg_abs(a - b) < epsilon);
    return result;
}

float rg_abs(float a) {
    float result = a >= 0.f ? a : -a;
    return result;
}

float rg_min(float a, float b) {
    float result = a < b ? a : b;
    return result;
}

float rg_max(float a, float b) {
    float result = a > b ? a : b;
    return result;
}

float rg_sin(float a) {
    float result = std::sin(a);
    return result;
}

float rg_asin(float a) {
    float result = std::asin(a);
    return result;
}

float rg_cos(float a) {
    float result = std::cos(a);
    return result;
}

float rg_acos(float a) {
    float result = std::acos(a);
    return result;
}

float rg_tan(float a) {
    float result = std::tan(a);
    return result;
}

float rg_atan(float a) {
    float result = std::atan(a);
    return result;
}

////////////////////////////////////////////////////////////////////
rg_vec2::rg_vec2() : x(0.f), y(0.f) {}
rg_vec2::rg_vec2(float s) : x(s), y(s) {}
rg_vec2::rg_vec2(float a, float b) : x(a), y(b) {}

float& rg_vec2::operator[](int i) {
    return _data[i];
}

float  rg_vec2::operator[](int i) const {
    return _data[i];
}

bool operator==(rg_vec2 a, rg_vec2 b) {
    const float eps = std::numeric_limits<float>::epsilon();
    bool result = (rg_equals(a.x, b.x, eps) && rg_equals(a.y, b.y, eps));
    return result;
}

bool operator!=(rg_vec2 a, rg_vec2 b) {
    bool result = !(a == b);
    return result;
}

#endif