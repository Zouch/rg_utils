#pragma once

//*****************************************************************************
bool rg_equals(float a, float b, float epsilon = 1e-5);

//*****************************************************************************
float rg_abs(float a);
float rg_min(float a, float b);
float rg_max(float a, float b);
float rg_sqrt(float a);

//*****************************************************************************
float rg_sin(float a);
float rg_asin(float a);
float rg_cos(float a);
float rg_acos(float a);
float rg_tan(float a);
float rg_atan(float a);

//*****************************************************************************
union rg_vec2 {
    struct { float u, v; };

    rg_vec2();
    rg_vec2(float s);
    rg_vec2(float a, float b);

    float* data();
    const float* data() const;

private:
    float _data[2];
};

//*****************************************************************************
union rg_vec3 {
    struct { float x, y, z; };
    struct { float r, g, b; };

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

bool operator==(rg_vec3 a, rg_vec3 b);
bool operator!=(rg_vec3 a, rg_vec3 b);

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

//*****************************************************************************
union rg_vec4 {
    struct { float x, y, z, w; };
    struct { float r, g, b, a; };
    struct { rg_vec3 xyz; float _w; };
    struct { rg_vec3 rgb; float _a; };

    rg_vec4();
    rg_vec4(float x);
    rg_vec4(float x, float y, float z, float w);

    float& operator[](int i);
    float  operator[](int i) const;

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

bool operator==(rg_vec4 a, rg_vec4 b);
bool operator!=(rg_vec4 a, rg_vec4 b);

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

//*****************************************************************************
struct rg_quat {
    float x, y, z, w;

    rg_quat();
    rg_quat(float x, float y, float z, float w);
};

//*****************************************************************************
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

//*****************************************************************************
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
    static rg_mat4 from_rotation(rg_quat r);
    static rg_mat4 from_scale(rg_vec3 s);
    static rg_mat4 from_trs(rg_vec3 t, rg_quat r, rg_vec3 s);

private:
    rg_vec4 _data[4];
};

rg_mat4 rg_look_at(rg_vec3 eye, rg_vec3 target, rg_vec3 up);
rg_mat4 rg_perspective(float fovy, float aspect, float znear, float zfar);

rg_mat4 operator*(const rg_mat4& a, const rg_mat4& b);

//*****************************************************************************
#define RG_MATH_IMPLEMENTATION
#if defined(RG_MATH_IMPLEMENTATION)

#include <cassert>
#include <cmath>
#include <limits>

//*****************************************************************************
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

float rg_sqrt(float a) {
    float result = std::sqrt(a);
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

//*****************************************************************************
rg_vec2::rg_vec2() : u(0.f), v(0.f) {}
rg_vec2::rg_vec2(float s) : u(s), v(s) {}
rg_vec2::rg_vec2(float a, float b) : u(a), v(b) {}

float* rg_vec2::data() { return _data; }
const float* rg_vec2::data() const { return _data; }

//*****************************************************************************
rg_vec3::rg_vec3() : x(0.f), y(0.f) {}
rg_vec3::rg_vec3(float s) : x(s), y(s), z(s) {}
rg_vec3::rg_vec3(float a, float b, float c) : x(a), y(b), z(c) {}

float& rg_vec3::operator[](int i) {
    return _data[i];
}

float  rg_vec3::operator[](int i) const {
    return _data[i];
}

rg_vec3& rg_vec3::operator+=(rg_vec3 v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

rg_vec3& rg_vec3::operator-=(rg_vec3 v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

rg_vec3& rg_vec3::operator*=(rg_vec3 v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    return *this;
}

rg_vec3& rg_vec3::operator/=(rg_vec3 v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    return *this;
}

rg_vec3& rg_vec3::operator+=(float s) {
    return operator+=(rg_vec3(s));
}

rg_vec3& rg_vec3::operator-=(float s) {
    return operator-=(rg_vec3(s));
}

rg_vec3& rg_vec3::operator*=(float s) {
    return operator*=(rg_vec3(s));
}

rg_vec3& rg_vec3::operator/=(float s) {
    return operator/=(rg_vec3(s));
}

float* rg_vec3::data() { return _data; }
const float* rg_vec3::data() const { return _data; }

rg_vec3 rg_vec3::axis_x() { return rg_vec3(1.f, 0.f, 0.f); }
rg_vec3 rg_vec3::axis_y() { return rg_vec3(0.f, 1.f, 0.f); }
rg_vec3 rg_vec3::axis_z() { return rg_vec3(0.f, 0.f, 1.f); }

bool operator==(rg_vec3 a, rg_vec3 b) {
    const float eps = std::numeric_limits<float>::epsilon();
    bool result = (rg_equals(a.x, b.x, eps) && rg_equals(a.y, b.y, eps));
    return result;
}

bool operator!=(rg_vec3 a, rg_vec3 b) {
    bool result = !(a == b);
    return result;
}

rg_vec3 operator+(rg_vec3 a, rg_vec3 b) {
    rg_vec3 result = a;
    result += b;
    return result;
}

rg_vec3 operator-(rg_vec3 a, rg_vec3 b) {
    rg_vec3 result = a;
    result -= b;
    return result;
}

rg_vec3 operator*(rg_vec3 a, rg_vec3 b) {
    rg_vec3 result = a;
    result *= b;
    return result;
}

rg_vec3 operator/(rg_vec3 a, rg_vec3 b) {
    rg_vec3 result = a;
    result /= b;
    return result;
}

rg_vec3 operator+(rg_vec3 v, float s) {
    rg_vec3 result = v;
    result += s;
    return result;
}

rg_vec3 operator-(rg_vec3 v, float s) {
    rg_vec3 result = v;
    result -= s;
    return result;
}

rg_vec3 operator*(rg_vec3 v, float s) {
    rg_vec3 result = v;
    result *= s;
    return result;
}

rg_vec3 operator/(rg_vec3 v, float s) {
    rg_vec3 result = v;
    result /= s;
    return result;
}

rg_vec3 operator+(float s, rg_vec3 v) {
    rg_vec3 result(s);
    result += v;
    return result;
}

rg_vec3 operator-(float s, rg_vec3 v) {
    rg_vec3 result(s);
    result -= v;
    return result;
}

rg_vec3 operator*(float s, rg_vec3 v) {
    rg_vec3 result(s);
    result *= v;
    return result;
}

rg_vec3 operator/(float s, rg_vec3 v) {
    rg_vec3 result(s);
    result /= v;
    return result;
}

float rg_dot(rg_vec3 a, rg_vec3 b) {
    float result = a.x * b.x + a.y * b.y + a.z * b.z;
    return result;
}

float rg_length_squared(rg_vec3 v) {
    float result = rg_dot(v, v);
    return result;
}

float rg_length(rg_vec3 v) {
    float result = rg_sqrt(rg_length_squared(v));
    return result;
}

float rg_inv_length(rg_vec3 v) {
    float result = 1.f / rg_length(v);
    return result;
}

rg_vec3 rg_normalize(rg_vec3 v) {
    rg_vec3 result = v * rg_inv_length(v);
    return result;
}

rg_vec3 rg_cross(rg_vec3 a, rg_vec3 b) {
    rg_vec3 result(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
    return result;
}

//*****************************************************************************

rg_vec4::rg_vec4() : x(0.f), y(0.f) {}
rg_vec4::rg_vec4(float s) : x(s), y(s), z(s) {}
rg_vec4::rg_vec4(float a, float b, float c, float d) : x(a), y(b), z(c), w(d) {}

float& rg_vec4::operator[](int i) {
    return _data[i];
}

float  rg_vec4::operator[](int i) const {
    return _data[i];
}

rg_vec4& rg_vec4::operator+=(rg_vec4 v) {
    x += v.x;
    y += v.y;
    z += v.z;
    w += v.w;
    return *this;
}

rg_vec4& rg_vec4::operator-=(rg_vec4 v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    w -= v.w;
    return *this;
}

rg_vec4& rg_vec4::operator*=(rg_vec4 v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;
    w *= v.w;
    return *this;
}

rg_vec4& rg_vec4::operator/=(rg_vec4 v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;
    w /= v.w;
    return *this;
}

rg_vec4& rg_vec4::operator+=(float s) {
    return operator+=(rg_vec4(s));
}

rg_vec4& rg_vec4::operator-=(float s) {
    return operator-=(rg_vec4(s));
}

rg_vec4& rg_vec4::operator*=(float s) {
    return operator*=(rg_vec4(s));
}

rg_vec4& rg_vec4::operator/=(float s) {
    return operator/=(rg_vec4(s));
}

float* rg_vec4::data() { return _data; }
const float* rg_vec4::data() const { return _data; }

rg_vec4 rg_vec4::axis_x() { return rg_vec4(1.f, 0.f, 0.f, 0.f); }
rg_vec4 rg_vec4::axis_y() { return rg_vec4(0.f, 1.f, 0.f, 0.f); }
rg_vec4 rg_vec4::axis_z() { return rg_vec4(0.f, 0.f, 1.f, 0.f); }
rg_vec4 rg_vec4::axis_w() { return rg_vec4(0.f, 0.f, 0.f, 1.f); }

bool operator==(rg_vec4 a, rg_vec4 b) {
    const float eps = std::numeric_limits<float>::epsilon();
    bool result = (rg_equals(a.x, b.x, eps) &&
                   rg_equals(a.y, b.y, eps) &&
                   rg_equals(a.z, b.z, eps) &&
                   rg_equals(a.w, b.w));
    return result;
}

bool operator!=(rg_vec4 a, rg_vec4 b) {
    bool result = !(a == b);
    return result;
}

rg_vec4 operator+(rg_vec4 a, rg_vec4 b) {
    rg_vec4 result = a;
    result += b;
    return result;
}

rg_vec4 operator-(rg_vec4 a, rg_vec4 b) {
    rg_vec4 result = a;
    result -= b;
    return result;
}

rg_vec4 operator*(rg_vec4 a, rg_vec4 b) {
    rg_vec4 result = a;
    result *= b;
    return result;
}

rg_vec4 operator/(rg_vec4 a, rg_vec4 b) {
    rg_vec4 result = a;
    result /= b;
    return result;
}

rg_vec4 operator+(rg_vec4 v, float s) {
    rg_vec4 result = v;
    result += s;
    return result;
}

rg_vec4 operator-(rg_vec4 v, float s) {
    rg_vec4 result = v;
    result -= s;
    return result;
}

rg_vec4 operator*(rg_vec4 v, float s) {
    rg_vec4 result = v;
    result *= s;
    return result;
}

rg_vec4 operator/(rg_vec4 v, float s) {
    rg_vec4 result = v;
    result /= s;
    return result;
}

rg_vec4 operator+(float s, rg_vec4 v) {
    rg_vec4 result(s);
    result += v;
    return result;
}

rg_vec4 operator-(float s, rg_vec4 v) {
    rg_vec4 result(s);
    result -= v;
    return result;
}

rg_vec4 operator*(float s, rg_vec4 v) {
    rg_vec4 result(s);
    result *= v;
    return result;
}

rg_vec4 operator/(float s, rg_vec4 v) {
    rg_vec4 result(s);
    result /= v;
    return result;
}

float rg_dot(rg_vec4 a, rg_vec4 b) {
    float result = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    return result;
}

float rg_length_squared(rg_vec4 v) {
    float result = rg_dot(v, v);
    return result;
}

float rg_length(rg_vec4 v) {
    float result = rg_sqrt(rg_length_squared(v));
    return result;
}

float rg_inv_length(rg_vec4 v) {
    float result = 1.f / rg_length(v);
    return result;
}

rg_vec4 rg_normalize(rg_vec4 v) {
    rg_vec4 result = v * rg_inv_length(v);
    return result;
}


//*****************************************************************************
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

//*****************************************************************************
rg_mat4::rg_mat4() {
    memset(_data, 0, 16 * sizeof(float));
}

rg_mat4::rg_mat4(float x) {
    memset(_data, 0, 16 * sizeof(float));
    for (int i = 0; i < 4; ++i) {
        _data[i][j] = x;
    }
}

rg_mat4::rg_mat4(rg_vec4 a, rg_vec4 b, rg_vec4 c, rg_vec4 d) {
    _data[0] = a;
    _data[1] = b;
    _data[2] = c;
    _data[3] = d;
}

rg_mat4::rg_mat4(float m11, float m21, float m31, float m41,
                 float m12, float m22, float m32, float m42,
                 float m13, float m23, float m33, float m43,
                 float m14, float m24, float m34, float m44) {
    _data[0][0] = m11; _data[1][0] = m21; _data[2][0] = m31; _data[3][0] = m41;
    _data[0][1] = m12; _data[1][1] = m22; _data[2][1] = m32; _data[3][1] = m42;
    _data[0][2] = m13; _data[1][2] = m23; _data[2][2] = m33; _data[3][2] = m43;
    _data[0][3] = m14; _data[1][3] = m24; _data[2][3] = m34; _data[3][3] = m44;
}

rg_vec4& rg_mat4::operator[](int i) {
    return _data[i];
}

rg_vec4  rg_mat4::operator[](int i) const {
    return _data[i];
}

float& rg_mat4::operator()(int i, int j) {
    return _data[i][j];
}

float  rg_mat4::operator()(int i, int j) const {
    return _data[i][j];
}

float* rg_mat4::data() { return _data[0].data(); }
const float* rg_mat4::data() const { return _data[0].data(); }

rg_mat4 rg_mat4::identity() {
    rg_mat4 result = rg_mat4(1.f);
    return result;
}

rg_mat4 rg_mat4::from_translation(rg_vec3 t) {
    rg_mat4 m = identity();
    m[3].xyz = t;
    return m;
}

rg_mat4 rg_mat4::from_rotation(rg_quat q) {
    rg_mat4 m1;
    m1[0] = rg_vec4( q.w, -q.z,  q.y, -q.x);
    m1[1] = rg_vec4( q.z,  q.w, -q.x, -q.y);
    m1[2] = rg_vec4(-q.y,  q.x,  q.w, -q.z);
    m1[3] = rg_vec4( q.x,  q.y,  q.z,  q.w);

    rg_mat4 m2;
    m2[0] = rg_vec4( q.w, -q.z,  q.y, q.x);
    m2[1] = rg_vec4( q.z,  q.w, -q.x, q.y);
    m2[2] = rg_vec4(-q.y,  q.x,  q.w, q.z);
    m2[3] = rg_vec4(-q.x, -q.y, -q.z, q.w);

    return m1 * m2;
}

rg_mat4 rg_mat4::from_scale(rg_vec3 s) {
    rg_mat4 m = identity();
    m[0][0] = s.x;
    m[1][1] = s.y;
    m[2][2] = s.z;
    return m;
}

rg_mat4 rg_mat4::from_trs(rg_vec3 t, rg_quat r, rg_vec3 s) {
    rg_mat4 T = from_translation(t);
    rg_mat4 R = from_rotation(r);
    rg_mat4 S = from_scale(s);

    return T * R * S;
}

rg_mat4 rg_look_at(rg_vec3 eye, rg_vec3 target, rg_vec3 up) {
    const rg_vec3 f = rg_normalize(target - eye);
    const rg_vec3 s = rg_normalize(rg_cross(f, up));
    const rg_vec3 u = rg_cross(s, f);

    rg_mat4 view = rg_mat4::identity();

    view[0][0] = s.x;
    view[1][0] = s.y;
    view[2][0] = s.z;
    view[0][1] = u.x;
    view[1][1] = u.y;
    view[2][1] = u.z;
    view[0][2] = -f.x;
    view[1][2] = -f.y;
    view[2][2] = -f.z;
    view[3][0] = -rg_dot(s, eye);
    view[3][1] = -rg_dot(u, eye);
    view[3][2] =  rg_dot(f, eye);

    return view;
}

rg_mat4 rg_perspective(float fovy, float aspect, float znear, float zfar) {
    rg_mat4 proj = rg_mat4::identity();

    const float tan_half_fov = rg_tan(fovy * 0.5f);
    const float fn = zfar - znear;

    proj[0][0] = 1.0f / (tan_half_fov * aspect);
    proj[1][1] = 1.0f / tan_half_fov;
    proj[2][2] = -zfar / fn;
    proj[3][2] = -(zfar * znear) / fn;
    proj[2][3] = -1.f;
    proj[3][3] = 0.0f;

    return proj;
}

inline rg_mat4 operator*(const rg_mat4& a, const rg_mat4& b)
{
    rg_mat4 result;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            float sum = 0;
            for (int k = 0; k < 4; ++k) {
                sum += a[k][i] * b[j][k];
            }
            result[i][j] = sum;
        }
    }

    return result;
}

//*****************************************************************************

#endif