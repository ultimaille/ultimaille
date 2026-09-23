#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cassert>

#include "vec.h"
#include "mat.h"

namespace UM {
    struct Quaternion;
    Quaternion operator*(const Quaternion &q, const double &val);
    Quaternion operator/(const Quaternion &q, const double &val);

    struct Quaternion {
//      Quaternion(const vec3 axis = {0, 0, 0}, const double alpha = 0) : v(axis*sin(alpha*.5)), w(cos(alpha*.5)) {}

        double norm2() const {
            return v.norm2() + w*w;
        }

        double norm() const {
            return std::sqrt(norm2());
        }

        Quaternion &normalize() { *this = (*this)/norm(); return *this; }
        Quaternion normalized() const { return *this/norm(); }

        Quaternion conjugate() const {
            return { -v, w };
        }

        vec3 rotate(const vec3& p) const;

        double &operator[](const int i)       { assert(i>=0 && i<4); return i<3 ? v[i] : w; }
        double  operator[](const int i) const { assert(i>=0 && i<4); return i<3 ? v[i] : w; }

        // convert to a 3x3 rotation matrix
        // https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
        mat<3,3> rotation_matrix() const {
            double n = norm2();
            double s = n>1e-20 ? 2./n : 0.;

            double wx = s * w * v.x;
            double wy = s * w * v.y;
            double wz = s * w * v.z;

            double xx = s * v.x * v.x;
            double xy = s * v.x * v.y;
            double xz = s * v.x * v.z;

            double yy = s * v.y * v.y;
            double yz = s * v.y * v.z;
            double zz = s * v.z * v.z;

            return {{{1.-yy-zz, xy-wz, xz+wy}, {xy+wz, 1.-xx-zz, yz-wx}, {xz-wy, yz+wx, 1.-xx-yy}}};
        }

        // z-y'-x" intrinsic Tait–Bryan angles: yaw, pitch and roll
        // to rotate a 3d column vector v, compute Rz(yaw) x Ry(pitch) x Rx(roll) x v,
        // where Rx, Ry and Rz are the basic rotation matrices following the right-hand rule
        vec3 euler_angles() const {
            return {
                atan2(2.*(w*v.x + v.y*v.z), 1.-2*(v.x*v.x + v.y*v.y)), // roll
                asin(2.*(w*v.y - v.z*v.x)),                            // pitch
                atan2(2*(w*v.z+v.x*v.y), 1.-2.*(v.y*v.y + v.z*v.z))    // yaw
                };
        }

        static Quaternion from_axis_angle(const vec3& axis, double angle) {
            assert(axis.norm2() > 1e-20);
            return { axis.normalized() * std::sin(angle/2), std::cos(angle/2) };
        }

	static Quaternion rotation_between(vec3 v0, vec3 v1) {
		const vec3 axis = cross(v0, v1);
		const double d = v0 * v1;
		if (d >  1 - 1e-12) return {}; // almost same vectors => no rotation
		if (d < -1 + 1e-12) {          // antipodal case: choose a deterministic perpendicular axis.
			axis = cross(v0, vec3{1, 0, 0});
			if (axis.norm2() < 1e-12)
				axis = cross(v0, vec3{0, 1, 0});
			return Quaternion::from_axis_angle(axis.normalized(), M_PI);
		}
		return Quaternion{axis, 1 + d}.normalized();
	}

        vec3 v = {0., 0., 0.};
        double w = {1.};
    };

    inline Quaternion operator*(const Quaternion &q, const double &val) {
        Quaternion ret = q;
        for (int i=4; i--; ret[i]*=val);
        return ret;
    }

    inline Quaternion operator/(const Quaternion &q, const double &val) {
        Quaternion ret = q;
        for (int i=4; i--; ret[i]/=val);
        return ret;
    }

    inline Quaternion operator*(const Quaternion &a, const Quaternion &b) {
        Quaternion res;
        res.w = a.w*b.w - a.v*b.v;
        res.v = a.w*b.v + b.w*a.v + cross(a.v, b.v);
        return res;
    }

    inline vec3 Quaternion::rotate(const vec3& p) const {
        Quaternion q = normalized();
        return (q * Quaternion{p, 0} * q.conjugate()).v;
    }


}
#endif //__QUATERNION_H__

