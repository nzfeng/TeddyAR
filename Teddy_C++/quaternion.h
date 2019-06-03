#pragma once

#include <string>
#include <vector>
#include <math.h>
#include "../Eigen/Dense"


/* A class for quaternions. Defines addition, subtraction, scalar 
 * multiplication, quaternion multiplication, and finding the conjugate, the 
 * norm, and inverse. Also has a function to return the corresponding 
 * rotation matrix.
 */
class Quaternion {

public:
	float s;
	Eigen::Vector3f v;

	// Initialize a new quaternion
	Quaternion() {
		s = 1.0;
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 0.0;
	}

	Quaternion(float r, float x, float y, float z) {
		s = r;
		v[0] = x;
		v[1] = y;
		v[2] = z;
	}

	Quaternion(float r, Eigen::Vector3f imag) {

		s = r;
		v = imag;
	}

	// Make a copy of an existing Quaternion object
	Quaternion(const Quaternion &q) {
		s = q.s;
		v = q.v;
	}

	// Quaternion assignment operator
	Quaternion & operator=(const Quaternion &q) {
		// if not self-assignment
		if (this != &q) {
			s = q.s;
			v = q.v;
		}
		return *this;
	}

	// Return true if two quaternions are equivalent.
	bool operator==(const Quaternion &q) const {
		return (s == q.s) && (v[0] == q.v[0]) 
						  && (v[1] == q.v[1]) 
						  && (v[2] == q.v[2]);
	}

	// Return conjugate.
	Quaternion conj() {
		return Quaternion(this->s, -this->v[0], -this->v[1], -this->v[2]);
	}

	// Return norm.
	float norm() {
		return sqrt(this->s*this->s + this->v[0]*this->v[0] \
			+ this->v[1]*this->v[1] + this->v[2]*this->v[2]);
	}

	// Quaternion addition.
	Quaternion operator+(const Quaternion &q) const {
		Quaternion sum;
		sum.s = this->s + q.s;
		sum.v[0] = this->v[0] + q.v[0];
		sum.v[1] = this->v[1] + q.v[1];
		sum.v[2] = this->v[2] + q.v[2];
		return sum;
	}

	// Quaternion subtraction.
	Quaternion operator-(const Quaternion &q) const {
		Quaternion diff;
		diff.s = this->s - q.s;
		diff.v[0] = this->v[0] - q.v[0];
		diff.v[1] = this->v[1] - q.v[1];
		diff.v[2] = this->v[2] - q.v[2];
		return diff;
	}

	// Scalar multiplication
	Quaternion operator*(float scalar) {
		Quaternion prod;
		prod.s *= scalar;
		prod.v *= scalar;
		return prod;
	}

	// Quaternion multiplication.
	Quaternion operator*(const Quaternion &q) const {
		Quaternion prod;
		prod.s = this->s * q.s - this->v.dot(q.v);
		prod.v = this->s * q.v + q.s * this->v + this->v.cross(q.v);
		return prod;
	}

	// Scalar division.
	Quaternion operator/(float scalar) {
		Quaternion quot;
		quot.s = this->s / scalar;
		quot.v = this->v / scalar;
		return quot;
	}

	// Return inverse.
	Quaternion inv() {
		Quaternion inverse = (*this).conj();
		inverse.s /= (*this).norm() * (*this).norm();
		inverse.v /= (*this).norm() * (*this).norm();
		return inverse;
	}

	// Return the corresponding rotation matrix.
	Eigen::Matrix4f toMatrix() {
		float qs = this->s;
		float qx = this->v[0];
		float qy = this->v[1];
		float qz = this->v[2];

		Eigen::Matrix4f R;
		R << 1.0 - 2.0*qy*qy - 2.0*qz*qz,
			 2.0 * (qx*qy - qz*qs),
			 2.0 * (qx*qz + qy*qs), 0.0,
			 2.0 * (qx*qy + qz*qs), 
			 1.0 - 2.0*qx*qx - 2.0*qz*qz,
			 2.0 * (qy*qz - qx*qs), 0.0,
			 2.0 * (qx*qz - qy*qs), 
			 2.0 * (qy*qz + qx*qs),
			 1.0 - 2.0*qx*qx - 2.0*qy*qy, 0.0,
			 0.0, 0.0, 0.0, 1.0;

		return R;
	}
};
