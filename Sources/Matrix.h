#ifndef MATRIX_H
#define MATRIX_H

#include "types.h"
#include "Atom.h"


// --- 三维向量处理类 ---
class RadVec3;
class Vec3;

template<typename T>
class Matrix2;


class RadVec3 {
public:
	double r;
	double theta;
	double varphi;
	RadVec3() :r(0), theta(0), varphi(0) {}
	RadVec3(double rval, double thetaval, double varphival) :r(rval), theta(thetaval), varphi(varphival) {}
	RadVec3(const RadVec3& vec) :r(vec.r), theta(vec.theta), varphi(vec.varphi) {}
	void Print(string sep = "\t")const
	{
		cout << r << sep << theta << sep << varphi << endl;
	}
	friend ostream& operator<<(ostream& out, const RadVec3& v)
	{
		out << v.r << "," << v.theta << "," << v.varphi;
		return out;
	}
};

class Vec3 {
private:
	double len; //长度缓存，-1表示未计算

public:
	double x;
	double y;
	double z;

	Vec3() :x(0), y(0), z(0), len(0) {}
	Vec3(double xval, double yval, double zval) :x(xval), y(yval), z(zval) {
		len = sqrt(x * x + y * y + z * z);
	}
	Vec3(const Vec3& vec) :x(vec.x), y(vec.y), z(vec.z) {
		len = sqrt(x * x + y * y + z * z);
	}
	Vec3(const RadVec3& rvec) {
		*this = SpheToRect(rvec);
		len = rvec.r;
	}
	Vec3(const Matrix2<double>& mat);
	double GetLen()const
	{
		return len;
	}
	friend double Cdot(const Vec3& a, const Vec3& b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	friend Vec3 Cross(const Vec3& a, const Vec3& b)
	{
		return Vec3(a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x);
	}
	friend Vec3 operator-(const Vec3& a, const Vec3& b)
	{
		return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
	}
	friend Vec3 operator+(const Vec3& a, const Vec3& b)
	{
		return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
	}
	friend double GetAngle_Rad(Vec3& a, Vec3& b)
	{
		double val = Cdot(a, b) / (a.GetLen() * b.GetLen());
		return acos(val);
	}
	bool Normalize()
	{
		len = sqrt(x * x + y * y + z * z);
		if (len == 0) return false;
		x /= len;
		y /= len;
		z /= len;
		len = 1.0;
		return true;
	}

	friend ostream& operator<<(ostream& out, const Vec3& v)
	{
		out << v.x << "," << v.y << "," << v.z << endl;
		return out;
	}
	friend Vec3 operator*(const double& k, const Vec3& vec)
	{
		return Vec3(k * vec.x, k * vec.y, k * vec.z);
	}

	Vec3 SpheToRect(const RadVec3& sv)
	{
		double x = sv.r * sin(sv.varphi) * cos(sv.theta);
		double y = sv.r * sin(sv.varphi) * sin(sv.theta);
		double z = sv.r * cos(sv.varphi);
		return Vec3(x, y, z);
	}
	void Print(string sep = "\t")const
	{
		cout<<fixed<<setprecision(4) << x << sep << y << sep << z << endl;
	}


	Vec3(double len, double x, double y, double z)
		: len(len), x(x), y(y), z(z)
	{
	}
};

template<typename T>
class Matrix2 {
private:
	int rsize;
	int csize;
	vector<vector <T>> data;
public:
	Matrix2() :rsize(0), csize(0), data() {}
	Matrix2(int row_size, int col_size) :rsize(row_size), csize(col_size), data(rsize, vector<T>(csize, 0)) {}
	Matrix2(const Matrix2& m) :rsize(m.rsize), csize(m.csize), data(m.data) {}
	Matrix2(vector<vector<T>> d) :rsize(d.size()), csize(d[0].size()), data(d) {}
	Matrix2(const Vec3& vec) :csize(3), rsize(1), data({ {vec.x, vec.y, vec.z} }) {}
    const vector<T>& operator[](int i) const {
		return data[i];
    };
	friend Matrix2<double> operator*(const Matrix2<double>&, const Matrix2<double>&);

	Matrix2<T> Trans();


	int GetRSize()const { return rsize; }
	int GetCSize()const { return csize; }

	void Print(string sep = "\t")const;

};

inline double DegToRad(const double& deg)
{
	return deg * PI / 180.0;
}

class CoordSys3 {
public:
	Vec3 ox;
	Vec3 oy;
	Vec3 oz;
	CoordSys3() :ox(), oy(), oz() {}
	CoordSys3(Vec3 vec_oz) :ox(), oy(), oz(vec_oz)
	{
		oz.Normalize();
		if (oz.y == 0 && oz.z == 0) ox.y = 1;
		else ox.x = 1;
		double c = Cdot(ox, oz);
		if (c != 0)	ox.z = -c / oz.z;
		ox.Normalize();
		oy = Cross(oz, ox);
		oy.Normalize();
	}
	Matrix2<double> GetMatrixR()
	{
		vector<vector<double>> R = {
			{ox.x, ox.y, ox.z},
			{oy.x, oy.y, oy.z},
			{oz.x, oz.y, oz.z}
		};
		return R;
	}
};

#endif 