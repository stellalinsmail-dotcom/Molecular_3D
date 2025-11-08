#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include "types.h"
//#include "Atom.h"


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
		out << v.x << "," << v.y << "," << v.z;
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
		//oz.Normalize();
		//if (oz.y == 0 && oz.z == 0) ox.y = 1;
		//else ox.x = 1;
		//double c = Cdot(ox, oz);
		//if (c != 0)	ox.z = -c / oz.z;
		//ox.Normalize();
		//oy = Cross(oz, ox);
		//oy.Normalize();

		oz.Normalize();
		if (oz.x == 0 && oz.y == 0)
		{
			if (oz.z > 0) ox.x = -1;
			else ox.x = 1;
		}
		else {
			ox.z = 1;
			double c = Cdot(ox, oz);
			if (oz.x != 0)	ox.x = -c / oz.x;
			else ox.y = -c / oz.y;
		}

		ox.Normalize();
		oy = Cross(oz, ox);
		oy.Normalize();
	}
	CoordSys3(Vec3 vec_oz,Vec3 ox_ref):ox(), oy(), oz(vec_oz)
	{
		oz.Normalize();
		oy = Cross(oz, ox_ref);
		oy.Normalize();
		ox = Cross(oy,oz);
		ox.Normalize();
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

// --- MMFF94参数表快表类定义 ---
template<typename TLine, typename TVal>
class FastMatrix1 {
private:
	vector<TVal> data;
public:
	FastMatrix1(const vector<TLine> table, vector<int> htb)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		data = vector<TVal>(dsize, TVal());
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			if (i < 0 || i > max_mtype) continue;
			int ti = htb[i];
			if (ti < 0) continue;
			data[ti] = line.GetVal();
		}
	}
	TVal operator[](const FAST_TABLE_INDEX& m)const
	{
		return data[m[0]];
	}
	int GetSize()const
	{
		return data.size();
	}
};

template<typename TLine, typename TVal>
class FastMatrix2 {
private:
	vector<vector<TVal>> data;
public:
	FastMatrix2(const vector<TLine> table, vector<int> htb)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		data = vector<vector<TVal>>(dsize, vector<TVal>(dsize, TVal()));
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			int j = line.GetMTypeJ();
			if (min(i, j) < 0 || max(i, j) > max_mtype) continue;
			int ti = htb[i];
			int tj = htb[j];
			if (min(ti, tj) < 0) continue;
			data[ti][tj] = line.GetVal();
			data[tj][ti] = line.GetVal();
		}
	}
	const TVal operator[](const FAST_TABLE_INDEX& m)const {
		return data[m[0]][m[1]];
	}
};

template<typename TLine, typename TVal>
class FastMatrix3 {
private:
	vector<vector<vector<TVal>>> data;
public:
	FastMatrix3(const vector<TLine> table, vector<int> htb)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		//cout <<"Max MType: " << max_mtype << endl;
		data = vector<vector<vector<TVal>>>(dsize, vector<vector<TVal>>(dsize, vector<TVal>(dsize, TVal())));
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			int j = line.GetMTypeJ();
			int k = line.GetMTypeK();
			if (Min3(i, j, k) < 0 || Max3(i, j, k) > max_mtype) continue;
			int ti = htb[i];
			int tj = htb[j];
			int tk = htb[k];
			if (Min3(ti, tj, tk) < 0) continue;
			data[ti][tj][tk] = line.GetVal();
			data[tk][tj][ti] = data[ti][tj][tk];
		}
	}
	const TVal operator[](const FAST_TABLE_INDEX& m)const {
		return data[m[0]][m[1]][m[2]];
	}
};
template<typename TLine, typename TVal>
class FastMatrix4 {
private:
	vector<vector<vector<vector<TVal>>>> data;
public:
	FastMatrix4(const vector<TLine> table, vector<int> htb, bool is_ij_solid = NO)
	{
		int dsize = htb.back() + 1;
		int max_mtype = htb.size() - 1;
		data = vector<vector<vector<vector<TVal>>>>(dsize, vector<vector<vector<TVal>>>(dsize, vector<vector<TVal>>(dsize, vector<TVal>(dsize, TVal()))));
		for (auto& line : table)
		{
			int i = line.GetMTypeI();
			int j = line.GetMTypeJ();
			int k = line.GetMTypeK();
			int l = line.GetMTypeL();
			if (Min4(i, j, k, l) < 0 || Max4(i, j, k, l) > max_mtype) continue;
			int ti = htb[i];
			int tj = htb[j];
			int tk = htb[k];
			int tl = htb[l];
			if (Min4(ti, tj, tk, tl) < 0) continue;
			if (is_ij_solid)
			{
				data[ti][tj][tk][tl] = line.GetVal();
				data[tk][tj][ti][tl] = data[ti][tj][tk][tl];
			}
			else
			{
				data[ti][tj][tk][tl] = line.GetVal();
				data[tl][tk][tj][ti] = line.GetVal();
			}
		}
	}
	const TVal operator[](const FAST_TABLE_INDEX& m)const {
		return data[m[0]][m[1]][m[2]][m[3]];
	}
};


#endif 