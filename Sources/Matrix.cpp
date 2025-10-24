#include "Matrix.h"

// --- Vec3 类实现 ---

Vec3::Vec3(const Matrix2<double>& mat) : Vec3() {
	if (mat.GetRSize() == 1 && mat.GetCSize() == 3)
	{
		x = mat[0][0];
		y = mat[0][1];
		z = mat[0][2];
	}
}


// --- Matrix2 类实现 --- 
Matrix2<double> operator*(const Matrix2<double>& a, const Matrix2<double>& b)
{
	if (a.GetCSize() != b.GetRSize()) { cout << "矩阵乘法：行列值不匹配！\n"; }
	Matrix2<double> result(a.rsize, b.csize);
	for (int i = 0; i < a.rsize; i++)
	{
		for (int j = 0; j < b.csize; j++)
		{
			double sum = 0;
			for (int k = 0; k < a.csize; k++)
			{
				sum += a.data[i][k] * b.data[k][j];
			}
			result.data[i][j] = sum;
		}
	}
	return result;
}

template<typename T>
void Matrix2<T>::Print(string sep)const
{
	for (int i = 0; i < rsize; i++)
	{
		for (int j = 0; j < csize; j++)
		{
			cout << data[i][j];
			if (j < csize - 1) cout << sep;
		}
		cout << endl;
	}
}  

template<typename T>
Matrix2<T> Matrix2<T>::Trans()
{
	Matrix2<T> result(csize, rsize);
	for (int i = 0; i < rsize; i++)
	{
		for (int j = 0; j < csize; j++)
		{
			result.data[j][i] = data[i][j];
		}
	}
	return result;
}

// 显式实例化：为 double 类型生成模板定义（修复模板实现放在 .cpp 无法被使用时的链接/编译错误）
template class Matrix2<double>;
