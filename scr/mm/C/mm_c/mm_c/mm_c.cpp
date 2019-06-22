// mm_c.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include <stdio.h>
#include <math.h>

#define DLLEXPORT extern "C" __declspec(dllexport) 

const double Pi = 3.14159265358979323846264338328;

DLLEXPORT int __stdcall show_matrix(double *matrix, int rows, int columns)
{
	int i, j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			printf("matrix[%d][%d] = %d\n", i, j, matrix[i*rows + j]);
		}
	}
	return 0;
}

double *cross(double *v1, double *v2)
{
	double ans[3];
	ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
	ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
	ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return ans;
};

double dot(double *v1, double *v2)
{
	double ans;
	ans = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	return ans;
};

double potential_dihedral(double *array1, double *array2, double *array3, double *array4, double *dihedraldata_1, int len)
{
	double delta1[3];
	for (int i = 0; i < 3; i++)delta1[i] = array2[1] - array1[i];
	double delta2[3];
	for (int i = 0; i < 3; i++)delta2[i] = array3[1] - array2[i];
	double delta3[3];
	for (int i = 0; i < 3; i++)delta2[i] = array4[1] - array3[i];
	
	*delta1 = *cross(delta1, delta2);
	*delta3 = *cross(delta2, delta3);

	double cosa = dot(delta1, delta3) / sqrt(dot(delta1, delta1) *dot(delta3, delta3));
	double phi;
	
	if (-1 < cosa && cosa < 1)
	{
		phi = acos(cosa);
	}
	else
	{
		if (cosa > 0) { phi = 0; }
		else { phi = Pi; };
	};

	double sum = 0;
	for (int i = 0; i < len; i++) sum += dihedraldata_1[3*i+2] * cos(dihedraldata_1[3*i+0] * phi - dihedraldata_1[3*i+1]);
	return sum;
};

DLLEXPORT void __stdcall fdihedral_c(double *array1, double *array2, double *array3, double *array4, double *dihedraldata, int len, double *ans)
{
	double array2_2[3] = { array2[0] + 0.00001, array2[1] + 0.00001, array2[2] + 0.00001 };
	for (int i = 0; i < 3; i++) {
		double array2_3[3] = { array2[0], array2[1], array2[2] };
		array2_3[i] = array2_2[i];
		ans[i] = -(potential_dihedral(array1, array2_3, array3, array4, dihedraldata, len) - potential_dihedral(array1, array2, array3, array4, dihedraldata, len)) / 0.00001;
	};
};

DLLEXPORT void __stdcall fdihedral_s(double *array1, double *array2, double *array3, double *array4, double *dihedraldata, int len, double *ans)
{
	double array1_2[3] = { array1[0] + 0.00001, array1[1] + 0.00001, array1[2] + 0.00001 };
	for (int i = 0; i < 3; i++) {
		double array1_3[3] = { array1[0], array1[1], array1[2] };
		array1_3[i] = array1_2[i];
		ans[i] = -(potential_dihedral(array1_3, array2, array3, array4, dihedraldata, len) - potential_dihedral(array1, array2, array3, array4, dihedraldata, len)) / 0.00001;
	};
};