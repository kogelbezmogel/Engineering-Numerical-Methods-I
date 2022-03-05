#ifndef __MAXMIN__
#define __MAXMIN__


#include <iostream>
#include <vector>

typedef std::vector< std::vector<double> > Matrix;
typedef std::vector<double> Vector;

double max_mat(const Matrix& mat ) {
	double m = mat[100][40];
	for(int i = 0; i < mat.size(); ++i)
		for(int j = 0; j < mat[0].size(); ++j)
			if(m < mat[i][j] && mat[i][j] != 0)
				m = mat[i][j];
	return m;
}


double min_mat(const Matrix& mat ) {
	double m = mat[100][40];
	for(int i = 0; i < mat.size(); ++i)
		for(int j = 0; j < mat[0].size(); ++j)
			if(m > mat[i][j] && mat[i][j] != 0 )
				m = mat[i][j];
	return m;
}

#endif