/*
 * cheby.h
 *
 *  Created on: Feb 16, 2020
 *      Author: ypmen
 */

#ifndef CHEBY_H_
#define CHEBY_H_

#include <iostream>
#include <string.h>

template <typename T>
T cheby1d_eval(T *coef, int n, T x)
{
	T b1 = 0.;
	T b2 = 0.;
	for (long int i=n-1; i>0; i--)
	{
		T b = 2.*x*b2 - b1 + coef[i];
		b1 = b2;
		b2 = b;
	}
	return x*b2 - b1 + 0.5*coef[0];
}

template <typename T>
T cheby2d_eval(T *coef, int m, int n, T x, T y)
{
	T b1 = 0.;
	T b2 = 0.;
	for (long int j=m-1; j>0; j--)
	{
		T b = 2.*y*b2 - b1 + cheby1d_eval(coef+j*n, n, x);
		b1 = b2;
		b2 = b;
	}
	return y*b2 - b1 + 0.5*cheby1d_eval(coef, n, x);
}

template <typename T>
void dcheby(T *coef, int m, int n)
{
	T *dcoef = new T [m*n];
	for (long int j=0; j<m; j++)
	{
		dcoef[j*n + n-1] = 0.;
		dcoef[j*n + n-2] = 2.*(n-1)*coef[j*n + n-1];
		for (long int i=n-3; i>=0; i--)
		{
			dcoef[j*n + i] = dcoef[j*n + i+2] + 2.*(i+1)*coef[j*n + i+1];
		}
	}
	memcpy(coef, dcoef, sizeof(T)*m*n);
	delete [] dcoef;
	dcoef = NULL;
}

#endif /* CHEBY_H_ */
