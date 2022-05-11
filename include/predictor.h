/*
 * predictor.h
 *
 *  Created on: Feb 16, 2020
 *      Author: ypmen
 */

#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include <string.h>
#include <vector>

using namespace std;

class Predictor
{
public:
	Predictor();
	Predictor(const Predictor &pred);
	Predictor & operator=(const Predictor &pred);
	~Predictor();
	void set_coef(long double *coefs, int m, int n);
	long double get_phase(long double mjd, long double freq);
	double get_pfold(long double mjd, long double freq);
	double get_fdfold(long double mjd, long double freq);
	long double get_phase_inverse(long double phase, long double freq);
private:
	void compute_dcoef();
	void compute_ddcoef();
public:
	long double dispersion_constant;
	long double mjd_low;
	long double mjd_up;
	long double freq_low;
	long double freq_up;
	string psrname;
	string sitename;
private:
	int ncoef_t;
	int ncoef_f;
	long double *coef;
	long double *dcoef;
	long double *ddcoef;
};

class Predictors
{
public:
	Predictors();
	Predictors(const Predictors &preds);
	Predictors & operator=(const Predictors &preds);
	Predictors(const string fname);
	~Predictors();
	void read_t2pred(const string fname);
	long double get_phase(long double mjd, long double freq);
	double get_fphase(long double mjd, long double freq);
	double get_pfold(long double mjd, long double freq);
	double get_fdfold(long double mjd, long double freq);
	long double get_phase_inverse(long double phase, long double freq);
	long double get_center_frequency();
	bool empty()
	{
		if (pred == NULL)
			return true;
		return false;
	}
public:
	string filename;
private:
	int get_nearest(long double mjd);
private:
	int npred;
	Predictor *pred;

};

#endif /* PREDICTOR_H_ */
