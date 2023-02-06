/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-07-12 11:45:14
 * @modify date 2022-07-12 11:45:14
 * @desc [function set for orbital solver]
 */

#include "avx2.h"

#include <iostream>
#include <iomanip>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

extern "C"
{
	void dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda,
				double *b, int *ldb, double *work, int *lwork, int *info);
}

py::array_t<double> lmfit(std::function<py::array_t<double> (py::array_t<double>, py::array_t<double>)> &func, py::array_t<double> &x, py::array_t<double> &y, py::array_t<double> &p0, py::array_t<double> &sigma, py::array_t<double> &diff_step, py::array_t<bool> &if_fit);

//#undef __AVX2__

#ifndef __AVX2__

#include <cmath>

py::array_t<double> compute_f0_f1(py::array_t<double> &vT, double f0, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_T = vT.unchecked<1>();
	size_t size = array_T.size();

	double *vf0_f1 = new double[size * 2];

	double sqrt_1_ee = std::sqrt(1. - e * e);

	for (size_t i = 0; i < size; i++)
	{
		double T = array_T(i);
		// calculate fractional orbit
		double orbit = (T - T0) / Pb;
		double norbit = std::floor(orbit);
		double forbit = orbit - norbit;

		// calculate eccentric anomaly
		double phase = 2. * M_PI * forbit;
		double E = phase + e * std::sin(phase) / std::sqrt(1. - 2. * e * std::cos(phase) + e * e);
		double dE = 1.;
		size_t nstep = 0;
		while (std::abs(dE) > 1e-14 && nstep++ < 8)
		{
			dE = (phase - (E - e * std::sin(E))) / (1. - e * std::cos(E));
			E += dE;
		}

		double sin_E = std::sin(E);
		double cos_E = std::cos(E);

		// calculate omega
		double omega = omega0 + omega_dot / (365.25 / Pb) * (norbit + forbit);

		// calculate v/c, a/c
		double sin_omega = std::sin(omega);
		double cos_omega = std::cos(omega);

		double temp = 2 * M_PI / (Pb * 86400. * (1. - e * cos_E));

		double doppler_factor1 = -asini_c * (-sin_omega * sin_E + sqrt_1_ee * cos_omega * cos_E) * temp;
		double doppler_factor2 = -asini_c * (-e * sin_omega + sin_omega * cos_E + sqrt_1_ee * cos_omega * sin_E) / (e * cos_E - 1.) * std::pow(temp, 2.);

		vf0_f1[i * 2 + 0] = f0 * (1. + doppler_factor1);
		vf0_f1[i * 2 + 1] = f0 * doppler_factor2;
	}

	py::capsule free_when_done(vf0_f1, [](void *f)
									 {
            double *vf0_f1 = reinterpret_cast<double *>(f);
            delete [] vf0_f1; });

	return py::array_t<double>(
		std::vector<size_t>{size, 2},
		std::vector<size_t>{2 * sizeof(double), sizeof(double)},
		vf0_f1,
		free_when_done);
}

py::array_t<double> compute_roemer_delay(py::array_t<double> &vT, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_T = vT.unchecked<1>();
	size_t size = array_T.size();

	double *vroemer_delay = new double[size];

	double sqrt_1_ee = std::sqrt(1. - e * e);

	for (size_t i = 0; i < size; i++)
	{
		double T = array_T(i);
		// calculate fractional orbit
		double orbit = (T - T0) / Pb;
		double norbit = std::floor(orbit);
		double forbit = orbit - norbit;

		// calculate eccentric anomaly
		double phase = 2. * M_PI * forbit;
		double E = phase + e * std::sin(phase) / std::sqrt(1. - 2. * e * std::cos(phase) + e * e);
		double dE = 1.;
		size_t nstep = 0;
		while (std::abs(dE) > 1e-14 && nstep++ < 8)
		{
			dE = (phase - (E - e * std::sin(E))) / (1. - e * std::cos(E));
			E += dE;
		}

		double sin_E = std::sin(E);
		double cos_E = std::cos(E);

		// calculate omega
		double omega = omega0 + omega_dot / (365.25 / Pb) * (norbit + forbit);

		// calculate v/c, a/c
		double sin_omega = std::sin(omega);
		double cos_omega = std::cos(omega);

		double temp = 2 * M_PI / (Pb * 86400. * (1. - e * cos_E));

		double roemer_delay = asini_c * ((cos_E - e) * sin_omega + sqrt_1_ee * sin_E * cos_omega);

		double doppler_factor1 = -asini_c * (-sin_omega * sin_E + sqrt_1_ee * cos_omega * cos_E) * temp;

		vroemer_delay[i] = roemer_delay * (1. + doppler_factor1);
	}

	py::capsule free_when_done(vroemer_delay, [](void *f)
									 {
            double *vroemer_delay = reinterpret_cast<double *>(f);
            delete [] vroemer_delay; });

	return py::array_t<double>(
		std::vector<size_t>{size},
		std::vector<size_t>{sizeof(double)},
		vroemer_delay,
		free_when_done);
}

double compute_roemer_delay2(double T, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	double sqrt_1_ee = std::sqrt(1. - e * e);

	// calculate fractional orbit
	double orbit = (T - T0) / Pb;
	double norbit = std::floor(orbit);
	double forbit = orbit - norbit;

	// calculate eccentric anomaly
	double phase = 2. * M_PI * forbit;
	double E = phase + e * std::sin(phase) / std::sqrt(1. - 2. * e * std::cos(phase) + e * e);
	double dE = 1.;
	size_t nstep = 0;
	while (std::abs(dE) > 1e-14 && nstep++ < 8)
	{
		dE = (phase - (E - e * std::sin(E))) / (1. - e * std::cos(E));
		E += dE;
	}

	double sin_E = std::sin(E);
	double cos_E = std::cos(E);

	// calculate omega
	double omega = omega0 + omega_dot / (365.25 / Pb) * (norbit + forbit);

	// calculate v/c, a/c
	double sin_omega = std::sin(omega);
	double cos_omega = std::cos(omega);

	double temp = 2 * M_PI / (Pb * 86400. * (1. - e * cos_E));

	double roemer_delay = asini_c * ((cos_E - e) * sin_omega + sqrt_1_ee * sin_E * cos_omega);

	double doppler_factor1 = -asini_c * (-sin_omega * sin_E + sqrt_1_ee * cos_omega * cos_E) * temp;

	return roemer_delay * (1. + doppler_factor1);
}

py::array_t<double> compute_Te(py::array_t<double> &vTb, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_Tb = vTb.unchecked<1>();
	size_t size = array_Tb.size();

	double *vTe = new double[size];

	for (size_t i = 0; i < size; i++)
	{
		double Tb = array_Tb[i];
		vTe[i] = Tb - compute_roemer_delay2(Tb, asini_c, Pb, T0, e, omega0, omega_dot) / 86400.;
	}

	py::capsule free_when_done(vTe, [](void *f)
									 {
            double *vTe = reinterpret_cast<double *>(f);
            delete [] vTe; });

	return py::array_t<double>(
		std::vector<size_t>{size},
		std::vector<size_t>{sizeof(double)},
		vTe,
		free_when_done);
}

py::array_t<double> timing(py::array_t<double> vTb, double f0, double f1, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_Tb = vTb.unchecked<1>();
	size_t size = array_Tb.size();

	double *vphase = new double[size];

	auto vTe = compute_Te(vTb, asini_c, Pb, T0, e, omega0, omega_dot);

	auto array_Te = vTe.unchecked<1>();

	double phsin = 0., phcos = 0.;
	for (size_t i = 0; i < size; i++)
	{
		double t = (array_Te(i) - T0) * 86400.;
		vphase[i] = f0 * t + 0.5 * f1 * t * t;

		phsin += std::sin(2. * M_PI * vphase[i]);
		phcos += std::cos(2. * M_PI * vphase[i]);
	}

	double phmean = std::atan2(phsin / size, phcos / size) / (2. * M_PI);

	for (size_t i = 0; i < size; i++)
	{
		vphase[i] -= phmean;
		vphase[i] -= std::round(vphase[i]);
	}

	py::capsule free_when_done(vphase, [](void *f)
									 {
            double *vphase = reinterpret_cast<double *>(f);
            delete [] vphase; });

	return py::array_t<double>(
		std::vector<size_t>{size},
		std::vector<size_t>{sizeof(double)},
		vphase,
		free_when_done);
}

#else

#include <stdlib.h>
#include <immintrin.h>

#include <math.h>
extern "C"
{
	__m256d _ZGVdN4v_cos(__m256d x);
	__m256d _ZGVdN4v_sin(__m256d x);
}

inline double haddd(__m256d a)
{
	__m256d b = _mm256_set_pd(0., 0., 0., 0.);
	a = _mm256_hadd_pd(a, b);

	return ((double *)&a)[0] + ((double *)&a)[2];
}

py::array_t<double> compute_f0_f1(py::array_t<double> &vT, double f0, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_T = vT.unchecked<1>();
	size_t size = array_T.size();

	double *vf0_f1 = new double[size * 2];

	__m256d avx_Pb = _mm256_set_pd(Pb, Pb, Pb, Pb);
	__m256d avx_T0 = _mm256_set_pd(T0, T0, T0, T0);
	double sqrt_1_ee = std::sqrt(1. - e * e);

	__m256d avx_threshold = _mm256_set_pd(1e-14, 1e-14, 1e-14, 1e-14);

	size_t aligned_size = (int)std::ceil(size / 4.) * 4;

	double *aligned_f0 = (double *)aligned_alloc(32, aligned_size * sizeof(double));
	double *aligned_f1 = (double *)aligned_alloc(32, aligned_size * sizeof(double));

	double *aligned_T = (double *)aligned_alloc(32, aligned_size * sizeof(double));
	for (size_t i = 0; i < size; i++)
	{
		aligned_T[i] = array_T(i);
	}

	for (size_t i = 0; i < size; i += 4)
	{
		__m256d avx_T = _mm256_load_pd(aligned_T + i);

		// calculate fractional orbit
		__m256d avx_orbit = _mm256_div_pd(_mm256_sub_pd(avx_T, avx_T0), avx_Pb);
		__m256d avx_norbit = _mm256_floor_pd(avx_orbit);
		__m256d avx_forbit = _mm256_sub_pd(avx_orbit, avx_norbit);

		// calculate eccentric anomaly
		__m256d avx_phase = 2. * M_PI * avx_forbit;
		__m256d avx_E = avx_phase + e * _ZGVdN4v_sin(avx_phase) / _mm256_sqrt_pd(1. - 2. * e * _ZGVdN4v_cos(avx_phase) + e * e);
		__m256d avx_dE = _mm256_set_pd(1., 1., 1., 1.);
		size_t nstep = 0;
		while (haddd(_mm256_or_pd(_mm256_cmp_pd(avx_dE, -avx_threshold, 1), _mm256_cmp_pd(avx_dE, avx_threshold, 14))) && nstep++ < 8)
		{
			avx_dE = (avx_phase - (avx_E - e * _ZGVdN4v_sin(avx_E))) / (1. - e * _ZGVdN4v_cos(avx_E));
			avx_E += avx_dE;
		}

		__m256d avx_sin_E = _ZGVdN4v_sin(avx_E);
		__m256d avx_cos_E = _ZGVdN4v_cos(avx_E);

		// calculate omega
		__m256d avx_omega = omega0 + omega_dot / (365.25 / Pb) * (avx_norbit + avx_forbit);

		// calculate v/c, a/c
		__m256d avx_sin_omega = _ZGVdN4v_sin(avx_omega);
		__m256d avx_cos_omega = _ZGVdN4v_cos(avx_omega);

		__m256d avx_temp = 2 * M_PI / (Pb * 86400. * (1. - e * avx_cos_E));

		__m256d avx_doppler_factor1 = -asini_c * (-avx_sin_omega * avx_sin_E + sqrt_1_ee * avx_cos_omega * avx_cos_E) * avx_temp;
		__m256d avx_doppler_factor2 = -asini_c * (-e * avx_sin_omega + avx_sin_omega * avx_cos_E + sqrt_1_ee * avx_cos_omega * avx_sin_E) / (e * avx_cos_E - 1.) * _mm256_mul_pd(avx_temp, avx_temp);

		__m256d avx_f0 = f0 * (1. + avx_doppler_factor1);
		__m256d avx_f1 = f0 * avx_doppler_factor2;

		_mm256_store_pd(aligned_f0 + i, avx_f0);
		_mm256_store_pd(aligned_f1 + i, avx_f1);
	}

	for (size_t i = 0; i < size; i++)
	{
		vf0_f1[i * 2 + 0] = aligned_f0[i];
		vf0_f1[i * 2 + 1] = aligned_f1[i];
	}

	free(aligned_T);
	free(aligned_f0);
	free(aligned_f1);

	py::capsule free_when_done(vf0_f1, [](void *f)
									 {
            double *vf0_f1 = reinterpret_cast<double *>(f);
            delete [] vf0_f1; });

	return py::array_t<double>(
		std::vector<size_t>{size, 2},
		std::vector<size_t>{2 * sizeof(double), sizeof(double)},
		vf0_f1,
		free_when_done);
}

py::array_t<double> compute_roemer_delay(py::array_t<double> &vT, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_T = vT.unchecked<1>();
	size_t size = array_T.size();

	double *vroemer_delay = new double[size];

	__m256d avx_Pb = _mm256_set_pd(Pb, Pb, Pb, Pb);
	__m256d avx_T0 = _mm256_set_pd(T0, T0, T0, T0);
	double sqrt_1_ee = std::sqrt(1. - e * e);

	__m256d avx_threshold = _mm256_set_pd(1e-14, 1e-14, 1e-14, 1e-14);

	size_t aligned_size = (int)std::ceil(size / 4.) * 4;

	double *aligned_roemer_delay = (double *)aligned_alloc(32, aligned_size * sizeof(double));

	double *aligned_T = (double *)aligned_alloc(32, aligned_size * sizeof(double));
	for (size_t i = 0; i < size; i++)
	{
		aligned_T[i] = array_T(i);
	}

	for (size_t i = 0; i < size; i += 4)
	{
		__m256d avx_T = _mm256_load_pd(aligned_T + i);

		// calculate fractional orbit
		__m256d avx_orbit = _mm256_div_pd(_mm256_sub_pd(avx_T, avx_T0), avx_Pb);
		__m256d avx_norbit = _mm256_floor_pd(avx_orbit);
		__m256d avx_forbit = _mm256_sub_pd(avx_orbit, avx_norbit);

		// calculate eccentric anomaly
		__m256d avx_phase = 2. * M_PI * avx_forbit;
		__m256d avx_E = avx_phase + e * _ZGVdN4v_sin(avx_phase) / _mm256_sqrt_pd(1. - 2. * e * _ZGVdN4v_cos(avx_phase) + e * e);
		__m256d avx_dE = _mm256_set_pd(1., 1., 1., 1.);
		size_t nstep = 0;
		while (haddd(_mm256_or_pd(_mm256_cmp_pd(avx_dE, -avx_threshold, 1), _mm256_cmp_pd(avx_dE, avx_threshold, 14))) && nstep++ < 8)
		{
			avx_dE = (avx_phase - (avx_E - e * _ZGVdN4v_sin(avx_E))) / (1. - e * _ZGVdN4v_cos(avx_E));
			avx_E += avx_dE;
		}

		__m256d avx_sin_E = _ZGVdN4v_sin(avx_E);
		__m256d avx_cos_E = _ZGVdN4v_cos(avx_E);

		// calculate omega
		__m256d avx_omega = omega0 + omega_dot / (365.25 / Pb) * (avx_norbit + avx_forbit);

		// calculate v/c, a/c
		__m256d avx_sin_omega = _ZGVdN4v_sin(avx_omega);
		__m256d avx_cos_omega = _ZGVdN4v_cos(avx_omega);

		__m256d avx_temp = 2 * M_PI / (Pb * 86400. * (1. - e * avx_cos_E));

		__m256d avx_roemer_delay = asini_c * ((avx_cos_E - e) * avx_sin_omega + sqrt_1_ee * avx_sin_E * avx_cos_omega);

		__m256d avx_doppler_factor1 = -asini_c * (-avx_sin_omega * avx_sin_E + sqrt_1_ee * avx_cos_omega * avx_cos_E) * avx_temp;

		avx_roemer_delay = avx_roemer_delay * (1. + avx_doppler_factor1);

		_mm256_store_pd(aligned_roemer_delay + i, avx_roemer_delay);
	}

	for (size_t i = 0; i < size; i++)
	{
		vroemer_delay[i] = aligned_roemer_delay[i];
	}

	free(aligned_T);
	free(aligned_roemer_delay);

	py::capsule free_when_done(vroemer_delay, [](void *f)
									 {
            double *vroemer_delay = reinterpret_cast<double *>(f);
            delete [] vroemer_delay; });

	return py::array_t<double>(
		std::vector<size_t>{size},
		std::vector<size_t>{sizeof(double)},
		vroemer_delay,
		free_when_done);
}

__m256d avx_compute_roemer_delay(__m256d avx_T, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	__m256d avx_Pb = _mm256_set_pd(Pb, Pb, Pb, Pb);
	__m256d avx_T0 = _mm256_set_pd(T0, T0, T0, T0);
	double sqrt_1_ee = std::sqrt(1. - e * e);

	__m256d avx_threshold = _mm256_set_pd(1e-14, 1e-14, 1e-14, 1e-14);

	// calculate fractional orbit
	__m256d avx_orbit = _mm256_div_pd(_mm256_sub_pd(avx_T, avx_T0), avx_Pb);
	__m256d avx_norbit = _mm256_floor_pd(avx_orbit);
	__m256d avx_forbit = _mm256_sub_pd(avx_orbit, avx_norbit);

	// calculate eccentric anomaly
	__m256d avx_phase = 2. * M_PI * avx_forbit;
	__m256d avx_E = avx_phase + e * _ZGVdN4v_sin(avx_phase) / _mm256_sqrt_pd(1. - 2. * e * _ZGVdN4v_cos(avx_phase) + e * e);
	__m256d avx_dE = _mm256_set_pd(1., 1., 1., 1.);
	size_t nstep = 0;
	while (haddd(_mm256_or_pd(_mm256_cmp_pd(avx_dE, -avx_threshold, 1), _mm256_cmp_pd(avx_dE, avx_threshold, 14))) && nstep++ < 8)
	{
		avx_dE = (avx_phase - (avx_E - e * _ZGVdN4v_sin(avx_E))) / (1. - e * _ZGVdN4v_cos(avx_E));
		avx_E += avx_dE;
	}

	__m256d avx_sin_E = _ZGVdN4v_sin(avx_E);
	__m256d avx_cos_E = _ZGVdN4v_cos(avx_E);

	// calculate omega
	__m256d avx_omega = omega0 + omega_dot / (365.25 / Pb) * (avx_norbit + avx_forbit);

	// calculate v/c, a/c
	__m256d avx_sin_omega = _ZGVdN4v_sin(avx_omega);
	__m256d avx_cos_omega = _ZGVdN4v_cos(avx_omega);

	__m256d avx_temp = 2 * M_PI / (Pb * 86400. * (1. - e * avx_cos_E));

	__m256d avx_roemer_delay = asini_c * ((avx_cos_E - e) * avx_sin_omega + sqrt_1_ee * avx_sin_E * avx_cos_omega);

	__m256d avx_doppler_factor1 = -asini_c * (-avx_sin_omega * avx_sin_E + sqrt_1_ee * avx_cos_omega * avx_cos_E) * avx_temp;

	return avx_roemer_delay * (1. + avx_doppler_factor1);
}

py::array_t<double> compute_Te(py::array_t<double> &vTb, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_Tb = vTb.unchecked<1>();
	size_t size = array_Tb.size();

	double *vTe = new double[size];

	size_t aligned_size = (int)std::ceil(size / 4.) * 4;

	double *aligned_Te = (double *)aligned_alloc(32, aligned_size * sizeof(double));

	double *aligned_Tb = (double *)aligned_alloc(32, aligned_size * sizeof(double));
	for (size_t i = 0; i < size; i++)
	{
		aligned_Tb[i] = array_Tb(i);
	}

	__m256d avx_threshold = _mm256_set_pd(1e-14, 1e-14, 1e-14, 1e-14);

	for (size_t i = 0; i < size; i += 4)
	{
		__m256d avx_Tb = _mm256_load_pd(aligned_Tb + i);

		__m256d avx_Te = avx_Tb - avx_compute_roemer_delay(avx_Tb, asini_c, Pb, T0, e, omega0, omega_dot) / 86400.;

		_mm256_store_pd(aligned_Te + i, avx_Te);
	}

	for (size_t i = 0; i < size; i++)
	{
		vTe[i] = aligned_Te[i];
	}

	free(aligned_Tb);
	free(aligned_Te);

	py::capsule free_when_done(vTe, [](void *f)
									 {
            double *vTe = reinterpret_cast<double *>(f);
            delete [] vTe; });

	return py::array_t<double>(
		std::vector<size_t>{size},
		std::vector<size_t>{sizeof(double)},
		vTe,
		free_when_done);
}

__m256d avx_compute_Te(__m256d avx_Tb, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	__m256d avx_Te = avx_Tb - avx_compute_roemer_delay(avx_Tb, asini_c, Pb, T0, e, omega0, omega_dot) / 86400.;

	return avx_Te;
}

py::array_t<double> timing(py::array_t<double> vTb, double f0, double f1, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot)
{
	auto array_Tb = vTb.cast<py::array_t<double>>().unchecked<1>();
	size_t size = array_Tb.size();

	double *vph = new double[size];

	size_t aligned_size = (int)std::ceil(size / 4.) * 4;

	double *aligned_ph = (double *)aligned_alloc(32, aligned_size * sizeof(double));

	double *aligned_Tb = (double *)aligned_alloc(32, aligned_size * sizeof(double));
	for (size_t i = 0; i < size; i++)
	{
		aligned_Tb[i] = array_Tb(i);
	}

	__m256d avx_2_PI = _mm256_set_pd(2. * M_PI, 2. * M_PI, 2. * M_PI, 2. * M_PI);

	__m256d avx_phsin = _mm256_set_pd(0., 0., 0., 0.);
	__m256d avx_phcos = _mm256_set_pd(0., 0., 0., 0.);
	for (size_t i = 0; i < (aligned_size - 4); i += 4)
	{
		__m256d avx_Tb = _mm256_load_pd(aligned_Tb + i);
		__m256d avx_Te = avx_compute_Te(avx_Tb, asini_c, Pb, T0, e, omega0, omega_dot);

		__m256d avx_t = (avx_Te - T0) * 86400.;
		__m256d avx_ph = f0 * avx_t + 0.5 * f1 * avx_t * avx_t;

		_mm256_store_pd(aligned_ph + i, avx_ph);

		avx_phsin += _ZGVdN4v_sin(avx_2_PI * avx_ph);
		avx_phcos += _ZGVdN4v_cos(avx_2_PI * avx_ph);
	}

	__m256d avx_Tb = _mm256_load_pd(aligned_Tb + (aligned_size - 4));
	__m256d avx_Te = avx_compute_Te(avx_Tb, asini_c, Pb, T0, e, omega0, omega_dot);
	__m256d avx_t = (avx_Te - T0) * 86400.;
	__m256d avx_ph = f0 * avx_t + 0.5 * f1 * avx_t * avx_t;
	_mm256_store_pd(aligned_ph + (aligned_size - 4), avx_ph);

	__m256d avx_phsin_tmp = _ZGVdN4v_sin(avx_2_PI * avx_ph);
	__m256d avx_phcos_tmp = _ZGVdN4v_cos(avx_2_PI * avx_ph);

	double cutoff = size + 4 - aligned_size;
	__m256d avx_cutoff = _mm256_set_pd(cutoff, cutoff, cutoff, cutoff);
	__m256d avx_range = _mm256_set_pd(3., 2., 1., 0.);
	__m256d avx_mask = _mm256_cmp_pd(avx_range, avx_cutoff, 1);

	avx_phsin_tmp = _mm256_and_pd(avx_phsin_tmp, avx_mask);
	avx_phcos_tmp = _mm256_and_pd(avx_phcos_tmp, avx_mask);

	avx_phsin += avx_phsin_tmp;
	avx_phcos += avx_phcos_tmp;

	double phsin = haddd(avx_phsin);
	double phcos = haddd(avx_phcos);

	double phmean = std::atan2(phsin / size, phcos / size) / (2. * M_PI);

	for (size_t i = 0; i < size; i += 4)
	{
		avx_ph = _mm256_load_pd(aligned_ph + i);
		avx_ph -= phmean;
		avx_ph -= _mm256_round_pd(avx_ph, (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));

		_mm256_store_pd(aligned_ph + i, avx_ph);
	}

	for (size_t i = 0; i < size; i++)
	{
		vph[i] = aligned_ph[i];
	}

	free(aligned_Tb);
	free(aligned_ph);

	py::capsule free_when_done(vph, [](void *f)
									 {
            double *vph = reinterpret_cast<double *>(f);
            delete [] vph; });

	return py::array_t<double>(
		std::vector<size_t>{size},
		std::vector<size_t>{sizeof(double)},
		vph,
		free_when_done);
}

#endif

double get_chi2(py::array_t<double> &data)
{
	py::buffer_info buf_data = data.request();

	size_t size = buf_data.shape[1];

	double *array_data = static_cast<double *>(buf_data.ptr);

	double *array_f = array_data + 0 * size;
	double *array_f_s = array_data + 1 * size;
	double *array_ferr = array_data + 2 * size;

	double chi2 = 0.;

	for (size_t i = 0; i < size; i++)
	{
		double tmp = (array_f[i] - array_f_s[i]) / array_ferr[i];
		chi2 += tmp * tmp;
	}

	return chi2;
}

py::array_t<double> correct_ra_dec(py::array_t<double> &vobs_ssb_arr, double a, double b, double da, double db)
{
	auto obs_ssb_arr = vobs_ssb_arr.unchecked<2>();

	if (obs_ssb_arr.shape(1) != 3)
		throw std::runtime_error("Input shape[-1] must be 3");

	size_t N = obs_ssb_arr.shape(0);

	std::vector<double> posPulsar = {std::cos(b) * std::cos(a), std::cos(b) * std::sin(a), std::sin(b)};
	std::vector<double> posPulsar_new = {std::cos(b + db) * std::cos(a + da), std::cos(b + db) * std::sin(a + da), std::sin(b + db)};

	std::vector<double> dpos = {posPulsar_new[0] - posPulsar[0], posPulsar_new[1] - posPulsar[1], posPulsar_new[2] - posPulsar[2]};

	std::vector<double> tcorr(N, 0.);
	for (size_t j = 0; j < N; j++)
	{
		tcorr[j] = obs_ssb_arr(j, 0) * dpos[0] + obs_ssb_arr(j, 1) * dpos[1] + obs_ssb_arr(j, 2) * dpos[2];
	}

	return py::array(py::cast(tcorr));
}

py::array_t<double> timing2(py::array_t<double> vTb, double da, double db, double f0, double f1, double asini_c, double Pb, double T0, double e, double omega0, double omega_dot, py::array_t<double> &vobs_ssb_arr, double a, double b)
{
	auto tcorr = correct_ra_dec(vobs_ssb_arr, a, b, da, db);

	std::vector<double> vTb_corr(vTb.size(), 0.);
	for (size_t i=0; i < vTb.size(); i++)
	{
		vTb_corr[i] = vTb.data()[i] + tcorr.data()[i] * 1.1574074074074073e-05;
	}

	return timing(py::array(py::cast(vTb_corr)), f0, f1, asini_c, Pb, T0, e, omega0, omega_dot);
}

py::array_t<double> lmfit2(py::array_t<double> &x, py::array_t<double> &y, py::array_t<double> &p0, py::array_t<double> &sigma, py::array_t<double> &diff_step, py::array_t<bool> &if_fit, py::array_t<double> &vobs_ssb_arr, double a, double b)
{
	std::function<py::array_t<double> (py::array_t<double>, py::array_t<double>)> func = [&vobs_ssb_arr, a, b](py::array_t<double> x, py::array_t<double> p0){
		double f0 = p0.data()[0];
		double f1 = p0.data()[1];
		double da = p0.data()[2];
		double db = p0.data()[3];
		double asini_c = p0.data()[4];
		double Pb = p0.data()[5];
		double T0 = p0.data()[6];
		double e = p0.data()[7];
		double omega0 = p0.data()[8];
		double omega_dot = p0.data()[9];

		return timing2(x, da, db, f0, f1, asini_c, Pb, T0, e, omega0, omega_dot, vobs_ssb_arr, a, b);
	};

	return lmfit(func, x, y, p0, sigma, diff_step, if_fit);
}

PYBIND11_MODULE(orbit_utils, m)
{
	m.def(
		"compute_f0_f1", &compute_f0_f1,
		"calculate f0,f1 with the orbital parameters");

	m.def(
		"compute_roemer_delay", &compute_roemer_delay,
		"calculate roemer_delay with the orbital parameters");

	m.def(
		"compute_Te", &compute_Te,
		"calculate Te with the orbital parameters");

	m.def(
		"timing", &timing,
		"perform timing");

	m.def(
		"timing2", &timing2,
		"perform timing with ra,dec");

	m.def(
		"get_chi2", &get_chi2,
		"calculate likelihood");

	m.def(
		"lmfit2", &lmfit2,
		"perform leastsq fitting");
}