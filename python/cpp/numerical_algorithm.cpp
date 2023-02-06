/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-12-14 11:35:37
 * @modify date 2022-12-14 11:35:37
 * @desc [description]
 */

#include <iostream>

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

py::array_t<double> lmfit(std::function<py::array_t<double> (py::array_t<double>, py::array_t<double>)> &func, py::array_t<double> &x, py::array_t<double> &y, py::array_t<double> &p0, py::array_t<double> &sigma, py::array_t<double> &diff_step, py::array_t<bool> &if_fit)
{
	size_t N_full = p0.size();
	size_t M = x.size();

	if (y.size() != M or sigma.size() != M or diff_step.size() != N_full or if_fit.size() != N_full)
	{
		throw std::runtime_error("Input shapes must match");
	}

	std::vector<double> p(N_full);

	for (size_t i = 0; i < N_full; i++)
	{
		p[i] = p0.data()[i];
	}

	size_t N = 0;
	for (size_t i = 0; i < N_full; i++)
	{
		if (if_fit.data()[i])
			N++;
	}

	double *A = new double[N * M];

	size_t k = 0;
	for (size_t i = 0; i < N_full; i++)
	{
		if (if_fit.data()[i])
		{
			p[i] = p0.data()[i] - diff_step.data()[i];
			auto f1 = func(x, pybind11::array(pybind11::cast(p)));

			p[i] = p0.data()[i] + diff_step.data()[i];
			auto f2 = func(x, pybind11::array(pybind11::cast(p)));

			for (size_t j = 0; j < M; j++)
			{
				A[k * M + j] = (f2.data()[j] - f1.data()[j]) / (2. * diff_step.data()[i]) / sigma.data()[j];
			}

			p[i] = p0.data()[i];

			k++;
		}
	}

	double *b = new double[M];

	auto f = func(x, p0);

	for (size_t j = 0; j < M; j++)
	{
		b[j] = (y.data()[j] - f.data()[j]) / sigma.data()[j];
	}

	int m = M, n = N, nrhs = 1, lda = M, ldb = M, info, lwork;

	double wkopt;

	double *work;

	lwork = -1;

	dgels_("No transpose", &m, &n, &nrhs, A, &lda, b, &ldb, &wkopt, &lwork, &info);

	lwork = (int)wkopt;

	work = new double[lwork];

	dgels_("No transpose", &m, &n, &nrhs, A, &lda, b, &ldb, work, &lwork, &info);

	if (info > 0)
	{
		std::cout << "The diagonal element %i of the triangular factor " << info << " of A is zero, so that A does not have full rank;" << std::endl;
		delete[] work;
		delete[] A;
		delete[] b;
		throw std::runtime_error("the least squares solution could not be computed.\n");
	}

	delete[] work;
	delete[] A;

	k = 0; 
	for (size_t i = 0; i < N_full; i++)
	{
		if (if_fit.data()[i])
		{
			p[i] = p0.data()[i] + b[k];
			k++;
		}
	}

	delete[] b;

	return pybind11::array(pybind11::cast(p));
}

PYBIND11_MODULE(numerical_algorithm, m)
{
	m.def(
		"lmfit", &lmfit,
		"linear least squares fit\nlmfit(func, x, y, p0, sigma, diff_step, if_fit)");
}