/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-15 13:52:50
 * @modify date 2022-08-15 13:52:50
 * @desc [description]
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include "databuffer.h"
#include "downsample.h"
#include "preprocesslite.h"
#include "equalize.h"
#include "rfi.h"

unsigned int num_threads = 1;

PYBIND11_MODULE(filtools, m)
{
	pybind11::class_<DataBuffer<float>>(m, "DataBuffer")
		.def(py::init<>())
		.def(py::init<long int, int>())
		.def("get_data",
			[](const DataBuffer<float> &m) {
				py::array_t<float> ret = py::cast(m.buffer);
				return ret.reshape({(size_t)m.nsamples, (size_t)m.nchans});
			}
		)
		.def("update_data",
			[](DataBuffer<float> &m, py::array_t<float> &data) {
				auto array_data = data.unchecked<2>();

				size_t nrow = array_data.shape(0);
				size_t ncol = array_data.shape(1);

				if (nrow != m.nsamples || ncol != m.nchans)
					throw std::runtime_error("Input shape must match");

				for (size_t j=0; j<nrow; j++)
					for (size_t i=0; i<ncol; i++)
						m.buffer[j * ncol + i] = array_data(j, i);
			}
		)
		.def_readonly("nsamples", &DataBuffer<float>::nsamples)
		.def_readonly("nchans", &DataBuffer<float>::nchans)
		.def_readwrite("tsamp", &DataBuffer<float>::tsamp)
		.def_readwrite("frequencies", &DataBuffer<float>::frequencies);

	pybind11::class_<Downsample, DataBuffer<float>>(m, "Downsample")
		.def(py::init<>())
		.def(py::init<int, int>())
		.def("prepare", &Downsample::prepare)
		.def("run", &Downsample::run)
		.def("get_data",
			[](const Downsample &m) {
				py::array_t<float> ret = py::cast(m.buffer);
				return ret.reshape({(size_t)m.nsamples, (size_t)m.nchans});
			}
		)
		.def("update_data",
			[](Downsample &m, py::array_t<float> &data) {
				auto array_data = data.unchecked<2>();

				size_t nrow = array_data.shape(0);
				size_t ncol = array_data.shape(1);

				if (nrow != m.nsamples || ncol != m.nchans)
					throw std::runtime_error("Input shape must match");

				for (size_t j=0; j<nrow; j++)
					for (size_t i=0; i<ncol; i++)
						m.buffer[j * ncol + i] = array_data(j, i);
			}
		)
		.def_readwrite("td", &Downsample::td)
		.def_readwrite("fd", &Downsample::fd)
		.def_readonly("nsamples", &Downsample::nsamples)
		.def_readonly("nchans", &Downsample::nchans)
		.def_readonly("tsamp", &Downsample::tsamp)
		.def_readonly("frequencies", &Downsample::frequencies);

	pybind11::class_<PreprocessLite, DataBuffer<float>>(m, "PreprocessLite")
		.def(py::init<>())
		.def("prepare", &PreprocessLite::prepare)
		.def("run", &PreprocessLite::run)
		.def("get_data",
			[](const PreprocessLite &m) {
				py::array_t<float> ret = py::cast(m.buffer);
				return ret.reshape({(size_t)m.nsamples, (size_t)m.nchans});
			}
		)
		.def("update_data",
			[](PreprocessLite &m, py::array_t<float> &data) {
				auto array_data = data.unchecked<2>();

				size_t nrow = array_data.shape(0);
				size_t ncol = array_data.shape(1);

				if (nrow != m.nsamples || ncol != m.nchans)
					throw std::runtime_error("Input shape must match");

				for (size_t j=0; j<nrow; j++)
					for (size_t i=0; i<ncol; i++)
						m.buffer[j * ncol + i] = array_data(j, i);
			}
		)
		.def_readwrite("threshold", &PreprocessLite::thresig)
		.def_readwrite("filltype", &PreprocessLite::filltype)
		.def_readonly("killrate", &PreprocessLite::killrate)
		.def_readwrite("td", &PreprocessLite::td)
		.def_readwrite("fd", &PreprocessLite::fd)
		.def_readonly("nsamples", &PreprocessLite::nsamples)
		.def_readonly("nchans", &PreprocessLite::nchans)
		.def_readonly("tsamp", &PreprocessLite::tsamp)
		.def_readonly("frequencies", &PreprocessLite::frequencies);

	pybind11::class_<Equalize, DataBuffer<float>>(m, "Equalize")
		.def(py::init<>())
		.def("prepare", &Equalize::prepare)
		.def("filter", &Equalize::filter)
		.def("run", &Equalize::run)
		.def("get_data",
			[](const Equalize &m) {
				py::array_t<float> ret = py::cast(m.buffer);
				return ret.reshape({(size_t)m.nsamples, (size_t)m.nchans});
			}
		)
		.def("update_data",
			[](Equalize &m, py::array_t<float> &data) {
				auto array_data = data.unchecked<2>();

				size_t nrow = array_data.shape(0);
				size_t ncol = array_data.shape(1);

				if (nrow != m.nsamples || ncol != m.nchans)
					throw std::runtime_error("Input shape must match");

				for (size_t j=0; j<nrow; j++)
					for (size_t i=0; i<ncol; i++)
						m.buffer[j * ncol + i] = array_data(j, i);
			}
		)
		.def_readonly("nsamples", &Equalize::nsamples)
		.def_readonly("nchans", &Equalize::nchans)
		.def_readonly("tsamp", &Equalize::tsamp)
		.def_readonly("frequencies", &Equalize::frequencies);

	pybind11::class_<RFI, DataBuffer<float>>(m, "RFI")
		.def(py::init<>())
		.def("prepare", &RFI::prepare)
		.def("zap", &RFI::zap)
		.def("zdot", &RFI::zdot)
		.def("zero", &RFI::zero)
		.def("kadaneF", &RFI::kadaneF)
		.def("get_data",
			[](const RFI &m) {
				py::array_t<float> ret = py::cast(m.buffer);
				return ret.reshape({(size_t)m.nsamples, (size_t)m.nchans});
			}
		)
		.def("update_data",
			[](RFI &m, py::array_t<float> &data) {
				auto array_data = data.unchecked<2>();

				size_t nrow = array_data.shape(0);
				size_t ncol = array_data.shape(1);

				if (nrow != m.nsamples || ncol != m.nchans)
					throw std::runtime_error("Input shape must match");

				for (size_t j=0; j<nrow; j++)
					for (size_t i=0; i<ncol; i++)
						m.buffer[j * ncol + i] = array_data(j, i);
			}
		)
		.def_readwrite("filltype", &RFI::filltype)
		.def_readonly("nsamples", &RFI::nsamples)
		.def_readonly("nchans", &RFI::nchans)
		.def_readonly("tsamp", &RFI::tsamp)
		.def_readonly("frequencies", &RFI::frequencies);
}