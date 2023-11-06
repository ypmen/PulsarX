/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-10-04 20:44:06
 * @modify date 2022-10-04 20:44:06
 * @desc [description]
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

#include "psrdatareader.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "filterbankwriter.h"

unsigned int num_threads = 1;

PYBIND11_MODULE(psrdata_reader, m)
{
	py::class_<MJD>(m, "MJD")
		.def(py::init<>())
		.def(py::init<long int, long int, double>())
		.def_readwrite("stt_imjd", &MJD::stt_imjd)
		.def_readwrite("stt_smjd", &MJD::stt_smjd)
		.def_readwrite("stt_offs", &MJD::stt_offs);

	py::class_<PSRDataReader>(m, "PSRDataReader")
		.def("get_fmin_fmax", [](PSRDataReader &m){
			double fmin = 0., fmax = 0.;
			m.get_fmin_fmax(fmin, fmax);
			py::tuple tup = py::make_tuple(fmin, fmax);
			return tup;
		})
		.def("get_fc_bw", [](PSRDataReader &m){
			double fc = 0., bw = 0.;
			m.get_fc_bw(fc, bw);
			py::tuple tup = py::make_tuple(fc, bw);
			return tup;
		})
		.def("get_fch1_foff", [](PSRDataReader &m){
			double fch1 = 0., foff = 0.;
			m.get_fch1_foff(fch1, foff);
			py::tuple tup = py::make_tuple(fch1, foff);
			return tup;
		})
		.def("get_chbw", &PSRDataReader::get_chbw)
		.def_readwrite("fnames", &PSRDataReader::fnames)
		.def_readwrite("sumif", &PSRDataReader::sumif)
		.def_readwrite("contiguous", &PSRDataReader::contiguous)
		.def_readwrite("verbose", &PSRDataReader::verbose)
		.def_readwrite("skip_start", &PSRDataReader::skip_start)
		.def_readwrite("skip_end", &PSRDataReader::skip_end)
		.def_readwrite("apply_zero_off", &PSRDataReader::apply_zero_off)
		.def_readwrite("apply_scloffs", &PSRDataReader::apply_scloffs)
		.def_readwrite("apply_wts", &PSRDataReader::apply_wts)
		.def_readwrite("telescope", &PSRDataReader::telescope)
		.def_readwrite("source_name", &PSRDataReader::source_name)
		.def_readwrite("ra", &PSRDataReader::ra)
		.def_readwrite("dec", &PSRDataReader::dec)
		.def_readwrite("beam", &PSRDataReader::beam)
		.def_readonly("start_mjd", &PSRDataReader::start_mjd)
		.def_readonly("mjd_starts", &PSRDataReader::mjd_starts)
		.def_readonly("mjd_ends", &PSRDataReader::mjd_ends)
		.def_readonly("idmap", &PSRDataReader::idmap)
		.def_readonly("nsamples", &PSRDataReader::nsamples)
		.def_readonly("nifs", &PSRDataReader::nifs)
		.def_readonly("nsblk", &PSRDataReader::nsblk)
		.def_readonly("nchans", &PSRDataReader::nchans)
		.def_readonly("tsamp", &PSRDataReader::tsamp)
		.def_readonly("frequencies", &PSRDataReader::frequencies)
		.def_readonly("is_end", &PSRDataReader::is_end);

	py::class_<FilterbankReader, PSRDataReader>(m, "FilterbankReader")
		.def(py::init<>())
		.def("check", &FilterbankReader::check)
		.def("read_header", &FilterbankReader::read_header)
		.def("read_data", static_cast<size_t (FilterbankReader::*)(DataBuffer<float> &, size_t, bool)>(&FilterbankReader::read_data));

	py::class_<PsrfitsReader, PSRDataReader>(m, "PsrfitsReader")
		.def("check", &PsrfitsReader::check)
		.def("read_header", &PsrfitsReader::read_header)
		.def("read_data", static_cast<size_t (PsrfitsReader::*)(DataBuffer<float> &, size_t, bool)>(&PsrfitsReader::read_data));

	py::class_<FilterbankWriter>(m, "FilterbankWriter")
		.def(py::init<>())
		.def("prepare", &FilterbankWriter::prepare)
		.def("run", &FilterbankWriter::run)
		.def_readwrite("outstd", &FilterbankWriter::outstd)
		.def_readwrite("outmean", &FilterbankWriter::outmean)
		.def_readwrite("fil", &FilterbankWriter::fil);
}