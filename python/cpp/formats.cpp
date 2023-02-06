/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-17 21:09:59
 * @modify date 2022-08-17 21:09:59
 * @desc [description]
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>

namespace py = pybind11;

#include <iostream>
#include "psrfits.h"
#include "filterbank.h"

PYBIND11_MODULE(formats, m)
{
	py::class_<Integration> integration(m, "Integration");
	integration.def(py::init<>())
		.def("get_data", [](Integration &m) -> py::object {
			if (m.mode == Integration::Mode::SEARCH)
			{
				switch (m.dtype)
				{
				case Integration::DataType::UINT8 :
					{
						auto data = py::array_t<unsigned char> (m.nsblk * m.npol * m.nchan);
						
						py::buffer_info buf = data.request();

						unsigned char *ptr = static_cast<unsigned char *>(buf.ptr);

						std::memcpy(ptr, m.data, sizeof(unsigned char) * m.nsblk * m.npol * m.nchan);

						return data.reshape({(size_t)m.nsblk, (size_t)m.npol, (size_t)m.nchan});
					};
					break;
				
				default:
					{
						std::cerr<<"Data type is not supported"<<std::endl;
						return py::cast<py::none>(Py_None);
					};
					break;
				}
			}
			else
			{
				switch (m.dtype)
				{
				default:
					{
						std::cerr<<"Data type is not supported"<<std::endl;
						return py::cast<py::none>(Py_None);
					};
					break;
				}
			}
		}, py::return_value_policy::move)
		.def("get_frequencies", [](Integration &m){
			auto frequencies = py::array_t<double> (m.nchan);
			py::buffer_info buf = frequencies.request();
			double *ptr = static_cast<double *>(buf.ptr);
			std::memcpy(ptr, m.frequencies, sizeof(double) * m.nchan);
			return frequencies;
		}, py::return_value_policy::move)
		.def("get_weights", [](Integration &m){
			auto weights = py::array_t<float> (m.nchan);
			py::buffer_info buf = weights.request();
			float *ptr = static_cast<float *>(buf.ptr);
			std::memcpy(ptr, m.weights, sizeof(float) * m.nchan);
			return weights;
		}, py::return_value_policy::move)
		.def("get_offsets", [](Integration &m){
			auto offsets = py::array_t<float> (m.npol * m.nchan);
			py::buffer_info buf = offsets.request();
			float *ptr = static_cast<float *>(buf.ptr);
			std::memcpy(ptr, m.offsets, sizeof(float) * m.npol * m.nchan);
			return offsets.reshape({(size_t)m.npol, (size_t)m.nchan});
		}, py::return_value_policy::move)
		.def("get_scales", [](Integration &m){
			auto scales = py::array_t<float> (m.npol * m.nchan);
			py::buffer_info buf = scales.request();
			float *ptr = static_cast<float *>(buf.ptr);
			std::memcpy(ptr, m.scales, sizeof(float) * m.npol * m.nchan);
			return scales.reshape({(size_t)m.npol, (size_t)m.nchan});
		}, py::return_value_policy::move)
		.def_readonly("mode", &Integration::mode)
		.def_readonly("dtype", &Integration::dtype)
		.def_readwrite("indexval", &Integration::indexval)
		.def_readwrite("folding_period", &Integration::folding_period)
		.def_readwrite("tsubint", &Integration::tsubint)
		.def_readwrite("offs_sub", &Integration::offs_sub)
		.def_readwrite("lst_sub", &Integration::lst_sub)
		.def_readwrite("ra_sub", &Integration::ra_sub)
		.def_readwrite("dec_sub", &Integration::dec_sub)
		.def_readwrite("glon_sub", &Integration::glon_sub)
		.def_readwrite("glat_sub", &Integration::glat_sub)
		.def_readwrite("fd_ang", &Integration::fd_ang)
		.def_readwrite("pos_ang", &Integration::pos_ang)
		.def_readwrite("par_ang", &Integration::par_ang)
		.def_readwrite("tel_az", &Integration::tel_az)
		.def_readwrite("tel_zen", &Integration::tel_zen)
		.def_readwrite("aux_dm", &Integration::aux_dm)
		.def_readwrite("aux_rm", &Integration::aux_rm)
		.def_readonly("npol", &Integration::npol)
		.def_readonly("nchan", &Integration::nchan)
		.def_readonly("nbin", &Integration::nbin)
		.def_readonly("nsblk", &Integration::nsblk)
		.def_readonly("nbits", &Integration::nbits);

	py::enum_<Integration::Mode>(integration, "Mode")
		.value("SEARCH", Integration::Mode::SEARCH)
		.value("FOLD", Integration::Mode::FOLD)
		.export_values();

	py::enum_<Integration::DataType>(integration, "DataType")
		.value("USHORT", Integration::DataType::USHORT)
		.value("SHORT", Integration::DataType::SHORT)
		.value("UINT1", Integration::DataType::UINT1)
		.value("UINT2", Integration::DataType::UINT2)
		.value("UINT4", Integration::DataType::UINT4)
		.value("UINT8", Integration::DataType::UINT8)
		.value("FLOAT", Integration::DataType::FLOAT)
		.export_values();

	// py::class_<MJD>(m, "MJD")
	// 	.def(py::init<>())
	// 	.def(py::init<long int, long int, double>())
	// 	.def(py::init<long double>())
	// 	.def_readonly("stt_imjd", &MJD::stt_imjd)
	// 	.def_readonly("stt_smjd", &MJD::stt_smjd)
	// 	.def_readonly("stt_offs", &MJD::stt_offs);

	py::class_<PrimaryHDU>(m, "PrimaryHDU")
		.def(py::init<>())
		.def("load", &PrimaryHDU::load)
		.def("unload", &PrimaryHDU::unload)
		.def_readwrite("start_mjd", &PrimaryHDU::start_mjd)
		.def("set_observer", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.observer, s.c_str());
		})
		.def("set_projid", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.projid, s.c_str());
		})
		.def("set_telesop", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.telesop, s.c_str());
		})
		.def("set_ibeam", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.ibeam, s.c_str());
		})
		.def("set_obs_mode", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.obs_mode, s.c_str());
		})
		.def_readwrite("chan_dm", &PrimaryHDU::chan_dm)
		.def("set_date_obs", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.date_obs, s.c_str());
		})
		.def("set_src_name", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.src_name, s.c_str());
		})
		.def("set_ra", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.ra, s.c_str());
		})
		.def("set_dec", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.dec, s.c_str());
		})
		.def("set_stt_crd1", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.stt_crd1, s.c_str());
		})
		.def("set_stt_crd2", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.stt_crd2, s.c_str());
		})
		.def("set_trk_mode", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.trk_mode, s.c_str());
		})
		.def("set_stp_crd1", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.stp_crd1, s.c_str());
		})
		.def("set_stp_crd2", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.stp_crd2, s.c_str());
		})
		.def("set_fd_mode", [](PrimaryHDU &m, const std::string &s){
			std::strcpy(m.fd_mode, s.c_str());
		})
		.def("get_observer", [](PrimaryHDU &m){
			return std::string(m.observer);
		})
		.def("get_projid", [](PrimaryHDU &m){
			return std::string(m.projid);
		})
		.def("get_telesop", [](PrimaryHDU &m){
			return std::string(m.telesop);
		})
		.def("get_ibeam", [](PrimaryHDU &m){
			return std::string(m.ibeam);
		})
		.def("get_obs_mode", [](PrimaryHDU &m){
			return std::string(m.obs_mode);
		})
		.def("get_date_obs", [](PrimaryHDU &m){
			return std::string(m.date_obs);
		})
		.def("get_src_name", [](PrimaryHDU &m){
			return std::string(m.src_name);
		})
		.def("get_ra", [](PrimaryHDU &m){
			return std::string(m.ra);
		})
		.def("get_dec", [](PrimaryHDU &m){
			return std::string(m.dec);
		})
		.def("get_stt_crd1", [](PrimaryHDU &m){
			return std::string(m.stt_crd1);
		})
		.def("get_stt_crd2", [](PrimaryHDU &m){
			return std::string(m.stt_crd2);
		})
		.def("get_trk_mode", [](PrimaryHDU &m){
			return std::string(m.trk_mode);
		})
		.def("get_stp_crd1", [](PrimaryHDU &m){
			return std::string(m.stp_crd1);
		})
		.def("get_stp_crd2", [](PrimaryHDU &m){
			return std::string(m.stp_crd2);
		})
		.def("get_fd_mode", [](PrimaryHDU &m){
			return std::string(m.fd_mode);
		});
	
	py::class_<PsrparamHDU>(m, "PsrparamHDU")
		.def(py::init<>())
		.def("load", static_cast<bool (PsrparamHDU::*)(fitsfile *)>(&PsrparamHDU::load))
		.def("unload", &PsrparamHDU::unload)
		.def_readwrite("text", &PsrparamHDU::text);

	py::class_<T2predictHDU>(m, "T2predictHDU")
		.def(py::init<>())
		.def("load", static_cast<bool (T2predictHDU::*)(fitsfile *)>(&T2predictHDU::load))
		.def("unload", &T2predictHDU::unload)
		.def_readwrite("text", &T2predictHDU::text);

	py::class_<SubintHDU>(m, "SubintHDU")
		.def(py::init<>())
		.def("resize", &SubintHDU::resize)
		.def("load_fold_data", [](SubintHDU &m, py::array_t<float> &data){
			py::buffer_info buf = data.request();

			if (buf.format != py::format_descriptor<float>::format())
				throw std::runtime_error("Data type should be float32");

			if (buf.ndim != 4)
        		throw std::runtime_error("Number of dimensions must be 4, (nsubint, npol, nchan, nbin)");

			int ns = buf.shape[0];
			int np = buf.shape[1];
			int nc = buf.shape[2];
			int nb = buf.shape[3];

			m.load_data(static_cast<float *> (buf.ptr), ns, np, nc, nb);
		})
		.def("load_search_data", [](SubintHDU &m, py::array_t<unsigned char> &data){
			py::buffer_info buf = data.request();

			if (buf.format != py::format_descriptor<unsigned char>::format())
				throw std::runtime_error("Data type should be uint8");

			if (buf.ndim != 4)
        		throw std::runtime_error("Number of dimensions must be 4, (nsubint, npol, nchan, nbin)");

			switch (m.dtype)
			{
				case Integration::UINT1:
				{
					if (m.nbits != 1)
						throw std::runtime_error("Data type is not set");
				};break;
				case Integration::UINT2:
				{
					if (m.nbits != 2)
						throw std::runtime_error("Data type is not set");
				};break;
				case Integration::UINT4:
				{
					if (m.nbits != 4)
						throw std::runtime_error("Data type is not set");
				};break;
				case Integration::UINT8:
				{
					if (m.nbits != 8)
						throw std::runtime_error("Data type is not set");
				};break;
				case Integration::FLOAT:
				{
					if (m.nbits != 32)
						throw std::runtime_error("Data type is not set");
				};break;
				default:
				{
					std::cerr<<"Data type is supported"<<std::endl;
					throw std::runtime_error("Data type is supported");
				};
			}

			int ns = buf.shape[0];
			int nb = buf.shape[1] * (8 / m.nbits);
			int np = buf.shape[2];
			int nc = buf.shape[3];

			m.load_data(static_cast<float *> (buf.ptr), ns, np, nc, nb);
		})
		.def("load_frequencies", [](SubintHDU &m, py::array_t<double> &frequencies){
			py::buffer_info buf = frequencies.request();
			if (buf.ndim != 1)
        		throw std::runtime_error("Number of dimensions must be 1");
			if (buf.format != py::format_descriptor<double>::format())
				throw std::runtime_error("Data type should be float64");
			if (buf.size != m.nchan)
				throw std::runtime_error("frequencies size is not matched");

			m.load_frequencies(static_cast<double *>(buf.ptr), buf.size);
		})
		.def("load_weights", [](SubintHDU &m, py::array_t<float> &weights){
			py::buffer_info buf = weights.request();
			if (buf.ndim != 1)
        		throw std::runtime_error("Number of dimensions must be 1");
			if (buf.format != py::format_descriptor<float>::format())
				throw std::runtime_error("Data type should be float32");
			if (buf.size != m.nchan)
				throw std::runtime_error("frequencies size is not matched");

			m.load_weights(static_cast<float *>(buf.ptr), buf.size);
		})
		.def("update_period", [](SubintHDU &m, double period){
			for (size_t k=0; k<m.nsubint; k++)
				m.integrations[k].folding_period = period;
		})
		.def("update_dm", [](SubintHDU &m, double dm){
			for (size_t k=0; k<m.nsubint; k++)
				m.integrations[k].aux_dm = dm;
		})
		.def("update_rm", [](SubintHDU &m, double rm){
			for (size_t k=0; k<m.nsubint; k++)
				m.integrations[k].aux_rm = rm;
		})
		.def("unload", &SubintHDU::unload)
		.def("unload_header", &SubintHDU::unload_header)
		.def("unload_data", &SubintHDU::unload_data)
		.def("load_header", &SubintHDU::load_header)
		.def("load_integration", static_cast<bool (SubintHDU::*)(fitsfile *, int, Integration &)>(&SubintHDU::load_integration))
		.def("get_integration", [](SubintHDU &m, int k) {
			return m.integrations[k];
		})
		.def_readonly("nbits", &SubintHDU::nbits)
		.def_readonly("nsubint", &SubintHDU::nsubint)
		.def_readonly("npol", &SubintHDU::npol)
		.def_readonly("nchan", &SubintHDU::nchan)
		.def_readonly("nbin", &SubintHDU::nbin)
		.def_readonly("nsblk", &SubintHDU::nsblk)
		.def_readonly("nsuboffs", &SubintHDU::nsuboffs)
		.def_readonly("nstot", &SubintHDU::nstot)
		.def_readonly("nsamples", &SubintHDU::nsamples)
		.def_readwrite("tbin", &SubintHDU::tbin)
		.def_readwrite("dm", &SubintHDU::dm)
		.def_readwrite("rm", &SubintHDU::rm)
		.def_readwrite("mode", &SubintHDU::mode)
		.def_readwrite("dtype", &SubintHDU::dtype);

	py::class_<fitsfile, std::shared_ptr<fitsfile>>(m, "Fitsfile");

	py::class_<Psrfits>(m, "Psrfits")
		.def(py::init<>())
		.def(py::init<const string, int>())
		.def("open", &Psrfits::open, py::arg("iomode")=0)
		.def("close", &Psrfits::close)
		.def("load_mode", &Psrfits::load_mode)
		.def("check_template", &Psrfits::check_template)
		.def("parse_template", &Psrfits::parse_template)
		.def("remove_hdu", &Psrfits::remove_hdu)
		.def_readwrite("filename", &Psrfits::filename)
		.def_readwrite("mode", &Psrfits::mode)
		.def_readwrite("primary", &Psrfits::primary)
		.def_readwrite("psrparam", &Psrfits::psrparam)
		.def_readwrite("t2predict", &Psrfits::t2predict)
		.def_readwrite("subint", &Psrfits::subint)
		.def_readonly("fptr", &Psrfits::fptr);

	py::class_<Filterbank>(m, "Filterbank")
		.def(py::init<>())
		.def(py::init<const string>())
		.def("free", &Filterbank::free)
		.def("close", &Filterbank::close)
		.def("read_header", &Filterbank::read_header)
		.def("read_data", static_cast<bool (Filterbank::*)()>(&Filterbank::read_data))
		.def("read_data", static_cast<bool (Filterbank::*)(long int, long int)>(&Filterbank::read_data))
		.def("read_data", static_cast<bool (Filterbank::*)(long int)>(&Filterbank::read_data))
		.def("update_rawdatafile", [](Filterbank &m, const std::string &rawdatafile){
			std::strcpy(m.rawdatafile, rawdatafile.c_str());
		})
		.def("update_source_name", [](Filterbank &m, const std::string &source_name){
			std::strcpy(m.source_name, source_name.c_str());
		})
		.def_readwrite("filename", &Filterbank::filename)
		.def_readwrite("header_size", &Filterbank::header_size)
		.def_readwrite("use_frequence_table", &Filterbank::use_frequence_table)
		.def_readwrite("telescope_id", &Filterbank::telescope_id)
		.def_readwrite("machine_id", &Filterbank::machine_id)
		.def_readwrite("data_type", &Filterbank::data_type)
		.def_readwrite("barycentric", &Filterbank::barycentric)
		.def_readwrite("pulsarcentric", &Filterbank::pulsarcentric)
		.def_readwrite("ibeam", &Filterbank::ibeam)
		.def_readwrite("nbeams", &Filterbank::nbeams)
		.def_readwrite("npuls", &Filterbank::npuls)
		.def_readwrite("nbins", &Filterbank::nbins)
		.def_readwrite("az_start", &Filterbank::az_start)
		.def_readwrite("za_start", &Filterbank::za_start)
		.def_readwrite("src_raj", &Filterbank::src_raj)
		.def_readwrite("src_dej", &Filterbank::src_dej)
		.def_readwrite("tstart", &Filterbank::tstart)
		.def_readwrite("tsamp", &Filterbank::tsamp)
		.def_readwrite("nbits", &Filterbank::nbits)
		.def_readwrite("nsamples", &Filterbank::nsamples)
		.def_readwrite("nifs", &Filterbank::nifs)
		.def_readwrite("nchans", &Filterbank::nchans)
		.def_readwrite("fch1", &Filterbank::fch1)
		.def_readwrite("foff", &Filterbank::foff)
		.def_readwrite("refdm", &Filterbank::refdm)
		.def_readwrite("period", &Filterbank::period)
		.def_readwrite("ndata", &Filterbank::ndata);
}
