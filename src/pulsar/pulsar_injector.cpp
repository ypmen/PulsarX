/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-08-16 07:36:47
 * @modify date 2021-08-16 07:36:47
 * @desc [description]
 */

#define EIGHTBIT 1

#define NSBLK 65536

#ifdef __AVX2__
#include "avx_mathfun.h"
#include <boost/align/aligned_allocator.hpp>
#endif

#include <vector>
#include <random>
#include <boost/program_options.hpp>
#include "logging.h"
#include "dedisperse.h"
#include "filterbank.h"
#include "filterbankwriter.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "mjd.h"
#include "utils.h"

using namespace boost::program_options;

unsigned int num_threads;

/**
 * @brief gaussian profile
 * 
 * @param phase 
 * @param phasewidth 
 * @param phasecentre 
 */
inline double gaussian_profile(double phase, double phasewidth, double phasecentre=0.5)
{
	phase -= phasecentre-0.5;
	phase -= floor(phase);
	return exp(-(phase-0.5)*(phase-0.5)/(phasewidth*phasewidth));
}

/**
 * @brief Get the phase using ppdot (phase=0 at t=0)
 * 
 * @param t 
 * @param p 
 * @param pdot 
 * @return double 
 */
inline double get_phase_ppdot(double t, double p, double pdot)
{
	return 1./p*t-0.5*pdot/(p*p)*t*t;
}

inline double dmdelay(double dm, double fh, double fl)
{
	//return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
	return 4149.38*dm*(1./(fl*fl)-1./(fh*fh));
}

/**
 * @brief generate dynamic spectra using p and pdot
 * 
 * @param data 
 * @param vt 
 * @param vf 
 * @param phasewidth 
 * @param dm 
 * @param f0 
 * @param f1 
 */

#ifdef __AVX2__
void spectra_ppdot(std::vector<float, boost::alignment::aligned_allocator<float, 32>> &data, long int nsamples, double tsamp, const std::vector<double> &vf, double phasewidth, MJD &start_epoch, MJD &ref_epoch, double fref, double dm, double f0, double f1, double phase0=0.)
#else
void spectra_ppdot(std::vector<float> &data, long int nsamples, double tsamp, const std::vector<double> &vf, double phasewidth, MJD &start_epoch, MJD &ref_epoch, double fref, double dm, double f0, double f1, double phase0=0.)
#endif
{
	long int nchans = vf.size();

	data.clear();
	data.resize(nsamples*nchans, 0.);
#ifdef __AVX2__
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> phase(nchans, 0.);
	std::vector<double, boost::alignment::aligned_allocator<double, 32>> dphase(nchans, 0.);
#else
	std::vector<double> phase(nchans, 0.);
	std::vector<double> dphase(nchans, 0);
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int j=0; j<nchans; j++)
	{
		double delay = dmdelay(-dm, fref, vf[j]);
		double t = (start_epoch-ref_epoch+0.5*tsamp).to_second();
		phase[j] = f0*t+0.5*f1*t*t + delay*f0 + phase0;
		phase[j] -= floor(phase[j]);
		double ffold = f0+f1*(start_epoch-ref_epoch+0.5*nsamples*tsamp).to_second();
		dphase[j] = tsamp*ffold;
	}

#ifdef __AVX2__
	if (nchans%8 == 0)
	{   
		__m256d avx_ww = _mm256_setr_pd(1./(phasewidth*phasewidth),1./(phasewidth*phasewidth),1./(phasewidth*phasewidth),1./(phasewidth*phasewidth));
		__m256d avx_zero = _mm256_setr_pd(0.,0.,0.,0.);
		__m256d avx_half = _mm256_setr_pd(0.5,0.5,0.5,0.5);
		for (long int i=0; i<nsamples; i++)
		{
			for (long int j=0; j<nchans/8; j++)
			{
				__m256d avx_phase0 = _mm256_load_pd(&phase[0]+j*8);
				__m256d avx_phase = _mm256_sub_pd(avx_phase0, _mm256_floor_pd(avx_phase0));
				avx_phase = _mm256_sub_pd(avx_phase, avx_half);
				avx_phase = _mm256_mul_pd(avx_phase, avx_phase);
				avx_phase = _mm256_sub_pd(avx_zero, avx_phase);
				avx_phase = _mm256_mul_pd(avx_phase, avx_ww);
				__m256 avx_low = _mm256_castps128_ps256(_mm256_cvtpd_ps(avx_phase));
				__m256d avx_dphase = _mm256_load_pd(&dphase[0]+j*8);
				_mm256_store_pd(&phase[0]+j*8, _mm256_add_pd(avx_phase0, avx_dphase));

				avx_phase0 = _mm256_load_pd(&phase[0]+j*8+4);
				avx_phase = _mm256_sub_pd(avx_phase0, _mm256_floor_pd(avx_phase0));
				avx_phase = _mm256_sub_pd(avx_phase, avx_half);
				avx_phase = _mm256_mul_pd(avx_phase, avx_phase);
				avx_phase = _mm256_sub_pd(avx_zero, avx_phase);
				avx_phase = _mm256_mul_pd(avx_phase, avx_ww);
				__m128 avx_high = _mm256_cvtpd_ps(avx_phase);
				avx_dphase = _mm256_load_pd(&dphase[0]+j*8+4);
				_mm256_store_pd(&phase[0]+j*8+4, _mm256_add_pd(avx_phase0, avx_dphase));

				__m256 avx_chi = _mm256_insertf128_ps(avx_low, avx_high, 1);
				avx_chi = exp256_ps(avx_chi);
				_mm256_store_ps(&data[0]+i*nchans+j*8, avx_chi);
			}
		}
	}
	else
	{
		for (long int i=0; i<nsamples; i++)
		{
			for (long int j=0; j<nchans; j++)
			{
				data[i*nchans+j] = gaussian_profile(phase[j], phasewidth);
				phase[j] += dphase[j];
			}
		}
	}
#else
	for (long int i=0; i<nsamples; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			data[i*nchans+j] = gaussian_profile(phase[j], phasewidth);
			phase[j] += dphase[j];
		}
	}
#endif
}

int main(const int argc, const char *argv[])
{
	init_logging();

	/* options */
	int verbose = 0;
	bool outbest = false;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("dm,d", value<double>()->default_value(0), "Dispersion measure (pc/cc)")
			("F0", value<double>()->default_value(1), "Frequency at centre epoch(s)")
			("F1", value<double>()->default_value(0), "Frequency derivative (s/s)")
			("width,w", value<double>()->default_value(0.1), "Pulse width in phase (0-1)")
			("snr", value<double>()->default_value(100), "S/N")
			("spi", value<double>()->default_value(0), "Spectra index")
			("psr", value<std::string>(), "Parameter file of injected pulsars [DM F0 F1 Width S/N]")
			("baseline", value<float>()->default_value(0.01), "The scale of baseline remove (s)")
			("zapthre", value<float>()->default_value(3), "Threshold in IQR for zapping channels")
			("seglen,l", value<float>()->default_value(1), "Time length per segment (s)")
			("rootname,o", value<std::string>()->default_value("J0000-00"), "Output rootname")
			("cont", "Input files are contiguous")
			("input,f", value<std::vector<std::string>>()->multitoken()->composing(), "Input files");

	positional_options_description pos_desc;
	pos_desc.add("input", -1);
	command_line_parser parser{argc, argv};
	parser.options(desc).style(command_line_style::default_style | command_line_style::allow_short);
	parser.options(desc).positional(pos_desc);
	parsed_options parsed_options = parser.run();

	variables_map vm;
	store(parsed_options, vm);
	notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << '\n';
		return 0;
	}

	if (vm.count("verbose"))
	{
		verbose = 1;
	}

	if (vm.count("input") == 0)
	{
		BOOST_LOG_TRIVIAL(error)<<"Error: no input file"<<std::endl;
		return -1;
	}

	bool contiguous = vm.count("cont");
	num_threads = vm["threads"].as<unsigned int>();
	std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();
	std::vector<double> jump = vm["jump"].as<std::vector<double>>();
	std::string rootname = vm["rootname"].as<std::string>();

	std::vector<double> vdm = {vm["dm"].as<double>()};
	std::vector<double> vF0 = {vm["F0"].as<double>()};
	std::vector<double> vF1 = {vm["F1"].as<double>()};
	std::vector<double> vphasewidth = {vm["width"].as<double>()/2.};
	std::vector<double> vsnr = {vm["snr"].as<double>()};
	std::vector<double> vspectraindex = {vm["spi"].as<double>()};

	if (vm.count("psr"))
	{
		vdm.clear();
		vF0.clear();
		vF1.clear();
		vphasewidth.clear();
		vsnr.clear();

		std::string filename = vm["psr"].as<std::string>();
		std::string line;
		std::ifstream psrfile(filename);
		while (getline(psrfile, line))
		{
			boost::trim(line);
			if(!std::isdigit(line[0])) continue;
			
			std::vector<std::string> parameters;
			boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);
			
			vdm.push_back(std::stod(parameters[0]));
			vF0.push_back(std::stod(parameters[1]));
			vF1.push_back(std::stod(parameters[2]));
			vphasewidth.push_back(std::stod(parameters[3])/2.);
			vsnr.push_back(std::stod(parameters[4]));
			if (parameters.size()>5)
				vspectraindex.push_back(std::stod(parameters[5]));
			else
				vspectraindex.push_back(vm["spi"].as<double>());
		}
		psrfile.close();
	}

	int npsr = vdm.size();

	long int nfil = fnames.size();
	Filterbank *fil = new Filterbank [nfil];
	for (long int i=0; i<nfil; i++)
	{
		fil[i].filename = fnames[i];
	}

	PSRDataReader * reader = new FilterbankReader;

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = false;
	reader->apply_wts = false;
	reader->apply_zero_off = false;
	reader->check();
	reader->read_header();

	int td = 1;
	int fd = 1;

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;
	long double tstart = reader->start_mjd.to_day();

	MJD start_epoch((long double)tstart);
	MJD ref_epoch((long double)tstart+ntotal*tsamp*0.5/86400.);

	short *buffer = new short [nchans];

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	BOOST_LOG_TRIVIAL(info)<<"writing filterbank header...";
	Filterbank tmpfil;
	reader->get_filterbank_template(tmpfil);
	for (long int k=0; k<npsr; k++)
	{
		tmpfil.filename = rootname + "_" + std::to_string(k+1) + ".fil";
		tmpfil.tstart = tstart + nstart*tsamp/86400.;
		tmpfil.write_header();
		tmpfil.close();
	}

	std::vector<std::ofstream> outfils;
	for (long int k=0; k<npsr; k++)
	{
		ofstream outfil;
		outfil.open(rootname + "_" + std::to_string(k+1) + ".fil", ios::binary|ios::app);
		outfils.push_back(std::move(outfil));
	}

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td*td;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	databuf.frequencies = reader->frequencies;

	double fref = 0.5*(*min_element(databuf.frequencies.begin(), databuf.frequencies.end())+*max_element(databuf.frequencies.begin(), databuf.frequencies.end()));

	std::vector<unsigned char> tempdata(ndump*nchans, 0);
#ifdef __AVX2__
	std::vector<float, boost::alignment::aligned_allocator<float, 32>> injectdata(ndump*nchans, 0.);
#else
	std::vector<float> injectdata(ndump*nchans, 0.);
#endif

	int sumif = nifs>2? 2:nifs;

	if (sumif != 1)
	{
		BOOST_LOG_TRIVIAL(error)<<"data not polarization summed"<<endl;
		exit(-1);
	}

	std::mt19937 generator;
	std::normal_distribution<double> distribution(0., 1.);

	BOOST_LOG_TRIVIAL(info)<<"injecting data...";

	size_t count = 0;
	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;
		databuf.counter += ndump;

		databuf.get_mean_rms();

		count = reader->get_count();
		
		double offset = (long double)(count-ndump)*(long double)tsamp;
		MJD tmp_epoch = start_epoch + offset;
		double factor1 = sqrt(4*0.989939)/(0.842732*sqrt(2.*M_PI));
		double sum2 = 0.;
		for (long int jj=0; jj<nchans; jj++)
		{
			sum2 += databuf.vars[jj];
		}
		double fmin = *std::min_element(databuf.frequencies.begin(), databuf.frequencies.end());
		for (long int ipsr=0; ipsr<npsr; ipsr++)
		{
			spectra_ppdot(injectdata, ndump, tsamp, databuf.frequencies, vphasewidth[ipsr], tmp_epoch, ref_epoch, fref, vdm[ipsr], vF0[ipsr], vF1[ipsr], 0.);
			std::vector<double> normf(nchans, 1.);
			double sum1 = 0.;
			for (long int jj=0; jj<nchans; jj++)
			{
				normf[jj] = std::pow(databuf.frequencies[jj]/fmin, -vspectraindex[ipsr]);
				sum1 += std::sqrt(databuf.vars[jj])*normf[jj];
			}
			double factor2 = sum1/std::sqrt(sum2);
			double alpha = vsnr[ipsr]/(factor1*factor2*std::sqrt(ntotal*vphasewidth[ipsr]));
			std::vector<double> Af(nchans, 0.);
			for (long int jj=0; jj<nchans; jj++)
			{
				Af[jj] = alpha*normf[jj]*(1.+std::sqrt(databuf.vars[jj]));
			}
			for (long int ii=0; ii<ndump; ii++)
			{
				for (long int jj=0; jj<nchans; jj++)
				{
					tempdata[ii*nchans+jj] = std::round(databuf.buffer[ii*nchans+jj] + (Af[jj]*injectdata[ii*nchans+jj])+distribution(generator));
				}
			}
			outfils[ipsr].write((char *)tempdata.data(), sizeof(unsigned char)*ndump*nchans);
		}

		std::fill(injectdata.begin(), injectdata.end(), 0);
	}

	if (verbose)
	{
		std::cerr<<"\r\rfinish "<<std::setprecision(2)<<std::fixed<<tsamp*count<<" seconds ";
		std::cerr<<"("<<100.*count/ntotal<<"%)"<<std::endl;
	}

	BOOST_LOG_TRIVIAL(info)<<"Done!";

	delete [] buffer;
	delete [] fil;

	return 0;
}