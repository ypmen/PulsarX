/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2025-10-13 09:30:11
 * @modify date 2025-10-13 09:30:11
 * @desc [description]
 */

 #include <vector>
#include <utility>
#include <string>

#include <boost/program_options.hpp>

#include "dedisperse.h"
#include "filterbank.h"
#include "mjd.h"
#include "logging.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"

using namespace boost::program_options;

unsigned int num_threads;

#define NSBLK 65536

int main(int argc, const char *argv[])
{
	init_logging();

	/* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("cont", "Input files are contiguous")
			("psrfits", "Input psrfits format data")
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
		cerr<<"Error: no input file"<<endl;
		return -1;
	}

	bool contiguous = vm.count("cont");
	num_threads = vm["threads"].as<unsigned int>();
	std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();

	bool apply_wts = false;
	bool apply_scloffs = false;
	bool apply_zero_off = false;

	if (vm.count("wts"))
		apply_wts = true;
	if (vm.count("scloffs"))
		apply_scloffs = true;
	if (vm.count("zero_off"))
		apply_zero_off = true;

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	string source_name = reader->source_name;
	string s_telescope = reader->telescope;
	int ibeam = reader->beam.empty() ? 0 : std::stoi(reader->beam);
	int telescope_id = get_telescope_id(s_telescope);

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;
	long int header_size = reader->header_size;

	char start_date[32];
	reader->start_mjd.to_date(start_date, 32, "%Y/%m/%d %H:%M:%S");

	std::cout<<std::left<<std::setw(32)<<"Header size (Byte)"<<" : "<<header_size<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Telescope"<<" : "<<s_telescope<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Source Name"<<" : "<<source_name<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Source RA (J2000)"<<" : "<<reader->ra<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Source DEC (J2000)"<<" : "<<reader->dec<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Central frequency (MHz)"<<" : "<<std::setprecision(18)<<0.5 * (reader->frequencies.front()+reader->frequencies.back())<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Bandwidth (MHz)"<<" : "<<std::setprecision(18)<<(reader->frequencies.back()-reader->frequencies.front()) / (nchans - 1) * nchans<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Number of channels"<<" : "<<nchans<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Number of IFs"<<" : "<<nifs<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Sample time (us)"<<" : "<<tsamp * 1E6<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Number of samples"<<" : "<<ntotal<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Observation length (s)"<<" : "<<ntotal * tsamp<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Number of bits per sample"<<" : "<<reader->nbits<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Time stamp of first sample (MJD)"<<" : "<<std::setprecision(18)<<reader->start_mjd.to_day()<<std::endl;
	std::cout<<std::left<<std::setw(32)<<"Gregorian date (YYYY/MM/DD)"<<" : "<<start_date<<std::endl;
	
	return 0;
}