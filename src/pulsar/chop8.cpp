/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2023-07-09 11:02:52
 * @modify date 2023-07-09 11:02:52
 * @desc [description]
 */

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
using namespace boost::program_options;

#include "logging.h"
#include "psrdatareader.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"

unsigned int num_threads;

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
			("skip,s", value<double>()->default_value(0), "Skip time length (s)")
			("read,r", value<double>()->default_value(1), "Read time length (s)")
			("cont", "Input files are contiguous")
			("psrfits", "Input psrfits format data")
			("output,o", value<std::string>()->default_value("chop"), "Output rootname")
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
		std::cerr<<"Error: no input file"<<std::endl;
		return -1;
	}

	bool contiguous = vm.count("cont");

	num_threads = vm["threads"].as<unsigned int>();

	std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	reader->fnames = fnames;
	reader->sumif = false;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = false;
	reader->apply_wts = false;
	reader->apply_zero_off = false;
	reader->check();
	reader->read_header();

	long double tstart = reader->start_mjd.to_day();
	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;
	double fch1 = 0.;
	double foff = 0.;
	std::string source_name = reader->source_name;
	reader->get_fch1_foff(fch1, foff);

	long int nstart = std::floor(vm["skip"].as<double>() / tsamp);
	long int nend = std::floor((vm["skip"].as<double>() + vm["read"].as<double>()) / tsamp);

	reader->skip_start = nstart;
	reader->skip_end = nend > ntotal ? 0 : ntotal - nend;

	reader->skip_head();

	std::string rootname = vm["output"].as<std::string>();

	int ndump = 1024;

	DataBuffer<unsigned char> databuf(ndump, nifs*nchans);

	Filterbank fil;
	fil.filename = rootname + ".fil";
	std::strcpy(fil.source_name, source_name.c_str());
	fil.tstart += tstart + nstart * tsamp / 86400.;
	fil.tsamp = tsamp;
	fil.nchans = nchans;
	fil.nbits = 8;
	fil.nifs = nifs;
	fil.fch1 = fch1;
	fil.foff = foff;

	if (!fil.write_header())
		BOOST_LOG_TRIVIAL(error)<<"Error: Can not write filterbank header";

	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		fwrite(databuf.buffer.data(), 1, sizeof(unsigned char) * ndump * nifs * nchans, fil.fptr);
	}

	fil.close();
}