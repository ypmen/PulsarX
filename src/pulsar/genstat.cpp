/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-04-10 20:59:51
 * @modify date 2024-04-10 20:59:51
 * @desc [description]
 */

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "stat2.h"
#include "logging.h"
#include "json.hpp"

using namespace boost::program_options;
using namespace std;

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
			("jump,j", value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("seglen,l", value<float>()->default_value(2), "Time length per segment (s)")
			("ra", value<string>()->default_value("00:00:00"), "RA (hhmmss.s)")
			("dec", value<string>()->default_value("00:00:00"), "DEC (ddmmss.s)")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("telescope", value<std::string>()->default_value("Fake"), "Telescope name")
			("source_name,s", value<std::string>()->default_value("J0000-00"), "Source name")
			("rootname,o", value<std::string>()->default_value("J0000-00"), "Output rootname")
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
	std::vector<double> jump = vm["jump"].as<std::vector<double>>();
	std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();
	std::string rootname = vm["rootname"].as<std::string>();

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

	if (!vm["ibeam"].defaulted())
		reader->beam = std::to_string(vm["ibeam"].as<int>());
	if (!vm["telescope"].defaulted())
		reader->telescope = vm["telescope"].as<std::string>();
	if (!vm["ra"].defaulted())
		reader->ra = vm["ra"].as<std::string>();
	if (!vm["dec"].defaulted())
		reader->dec = vm["dec"].as<std::string>();
	if (!vm["srcname"].defaulted())
		reader->source_name = vm["source_name"].as<std::string>();

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;
	reader->skip_start = nstart;
	reader->skip_end = ntotal - nend;
	reader->skip_head();

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp);

	DataBuffer<unsigned char> databuf(ndump, nchans);
	databuf.closable = false;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader->frequencies.data(), sizeof(double)*nchans);

	Stat2 stat;
	stat.prepare(databuf);

	std::ofstream outdat;
	outdat.open(rootname + "_hist.dat", std::ios::binary);

	size_t nsubint = 0;
	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		databuf.counter += ndump;

		stat.run(databuf);

		outdat.write((char *)(stat.hists.data()), stat.hists.size() * sizeof(int));

		nsubint++;
	}

	outdat.close();

	nlohmann::json info = {
		{"dtype", "int32"},
		{"nbin", 256},
		{"nchans", nchans},
		{"nsubint", nsubint},
		{"tsubint", ndump * tsamp},
		{"fch1", stat.frequencies[0]},
		{"foff", stat.frequencies[1] - stat.frequencies[0]} 
	};

	std::string s_info = info.dump();

	std::ofstream outheader(rootname + "_hist.inf");
    outheader << s_info;
    outheader.close();

	return 0;
}