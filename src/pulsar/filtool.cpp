/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-05-27 20:57:23
 * @modify date 2021-05-27 20:57:23
 * @desc [description]
 */

#include <vector>
#include <utility>
#include <string>

#include <boost/program_options.hpp>

#include "dedisperse.h"
#include "databuffer.h"
#include "patch.h"
#include "preprocesslite.h"
#include "filmaker.h"
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
			("jump,j", value<std::vector<double>>()->multitoken()->default_value(std::vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("td", value<int>()->default_value(1), "Time downsample for preprocessing")
			("fd", value<int>()->default_value(1), "Frequency downsample for preprocessing")
			("nbits", value<int>()->default_value(8), "Nbits of output filterbank")
			("mean", value<float>()->default_value(128), "Mean value of output filterbank")
			("std", value<float>()->default_value(6), "Standard deviation of output filterbank")
			("filplan", value<std::string>(), "Input filterbank plan file")
			("seglen,l", value<float>()->default_value(2), "Time length per segment (s)")
			("zapthre", value<float>()->default_value(4), "Threshold in IQR for zapping channels")
			("ra", value<string>()->default_value("00:00:00"), "RA (hhmmss.s)")
			("dec", value<string>()->default_value("00:00:00"), "DEC (ddmmss.s)")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("telescope", value<std::string>()->default_value("Fake"), "Telescope name")
			("outref", value<std::string>(), "Reference baseline data")
			("baseline", value<std::vector<float>>()->multitoken()->default_value(std::vector<float>{0.0, 0.0}, "0.0, 0.0"), "The scale of baseline remove (s)")
			("rfi,z", value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("threMask", value<float>()->default_value(10), "S/N threshold of Mask")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("threPatch", value<float>()->default_value(10), "IQR threshold of patch for bad data")
			("widthPatch", value<float>()->default_value(0.4), "Width threshold (s) of patch for bad data")
			("fillPatch", value<std::string>()->default_value("none"), "Fill the bad data by [none, mean, rand] in patch")
			("fill", value<string>()->default_value("mean"), "Fill the zapped samples by [mean, rand]")
			("source_name,s", value<std::string>()->default_value("J0000-00"), "Source name")
			("rootname,o", value<std::string>()->default_value("J0000-00"), "Output rootname")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("mask", "save rfi mask")
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

	long double tstart = reader->start_mjd.to_day();

	string source_name = reader->source_name;
	string s_telescope = reader->telescope;
	int ibeam = reader->beam.empty() ? 0 : std::stoi(reader->beam);
	int telescope_id = get_telescope_id(s_telescope);
	double src_raj = 0., src_dej;
	if (!reader->ra.empty())
	{
		string ra = reader->ra;
		ra.erase(remove(ra.begin(), ra.end(), ':'), ra.end());
		src_raj = stod(ra);
	}

	if (!reader->dec.empty())
	{
		string dec = reader->dec;
		dec.erase(remove(dec.begin(), dec.end(), ':'), dec.end());
		src_dej = stod(dec);
	}

	std::vector<FilMaker> filmakers;
	plan(vm, filmakers);

	/** pre-downsample */
	int td = vm["td"].as<int>();
	int fd = vm["fd"].as<int>();

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	vector<int> tds;
	for (auto fm=filmakers.begin(); fm!=filmakers.end(); ++fm)
	{
		tds.push_back((*fm).td*td);
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm/2*td_lcm*2;

	// read zero dm data
	std::vector<float> outref;
	long double tstart_zero = 0.;
	double tsamp_zero = 0.;
	if (vm.count("outref"))
	{
		Filterbank fil(vm["outref"].as<std::string>());
		fil.read_header();
		fil.read_data();

		assert(fil.nchans == 1);

		tsamp_zero = fil.tsamp;
		tstart_zero = fil.tstart;
		outref.resize(fil.nsamples, 0.);
		for (size_t i=0; i<outref.size(); i++)
		{
			outref[i] = ((float *)fil.data)[i];
		}

		fil.close();

		//assert(tsamp_zero == tsamp);
	}

	DataBuffer<float> databuf(ndump, nchans);
	databuf.closable = false;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader->frequencies.data(), sizeof(double)*nchans);

	Patch patch;
	patch.filltype = vm["fillPatch"].as<string>();
	patch.width = vm["widthPatch"].as<float>();
	patch.threshold = vm["threPatch"].as<float>();
	patch.prepare(databuf);
	patch.close();

	PreprocessLite prep;
	prep.td = td;
	prep.fd = fd;
	prep.thresig = vm["zapthre"].as<float>();
	prep.filltype = vm["fill"].as<string>();
	prep.prepare(databuf);

	// handle input data
	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;
	reader->skip_start = nstart;
	reader->skip_end = ntotal - nend;
	reader->skip_head();

	float ndiv = reader->tsamp / tsamp_zero;

	int td_ref = std::round(prep.td * reader->tsamp / tsamp_zero);

	std::vector<float> outref_ds(outref.size()/td_ref, 0.); 
	if (vm.count("outref"))
	{
		// align start sample
		if (tstart_zero > tstart + nstart * tsamp / 86400.)
		{
			reader->skip_start = std::round((tstart_zero - tstart) * 86400. / tsamp);
			reader->skip_head();
		}
		else
		{
			size_t zero_offset = std::round((tstart + nstart * tsamp / 86400. - tstart_zero) * 86400. / tsamp_zero);
			for (size_t i=zero_offset; i<outref.size(); i++)
			{
				outref[i-zero_offset] = outref[i];
			}
		}

		// align end sample
		if (tstart_zero + outref.size() * tsamp_zero / 86400. < tstart + ntotal * tsamp / 86400.)
		{
			reader->skip_end = std::round(((tstart + ntotal * tsamp / 86400.) - (tstart_zero + outref.size() * tsamp_zero / 86400.)) * 86400. / tsamp);
		}

		for (size_t k=0; k<td_ref; k++)
		{
			for (size_t i=0; i<outref.size()/td_ref; i++)
			{
				outref_ds[i] += outref[i * td_ref + k];
			}
		}
	}

	std::vector<float> bswidth = vm["baseline"].as<std::vector<float>>();

	long int noutfil = filmakers.size();
	for (long int k=0; k<noutfil; k++)
	{
		reader->get_filterbank_template(filmakers[k].filwriter.fil);
		filmakers[k].ibeam = ibeam;
		filmakers[k].rootname = rootname;
		filmakers[k].source_name = source_name;
		filmakers[k].src_raj = src_raj;
		filmakers[k].src_dej = src_dej;
		filmakers[k].telescope_id = telescope_id;
		filmakers[k].filltype = vm["fill"].as<string>();

		filmakers[k].bswidth = bswidth[1];

		if (vm.count("outref"))
		{
			filmakers[k].outref = outref_ds;
		}

		filmakers[k].prepare(prep);
	}

	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		databuf.counter += ndump;

		DataBuffer<float> *data = patch.filter(databuf);
		data = prep.run(*data);

		for (long int ioutfil=0; ioutfil<noutfil; ioutfil++)
		{
			data->isbusy = true;
			filmakers[ioutfil].run(*data);
		}

		//databuf.open();
	}

	databuf.close();

	if (vm.count("mask"))
	{
		Filterbank maskfil;
		reader->get_filterbank_template(maskfil);
		maskfil.filename = rootname + "_mask.fil";
		std::strcpy(maskfil.source_name, source_name.c_str());
		maskfil.telescope_id = telescope_id;
		maskfil.src_raj = src_raj;
		maskfil.src_dej = src_dej;
		maskfil.ibeam = ibeam;
		maskfil.nbits = 8;
		maskfil.tsamp = ndump * tsamp;

		maskfil.write_header();
		for (auto tmpmask=prep.mask.begin(); tmpmask!=prep.mask.end(); ++tmpmask)
		{
			fwrite(tmpmask->data(), 1, sizeof(unsigned char) * tmpmask->size(), maskfil.fptr);
		}

		maskfil.close();
	}

	return 0;
}