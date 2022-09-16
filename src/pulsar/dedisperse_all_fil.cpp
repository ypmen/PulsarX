/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-11-04 20:47:29
 * @modify date 2020-11-04 20:47:29
 * @desc [description]
 */

#define FAST 1
#define NSBLK 8192

#include <iostream>
#include <iomanip>
#include <string.h>
#include <boost/program_options.hpp>
#include <sys/time.h>
#include <sys/resource.h>

#include "pulsarsearch.h"
#include "subdedispersion.h"
#include "dedisperse.h"
#include "filterbank.h"
#include "utils.h"
#include "mjd.h"
#include "patch.h"
#include "preprocesslite.h"
#include "logging.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
unsigned int dbscan_radius;
unsigned int dbscan_k;

bool dumptim=true;

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
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{28800, 0}, "7200, 0"), "Segment length and following jump length (s)")
			("td", value<int>()->default_value(1), "Time downsample for preprocessing")
			("fd", value<int>()->default_value(1), "Frequency downsample for preprocessing")
			("zapthre", value<float>()->default_value(4), "Threshold in IQR for zapping channels")
			("dms", value<double>()->default_value(0), "DM start")
			("ddm", value<double>()->default_value(1), "DM step")
			("ndm", value<int>()->default_value(200), "Number of DM")
			("ddplan", value<string>(), "Input ddplan file")
			("seglen,l", value<float>()->default_value(2), "Time length per segment (s)")
			("mean", value<float>()->default_value(0), "Mean of dedispersed time series")
			("std", value<float>()->default_value(3), "Standard deviation of dedispersed time series")
			("nbits", value<int>()->default_value(8), "Data type of dedispersed time series")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("baseline", value<vector<float>>()->multitoken()->default_value(vector<float>{0.0, 0.0}, "0.0, 0.0"), "The scale of baseline remove (s)")
			("rfi,z", value<vector<string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("threMask", value<float>()->default_value(3), "S/N threshold of Mask")
			("threPatch", value<float>()->default_value(10), "IQR threshold of patch for bad data")
			("widthPatch", value<float>()->default_value(0.4), "Width threshold (s) of patch for bad data")
			("fillPatch", value<std::string>()->default_value("none"), "Fill the bad data by [none, mean, rand] in patch")
			("fill", value<string>()->default_value("mean"), "Fill the zapped samples by [mean, rand]")
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("format", value<string>()->default_value("pulsarx"), "Output format of dedispersed data [pulsarx(default),sigproc,presto]")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("cont", "Input files are contiguous")
			("psrfits", "Input psrfits format data")
			("input,f", value<vector<string>>()->multitoken()->composing(), "Input files");

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
		BOOST_LOG_TRIVIAL(error)<<"Error: no input file"<<endl;
		return -1;
	}

	bool contiguous = vm.count("cont");
	string rootname = vm["rootname"].as<string>();
	num_threads = vm["threads"].as<unsigned int>();
	vector<double> jump = vm["jump"].as<vector<double>>();
	vector<string> fnames = vm["input"].as<vector<string>>();
	float outmean = vm["mean"].as<float>();
	float outstd = vm["std"].as<float>();
	int outnbits = vm["nbits"].as<int>();

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

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();
	reader->nsblk = 8192;

	int ibeam = reader->beam.empty() ? 0 : std::stoi(reader->beam);

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	vector<PulsarSearch> search;
	plan(vm, search);

	size_t ndm_total = 0;

	vector<int> tds;
	for (auto sp=search.begin(); sp!=search.end(); ++sp)
	{
		tds.push_back((*sp).td*vm["td"].as<int>());

		ndm_total += sp->ndm;
	}

	if (vm["format"].as<string>() == "sigproc" or vm["format"].as<string>() == "presto")
	{
		struct rlimit rlim;
		int status = getrlimit(RLIMIT_NOFILE, &rlim);
		if (status)
			BOOST_LOG_TRIVIAL(warning)<<"Can't get the maximum file descriptor number";
		int rlimit_nofile = rlim.rlim_cur;
		if (ndm_total > rlimit_nofile)
		{
			BOOST_LOG_TRIVIAL(error)<<ndm_total<<" dedisersed files is larger than the maximum file descriptor number, please check ulimit -a";
			exit(-1);
		}

		BOOST_LOG_TRIVIAL(info)<<"create "<<ndm_total<<" "<<vm["format"].as<string>()<<" dedispersed files";
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm/2*td_lcm*2;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.closable = true;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader->frequencies.data(), sizeof(double)*nchans);

	Patch patch;
	patch.filltype = vm["fillPatch"].as<string>();
	patch.width = vm["widthPatch"].as<float>();
	patch.threshold = vm["threPatch"].as<float>();
	patch.prepare(databuf);
	patch.close();

	PreprocessLite prep;
	prep.td = vm["td"].as<int>();
	prep.fd = vm["fd"].as<int>();
	prep.thresig = vm["zapthre"].as<float>();
	prep.filltype = vm["fill"].as<string>();
	prep.prepare(databuf);

	long int nseg = std::ceil(jump[0]/tsamp / ndump) * ndump;
	long int njmp = std::ceil(jump[1]/tsamp / ndump) * ndump;

	long int ncover = 0;

	stringstream ss_ibeam;
	if (vm.count("incoherent"))
		ss_ibeam << "ifbf" << setw(5) << setfill('0') << ibeam;
	else
		ss_ibeam << "cfbf" << setw(5) << setfill('0') << ibeam;
	string s_ibeam = ss_ibeam.str();

	ncover++;
	long int nsearch = search.size();
	for (long int k=0; k<nsearch; k++)
	{
		search[k].ibeam = ibeam;
		search[k].rootname = rootname + "_" + s_ibeam + "_Plan" + to_string(k+1) + "_" + to_string(ncover);
		reader->get_filterbank_template(search[k].fildedisp);
		search[k].fildedisp.fch1 = (databuf.frequencies.front()+databuf.frequencies.back())/2.;
		search[k].fildedisp.foff = databuf.frequencies.back()-databuf.frequencies.front();
		search[k].fildedisp.tstart = reader->start_mjd.to_day();
		search[k].filltype = vm["fill"].as<string>();
		search[k].prepare(prep);
	}

	long int ntot = 0;
	long int jmpcont = 0;
	long int count = 0;
	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		if (ntot == nseg)
		{
			if (jmpcont < njmp)
			{
				jmpcont += ndump;
				continue;
			}
			
			ntot = 0;
			jmpcont = 0;

			ncover++;
			for (long int k=0; k<nsearch; k++)
			{
				if (vm["format"].as<string>() == "presto")
				{
					search[k].dedisp.makeinf(search[k].fildedisp);
				}
				else if (vm["format"].as<string>() != "sigproc")
				{
					search[k].dedisp.modifynblock();
				}
				search[k].dedisp.rootname = rootname + "_" + s_ibeam + "_Plan" + to_string(k+1) + "_" + to_string(ncover);
				search[k].dedisp.prepare(search[k].rfi);
				search[k].fildedisp.tstart = (reader->start_mjd + count*tsamp/86400.).to_day();
				search[k].dedisp.preparedump(search[k].fildedisp, outnbits, vm["format"].as<string>());
			}
		}

		ntot += ndump;
		count += ndump;

		patch.filter(databuf);
		prep.run(databuf);
		for (auto sp=search.begin(); sp!=search.end(); ++sp)
		{
			prep.isbusy = true;
			(*sp).run(prep);
		}

		databuf.open();
	}

	if (vm["format"].as<string>() == "presto")
	{
		for (auto sp=search.begin(); sp!=search.end(); ++sp)
		{
			(*sp).dedisp.makeinf((*sp).fildedisp);
		}
	}
	else if (vm["format"].as<string>() != "sigproc")
	{
		for (auto sp=search.begin(); sp!=search.end(); ++sp)
		{
			(*sp).dedisp.modifynblock();
		}
	}

	return 0;
}
