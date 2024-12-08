/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-10 21:48:03
 * @modify date 2022-08-10 21:48:03
 * @desc [description]
 */

#define EIGHTBIT 1
#define NSBLK 65536

#include "config.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <utility>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "logging.h"
#include "dedisperse.h"
#include "archivewriter.h"

#include "pulsarplot.h"

#include "gridsearch.h"
#include "archivelite.h"
#include "dedispersionX.h"
#include "databuffer.h"
#include "downsample.h"
#include "rfi.h"
#include "equalize.h"
#include "filterbank.h"
#include "mjd.h"
#include "utils.h"
#include "constants.h"
#include "preprocesslite.h"
#include "baseline.h"
#include "patch.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
bool dumptim=false;

void produce(variables_map &vm, std::list<double> &dmlist, vector<Pulsar::ArchiveLite> &folder);

int main(int argc, const char *argv[])
{
	init_logging();

	/* options */
	int verbose = 0;
	bool outbest = false;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("save_memory", "Trade off memory with speed")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("frac", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 1}, "0, 1"), "Reading data between start and end fraction")
			("td", value<int>()->default_value(1), "Time downsample")
			("fd", value<int>()->default_value(1), "Frequency downsample")
			("zapthre", value<float>()->default_value(4), "Threshold in IQR for zapping channels")
			("ddminit", value<double>()->default_value(0.), "Manually set initial dm step")
			("dmboost", value<double>()->default_value(200.), "Manually set boost dm after which optimal dm step is used")
			("dm", value<double>()->default_value(0), "DM (pc/cc)")
			("f0", value<double>()->default_value(0), "F0 (Hz)")
			("f1", value<double>()->default_value(0), "F1 (Hz/s)")
			("acc", value<double>()->default_value(0), "Acceleration (m/s/s)")
			("pepoch", value<double>(), "F0/F1/acc epoch (MJD)")
			("scale", value<int>()->default_value(1), "F0,F1,dm search range scale in phase")
			("nosearch", "Do not search dm,f0,f1")
			("nodmsearch", "Do not search dm")
			("nof0search", "Do not search f0")
			("nof1search", "Do not search f1")
			("noplot", "Do not generate figures")
			("noarch", "Do not generate archives")
			("parfile", value<vector<string>>(), "Input pulsar par files")
			("t2pred", value<vector<string>>(), "Input t2pred files")
			("candfile", value<string>(), "Input cand file")
			("template", value<string>(), "Input fold template file")
			("nbin,b", value<int>()->default_value(32), "Number of bins per period")
			("nbinplan", value<vector<float>>()->multitoken()->default_value(vector<float>{0.1, 128}, "0.1, 128"), "Using nbins of nbin,nbin1,nbin2,... for 0,p1,p2,... (s)")
			("tsubint,L", value<double>()->default_value(10), "Time length per integration (s)")
			("nsubband,n", value<int>()->default_value(32), "Number of subband")
			("blocksize", value<int>()->default_value(2), "Size of data block (s)")
			("srcname", value<string>()->default_value("PSRJ0000+00"), "Souce name")
			("telescope", value<string>()->default_value("Fake"), "Telescope name")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("ra", value<string>()->default_value("00:00:00"), "RA (hhmmss.s)")
			("dec", value<string>()->default_value("00:00:00"), "DEC (ddmmss.s)")
			("baseline", value<vector<float>>()->multitoken()->default_value(vector<float>{0.0, 0.0}, "0.0, 0.0"), "The scale of baseline remove (s)")
			("clfd", value<double>()->default_value(-1), "CLFD q value, if q<=0, CLFD will not be applied")
			("bandcorr", value<double>()->default_value(10), "CLFD correlation bandwidth threshold to differentiate RFI and real signal (MHz)")
			("rfi,z", value<vector<string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("tdRFI", value<int>()->default_value(1), "Time downsample of RFI")
			("fdRFI", value<int>()->default_value(1), "Frequency downsample of RFI")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("threMask", value<float>()->default_value(10), "S/N threshold of Mask")
			("threPatch", value<float>()->default_value(10), "IQR threshold of patch for bad data")
			("widthPatch", value<float>()->default_value(0.4), "Width threshold (s) of patch for bad data")
			("fillPatch", value<std::string>()->default_value("none"), "Fill the bad data by [none, mean, rand] in patch")
			("fill", value<string>()->default_value("rand"), "Fill the zapped samples by [mean, rand]")
			("render", "Using new folding algorithm (default)")
			("dspsr", "Using dspsr folding algorithm")
			("presto", "Using presto folding algorithm (in test phase, will be default)")
			("plotx", "Using PlotX for plotting (default; not used any more)")
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
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
	if (vm.count("outbest"))
	{
		outbest = true;
	}
	if (vm.count("input") == 0)
	{
		BOOST_LOG_TRIVIAL(error)<<"Error: no input file"<<endl;
		return -1;
	}
	if (vm.count("template") == 0)
	{
		BOOST_LOG_TRIVIAL(error)<<"Error: no template file"<<endl;
		return -1;
	}
	else
	{
		std::ifstream tmp(vm["template"].as<std::string>().c_str());
		if (!tmp.good())
		{
			BOOST_LOG_TRIVIAL(error)<<"Error: can not open file "<<vm["template"].as<std::string>()<<endl;
			return -1;
		}
	}

	int scale = vm["scale"].as<int>();
	bool nosearch = vm.count("nosearch");
	bool nodmsearch = vm.count("nodmsearch");
	bool nof0search = vm.count("nof0search");
	bool nof1search = vm.count("nof1search");
	if (nodmsearch && nof0search && nof1search) nosearch = true;
	if (nosearch)
	{
		nodmsearch = true;
		nof0search = true;
		nof1search = true;
	}
	bool noplot = vm.count("noplot");
	bool noarch = vm.count("noarch");
	bool contiguous = vm.count("cont");
	string rootname = vm["rootname"].as<string>();
	num_threads = vm["threads"].as<unsigned int>();
	vector<double> jump = vm["jump"].as<vector<double>>();
	vector<string> fnames = vm["input"].as<vector<string>>();

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
		reader->source_name = vm["srcname"].as<std::string>();

	reader->fnames = fnames;
	reader->sumif = true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	string src_name = reader->source_name;
	string s_telescope = reader->telescope;
	int ibeam = reader->beam.empty() ? 0 : std::stoi(reader->beam);

	/** downsample */
	int td = vm["td"].as<int>();
	int fd = vm["fd"].as<int>();

	/** rfi */
	vector<pair<double, double>> zaplist;
	vector<vector<string>> rfilist;
	vector<string> rfi_opts;
	if (vm.count("rfi"))
	{
		rfi_opts = vm["rfi"].as<vector<string>>();
		for (auto opt=rfi_opts.begin(); opt!=rfi_opts.end(); ++opt)
		{
			if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
			{
				vector<string> temp{*opt, *(opt+1), *(opt+2)};       
				rfilist.push_back(temp);
				advance(opt, 2);
			}
			else if (*opt == "zap")
			{
				zaplist.push_back(pair<double, double>(stod(*(opt+1)), stod(*(opt+2))));
				advance(opt, 2);
			}
			else if (*opt=="zero" or *opt=="zdot")
			{
				vector<string> temp{*opt};
				rfilist.push_back(temp);
			}
		}
	}

	double bandlimit = vm["bandlimit"].as<double>();
	double bandlimitKT = vm["bandlimitKT"].as<double>();
	double widthlimit = vm["widthlimit"].as<double>();
	int tdRFI = vm["tdRFI"].as<int>();
	int fdRFI = vm["fdRFI"].as<int>();
	float threKadaneF = vm["threKadaneF"].as<float>();
	float threKadaneT = vm["threKadaneT"].as<float>();
	float threMask = vm["threMask"].as<float>();

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp;
	int nifs = reader->nifs;
	long int ntotal = reader->nsamples;

	vector<Pulsar::ArchiveLite> folder;

	std::list<double> dmlist;
	produce(vm, dmlist, folder);

	double maxdm = *std::max_element(dmlist.begin(), dmlist.end());
	double fmin = 1e6;
	double fmax = 0.;
	reader->get_fmin_fmax(fmin, fmax);

	double ddm_init = vm["ddminit"].as<double>();
	double dm_boost = vm["dmboost"].as<double>();
	size_t maxds = Pulsar::DedispersionX::get_maxds(maxdm, tsamp * td, fmax, fmin, nchans / fd, ddm_init, dm_boost);

	int blocksize = vm["blocksize"].as<int>();
	long int ndump = ceil((blocksize/tsamp)/(td*maxds))*(td*maxds);
	int nblock = ceil(vm["tsubint"].as<double>())/blocksize;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.closable = true;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader->frequencies.data(), sizeof(double)*nchans);

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	double startfrac = vm["frac"].as<vector<double>>()[0];
	double endfrac = vm["frac"].as<vector<double>>()[1];
	if (startfrac != 0 or endfrac != 1)
	{
		nstart = startfrac*ntotal;
		nend = endfrac*ntotal;
	}

	reader->skip_start = nstart;
	reader->skip_end = ntotal - nend;
	reader->skip_head();

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

	Equalize equalize;
	equalize.prepare(prep);
	equalize.close();
	equalize.closable = true;

	BaseLine baseline;
	baseline.width = vm["baseline"].as<vector<float>>().back();
	baseline.prepare(equalize);
	baseline.close();
	baseline.closable = true;

	RFI rfi;
	rfi.filltype = vm["fill"].as<string>();
	rfi.prepare(baseline);
	rfi.close();
	rfi.closable = true;

	std::string rfi_flags;
	for (auto irfi = rfilist.begin(); irfi!=rfilist.end(); ++irfi)
	{
		for (auto r = irfi->begin(); r!=irfi->end(); ++r)
		rfi_flags += *r + " ";
	}
	std::vector<std::pair<std::string, std::string>> meta = {
		{"nsamples", std::to_string(rfi.nsamples)},
		{"nchans", std::to_string(rfi.nchans)},
		{"tsamp", std::to_string(rfi.tsamp)},
		{"rfi_flags", rfi_flags},
		{"kadaneF_snr_thre", std::to_string(threKadaneF)},
		{"kadaneF_width_thre", std::to_string(widthlimit)},
		{"filltype", vm["fill"].as<string>()}
	};
	format_logging("RFI Mitigation Info", meta);

	/** generate tempo2 predictor file */
	if (vm.count("t2pred"))
	{
		BOOST_LOG_TRIVIAL(info)<<"read tempo2 predictor file";

		std::vector<std::string> t2preds = vm["t2pred"].as<std::vector<std::string>>();
		
		for (long int k=0; k<folder.size(); k++)
		{
			folder[k].use_t2pred = true;
			folder[k].pred.read_t2pred(t2preds[k]);
		}
	}
	else if (vm.count("parfile"))
	{
		BOOST_LOG_TRIVIAL(info)<<"generate tempo2 predictor file";

		std::vector<std::string> parfiles = vm["parfile"].as<std::vector<std::string>>();
		
		/**
		 * @brief read DM from parfile
		 * 
		 */
		for (long int k=0; k<folder.size(); k++)
		{
			string filename = parfiles[k];

			double mjd_start = (reader->mjd_starts[reader->idmap.front()].to_day())-0.01;
			double mjd_end = (reader->mjd_ends[reader->idmap.back()].to_day())+1;
			double fmin = databuf.frequencies[0];
			double fmax = databuf.frequencies[0];
			for (long int j=0; j<nchans; j++)
			{
				fmin = fmin<databuf.frequencies[j] ? fmin:databuf.frequencies[j];
				fmax = fmax>databuf.frequencies[j] ? fmax:databuf.frequencies[j];
			}
			double chanwidth = abs(databuf.frequencies[1]-databuf.frequencies[0]);
			fmin -= 0.5*chanwidth;
			fmax += 0.5*chanwidth;

			string cmd = "tempo2 -npsr 1 -f " + filename + " -pred ";
			cmd += "'";
			if (s_telescope.empty())
				cmd += "fake";
			else
				cmd += s_telescope;
			cmd += " ";
			cmd += to_string(mjd_start);
			cmd += " ";
			cmd += to_string(mjd_end);
			cmd += " ";
			cmd += to_string(fmin);
			cmd += " ";
			cmd += to_string(fmax);
			cmd += " ";
			cmd += "12 2";
			cmd += " ";
			cmd += "3600";
			cmd += "'";
			cmd += " > /tmp/stdout.txt 2> /tmp/stderr.txt";
			BOOST_LOG_TRIVIAL(info)<<cmd;
			system(cmd.c_str());

			folder[k].use_t2pred = true;
			folder[k].pred.read_t2pred("t2pred.dat");

			system("rm pred.tim t2pred.dat");
		}
	}

	Pulsar::DedispersionX dedisp;
	dedisp.ddm_init = ddm_init;
	dedisp.dm_boost = dm_boost;

	dedisp.nsubband = vm["nsubband"].as<int>();
	dedisp.dmlist = dmlist;
	dedisp.prepare(rfi);

	size_t ndmseg = dedisp.get_nsegments();

	if (vm.count("save_memory") == 0)
	{
		for (size_t l=0; l<ndmseg; l++)
		{
			dedisp.alloc(l);
		}
	}

	// calculate dm batch
	std::vector<std::vector<size_t>> batchs;
	for (size_t l=0; l<ndmseg; l++)
	{
		double dms_tmp = dedisp.get_dms(l);
		double ddm_tmp = dedisp.get_ddm(l);

		std::vector<size_t> batch;
		for (long int k=0; k<folder.size(); k++)
		{
			double dm_tmp = folder[k].dm;
			if ((dm_tmp >= dms_tmp) && (dm_tmp < dms_tmp + rfi.nchans * ddm_tmp))
			{
				batch.push_back(k);
			}
		}
		batchs.push_back(batch);
	}

	long int ncand = folder.size();
	if (ncand <= 0)
	{
		BOOST_LOG_TRIVIAL(warning)<<"no candidate to fold"<<std::endl;
		exit(0);
	}

	std::vector<size_t> offsets(ncand, 0);
	std::vector<size_t> ndumps(ncand, 0);
	std::vector<double> tsamps(ncand, 0.);
	std::vector<size_t> ids(ncand, 0.);

	for (long int k=0; k<ncand; k++)
	{
		DataBuffer<float> tmpsubdata;
		dedisp.get_subdata_tem(folder[k].dm, tmpsubdata);

		offsets[k] = dedisp.get_offset(folder[k].dm);
		ndumps[k] = dedisp.get_ndump(folder[k].dm);
		tsamps[k] = dedisp.get_tsamp(folder[k].dm);
		ids[k] = dedisp.get_id(folder[k].dm);

		folder[k].start_mjd = reader->start_mjd + nstart * tsamp + (ceil(1. * offsets[k] / ndumps[k]) * ndumps[k] - offsets[k]) * tsamps[k];
		if (vm.count("pepoch"))
			folder[k].ref_epoch = MJD(vm["pepoch"].as<double>());
		else
			folder[k].ref_epoch = reader->start_mjd + (ntotal*tsamp/2.);
		folder[k].resize(1, tmpsubdata.nchans, folder[k].nbin);
		folder[k].nblock = nblock;
		folder[k].prepare(tmpsubdata);
	}

	BOOST_LOG_TRIVIAL(info)<<"start folding...";
	
	while (!reader->is_end)
	{
		if (reader->read_data(databuf, ndump) != ndump) break;

		databuf.counter += ndump;

		DataBuffer<float> *data = patch.filter(databuf);
		
		data = prep.run(*data);

		data = equalize.filter(*data);

		data = baseline.filter(*data);

		data = rfi.zap(*data, zaplist);
		if (rfi.isbusy) rfi.closable = false;

		for (auto irfi = rfilist.begin(); irfi!=rfilist.end(); ++irfi)
		{
			if ((*irfi)[0] == "mask")
			{
				data = rfi.mask(*data, threMask, stoi((*irfi)[1]), stoi((*irfi)[2]));
				if (rfi.isbusy) rfi.closable = false;
			}
			else if ((*irfi)[0] == "kadaneF")
			{
				data = rfi.kadaneF(*data, threKadaneF*threKadaneF, widthlimit, stoi((*irfi)[1]), stoi((*irfi)[2]));
				if (rfi.isbusy) rfi.closable = false;
			}
			else if ((*irfi)[0] == "kadaneT")
			{
				data = rfi.kadaneT(*data, threKadaneT*threKadaneT, bandlimitKT, stoi((*irfi)[1]), stoi((*irfi)[2]));
				if (rfi.isbusy) rfi.closable = false;
			}
			else if ((*irfi)[0] == "zdot")
			{
				data = rfi.zdot(*data);
				if (rfi.isbusy) rfi.closable = false;
			}
			else if ((*irfi)[0] == "zero")
			{
				data = rfi.zero(*data);
				if (rfi.isbusy) rfi.closable = false;
			}
		}

		data->closable = true;

		dedisp.prerun(*data);

		for (size_t l=0; l<ndmseg; l++)
		{
			if (vm.count("save_memory"))
				dedisp.alloc(l);
			dedisp.run(l);

#ifndef _OPENMP
			DataBuffer<float> subdata;
			for (auto b=batchs[l].begin(); b!=batchs[l].end(); ++b)
			{
				size_t k = *b;

				dedisp.get_subdata(folder[k].dm, subdata);
				
				if (dedisp.get_counter(folder[k].dm) >= offsets[k] + ndumps[k])
				{
					BOOST_LOG_TRIVIAL(debug)<<"fold cand "<<k;

					if (vm.count("presto"))
						folder[k].runPresto(subdata);
					else if (vm.count("dspsr"))
						folder[k].runDspsr(subdata);
					else
						folder[k].runTRLSM(subdata);	
				}
			}
#else
			std::vector<DataBuffer<float>> subdatas(num_threads);
#pragma omp parallel for num_threads(num_threads)
			for (size_t m=0; m<batchs[l].size(); m++)
			{
				size_t k = batchs[l][m];

				int thread_id = omp_get_thread_num();

				dedisp.get_subdata(folder[k].dm, subdatas[thread_id]);
				
				if (dedisp.get_counter(folder[k].dm) >= offsets[k] + ndumps[k])
				{
					BOOST_LOG_TRIVIAL(debug)<<"fold cand "<<k;

					if (vm.count("presto"))
						folder[k].runPresto(subdatas[thread_id]);
					else if (vm.count("dspsr"))
						folder[k].runDspsr(subdatas[thread_id]);
					else
						folder[k].runTRLSM(subdatas[thread_id]);
				}
			}
#endif
			if (vm.count("save_memory"))
				dedisp.free(l);
		}

		dedisp.postrun(*data);

		databuf.open();
	}
	databuf.close();

	/**
	 * @brief flush the end data
	 * 
	 */

	std::vector<int> nlefts(ncand, 0);
	for (long int k=0; k<ncand; k++)
	{
		nlefts[k] = dedisp.get_offset(folder[k].dm)/dedisp.get_ndump(folder[k].dm);
	}

	int nleft_max = *std::max_element(nlefts.begin(), nlefts.end());

	for (long int l=0; l<nleft_max; l++)
	{
		rfi.open();
		
		dedisp.prerun(rfi);

		DataBuffer<float> subdata;
		
		for (size_t m=0; m<ndmseg; m++)
		{
			if (vm.count("save_memory"))
				dedisp.alloc(m);
			dedisp.run(m);
			for (auto b=batchs[m].begin(); b!=batchs[m].end(); ++b)
			{
				size_t k = *b;

				dedisp.get_subdata(folder[k].dm, subdata);

				if (nlefts[k]-- > 0)
				{
					if (vm.count("presto"))
						folder[k].runPresto(subdata);
					else if (vm.count("dspsr"))
						folder[k].runDspsr(subdata);
					else
						folder[k].runTRLSM(subdata);
				}
			}
			if (vm.count("save_memory"))
				dedisp.free(m);
		}
		dedisp.postrun(rfi);
	}
	
	dedisp.close();

	double tint = ntotal*tsamp;

	vector<Pulsar::GridSearch> gridsearch;
	for (long int k=0; k<ncand; k++)
	{
		Pulsar::GridSearch gs;
		gridsearch.push_back(gs);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int k=0; k<ncand; k++)
	{
		double dm = folder[k].dm;
		double f0 = folder[k].f0;
		double f1 = folder[k].f1;

		if (!nodmsearch)
		{
			gridsearch[k].ddmstart = -scale*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin);
			gridsearch[k].ddmstep = 1./scale*abs(gridsearch[k].ddmstart/folder[k].nbin);
			gridsearch[k].nddm = 2*abs(scale)*folder[k].nbin;
		}
		else
		{
			gridsearch[k].ddmstart = 0.;
			gridsearch[k].ddmstep = 0.;
			gridsearch[k].nddm = 1;
		}

		if (!nof0search)
		{
			gridsearch[k].df0start = -scale*1./tint;
			gridsearch[k].df0step = 1./scale*abs(gridsearch[k].df0start/folder[k].nbin);
			gridsearch[k].ndf0 = 2*abs(scale)*folder[k].nbin;
		}
		else
		{
			gridsearch[k].df0start = 0.;
			gridsearch[k].df0step = 0.;
			gridsearch[k].ndf0 = 1;
		}
		
		if (!nof1search)
		{
			gridsearch[k].df1start = -scale*2./(tint*tint);
			gridsearch[k].df1step = 1./scale*abs(gridsearch[k].df1start/folder[k].nbin);
			gridsearch[k].ndf1 = 2*abs(scale)*folder[k].nbin;
		}
		else
		{
			gridsearch[k].df1start = 0.;
			gridsearch[k].df1step = 0.;
			gridsearch[k].ndf1 = 1;
		}

		gridsearch[k].clfd_q = vm["clfd"].as<double>();
		gridsearch[k].bandcorr = vm["bandcorr"].as<double>();

		gridsearch[k].prepare(folder[k]);
		folder[k].close();

		BOOST_LOG_TRIVIAL(info)<<"dedisperse with the inital dm firstly "<<folder[k].dm;

		gridsearch[k].dm = 0.;
		gridsearch[k].bestddm = folder[k].dm;
		gridsearch[k].dmsearch = true;
		gridsearch[k].correct_delays();
		gridsearch[k].dmsearch = false;

		if (!nosearch)
		{
			BOOST_LOG_TRIVIAL(info)<<"cand "<<k<<": initial dm(pc/cc)="<<gridsearch[k].dm<<", f0(Hz)="<<gridsearch[k].f0<<", f1(Hz/s)="<<gridsearch[k].f1;
			BOOST_LOG_TRIVIAL(info)<<"search scale in phase is "<<scale<<std::endl
			<<"dm search range: delta dm start="<<gridsearch[k].ddmstart<<", dm step="<<gridsearch[k].ddmstep<<", number of dm="<<gridsearch[k].nddm<<std::endl
			<<"f0 search range: delta f0 start="<<gridsearch[k].df0start<<", f0 step="<<gridsearch[k].df0step<<", number of f0="<<gridsearch[k].ndf0<<std::endl
			<<"f1 search range: delta f1 start="<<gridsearch[k].df1start<<", f1 step="<<gridsearch[k].df1step<<", number of f1="<<gridsearch[k].ndf1<<std::endl;

			double dm0 = gridsearch[k].dm;
			double f00 = gridsearch[k].f0;
			double f10 = gridsearch[k].f1;
			double dm1 = gridsearch[k].dm + 2*gridsearch[k].ddmstep;
			double f01 = gridsearch[k].f0 + 2*gridsearch[k].df0step;
			double f11 = gridsearch[k].f1 + 2*gridsearch[k].df1step;
			int cont = 0;
			while ((abs(dm0-dm1)>gridsearch[k].ddmstep or abs(f00-f01)>gridsearch[k].df0step or abs(f10-f11)>gridsearch[k].df1step) and cont<8)
			{
				dm0 = gridsearch[k].dm;
				f00 = gridsearch[k].f0;
				f10 = gridsearch[k].f1;

				if (!nof0search || !nof1search)
				{
					gridsearch[k].runFFdot();
					gridsearch[k].bestprofiles();
				}
				if (!nodmsearch)
				{
					gridsearch[k].runDM();
					gridsearch[k].bestprofiles();
				}

				dm1 = gridsearch[k].dm;
				f01 = gridsearch[k].f0;
				f11 = gridsearch[k].f1;

				cont++;
				BOOST_LOG_TRIVIAL(debug)<<"iteration "<<cont<<": dm(pc/cc)="<<gridsearch[k].dm<<", f0(Hz)="<<gridsearch[k].f0<<", f1(Hz/s)="<<gridsearch[k].f1;
			}
		}

		/**
		 * @brief recalculate the mxsbr_ffdot and vsnr_dm with scale = 3 for plotting
		 * 
		 */

		dm = gridsearch[k].dm;
		f0 = gridsearch[k].f0;
		f1 = gridsearch[k].f1;

		if (!nodmsearch)
		{
			double ddmstart = -3*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin);
			gridsearch[k].ddmstep = 1./3*abs(ddmstart/folder[k].nbin);
			int nddm = 2*3*folder[k].nbin;

			gridsearch[k].ddmstart = std::max(-3*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin), -dm);
			gridsearch[k].nddm = (ddmstart+gridsearch[k].ddmstep*nddm-gridsearch[k].ddmstart)/gridsearch[k].ddmstep;
		}
		else
		{
			gridsearch[k].ddmstart = 0.;
			gridsearch[k].ddmstep = 0.;
			gridsearch[k].nddm = 1;
		}

		if (!nof0search)
		{
			gridsearch[k].df0start = -3*1./tint;
			gridsearch[k].df0step = 1./3*abs(gridsearch[k].df0start/folder[k].nbin);
			gridsearch[k].ndf0 = 2*3*folder[k].nbin;
		}
		else
		{
			gridsearch[k].df0start = 0.;
			gridsearch[k].df0step = 0.;
			gridsearch[k].ndf0 = 1;
		}

		if (!nof1search)
		{
			gridsearch[k].df1start = -3*2./(tint*tint);
			gridsearch[k].df1step = 1./3*abs(gridsearch[k].df1start/folder[k].nbin);
			gridsearch[k].ndf1 = 2*3*folder[k].nbin;
		}
		else
		{
			gridsearch[k].df1start = 0.;
			gridsearch[k].df1step = 0.;
			gridsearch[k].ndf1 = 1.;
		}

		if (!nosearch)
		{
			if (!nof0search || !nof1search)
				gridsearch[k].runFFdot();
			
			if (!nodmsearch)
				gridsearch[k].runDM();
		}
	}

	/** form obsinfo*/
	std::map<std::string, std::string> obsinfo;
	//source name
	obsinfo["Source_name"] = src_name;
	//start mjd
	stringstream ss_mjd;
	ss_mjd << setprecision(10) << fixed << reader->start_mjd.to_day();
	string s_mjd = ss_mjd.str();
	obsinfo["Date"] = s_mjd;
	//ra dec string
	string s_ra = reader->ra;
	string s_dec = reader->dec;
	obsinfo["RA"] = s_ra;
	obsinfo["DEC"] = s_dec;
	//telescope
	obsinfo["Telescope"] = s_telescope;
	//beam
	stringstream ss_ibeam;
	if (vm.count("incoherent"))
		ss_ibeam << "ifbf" << setw(5) << setfill('0') << ibeam;
	else
		ss_ibeam << "cfbf" << setw(5) << setfill('0') << ibeam;
	string s_ibeam = ss_ibeam.str();
	obsinfo["Beam"] = s_ibeam;	
	//data filename
	obsinfo["Filename"] = fnames[reader->idmap[0]];
	//observation length
	obsinfo["Obslen"] = to_string(tint);
	//observation frequency
	obsinfo["Fcentre"] = to_string(0.5*(databuf.frequencies.front()+databuf.frequencies.back()));
	obsinfo["Bandwidth"] = to_string((databuf.frequencies.back()-databuf.frequencies.front())/(databuf.frequencies.size()-1)*databuf.frequencies.size());
	obsinfo["Nchan"] = to_string(databuf.frequencies.size());

	double gl = 0., gb = 0.;
#ifdef HAVE_SOFA
	get_gl_gb(gl, gb, s_ra, s_dec);
#endif
	obsinfo["GL"] = to_string(gl);
	obsinfo["GB"] = to_string(gb);

	double ymw16_maxdm = 0.;
	ymw16_maxdm = get_maxdm_ymw16(gl, gb);
	
	obsinfo["MaxDM_YMW16"] = to_string(ymw16_maxdm);

	//pepoch
	stringstream ss_pepoch;
	ss_pepoch << setprecision(9) << fixed << folder[0].ref_epoch.to_day();
	string s_pepoch = ss_pepoch.str();
	obsinfo["Pepoch"] = s_pepoch;

	stringstream ss_pepoch_in_candfile;
	ss_pepoch_in_candfile << setprecision(13) << fixed << folder[0].ref_epoch.to_day();
	string s_pepoch_in_candfile = ss_pepoch_in_candfile.str();

	BOOST_LOG_TRIVIAL(info)<<"write information to cand file...";

	std::ofstream outfile;
	outfile.open(rootname + "_" + obsinfo["Date"] + "_" + s_ibeam + ".cands");

	outfile<<"#Filename "<<obsinfo["Filename"]<<endl;
	outfile<<"#Telescope "<<obsinfo["Telescope"]<<endl;
	outfile<<"#Source_name "<<obsinfo["Source_name"]<<endl;
	outfile<<"#Beam "<<obsinfo["Beam"]<<endl;
	outfile<<"#Date "<<obsinfo["Date"]<<endl;
	outfile<<"#RA "<<obsinfo["RA"]<<endl;
	outfile<<"#DEC "<<obsinfo["DEC"]<<endl;
	outfile<<"#GL "<<obsinfo["GL"]<<endl;
	outfile<<"#GB "<<obsinfo["GB"]<<endl;
	outfile<<"#MaxDM_YMW16 "<<obsinfo["MaxDM_YMW16"]<<endl;
	outfile<<"#Pepoch "<<s_pepoch_in_candfile<<endl;
	outfile<<"#id       dm_old      dm_new      dm_err		dist_ymw16     f0_old     f0_new        f0_err      f1_old     f1_new       f1_err      acc_old        acc_new      acc_err      S/N        S/N_new"<<endl;

	for (long int k=0; k<ncand; k++)
	{
		stringstream ss_id;
		ss_id << setw(5) << setfill('0') << k+1;
		string s_id = ss_id.str();

		if (!noarch)
		{
			BOOST_LOG_TRIVIAL(info)<<"generate archive for cand "<<k<<"...";

			ArchiveWriter writer;
			writer.template_file = vm["template"].as<string>();
			writer.mode = Integration::FOLD;
			writer.ibeam = 1;
			writer.src_name = src_name;
			writer.ra = s_ra;
			writer.dec = s_dec;
			writer.rootname = rootname + "_" + obsinfo["Date"] + "_" + s_ibeam + "_" + s_id;

			writer.prepare(folder[k], gridsearch[k]);
			writer.run(folder[k], gridsearch[k]);
		}

		double f1_in_candfile = gridsearch[k].f1;

		double c = 1;
		double a = 0.96, b = 1.806;
		
		gridsearch[k].get_snr_width(c);
		gridsearch[k].get_error(obsinfo);

		double ymw16_dist = get_dist_ymw16(gl, gb, gridsearch[k].dm);

		/**
		 * @brief output best and old parameters to bestpar file
		 * 
		 */
		/**
		 * @brief id    dm_old  dm_new  f0_old  f0_new  f1_old  f1_new  acc_old acc_new S/N_old S/N_new
		 * 
		 */

		outfile<<k+1<<"\t\t";
		outfile<<fixed<<setprecision(8)<<folder[k].dm<<"\t\t";
		outfile<<fixed<<setprecision(8)<<gridsearch[k].dm<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_dm<<"\t\t";
		outfile<<fixed<<setprecision(1)<<ymw16_dist<<"\t\t";
		outfile<<setprecision(15)<<folder[k].f0<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].f0<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_f0<<"\t\t";
		outfile<<setprecision(15)<<folder[k].f1<<"\t\t";
		outfile<<setprecision(15)<<f1_in_candfile<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_f1<<"\t\t";
		outfile<<setprecision(15)<<-folder[k].f1/folder[k].f0*CONST_C<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].acc<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_acc<<"\t\t";
		outfile<<fixed<<setprecision(5)<<folder[k].snr<<"\t\t";
		outfile<<fixed<<setprecision(5)<<gridsearch[k].snr<<endl;

		if (!noplot)
		{
			BOOST_LOG_TRIVIAL(info)<<"generate png for cand "<<k<<"...";

			obsinfo["Dist_YMW16"] = to_string(ymw16_dist);
			Pulsar::PulsarPlot psrplot;
			psrplot.plot(folder[k], gridsearch[k], obsinfo, k+1, rootname, true);
		}
	}

	outfile.close();

	BOOST_LOG_TRIVIAL(info)<<"done!";

	return 0;
}

void produce(variables_map &vm, std::list<double> &dmlist, vector<Pulsar::ArchiveLite> &folder)
{
	Pulsar::ArchiveLite fdr;

	/** archive */
	fdr.dm = vm["dm"].as<double>();
	fdr.f0 = vm["f0"].as<double>();
	fdr.f1 = vm["f1"].as<double>();
	fdr.acc = vm["acc"].as<double>();
	fdr.nbin = vm["nbin"].as<int>();

	std::vector<float> vp0;
	std::vector<int> vnbin;
	std::vector<float> nbinplan = vm["nbinplan"].as<std::vector<float>>();
	for (auto pb=nbinplan.begin(); pb!=nbinplan.end(); ++pb)
	{
		vp0.push_back(*pb);
		vnbin.push_back(*(++pb));
	}
	std::vector<size_t> idx = argsort2<float>(vp0);

	dmlist.push_back(vm["dm"].as<double>());

	if (vm.count("candfile"))
	{
		string filename = vm["candfile"].as<string>();
		string line;
		ifstream candfile(filename);
		
		dmlist.clear();
		while (getline(candfile, line))
		{
			boost::trim(line);
			if(!isdigit(line[0])) continue;
			
			vector<string> parameters;
			boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);

			dmlist.push_back(stod(parameters[1]));

			fdr.dm = stod(parameters[1]);
			fdr.acc = stod(parameters[2]);
			fdr.f0 = stod(parameters[3]);
			fdr.f1 = stod(parameters[4]);
			fdr.snr = stod(parameters[5]);

			for (long int k=0; k<vp0.size(); k++)
			{
				if (fdr.f0 <= 1./vp0[idx[k]])
				{
					fdr.nbin = vnbin[idx[k]];
					break;
				}
			}

			folder.push_back(fdr);
		}

	}
	else if (vm.count("parfile"))
	{
		dmlist.clear();

		std::vector<std::string> parfiles = vm["parfile"].as<std::vector<std::string>>();
		for (long int k=0; k<parfiles.size(); k++)
		{
			string filename = parfiles[k];
			string line;
			ifstream parfile(filename);
			int id = 0;
			double f0 = 0.;
			while (getline(parfile, line))
			{
				vector<string> items;
				boost::split(items, line, boost::is_any_of("\t "), boost::token_compress_on);
				if (items[0] == "DM")
				{
					fdr.dm = stod(items[1]);
				}
				else if (items[0] == "F0")
				{
					f0 = stod(items[1]);
				}
			}

			dmlist.push_back(fdr.dm);

			for (long int k=0; k<vp0.size(); k++)
			{
				if (f0 <= 1./vp0[idx[k]])
				{
					fdr.nbin = vnbin[idx[k]];
					break;
				}
			}

			folder.push_back(fdr);
		}
	}
	else
	{
		folder.push_back(fdr);
	}
}
