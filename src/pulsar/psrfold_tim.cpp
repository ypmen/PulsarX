/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2025-02-13 09:45:30
 * @modify date 2025-02-13 09:45:30
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
#include "presto.h"
#include "databuffer.h"
#include "mjd.h"
#include "utils.h"
#include "constants.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
bool dumptim=false;

void produce(variables_map &vm, vector<Pulsar::ArchiveLite> &folder);

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
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("frac", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 1}, "0, 1"), "Reading data between start and end fraction")
			("f0", value<double>()->default_value(0), "F0 (Hz)")
			("f1", value<double>()->default_value(0), "F1 (Hz/s)")
			("f2", value<double>()->default_value(0), "F2 (Hz/s/s)")
			("acc", value<double>()->default_value(0), "Acceleration (m/s/s)")
			("binary", "Fold binary pulsar")
			("Pb", value<double>(), "The orbital period (day)")
			("A1", value<double>(), "The projected orbital semi-major axis (lt-sec)")
			("T0", value<double>(), "The time of periastron passage (day)")
			("ECC", value<double>()->default_value(0.), "The orbital eccentricity")
			("OM", value<double>()->default_value(0.), "Longitude of periastron (deg)")
			("pepoch", value<double>(), "F0/F1/F2/acc epoch (MJD)")
			("scale", value<int>()->default_value(1), "F0,F1,F2,dm search range scale in phase")
			("nosearch", "Do not search dm,f0,f1,f2")
			("nof0search", "Do not search f0")
			("nof1search", "Do not search f1")
			("f2search", "search f2")
			("noplot", "Do not generate figures")
			("noarch", "Do not generate archives")
			("parfile", value<vector<string>>(), "Input pulsar par files")
			("t2pred", value<vector<string>>(), "Input t2pred files")
			("candfile", value<string>(), "Input cand file")
			("template", value<string>(), "Input fold template file")
			("nbin,b", value<int>()->default_value(32), "Number of bins per period")
			("nbinplan", value<vector<float>>()->multitoken()->default_value(vector<float>{0.1, 128}, "0.1, 128"), "Using nbins of nbin,nbin1,nbin2,... for 0,p1,p2,... (s)")
			("tsubint,L", value<double>()->default_value(10), "Time length per integration (s)")
			("blocksize", value<int>()->default_value(2), "Size of data block (s)")
			("srcname", value<string>()->default_value("PSRJ0000+00"), "Souce name")
			("telescope", value<string>()->default_value("Fake"), "Telescope name")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("ra", value<string>()->default_value("00:00:00"), "RA (hhmmss.s)")
			("dec", value<string>()->default_value("00:00:00"), "DEC (ddmmss.s)")
			("cdm", value<double>()->default_value(0), "Coherent dedispersion DM (pc/cc)")
			("clfd", value<double>()->default_value(-1), "CLFD q value, if q<=0, CLFD will not be applied")
			("bandcorr", value<double>()->default_value(10), "CLFD correlation bandwidth threshold to differentiate RFI and real signal (MHz)")
			("render", "Using new folding algorithm")
			("dspsr", "Using dspsr folding algorithm")
			("presto", "Using presto folding algorithm (default)")
			("plotx", "Using PlotX for plotting (default; not used any more)")
			("output_width", "Output boxcar width (s)")
			("saveimage", "Save images to fits")
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("input,f", value<string>(), "Input timeseries file");

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
	bool nodmsearch = true;
	bool nof0search = vm.count("nof0search");
	bool nof1search = vm.count("nof1search");
	bool f2search = vm.count("f2search");
	if (nodmsearch && nof0search && nof1search && (!f2search)) nosearch = true;
	if (nosearch)
	{
		nodmsearch = true;
		nof0search = true;
		nof1search = true;
		f2search = false;
	}
	bool noplot = vm.count("noplot");
	bool noarch = vm.count("noarch");
	string rootname = vm["rootname"].as<string>();
	num_threads = vm["threads"].as<unsigned int>();
	vector<double> jump = vm["jump"].as<vector<double>>();

	std::string fname = vm["input"].as<string>();
	std::string fname_root = fname.substr(0, fname.find_last_of("."));

	PRESTO::Info reader;
	reader.verbose = verbose;
	reader.open(fname);
	reader.read(fname_root + ".inf");

	if (!vm["ibeam"].defaulted())
		reader.beam = std::to_string(vm["ibeam"].as<int>());
	if (!vm["telescope"].defaulted())
		reader.telescope = vm["telescope"].as<std::string>();
	if (!vm["ra"].defaulted())
		reader.ra = vm["ra"].as<std::string>();
	if (!vm["dec"].defaulted())
		reader.dec = vm["dec"].as<std::string>();
	if (!vm["srcname"].defaulted())
		reader.source_name = vm["srcname"].as<std::string>();

	string src_name = reader.source_name;
	string s_telescope = reader.telescope;
	int ibeam = reader.beam.empty() ? 0 : std::stoi(reader.beam);

	double tsamp = reader.tsamp;
	int nifs = 1;
	long int ntotal = reader.nsamples;

	vector<Pulsar::ArchiveLite> folder;

	produce(vm, folder);

	int blocksize = vm["blocksize"].as<int>();
	long int ndump = ceil((blocksize/tsamp));
	int nblock = ceil(vm["tsubint"].as<double>())/blocksize;

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	double startfrac = vm["frac"].as<vector<double>>()[0];
	double endfrac = vm["frac"].as<vector<double>>()[1];
	if (startfrac != 0 or endfrac != 1)
	{
		nstart = startfrac*ntotal;
		nend = endfrac*ntotal;
	}

	reader.skip_start = nstart;
	reader.skip_end = ntotal-nend;
	reader.skip_head();

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

			double mjd_start = (reader.epoch)-0.01;
			double mjd_end = (reader.epoch)+1;
			double fmin = reader.fch1;
			double fmax = reader.fch1 + (reader.nchans-1) * reader.foff;
			fmin -= 0.5*reader.foff;
			fmax += 0.5*reader.foff;

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

	long int ncand = folder.size();
	if (ncand <= 0)
	{
		BOOST_LOG_TRIVIAL(warning)<<"no candidate to fold"<<std::endl;
		exit(0);
	}

	DataBuffer<float> subdata;

	for (long int k=0; k<ncand; k++)
	{
		subdata.resize(ndump, 1);
		subdata.tsamp = tsamp;
		subdata.frequencies = std::vector<double>{reader.fch1+(reader.nchans-1)*reader.foff};

		folder[k].start_mjd = MJD(reader.epoch) + nstart * tsamp;
		if (vm.count("pepoch"))
			folder[k].ref_epoch = MJD(vm["pepoch"].as<double>());
		else
			folder[k].ref_epoch = MJD(reader.epoch) + (ntotal*tsamp/2.);
		folder[k].resize(1, subdata.nchans, folder[k].nbin);
		folder[k].nblock = nblock;
		folder[k].prepare(subdata);
		folder[k].dedispersed = true;
	}

	BOOST_LOG_TRIVIAL(info)<<"start folding...";
	
	while (!reader.is_end)
	{
		if (reader.read_data(subdata.buffer, ndump) != ndump) break;

		subdata.counter += ndump;

#pragma omp parallel for num_threads(num_threads)
		for (int k=0; k<ncand; k++)
		{
			BOOST_LOG_TRIVIAL(debug)<<"fold cand "<<k;

			if (vm.count("render"))
				folder[k].runTRLSM(subdata);
			else if (vm.count("dspsr"))
				folder[k].runDspsr(subdata);
			else
				folder[k].runPresto(subdata);
		}
	}

	reader.close();

	/**
	 * @brief flush the end data
	 * 
	 */

	double tint = ntotal*tsamp;

	vector<Pulsar::GridSearch> gridsearch;
	for (long int k=0; k<ncand; k++)
	{
		Pulsar::GridSearch gs;
		gs.scale = scale;
		gs.tint = tint;
		gs.nodmsearch = nodmsearch;
		gs.nof0search = nof0search;
		gs.nof1search = nof1search;
		gs.f2search = f2search;
		gs.nosearch = nosearch;
		gs.clfd_q = vm["clfd"].as<double>();
		gs.bandcorr = vm["bandcorr"].as<double>();
		gridsearch.push_back(gs);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (long int k=0; k<ncand; k++)
	{
		gridsearch[k].prepare(folder[k]);

		folder[k].close();

		gridsearch[k].run(k);
	}

	/** form obsinfo*/
	std::map<std::string, std::string> obsinfo;
	//source name
	obsinfo["Source_name"] = src_name;
	//start mjd
	stringstream ss_mjd;
	ss_mjd << setprecision(10) << fixed << reader.epoch;
	string s_mjd = ss_mjd.str();
	obsinfo["Date"] = s_mjd;
	//ra dec string
	string s_ra = reader.ra;
	string s_dec = reader.dec;
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
	obsinfo["Filename"] = fname;
	//observation length
	obsinfo["Obslen"] = to_string(tint);
	//observation frequency
	obsinfo["Fcentre"] = to_string(subdata.frequencies[0]);
	obsinfo["Bandwidth"] = to_string(reader.foff);
	obsinfo["Nchan"] = to_string(1);

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

	outfile<<"#id	dm_old	dm_new	dm_err	dist_ymw16	f0_old	f0_new	f0_err	f1_old	f1_new	f1_err";

	if (f2search)
		outfile<<" f2_old	f2_new	f2_err";

	outfile<<" acc_old	acc_new	acc_err	S/N	S/N_new";

	if (vm.count("output_width"))
		outfile<<" boxcar_width";
	
	outfile<<std::endl;

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
		double f2_in_candfile = gridsearch[k].f2;

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
		if (f2search)
		{
			outfile<<setprecision(15)<<folder[k].f2<<"\t\t";
			outfile<<setprecision(15)<<f2_in_candfile<<"\t\t";
			outfile<<setprecision(15)<<gridsearch[k].err_f2<<"\t\t";
		}
		outfile<<setprecision(15)<<-folder[k].f1/folder[k].f0*CONST_C<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].acc<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_acc<<"\t\t";
		outfile<<fixed<<setprecision(5)<<folder[k].snr<<"\t\t";
		outfile<<fixed<<setprecision(5)<<gridsearch[k].snr;
		if (vm.count("output_width"))
			outfile<<"\t\t"<<fixed<<setprecision(5)<<gridsearch[k].width<<endl;

		outfile<<std::endl;

		if (!noplot)
		{
			BOOST_LOG_TRIVIAL(info)<<"generate png for cand "<<k<<"...";

			obsinfo["Dist_YMW16"] = to_string(ymw16_dist);
			obsinfo["Coherent_DM"] = to_string(0);
			Pulsar::PulsarPlot psrplot;
			psrplot.plot(folder[k], gridsearch[k], obsinfo, k+1, rootname, true, vm.count("saveimage"));
		}
	}

	outfile.close();

	BOOST_LOG_TRIVIAL(info)<<"done!";

	return 0;
}

void produce(variables_map &vm, vector<Pulsar::ArchiveLite> &folder)
{
	Pulsar::ArchiveLite fdr;

	/** archive */
	fdr.f0 = vm["f0"].as<double>();
	fdr.f1 = vm["f1"].as<double>();
	fdr.f2 = vm["f2"].as<double>();
	fdr.acc = vm["acc"].as<double>();
	fdr.nbin = vm["nbin"].as<int>();

	if (vm.count("binary"))
	{
		std::vector<double> kepler_params;
		kepler_params.push_back(vm["f0"].as<double>());
		kepler_params.push_back(vm["Pb"].as<double>());
		kepler_params.push_back(vm["A1"].as<double>());
		kepler_params.push_back(vm["T0"].as<double>());
		kepler_params.push_back(vm["OM"].as<double>()/180.*M_PI);
		kepler_params.push_back(vm["ECC"].as<double>());

		fdr.orb = Pulsar::Kepler(kepler_params);
		fdr.use_kepler = true;
	}

	std::vector<float> vp0;
	std::vector<int> vnbin;
	std::vector<float> nbinplan = vm["nbinplan"].as<std::vector<float>>();
	for (auto pb=nbinplan.begin(); pb!=nbinplan.end(); ++pb)
	{
		vp0.push_back(*pb);
		vnbin.push_back(*(++pb));
	}
	std::vector<size_t> idx = argsort2<float>(vp0);

	if (vm.count("candfile"))
	{
		string filename = vm["candfile"].as<string>();
		string line;
		ifstream candfile(filename);
		
		while (getline(candfile, line))
		{
			boost::trim(line);
			if(!isdigit(line[0])) continue;
			
			vector<string> parameters;
			boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);

			if (parameters.size() == 6)
			{
				fdr.dm = stod(parameters[1]);
				fdr.acc = stod(parameters[2]);
				fdr.f0 = stod(parameters[3]);
				fdr.f1 = stod(parameters[4]);
				fdr.snr = stod(parameters[5]);
			}
			else if (parameters.size() == 7)
			{
				fdr.dm = stod(parameters[1]);
				fdr.acc = stod(parameters[2]);
				fdr.f0 = stod(parameters[3]);
				fdr.f1 = stod(parameters[4]);
				fdr.f2 = stod(parameters[5]);
				fdr.snr = stod(parameters[6]);
			}
			else if (parameters.size() == 12)
			{
				fdr.dm = stod(parameters[1]);
				fdr.acc = stod(parameters[2]);
				fdr.f1 = stod(parameters[4]);
				fdr.f2 = stod(parameters[5]);
				
				std::vector<double> kepler_params{
					stod(parameters[3]),
					stod(parameters[6]),
					stod(parameters[7]),
					stod(parameters[8]),
					stod(parameters[9]),
					stod(parameters[10])
				};

				fdr.orb = Pulsar::Kepler(kepler_params);
				fdr.use_kepler = true;

				fdr.snr = stod(parameters[11]);
			}
			
			fdr.nbin = vm["nbin"].as<int>();

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

			fdr.nbin = vm["nbin"].as<int>();

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
