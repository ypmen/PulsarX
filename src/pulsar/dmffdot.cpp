/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-04-27 16:38:32
 * @modify date 2021-04-27 16:38:32
 * @desc [optimize DM,F0,F1]
 */

#include "config.h"

#include <boost/program_options.hpp>

#include "logging.h"
#include "psrfits.h"
#include "dedisperse.h"
#include "archivelite.h"
#include "gridsearch.h"
#include "pulsarplot.h"
#include "archivewriter.h"
#include "constants.h"

using namespace boost::program_options;

unsigned int num_threads = 1;

int main(int argc, char *argv[])
{
	init_logging();

    /* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("scale", value<int>()->default_value(1), "F0,F1,dm search range scale in phase")
			("scales", value<std::vector<int>>()->multitoken(), "F0,F1,dm search range scales in phase for each archive")
			("nosearch", "Do not search dm,f0,f1")
			("noplot", "Do not generate figures")
			("candfile", value<std::vector<std::string>>()->multitoken()->zero_tokens(), "Input candfile and generate a new one")
			("correct", "Correct archive to the original fold based on candfile")
			("update", "Update archive")
			("template", value<std::string>(), "Input fold template file")
			("srcname", value<std::string>()->default_value("PSRJ0000+00"), "Souce name")
			("telescope", value<std::string>()->default_value("Fake"), "Telescope name")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("ra", value<double>()->default_value(0), "RA (hhmmss.s)")
			("dec", value<double>()->default_value(0), "DEC (ddmmss.s)")
			("clfd", value<double>()->default_value(-1), "CLFD q value, if q<=0, CLFD will not be applied")
#ifdef HAVE_PLOTX
			("plotx", "Using PlotX for plotting")
#endif
			("rootname,o", value<std::string>()->default_value("J0000-00"), "Output rootname")
			("input,f", value<std::vector<std::string>>()->multitoken()->composing(), "Input files");

    command_line_parser parser{argc, argv};
    parser.options(desc).style(command_line_style::default_style | command_line_style::allow_short);
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
		BOOST_LOG_TRIVIAL(error)<<"Error: no input file";
		return -1;
	}
	if (vm.count("correct"))
	{
		if (vm["candfile"].as<std::vector<std::string>>().size() == 0)
		{
			BOOST_LOG_TRIVIAL(error)<<"Error: no cand file";
			return -1;
		}
	}
	if (vm.count("update"))
	{
		if (vm.count("template") == 0)
		{
			BOOST_LOG_TRIVIAL(error)<<"Error: no template file";
			return -1;
		}
	}

	int scale = vm["scale"].as<int>();
	bool nosearch = vm.count("nosearch");
	bool noplot = vm.count("noplot");
	bool noarch = vm.count("noarch");
    string rootname = vm["rootname"].as<std::string>();
	string src_name = vm["srcname"].as<std::string>();
	string s_telescope = vm["telescope"].as<std::string>();
    num_threads = vm["threads"].as<unsigned int>();
    std::vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();
	double src_raj = vm["ra"].as<double>();
	double src_dej = vm["dec"].as<double>();

	Psrfits psf;
	psf.filename = fnames[0];
    psf.open();
    psf.primary.load(psf.fptr);
    psf.load_mode();
    psf.subint.load_header(psf.fptr);
    
    int ibeam = vm["ibeam"].as<int>();
    if (vm["ibeam"].defaulted())
	{
		if (strcmp(psf.primary.ibeam, "") != 0)
			ibeam = stoi(psf.primary.ibeam);
	}

	if (vm["srcname"].defaulted())
	{
		if (strcmp(psf.primary.src_name, "") != 0)
			src_name = psf.primary.src_name;
	}

	if (vm["telescope"].defaulted())
	{
		if (strcmp(psf.primary.telesop, "") != 0)
			s_telescope = psf.primary.telesop;
	}

	if (vm["ra"].defaulted())
	{
		if (strcmp(psf.primary.ra, "") != 0)
		{
			string ra = psf.primary.ra;
			ra.erase(remove(ra.begin(), ra.end(), ':'), ra.end());
			src_raj = stod(ra);
		}
	}
	if (vm["dec"].defaulted())
	{
		if (strcmp(psf.primary.dec, "") != 0)
		{
			string dec = psf.primary.dec;
			dec.erase(remove(dec.begin(), dec.end(), ':'), dec.end());
			src_dej = stod(dec);
		}
	}

    psf.close();

	/** form obsinfo*/
	std::map<std::string, std::string> obsinfo;
	//source name
	obsinfo["Source_name"] = src_name;
	//ra dec string
	string s_ra, s_dec;
	get_s_radec(src_raj, src_dej, s_ra, s_dec);
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
	obsinfo["Nchan"] = to_string(1000000);

	double gl = 0., gb = 0.;
#ifdef HAVE_SOFA
	get_gl_gb(gl, gb, s_ra, s_dec);
#endif
	obsinfo["GL"] = to_string(gl);
	obsinfo["GB"] = to_string(gb);

	double ymw16_maxdm = 0.;
	ymw16_maxdm = get_maxdm_ymw16(gl, gb);

	obsinfo["MaxDM_YMW16"] = to_string(ymw16_maxdm);

	BOOST_LOG_TRIVIAL(info)<<std::endl
	<<"Source name: "<<obsinfo["Source_name"]<<std::endl
	<<"RA: "<<obsinfo["RA"]<<std::endl
	<<"DEC: "<<obsinfo["DEC"]<<std::endl
	<<"Telescope: "<<obsinfo["Telescope"]<<std::endl
	<<"Beam: "<<obsinfo["Beam"];

	std::vector<std::vector<double>> vdmfa_dmffdot;

	std::ofstream outfile;
	if (vm.count("candfile"))
	{
		if (vm["candfile"].as<std::vector<std::string>>().size() != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read header of input candfile...";

			std::ifstream infile(vm["candfile"].as<std::vector<std::string>>()[0]);

			std::vector<std::string> header;
			std::string line;
			while (getline(infile, line))
			{
				boost::trim(line);
				if (line.rfind("#", 0) == 0)
				{
					header.push_back(line);
				}

				if (vm.count("correct") && isdigit(line[0]))
				{
					std::vector<std::string> parameters;
            		boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);

					std::vector<double> dmfa_dmffdot(6);
					dmfa_dmffdot[0] = std::stod(parameters[1]);
					dmfa_dmffdot[1] = std::stod(parameters[5]);
					dmfa_dmffdot[2] = std::stod(parameters[11]);
					dmfa_dmffdot[3] = std::stod(parameters[2]);
					dmfa_dmffdot[4] = std::stod(parameters[6]);
					dmfa_dmffdot[5] = std::stod(parameters[9]);
					vdmfa_dmffdot.push_back(dmfa_dmffdot);
				}
			}

			infile.close();

			BOOST_LOG_TRIVIAL(info)<<"write header of ouput candfile...";
			outfile.open(vm["candfile"].as<std::vector<std::string>>()[0]);
			for (auto line=header.begin(); line!=header.end(); ++line)
			{
				outfile<<*line<<std::endl;
			}
		}
		else
		{
			BOOST_LOG_TRIVIAL(info)<<"write header of ouput candfile...";

			stringstream ss_mjd;
			ss_mjd << setprecision(10) << fixed << psf.primary.start_mjd.to_day();
			string s_mjd = ss_mjd.str();

			outfile.open(rootname + "_" + s_mjd + "_" + s_ibeam + ".cands");

			outfile<<"#Filename "<<"--"<<endl;
			outfile<<"#Telescope "<<obsinfo["Telescope"]<<endl;
			outfile<<"#Source_name "<<obsinfo["Source_name"]<<endl;
			outfile<<"#Beam "<<obsinfo["Beam"]<<endl;
			outfile<<"#Date "<<s_mjd<<endl;
			outfile<<"#RA "<<obsinfo["RA"]<<endl;
			outfile<<"#DEC "<<obsinfo["DEC"]<<endl;
			outfile<<"#GL "<<obsinfo["GL"]<<endl;
			outfile<<"#GB "<<obsinfo["GB"]<<endl;
			outfile<<"#MaxDM_YMW16 "<<obsinfo["MaxDM_YMW16"]<<endl;
			outfile<<"#Pepoch "<<"--"<<endl;
			outfile<<"#id       dm_old      dm_new      dm_err		dist_ymw16     f0_old     f0_new        f0_err      f1_old     f1_new       f1_err      acc_old        acc_new      acc_err      S/N        S/N_new"<<endl;
		}
	}

	long int narch = fnames.size();

	std::vector<int> scales;
	if (vm.count("scales"))
	{
		scales = vm["scales"].as<std::vector<int>>();
		for (long int k=0; k<narch-scales.size(); k++) scales.push_back(scales.back());
	}

	BOOST_LOG_TRIVIAL(info)<<"read "<<narch<<" archives...";

    /* read archive */
	for (long int k=0; k<narch; k++)
	{
		if (!scales.empty()) scale = scales[k];

		std::string fname = fnames[k];

		BOOST_LOG_TRIVIAL(info)<<"read "<<fname;

		Pulsar::ArchiveLite arch;

		arch.read_archive(fname);

		double tint = 0.;
		for (auto it=arch.profiles.begin(); it!=arch.profiles.end(); ++it)
		{
			tint += it->tsubint;
		}

		double dm = arch.dm;
		double f0 = arch.f0;
		double f1 = arch.f1;
		double fcentre = 0.5*(arch.frequencies.front()+arch.frequencies.back());
		double bandwidth = (arch.frequencies.back()-arch.frequencies.front())/(arch.frequencies.size()-1)*arch.frequencies.size();
		double fmax = std::max(fcentre-0.5*bandwidth, fcentre+0.5*bandwidth);
		double fmin = std::min(fcentre-0.5*bandwidth, fcentre+0.5*bandwidth);

		BOOST_LOG_TRIVIAL(info)<<"inital DM (pc/cc)="<<dm<<", f0 (Hz)="<<f0<<" ,f1 (Hz/s)="<<f1;

		/* optimize */
		Pulsar::GridSearch gridsearch;
		gridsearch.ddmstart = -scale*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin);
		gridsearch.ddmstep = 1./scale*abs(gridsearch.ddmstart/arch.nbin);
		gridsearch.nddm = 2*abs(scale)*arch.nbin;
		gridsearch.df0start = -scale*1./tint;
		gridsearch.df0step = 1./scale*abs(gridsearch.df0start/arch.nbin);
		gridsearch.ndf0 = 2*abs(scale)*arch.nbin;
		gridsearch.df1start = -scale*2./(tint*tint);
		gridsearch.df1step = 1./scale*abs(gridsearch.df1start/arch.nbin);
		gridsearch.ndf1 = 2*abs(scale)*arch.nbin;

		gridsearch.clfd_q = vm["clfd"].as<double>();

		gridsearch.prepare(arch);

		if (vm.count("correct"))
		{
			double dm_new = vdmfa_dmffdot[k][3];
			double f0_new = vdmfa_dmffdot[k][4];
			double f1_new = vdmfa_dmffdot[k][5];
			double acc_new = fdot2acc(f1_new, f0_new);

			if (std::abs(dm-dm_new) > 0.001)
			{
				BOOST_LOG_TRIVIAL(error)<<"order in candfile and input archives are not consist";
				return -1;
			}

			double dm_old = vdmfa_dmffdot[k][0];
			double f0_old = vdmfa_dmffdot[k][1];
			double acc_old = vdmfa_dmffdot[k][2];
			double f1_old = acc2fdot(acc_old, f0_old);

			BOOST_LOG_TRIVIAL(info)<<"correct DM/F0/acc to input values "<<std::endl
			<<"from DM (pc/cc)="<<dm_new<<", f0 (Hz)="<<f0_new<<", f1 (Hz/s)="<<f1_new<<" ,acc (m/s/s)="<<acc_new<<std::endl
			<<"to   DM (pc/cc)="<<dm_old<<", f0 (Hz)="<<f0_old<<", f1 (Hz/s)="<<f1_old<<" ,acc (m/s/s)="<<acc_old;

			gridsearch.bestddm = dm_old-dm;
			gridsearch.bestdf0 = f0_old-f0;
			gridsearch.bestdf1 = f1_old-f1;

			gridsearch.dmsearch = true;
			gridsearch.ffdotsearch = true;
			gridsearch.bestprofiles();
			gridsearch.dmsearch = false;
			gridsearch.ffdotsearch = false;
		}

		BOOST_LOG_TRIVIAL(info)<<"optimization...";

		if (!nosearch)
		{
			double dm0 = gridsearch.dm;
			double f00 = gridsearch.f0;
			double f10 = gridsearch.f1;
			double dm1 = gridsearch.dm + 2*gridsearch.ddmstep;
			double f01 = gridsearch.f0 + 2*gridsearch.df0step;
			double f11 = gridsearch.f1 + 2*gridsearch.df1step;
			int cont = 0;
			BOOST_LOG_TRIVIAL(debug)<<gridsearch.dm<<" "<<gridsearch.f0<<" "<<gridsearch.f1;
			while ((abs(dm0-dm1)>gridsearch.ddmstep or abs(f00-f01)>gridsearch.df0step or abs(f10-f11)>gridsearch.df1step) and cont<8)
			{
				dm0 = gridsearch.dm;
				f00 = gridsearch.f0;
				f10 = gridsearch.f1;

				gridsearch.runFFdot();
				gridsearch.bestprofiles();
				gridsearch.runDM();
				gridsearch.bestprofiles();

				dm1 = gridsearch.dm;
				f01 = gridsearch.f0;
				f11 = gridsearch.f1;

				cont++;
				BOOST_LOG_TRIVIAL(debug)<<gridsearch.dm<<" "<<gridsearch.f0<<" "<<gridsearch.f1;
			}
		}

		/**
		 * @brief recalculate the mxsbr_ffdot and vsnr_dm with scale = 3 for plotting
		 * 
		 */

		dm = gridsearch.dm;
		f0 = gridsearch.f0;
		f1 = gridsearch.f1;

		double ddmstart = -3*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin);
		gridsearch.ddmstep = 1./3*abs(ddmstart/arch.nbin);
		int nddm = 2*3*arch.nbin;

		gridsearch.ddmstart = std::max(-3*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin), -dm);
		gridsearch.nddm = (ddmstart+gridsearch.ddmstep*nddm-gridsearch.ddmstart)/gridsearch.ddmstep;

		gridsearch.df0start = -3*1./tint;
		gridsearch.df0step = 1./3*abs(gridsearch.df0start/arch.nbin);
		gridsearch.ndf0 = 2*3*arch.nbin;
		gridsearch.df1start = -3*2./(tint*tint);
		gridsearch.df1step = 1./3*abs(gridsearch.df1start/arch.nbin);
		gridsearch.ndf1 = 2*3*arch.nbin;

		if (!nosearch)
		{
			gridsearch.runFFdot();
			gridsearch.runDM();
		}

		//start mjd
		stringstream ss_mjd;
		ss_mjd << setprecision(10) << fixed << arch.start_mjd.to_day();
		string s_mjd = ss_mjd.str();
		obsinfo["Date"] = s_mjd;

		//data filename
		obsinfo["Filename"] = fname;
		//observation length
		obsinfo["Obslen"] = to_string(tint);
		//frequency
		obsinfo["Fcentre"] = to_string(fcentre);
		obsinfo["Bandwidth"] = to_string(bandwidth);

		double pepoch_offset = arch.profiles[arch.profiles.size()/2].offs_sub;

		stringstream ss_pepoch;
		ss_pepoch << setprecision(9) << fixed << (arch.start_mjd+pepoch_offset).to_day();
		string s_pepoch = ss_pepoch.str();
		obsinfo["Pepoch"] = s_pepoch;

		if (vm.count("update"))
		{
			BOOST_LOG_TRIVIAL(info)<<"update archive...";
			size_t lastindex = fname.find_last_of("."); 
			std::string basename = fname.substr(0, lastindex);

			ArchiveWriter writer;
			writer.template_file = vm["template"].as<string>();
			writer.mode = Integration::FOLD;
			writer.ibeam = 1;
			writer.src_name = src_name;
			writer.ra = s_ra;
			writer.dec = s_dec;
			writer.rootname = basename;

			writer.prepare(arch, gridsearch);
			writer.run(arch, gridsearch);
		}

		gridsearch.get_snr_width();
		gridsearch.get_error(obsinfo);

		double ymw16_dist = get_dist_ymw16(gl, gb, gridsearch.dm);

		if (vm.count("candfile"))
		{
			outfile<<k+1<<"\t\t";
			if (vm.count("correct"))
				outfile<<fixed<<setprecision(8)<<vdmfa_dmffdot[k][0]<<"\t\t";
			else
				outfile<<fixed<<setprecision(8)<<arch.dm<<"\t\t";
			outfile<<fixed<<setprecision(8)<<gridsearch.dm<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.err_dm<<"\t\t";
			outfile<<fixed<<setprecision(1)<<ymw16_dist<<"\t\t";
			if (vm.count("correct"))
				outfile<<setprecision(15)<<vdmfa_dmffdot[k][1]<<"\t\t";
			else
				outfile<<setprecision(15)<<arch.f0<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.f0<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.err_f0<<"\t\t";
			if (vm.count("correct"))
				outfile<<setprecision(15)<<acc2fdot(vdmfa_dmffdot[k][2], vdmfa_dmffdot[k][1])<<"\t\t";
			else
				outfile<<setprecision(15)<<arch.f1<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.f1<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.err_f1<<"\t\t";
			if (vm.count("correct"))
				outfile<<setprecision(15)<<vdmfa_dmffdot[k][2]<<"\t\t";
			else
				outfile<<setprecision(15)<<fdot2acc(arch.f1, arch.f0)<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.acc<<"\t\t";
			outfile<<setprecision(15)<<gridsearch.err_acc<<"\t\t";
			outfile<<fixed<<setprecision(5)<<arch.snr<<"\t\t";
			outfile<<fixed<<setprecision(5)<<gridsearch.snr<<endl;
		}

		if (!noplot)
		{
			BOOST_LOG_TRIVIAL(info)<<"plotting...";
			obsinfo["Dist_YMW16"] = to_string(ymw16_dist);

			size_t lastindex = fname.find_last_of("."); 
			std::string basename = fname.substr(0, lastindex);
			obsinfo["Basename"] = basename;

			Pulsar::PulsarPlot psrplot;
			psrplot.plot(arch, gridsearch, obsinfo, 1, rootname, vm.count("plotx"));
		}
	}

	outfile.close();

	BOOST_LOG_TRIVIAL(info)<<"done!";

    return 0;
}