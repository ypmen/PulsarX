/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-04-27 16:38:32
 * @modify date 2021-04-27 16:38:32
 * @desc [optimize DM,F0,F1]
 */

#include "config.h"

#include <boost/program_options.hpp>

#include "psrfits.h"
#include "dedisperse.h"
#include "archivelite.h"
#include "gridsearch.h"
#include "pulsarplot.h"

using namespace boost::program_options;

unsigned int num_threads = 1;

int main(int argc, char *argv[])
{
    /* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("scale", value<int>()->default_value(1), "F0,F1,dm search range scale in phase")
			("nosearch", "Do not search dm,f0,f1")
			("noplot", "Do not generate figures")
			("srcname", value<std::string>()->default_value("PSRJ0000+00"), "Souce name")
			("telescope", value<std::string>()->default_value("Fake"), "Telescope name")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("ra", value<double>()->default_value(0), "RA (hhmmss.s)")
			("dec", value<double>()->default_value(0), "DEC (ddmmss.s)")
			("clfd", value<double>()->default_value(-1), "CLFD q value, if q<=0, CLFD will not be applied")
			("rfi,z", value<vector<std::string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [zap fl fh]")
#ifdef HAVE_PLOTX
			("plotx", "Using PlotX for plotting")
#endif
			("rootname,o", value<std::string>()->default_value("J0000-00"), "Output rootname")
			("input,f", value<std::string>()->multitoken()->composing(), "Input files");

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

	int scale = vm["scale"].as<int>();
	bool nosearch = vm.count("nosearch");
	bool noplot = vm.count("noplot");
	bool noarch = vm.count("noarch");
    string rootname = vm["rootname"].as<std::string>();
	string src_name = vm["srcname"].as<std::string>();
	string s_telescope = vm["telescope"].as<std::string>();
    num_threads = vm["threads"].as<unsigned int>();
    std::string fname = vm["input"].as<std::string>();
	double src_raj = vm["ra"].as<double>();
	double src_dej = vm["dec"].as<double>();

	Psrfits psf;
	psf.filename = fname;
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

    /* read archive */
    Pulsar::ArchiveLite arch;
    arch.read_archive(fname);

    double tint = 0;
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

    if (!nosearch)
    {
        double dm0 = gridsearch.dm;
        double f00 = gridsearch.f0;
        double f10 = gridsearch.f1;
        double dm1 = gridsearch.dm + 2*gridsearch.ddmstep;
        double f01 = gridsearch.f0 + 2*gridsearch.df0step;
        double f11 = gridsearch.f1 + 2*gridsearch.df1step;
        int cont = 0;
        //cout<<gridsearch.dm<<" "<<gridsearch.f0<<" "<<gridsearch.f1<<endl;
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
            //cout<<gridsearch.dm<<" "<<gridsearch.f0<<" "<<gridsearch.f1<<endl;
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


    /** form obsinfo*/
	std::map<std::string, std::string> obsinfo;
	//source name
	obsinfo["Source_name"] = src_name;
	//start mjd
	stringstream ss_mjd;
    ss_mjd << setprecision(10) << fixed << arch.start_mjd.to_day();
    string s_mjd = ss_mjd.str();
	obsinfo["Date"] = s_mjd;
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
	//data filename
	obsinfo["Filename"] = fname;
	//observation length
	obsinfo["Obslen"] = to_string(tint);
    //frequency
    obsinfo["Fcentre"] = to_string(fcentre);
    obsinfo["Bandwidth"] = to_string(bandwidth);
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

	double pepoch_offset = arch.profiles[arch.profiles.size()/2].offs_sub;
	
	//pepoch
	stringstream ss_pepoch;
    ss_pepoch << setprecision(9) << fixed << (arch.start_mjd+pepoch_offset).to_day();
    string s_pepoch = ss_pepoch.str();
	obsinfo["Pepoch"] = s_pepoch;

    gridsearch.get_snr_width();
	gridsearch.get_error(obsinfo);

    double ymw16_dist = get_dist_ymw16(gl, gb, gridsearch.dm);

    if (!noplot)
    {
        obsinfo["Dist_YMW16"] = to_string(ymw16_dist);

		size_t lastindex = fname.find_last_of("."); 
		std::string basename = fname.substr(0, lastindex);
		obsinfo["Basename"] = basename;

        Pulsar::PulsarPlot psrplot;
        psrplot.plot(arch, gridsearch, obsinfo, 1, rootname, vm.count("plotx"));
    }

    return 0;
}