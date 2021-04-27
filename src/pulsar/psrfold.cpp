/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 15:45:25
 * @modify date 2020-06-07 15:45:25
 * @desc [fold candidates]
 */

#define FAST 1
#define GROUPSIZE 64

#include "config.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <utility>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp> 

#include "dedisperse.h"
#include "archivewriter.h"

#include "pulsarplot.h"

#include "gridsearch.h"
#include "archivelite.h"
#include "dedispersionliteU.h"
#include "databuffer.h"
#include "downsample.h"
#include "rfi.h"
#include "equalize.h"
#include "psrfits.h"
#include "mjd.h"
#include "utils.h"
#include "constants.h"
#include "preprocess.h"
#include "baseline.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
bool dumptim=false;

void produce(variables_map &vm, std::vector<std::vector<double>> &dmsegs, vector<Pulsar::ArchiveLite> &folder);

int main(int argc, const char *argv[])
{
    /* options */
	int verbose = 0;
	bool outbest = false;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{0, 0}, "0, 0"), "Time jump at the beginning and end (s)")
			("td", value<int>()->default_value(1), "Time downsample")
			("fd", value<int>()->default_value(1), "Frequency downsample")
			("zapthre", value<float>()->default_value(4), "Threshold in IQR for zapping channels")
			("dm", value<double>()->default_value(0), "DM (pc/cc)")
            ("f0", value<double>()->default_value(0), "F0 (Hz)")
            ("f1", value<double>()->default_value(0), "F1 (Hz/s)")
			("acc", value<double>()->default_value(0), "Acceleration (m/s/s)")
			("scale", value<int>()->default_value(1), "F0,F1,dm search range scale in phase")
			("nosearch", "Do not search dm,f0,f1")
			("noplot", "Do not generate figures")
			("noarch", "Do not generate archives")
			("candfile", value<string>(), "Input cand file")
			("template", value<string>(), "Input fold template file")
			("nbin,b", value<int>()->default_value(32), "Number of bins per period")
			("nbinplan", value<vector<float>>()->multitoken()->default_value(vector<float>{0.1, 128}, "0.1, 128"), "Using nbins of nbin,nbin1,nbin2,... for 0,p1,p2,... (s)")
			("tsubint,L", value<double>()->default_value(1), "Time length per integration (s)")
			("nsubband,n", value<int>()->default_value(32), "Number of subband")
			("srcname", value<string>()->default_value("PSRJ0000+00"), "Souce name")
			("telescope", value<string>()->default_value("Fake"), "Telescope name")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("ra", value<double>()->default_value(0), "RA (hhmmss.s)")
			("dec", value<double>()->default_value(0), "DEC (ddmmss.s)")
			("baseline", value<vector<float>>()->multitoken()->default_value(vector<float>{0.0, 0.0}, "0.0, 0.0"), "The scale of baseline remove (s)")
			("clfd", value<double>()->default_value(-1), "CLFD q value, if q<=0, CLFD will not be applied")
			("rfi,z", value<vector<string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("tdRFI", value<int>()->default_value(1), "Time downsample of RFI")
			("fdRFI", value<int>()->default_value(1), "Frequency downsample of RFI")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("threMask", value<float>()->default_value(10), "S/N threshold of Mask")
            ("render", "Using new folding algorithm (deprecated, used by default)")
			("dspsr", "Using dspsr folding algorithm")
#ifdef HAVE_PLOTX
			("plotx", "Using PlotX for plotting")
#endif
			("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("cont", "Input files are contiguous")
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
		cerr<<"Error: no input file"<<endl;
		return -1;
	}

	int scale = vm["scale"].as<int>();
	bool nosearch = vm.count("nosearch");
	bool noplot = vm.count("noplot");
	bool noarch = vm.count("noarch");
	bool contiguous = vm.count("cont");
    string rootname = vm["rootname"].as<string>();
	string src_name = vm["srcname"].as<string>();
	string s_telescope = vm["telescope"].as<string>();
    num_threads = vm["threads"].as<unsigned int>();
    vector<double> jump = vm["jump"].as<vector<double>>();
    vector<string> fnames = vm["input"].as<vector<string>>();
	double src_raj = vm["ra"].as<double>();
	double src_dej = vm["dec"].as<double>();

	long int npsf = fnames.size();
	Psrfits *psf = new Psrfits [npsf];
	for (long int i=0; i<npsf; i++)
	{
		psf[i].filename = fnames[i];
	}

	vector<MJD> tstarts;
	vector<MJD> tends;
	long int ntotal = 0;
	for (long int i=0; i<npsf; i++)
	{
		psf[i].open();
		psf[i].primary.load(psf[i].fptr);
		psf[i].load_mode();
		psf[i].subint.load_header(psf[i].fptr);
		ntotal += psf[i].subint.nsamples;
		tstarts.push_back(psf[i].primary.start_mjd);
		tends.push_back(psf[i].primary.start_mjd+psf[i].subint.nsamples*psf[i].subint.tbin);
		psf[i].close();
	}
	vector<size_t> idx = argsort(tstarts);
	for (long int i=0; i<npsf-1; i++)
	{
		if (abs((tends[idx[i]]-tstarts[idx[i+1]]).to_second())>0.5*psf[idx[i]].subint.tbin)
		{
			if (contiguous)
			{
				cerr<<"Warning: time not contiguous"<<endl;
			}
			else
			{
				cerr<<"Error: time not contiguous"<<endl;
				exit(-1);
			}
		}
	}

    psf[0].open();
	psf[0].primary.load(psf[0].fptr);
	psf[0].load_mode();
	psf[0].subint.load_header(psf[0].fptr);

	if (psf[0].mode != Integration::SEARCH)
	{
		cerr<<"Error: mode is not SEARCH"<<endl;
		exit(-1);
	}

	int ibeam = vm["ibeam"].as<int>();

	if (vm["ibeam"].defaulted())
	{
		if (strcmp(psf[0].primary.ibeam, "") != 0)
			ibeam = stoi(psf[0].primary.ibeam);
	}

	if (vm["srcname"].defaulted())
	{
		if (strcmp(psf[0].primary.src_name, "") != 0)
			src_name = psf[0].primary.src_name;
	}

	if (vm["telescope"].defaulted())
	{
		if (strcmp(psf[0].primary.telesop, "") != 0)
			s_telescope = psf[0].primary.telesop;
	}

	if (vm["ra"].defaulted())
	{
		if (strcmp(psf[0].primary.ra, "") != 0)
		{
			string ra = psf[0].primary.ra;
			ra.erase(remove(ra.begin(), ra.end(), ':'), ra.end());
			src_raj = stod(ra);
		}
	}
	if (vm["dec"].defaulted())
	{
		if (strcmp(psf[0].primary.dec, "") != 0)
		{
			string dec = psf[0].primary.dec;
			dec.erase(remove(dec.begin(), dec.end(), ':'), dec.end());
			src_dej = stod(dec);
		}
	}

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

	Integration it;
	psf[0].subint.load_integration(psf[0].fptr, 0, it);

	long int nchans = it.nchan;
    double tsamp = psf[0].subint.tbin;
    int nifs = it.npol;

	short *buffer = new short [nchans];

	long int ndump = ceil((1./tsamp)/td)*td;
	int nblock = ceil(vm["tsubint"].as<double>());

	DataBuffer<short> databuf(ndump, nchans);
	databuf.closable = true;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], it.frequencies, sizeof(double)*nchans);

	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	Preprocess prep;
	prep.td = vm["td"].as<int>();
	prep.fd = vm["fd"].as<int>();
	prep.thresig = vm["zapthre"].as<float>();
	prep.width = vm["baseline"].as<vector<float>>().front();
	prep.prepare(databuf);

	Downsample downsample;
	downsample.td = td;
    downsample.fd = fd;
    downsample.prepare(prep);
	downsample.close();
	downsample.closable = true;

	Equalize equalize;
	equalize.prepare(downsample);
	equalize.close();
	equalize.closable = true;

	BaseLine baseline;
	baseline.width = vm["baseline"].as<vector<float>>().back();
	baseline.prepare(equalize);
	baseline.close();
	baseline.closable = true;

    RFI rfi;
	rfi.prepare(baseline);
	rfi.close();
	rfi.closable = true;
    
    Pulsar::DedispersionLiteU dedisp;
	vector<Pulsar::ArchiveLite> folder;

	std::vector<std::vector<double>> dmsegs;
	produce(vm, dmsegs, folder);
	dedisp.nsubband = vm["nsubband"].as<int>();
	double maxdm = 0;
	for (auto dmseg=dmsegs.begin(); dmseg!=dmsegs.end(); ++dmseg)
	{
		for (auto dm=(*dmseg).begin(); dm!=(*dmseg).end(); ++dm)
		{
			if (*dm > maxdm) maxdm = *dm;
		}
	}
	dedisp.maxdm = maxdm;
	dedisp.groupsize = GROUPSIZE;
    dedisp.prepare(rfi);

    DataBuffer<float> subdata;
    dedisp.get_subdata(subdata, 0);

	long int ncand = folder.size();

    for (long int k=0; k<ncand; k++)
	{
		folder[k].start_mjd = tstarts[idx[0]]+(ceil(1.*dedisp.offset/ndump)*ndump-dedisp.offset)*tsamp*td;
		folder[k].ref_epoch = tstarts[idx[0]]+(ntotal*tsamp/2.);
        folder[k].resize(1, subdata.nchans, folder[k].nbin);
		folder[k].nblock = nblock;
		folder[k].prepare(subdata);
	}

    psf[0].close();

    int sumif = nifs>2? 2:nifs;
	
	long int ntot = 0;
	long int ntot2 = 0;
	long int count = 0;
    long int bcnt1 = 0;
	for (long int idxn=0; idxn<npsf; idxn++)
	{
		long int n = idx[idxn];

		psf[n].open();
		psf[n].primary.load(psf[n].fptr);
		psf[n].load_mode();
		psf[n].subint.load_header(psf[n].fptr);

		for (long int s=0; s<psf[n].subint.nsubint; s++)
		{
			if (verbose)
			{
				cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
				cerr<<"("<<100.*count/ntotal<<"%)";
			}

			psf[n].subint.load_integration_data(psf[n].fptr, s, it);
#ifdef FAST
			unsigned char *pcur = (unsigned char *)(it.data);
#endif
			for (long int i=0; i<it.nsblk; i++)
			{
				count++;
				if (count-1<nstart or count-1>nend)
				{
					pcur += it.npol*it.nchan;
					continue;
				}

				memset(buffer, 0, sizeof(short)*nchans);
				long int m = 0;
				for (long int k=0; k<sumif; k++)
				{
					for (long int j=0; j<nchans; j++)
					{
						buffer[j] +=  pcur[m++];
					}
				}

                memcpy(&databuf.buffer[0]+bcnt1*nchans, buffer, sizeof(short)*1*nchans);
				databuf.counter++;
                bcnt1++;
				ntot++;

				if (ntot%ndump == 0)
				{
					DataBuffer<float> *data = prep.run(databuf);

    				data = downsample.run(*data);

    				data = equalize.run(*data);

					data = baseline.run(*data);

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

					for (long int k=0; k<ncand; k++)
					{
						if (k%GROUPSIZE == 0)
						{
							dedisp.updatedm(dmsegs[k/GROUPSIZE]);
							dedisp.run();
						}
						dedisp.get_subdata(subdata, k%GROUPSIZE);
                        
						if (dedisp.counter >= dedisp.offset+dedisp.ndump)
						{
							if (vm.count("dspsr"))
								folder[k].runDspsr(subdata);
							else
								folder[k].runTRLSM(subdata);				
						}
					}

					dedisp.postrun();

                    bcnt1 = 0;
					databuf.open();
				}

				pcur += it.npol*it.nchan;
			}
		}
		psf[n].close();
	}
	databuf.close();

	/**
	 * @brief flush the end data
	 * 
	 */
	int ngroup = dmsegs.size();

	int nleft = dedisp.offset/ndump;
	for (long int l=0; l<nleft; l++)
	{
		rfi.open();
		dedisp.prerun(rfi);
		for (long int k=0; k<ncand; k++)
		{
			dedisp.updatedm(dmsegs[k/GROUPSIZE]);
			dedisp.run();
			dedisp.get_subdata(subdata, k%GROUPSIZE);

			if (vm.count("dspsr"))
				folder[k].runDspsr(subdata);
			else
				folder[k].runTRLSM(subdata);
		}
		dedisp.postrun();
	}
	
	dedisp.close();

	double fmin = 1e6;
    double fmax = 0.;
    for (long int j=0; j<databuf.nchans; j++)
    {
        fmax = databuf.frequencies[j]>fmax? databuf.frequencies[j]:fmax;
        fmin = databuf.frequencies[j]<fmin? databuf.frequencies[j]:fmin;
    }

	double tint = ntotal*tsamp;

	vector<Pulsar::GridSearch> gridsearch(ncand);
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

		gridsearch[k].ddmstart = -scale*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin);
		gridsearch[k].ddmstep = 1./scale*abs(gridsearch[k].ddmstart/folder[k].nbin);
		gridsearch[k].nddm = 2*abs(scale)*folder[k].nbin;
		gridsearch[k].df0start = -scale*1./tint;
		gridsearch[k].df0step = 1./scale*abs(gridsearch[k].df0start/folder[k].nbin);
		gridsearch[k].ndf0 = 2*abs(scale)*folder[k].nbin;
		gridsearch[k].df1start = -scale*2./(tint*tint);
		gridsearch[k].df1step = 1./scale*abs(gridsearch[k].df1start/folder[k].nbin);
		gridsearch[k].ndf1 = 2*abs(scale)*folder[k].nbin;

		gridsearch[k].clfd_q = vm["clfd"].as<double>();

		gridsearch[k].prepare(folder[k]);
		folder[k].close();
		
		if (!nosearch)
		{
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

				gridsearch[k].runFFdot();
				gridsearch[k].bestprofiles();
				gridsearch[k].runDM();
				gridsearch[k].bestprofiles();

				dm1 = gridsearch[k].dm;
				f01 = gridsearch[k].f0;
				f11 = gridsearch[k].f1;

				cont++;
				//cout<<gridsearch[k].dm<<" "<<gridsearch[k].f0<<" "<<gridsearch[k].f1<<endl;
			}
		}

		/**
		 * @brief recalculate the mxsbr_ffdot and vsnr_dm with scale = 3 for plotting
		 * 
		 */

		dm = gridsearch[k].dm;
		f0 = gridsearch[k].f0;
		f1 = gridsearch[k].f1;

		double ddmstart = -3*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin);
		gridsearch[k].ddmstep = 1./3*abs(ddmstart/folder[k].nbin);
		int nddm = 2*3*folder[k].nbin;

		gridsearch[k].ddmstart = std::max(-3*1./f0/Pulsar::DedispersionLite::dmdelay(1, fmax, fmin), -dm);
		gridsearch[k].nddm = (ddmstart+gridsearch[k].ddmstep*nddm-gridsearch[k].ddmstart)/gridsearch[k].ddmstep;

		gridsearch[k].df0start = -3*1./tint;
		gridsearch[k].df0step = 1./3*abs(gridsearch[k].df0start/folder[k].nbin);
		gridsearch[k].ndf0 = 2*3*folder[k].nbin;
		gridsearch[k].df1start = -3*2./(tint*tint);
		gridsearch[k].df1step = 1./3*abs(gridsearch[k].df1start/folder[k].nbin);
		gridsearch[k].ndf1 = 2*3*folder[k].nbin;

		if (!nosearch)
		{
			gridsearch[k].runFFdot();
			gridsearch[k].runDM();
		}
	}

	/** form obsinfo*/
	std::map<std::string, std::string> obsinfo;
	//source name
	obsinfo["Source_name"] = src_name;
	//start mjd
	stringstream ss_mjd;
    ss_mjd << setprecision(10) << fixed << tstarts[idx[0]].to_day();
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
	obsinfo["Filename"] = fnames[idx[0]];
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

	double pepoch_offset = ntotal*tsamp/2.;
	
	//pepoch
	stringstream ss_pepoch;
    ss_pepoch << setprecision(9) << fixed << (tstarts[idx[0]]+pepoch_offset).to_day();
    string s_pepoch = ss_pepoch.str();
	obsinfo["Pepoch"] = s_pepoch;

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
	outfile<<"#Pepoch "<<obsinfo["Pepoch"]<<endl;
	outfile<<"#id       dm_old      dm_new      dm_err		dist_ymw16     f0_old     f0_new        f0_err      f1_old     f1_new       f1_err      acc_old        acc_new      acc_err      S/N        S/N_new"<<endl;

	for (long int k=0; k<ncand; k++)
	{
		stringstream ss_id;
		ss_id << setw(5) << setfill('0') << k+1;
		string s_id = ss_id.str();

		if (!noarch)
		{
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

		gridsearch[k].get_snr_width();
		gridsearch[k].get_error(obsinfo);

		double ymw16_dist = get_dist_ymw16(gl, gb, gridsearch[k].dm);

		/**
		 * @brief output best and old parameters to bestpar file
		 * 
		 */
		/**
		 * @brief id    dm_old  dm_new  dist_ymw16 f0_old  f0_new  f1_old  f1_new  acc_old acc_new S/N_old S/N_new
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
		outfile<<setprecision(15)<<gridsearch[k].f1<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_f1<<"\t\t";
		outfile<<setprecision(15)<<folder[k].f1/folder[k].f0*CONST_C<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].acc<<"\t\t";
		outfile<<setprecision(15)<<gridsearch[k].err_acc<<"\t\t";
		outfile<<fixed<<setprecision(5)<<folder[k].snr<<"\t\t";
		outfile<<fixed<<setprecision(5)<<gridsearch[k].snr<<endl;

		if (!noplot)
		{
			obsinfo["Dist_YMW16"] = to_string(ymw16_dist);
			Pulsar::PulsarPlot psrplot;
			psrplot.plot(folder[k], gridsearch[k], obsinfo, k+1, rootname, vm.count("plotx"));
		}
	}

	outfile.close();

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/ntotal<<"%)"<<endl;
	}

	delete [] buffer;
	delete [] psf;

    return 0;
}

void produce(variables_map &vm, std::vector<std::vector<double>> &dmsegs, vector<Pulsar::ArchiveLite> &folder)
{
	std::vector<double> dmseg;

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

    dmseg.push_back(vm["dm"].as<double>());

    if (vm.count("candfile"))
    {
        string filename = vm["candfile"].as<string>();
        string line;
        ifstream candfile(filename);
        
        dmseg.clear();
		int gcnt = 0;
        while (getline(candfile, line))
        {
			boost::trim(line);
			if(!isdigit(line[0])) continue;
			
            vector<string> parameters;
            boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);

            dmseg.push_back(stod(parameters[1]));

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

			if (++gcnt == GROUPSIZE)
			{
				dmsegs.push_back(dmseg);
				dmseg.clear();
				gcnt = 0;
			}
            folder.push_back(fdr);
        }

		if (gcnt)
		{
			dmsegs.push_back(dmseg);
			dmseg.clear();
			gcnt = 0;
		}
    }
    else
    {
		dmsegs.push_back(dmseg);
        folder.push_back(fdr);
    }
}
