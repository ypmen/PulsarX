/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-05-19 08:58:11
 * @modify date 2020-05-19 08:58:11
 * @desc [dedisperse for pulsar searching]
 */

#define FAST 1

#include <iostream>
#include <iomanip>
#include <string.h>
#include <boost/program_options.hpp>

#include "pulsarsearch.h"
#include "subdedispersion.h"
#include "dedisperse.h"
#include "psrfits.h"
#include "utils.h"
#include "preprocesslite.h"

using namespace std;
using namespace boost::program_options;

unsigned int num_threads;
unsigned int dbscan_radius;
unsigned int dbscan_k;

bool dumptim=true;

int main(int argc, const char *argv[])
{
    /* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("jump,j", value<vector<double>>()->multitoken()->default_value(vector<double>{7200, 0}, "7200, 0"), "Segment length and following jump length (s)")
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
			("baseline", value<vector<float>>()->multitoken()->default_value(vector<float>{0.0, 0.0}, "0.0, 0.1"), "The scale of baseline remove (s)")
			("rfi,z", value<vector<string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("threMask", value<float>()->default_value(3), "S/N threshold of Mask")
            ("rootname,o", value<string>()->default_value("J0000-00"), "Output rootname")
			("format", value<string>()->default_value("pulsarx"), "Output format of dedispersed data [pulsarx(default),sigproc,presto]")
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
	if (vm.count("input") == 0)
	{
		cerr<<"Error: no input file"<<endl;
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

	Integration it;
	psf[0].subint.load_integration(psf[0].fptr, 0, it);

	long int nchans = it.nchan;
    double tsamp = psf[0].subint.tbin;
    int nifs = it.npol;

	float *buffer = new float [nchans];

	vector<PulsarSearch> search;
	plan(vm, search);

	vector<int> tds;
	for (auto sp=search.begin(); sp!=search.end(); ++sp)
	{
		tds.push_back((*sp).td*vm["td"].as<int>());
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm/2*td_lcm*2;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.closable = true;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], it.frequencies, sizeof(double)*nchans);

	PreprocessLite prep;
	prep.td = vm["td"].as<int>();
	prep.fd = vm["fd"].as<int>();
	prep.thresig = vm["zapthre"].as<float>();
	prep.prepare(databuf);

	long int nseg = jump[0]/tsamp;
	long int njmp = jump[1]/tsamp;

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
		search[k].fildedisp.tstart = tstarts[idx[0]].to_day();
		search[k].fildedisp.ibeam = ibeam;
		search[k].fildedisp.fch1 = (databuf.frequencies.front()+databuf.frequencies.back())/2.;
		search[k].fildedisp.foff = databuf.frequencies.back()-databuf.frequencies.front();
		search[k].fildedisp.nchans = databuf.frequencies.size();
		search[k].prepare(prep);
	}

	psf[0].close();

	int sumif = nifs>2? 2:nifs;
	
    long int jmpcont = 0;
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
                if (ntot == nseg)
                {
                    if (jmpcont++ < njmp)
                    {
                        pcur += it.npol*it.nchan;
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
                        search[k].fildedisp.tstart = (tstarts[idx[0]] + count*tsamp/86400.).to_day();
						search[k].dedisp.preparedump(search[k].fildedisp, outnbits, vm["format"].as<string>());
	                }
                }

				memset(buffer, 0, sizeof(float)*nchans);
				long int m = 0;
				for (long int k=0; k<sumif; k++)
				{
					for (long int j=0; j<nchans; j++)
					{
						buffer[j] +=  pcur[m++];
					}
				}

				memcpy(&databuf.buffer[0]+bcnt1*nchans, buffer, sizeof(float)*1*nchans);
                databuf.counter++;
				bcnt1++;
                ntot++;

				if (ntot%ndump == 0)
				{
					prep.run(databuf);
					for (auto sp=search.begin(); sp!=search.end(); ++sp)
					{
						prep.isbusy = true;
						(*sp).run(prep);
					}
                    bcnt1 = 0;
					databuf.open();
				}

				pcur += it.npol*it.nchan;
			}
		}
		psf[n].close();
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

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/ntotal<<"%)"<<endl;
	}

	delete [] buffer;
	delete [] psf;

    return 0;
}
