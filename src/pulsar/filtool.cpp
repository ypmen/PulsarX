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
#include "preprocesslite.h"
#include "filmaker.h"
#include "filterbank.h"
#include "mjd.h"
#include "logging.h"

using namespace boost::program_options;

unsigned int num_threads;

#define NSBLK 65536
#define FAST 1

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
			("ra", value<double>()->default_value(0), "RA (hhmmss.s)")
			("dec", value<double>()->default_value(0), "DEC (ddmmss.s)")
			("ibeam,i", value<int>()->default_value(1), "Beam number")
			("telescope", value<std::string>()->default_value("Fake"), "Telescope name")
			("baseline", value<std::vector<float>>()->multitoken()->default_value(std::vector<float>{0.0, 0.0}, "0.0, 0.0"), "The scale of baseline remove (s)")
			("rfi,z", value<std::vector<std::string>>()->multitoken()->zero_tokens()->composing(), "RFI mitigation [[mask tdRFI fdRFI] [kadaneF tdRFI fdRFI] [kadaneT tdRFI fdRFI] [zap fl fh] [zdot] [zero]]")
			("bandlimit", value<double>()->default_value(10), "Band limit of RFI mask (MHz)")
			("bandlimitKT", value<double>()->default_value(10), "Band limit of RFI kadaneT (MHz)")
			("widthlimit", value<double>()->default_value(50e-3), "Width limit of RFI kadaneF (s)")
			("threMask", value<float>()->default_value(10), "S/N threshold of Mask")
			("threKadaneF", value<float>()->default_value(7), "S/N threshold of KadaneF")
			("threKadaneT", value<float>()->default_value(7), "S/N threshold of KadaneT")
			("fill", value<string>()->default_value("mean"), "Fill the zapped samples by [mean, rand]")
			("source_name,s", value<std::string>()->default_value("J0000-00"), "Source name")
			("rootname,o", value<std::string>()->default_value("J0000-00"), "Output rootname")
			("cont", "Input files are contiguous")
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
	std::string source_name = vm["source_name"].as<std::string>();
	std::string rootname = vm["rootname"].as<std::string>();
	std::string s_telescope = vm["telescope"].as<std::string>();
	int telescope_id = get_telescope_id(s_telescope);
	int ibeam = vm["ibeam"].as<int>();
	double src_raj = vm["ra"].as<double>();
	double src_dej = vm["dec"].as<double>();

	// sort fils and get nsamples
	long int nfil = fnames.size();
	Filterbank *fil = new Filterbank [nfil];
	for (long int i=0; i<nfil; i++)
	{
		fil[i].filename = fnames[i];
	}

	vector<MJD> tstarts;
	vector<MJD> tends;
	long int ntotal = 0;
	for (long int i=0; i<nfil; i++)
	{
		fil[i].read_header();
		ntotal += fil[i].nsamples;
		MJD tstart(fil[i].tstart);
		tstarts.push_back(tstart);
		tends.push_back(tstart+fil[i].nsamples*fil[i].tsamp);
	}
	vector<size_t> idx = argsort(tstarts);
	for (long int i=0; i<nfil-1; i++)
	{
		if (abs((tends[idx[i]]-tstarts[idx[i+1]]).to_second())>0.5*fil[idx[i]].tsamp)
		{
			if (contiguous)
			{
				std:cerr<<"Warning: time not contiguous"<<endl;
			}
			else
			{
				std::cerr<<"Error: time not contiguous"<<endl;
				exit(-1);
			}
		}
	}

	// update header
	if (vm["source_name"].defaulted())
	{
		if (strcmp(fil[0].source_name, "") != 0)
			source_name = fil[0].source_name;
	}

	if (vm["telescope"].defaulted())
	{
		telescope_id =  fil[0].telescope_id;
	}

	if (vm["ibeam"].defaulted())
	{
		if (fil[0].ibeam != 0)
			ibeam = fil[0].ibeam;
	}

	if (vm["ra"].defaulted())
	{
		if (fil[0].src_raj != 0.)
		{
			src_raj = fil[0].src_raj;
		}
	}
	if (vm["dec"].defaulted())
	{
		if (fil[0].src_dej != 0.)
		{
			src_dej = fil[0].src_dej;
		}
	}

	std::vector<FilMaker> filmakers;
	plan(vm, filmakers);

	/** pre-downsample */
	int td = vm["td"].as<int>();
	int fd = vm["fd"].as<int>();

	long int nchans = fil[0].nchans;
	double tsamp = fil[0].tsamp;
	int nifs = fil[0].nifs;

	float *buffer = new float [nchans];

	vector<int> tds;
	for (auto fm=filmakers.begin(); fm!=filmakers.end(); ++fm)
	{
		tds.push_back((*fm).td*td);
	}

	long int td_lcm = findlcm(&tds[0], tds.size());

	long int ndump = (int)(vm["seglen"].as<float>()/tsamp)/td_lcm/2*td_lcm*2;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.closable = true;
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], fil[0].frequency_table, sizeof(double)*nchans);

	PreprocessLite prep;
	prep.td = td;
	prep.fd = fd;
	prep.thresig = vm["zapthre"].as<float>();
	prep.filltype = vm["fill"].as<string>();
	prep.prepare(databuf);

	long int noutfil = filmakers.size();
	for (long int k=0; k<noutfil; k++)
	{
		filmakers[k].filwriter.fil = fil[idx[0]];
		filmakers[k].ibeam = ibeam;
		filmakers[k].rootname = rootname;
		filmakers[k].source_name = source_name;
		filmakers[k].src_raj = src_raj;
		filmakers[k].src_dej = src_dej;
		filmakers[k].telescope_id = telescope_id;
		filmakers[k].filltype = vm["fill"].as<string>();
		filmakers[k].prepare(prep);
	}

	// handle input data
	long int nstart = jump[0]/tsamp;
	long int nend = ntotal-jump[1]/tsamp;

	int sumif = nifs>2? 2:nifs;
	
	long int ntot = 0;
	long int ntot2 = 0;
	long int count = 0;
	long int bcnt1 = 0;
	for (long int idxn=0; idxn<nfil; idxn++)
	{
		long int n = idx[idxn];
		long int nseg = ceil(1.*fil[n].nsamples/NSBLK);
		long int ns_filn = 0;

		for (long int s=0; s<nseg; s++)
		{
			if (verbose)
			{
				std::cerr<<"\r\rfinish "<<std::setprecision(2)<<std::fixed<<tsamp*count<<" seconds ";
				std::cerr<<"("<<100.*count/ntotal<<"%)";
			}

			fil[n].read_data(NSBLK);
			assert(fil[n].ndata != 0);
#ifdef FAST
			unsigned char *pcur = (unsigned char *)(fil[n].data);
#endif
			for (long int i=0; i<NSBLK; i++)
			{
				count++;

				if (count-1<nstart or count-1>nend)
				{
					if (++ns_filn == fil[n].nsamples)
					{
						goto next;
					}
					pcur += nifs*nchans;
					continue;
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
					DataBuffer<float> *data = prep.run(databuf);

					for (long int ioutfil=0; ioutfil<noutfil; ioutfil++)
					{
						data->isbusy = true;
						filmakers[ioutfil].run(*data);
					}

					bcnt1 = 0;
					databuf.open();
				}

				if (++ns_filn == fil[n].nsamples)
				{
					goto next;
				}
				pcur += nifs*nchans;
			}
		}
		next:
		fil[n].close();
	}
	databuf.close();

	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/ntotal<<"%)"<<endl;
	}

	delete [] buffer;
	delete [] fil;

	return 0;
}