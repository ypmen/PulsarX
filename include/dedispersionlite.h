/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 17:01:10
 * @modify date 2020-06-07 17:01:10
 * @desc [description]
 */

#ifndef DEDISPERSIONLITE
#define DEDISPERSIONLITE

#include <fstream>

#include "databuffer.h"

using namespace std;

namespace Pulsar
{
	class DedispersionLite
	{
	public:
		DedispersionLite();
		~DedispersionLite();
		void close()
		{
			delayn.clear();
			delayn.shrink_to_fit();

			buffer.clear();
			buffer.shrink_to_fit();

			buffersub.clear();
			buffersub.shrink_to_fit();

			buffertim.clear();
			buffertim.shrink_to_fit();
		}
		void prepare(DataBuffer<float> &databuffer);
		void run(DataBuffer<float> &databuffer);
		void dumpsubdata(const string &rootname, int idm) const
		{
			ofstream outfile;
			outfile.open(rootname+".sub", ios::binary|ios::app);

			outfile.write((char *)(&buffersub[0]+idm*ndump*nsubband), sizeof(float)*ndump*nsubband);
		
			outfile.close();
		}

		void dumptimdata(const string &rootname, int idm) const
		{
			ofstream outfile;
			outfile.open(rootname+".tim", ios::binary|ios::app);

			outfile.write((char *)(&buffertim[0]+idm*ndump), sizeof(float)*ndump);
		
			outfile.close();
		}

		void get_subdata(DataBuffer<float> &databuffer, int idm)
		{
			databuffer.resize(ndump, nsubband);
			databuffer.tsamp = tsamp;
			databuffer.frequencies = frequencies_sub;

			for (long int i=0; i<ndump; i++)
			{
				for (long int j=0; j<nsubband; j++)
				{
					databuffer.buffer[i*nsubband+j] = buffersub[idm*ndump*nsubband+i*nsubband+j];
				}
			}

			databuffer.counter += ndump;
		}
	public:
		long int nsubband;
		vector<double> vdm;
	public:
		long int counter;
		long int offset;
		long int nchans;
		double tsamp;
		long int ndump;
		long int nsamples;
		vector<double> frequencies;
		vector<long int> delayn;
		vector<float> buffer;
		vector<double> frequencies_sub;
		vector<float> buffertim;
		vector<float> buffersub;
	public:
		static double dmdelay(double dm, double fh, double fl)
		{
			return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
		}
	};
}

#endif /* DEDISPERSIONLITE */
