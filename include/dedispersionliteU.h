/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 17:01:10
 * @modify date 2020-06-07 17:01:10
 * @desc [description]
 */

#ifndef DEDISPERSIONLITEU_H
#define DEDISPERSIONLITEU_H

#include <fstream>

#include "databuffer.h"

using namespace std;

namespace Pulsar
{
	class DedispersionLiteU
	{
	public:
		DedispersionLiteU();
		~DedispersionLiteU();
		void close()
		{
			delayn.clear();
			delayn.shrink_to_fit();

			buffer.clear();
			buffer.shrink_to_fit();

			bufferT.clear();
			bufferT.shrink_to_fit();

			buffersub.clear();
			buffersub.shrink_to_fit();
		}
		void prepare(DataBuffer<float> &databuffer);
		void updatedm(const std::vector<double> &dms);
		void prerun(DataBuffer<float> &databuffer);
		void run();
		void postrun();
		void dumpsubdata(const string &rootname, int idm) const
		{
			ofstream outfile;
			outfile.open(rootname+".sub", ios::binary|ios::app);

			outfile.write((char *)(&buffersub[0]+idm*ndump*nsubband), sizeof(float)*ndump*nsubband);
		
			outfile.close();
		}

		void get_subdata(DataBuffer<float> &databuffer, int idm)
		{
			databuffer.resize(ndump, nsubband);
			databuffer.tsamp = tsamp;
			databuffer.frequencies = frequencies_sub;

			if (!buffersub.empty())
			{
				for (long int i=0; i<ndump; i++)
				{
					for (long int j=0; j<nsubband; j++)
					{
						databuffer.buffer[i*nsubband+j] = buffersub[idm*ndump*nsubband+i*nsubband+j];
					}
				}
				
				databuffer.means = means;
				databuffer.vars = vars;
				databuffer.mean_var_ready = mean_var_ready;

				databuffer.counter += ndump;
			}
		}
	public:
		int nsubband;
		vector<double> vdm;
		double maxdm;
		int groupsize;
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
		vector<float> bufferT;
		vector<double> frequencies_sub;
		vector<float> buffersub;
		std::vector<double> means;
		std::vector<double> vars;
		bool mean_var_ready;
	public:
		static double dmdelay(double dm, double fh, double fl)
		{
			return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
		}
	};
}

#endif /* DEDISPERSIONLITEU_H */
