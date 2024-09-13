/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-10-13 21:17:40
 * @modify date 2022-10-13 21:17:40
 * @desc [description]
 */

#ifndef DERMDM_H
#define DERMDM_H

#include <list>

#include "databuffer.h"
#include "cdownsample.h"
#include <complex>
#include "constants.h"
#include "logging.h"

namespace Pulsar
{
	class DeRMDM
	{
	public:
		DeRMDM();
		~DeRMDM();
		void close();
		void prepare(DataBuffer<std::complex<float>> &databuffer);
		void run(DataBuffer<std::complex<float>> &databuffer);
		void run();
		void get_subdata(double dm, double rm, DataBuffer<std::complex<float>> &subdata);

		template<typename Iterator>
		void update_dmhit(Iterator begin, Iterator end)
		{
			std::vector<bool> hit0(nchans, false);

			for (auto dm=begin; dm!=end; ++dm)
			{
				size_t dmid = (*dm - dms) / ddm;
				hit0[mapdm[dmid]] = true;
			}

			hitdm.resize((maxdepth + 1) * nchans, false);

			for (size_t i=0; i<hit0.size(); i++)
			{
				if (hit0[i])
				{
					for (size_t k=0; k<maxdepth+1; k++)
					{
						for (size_t j=0; j<nchans; j+=(nchans>>k))
						{
							hitdm[k * nchans + (i + j) % nchans] = true;
						}
					}
				}
			}
		}

		void dmhit_all()
		{
			hitdm.resize((maxdepth + 1) * nchans, true);
		}

		template<typename Iterator>
		void update_rmhit(Iterator begin, Iterator end)
		{
			std::vector<bool> hit0(numrm, false);

			for (auto rm=begin; rm!=end; ++rm)
			{
				size_t rmid = (*rm - rms) / drm;
				hit0[maprm[rmid]] = true;
			}

			hitrm.resize(depthrm * numrm, false);

			for (size_t i=0; i<hit0.size(); i++)
			{
				if (hit0[i])
				{
					for (size_t k=0; k<depthrm; k++)
					{
						for (size_t j=0; j<numrm; j+=(numrm>>k))
						{
							hitrm[k * numrm + (i + j) % numrm] = true;
						}
					}
				}
			}
		}

		void rmhit_all()
		{
			hitrm.resize(depthrm * numrm, true);
		}

		void resize_cache()
		{
			cache00.resize(nsamples, 0.);
			cache01.resize(nsamples, 0.);
			cache10.resize(nsamples, 0.);
			cache11.resize(nsamples, 0.);
		}

		size_t get_nchans(){return nchans;}
		size_t get_nsamples(){return nsamples;}
		size_t get_numrmsub(){return numrmsub;}
		size_t get_numrm(){return numrm;}
		size_t get_ndump(){return ndump;}
		bool is_ready(){return ready;}

		void update_nsamples(size_t n){nsamples = n;offset = nsamples-ndump;}

	public:
		void transform(double &freq, size_t depth, size_t ichan);

		void update_dmdelay();
		void update_dmdelay_rec(std::vector<size_t> &vdmid, double &freq, size_t depth, size_t ichan);
		void update_dmmap();

		void update_rmdelay();
		void update_rmdelay_rec(std::vector<size_t> &vrmid, double &freq, size_t depth, size_t ichan);
		void update_rmmap();

		void update_freq_rec(double &freq, size_t depth, size_t ichan);

	public:
		double dms;
		double ddm;
		double rms;
		double drm;
		size_t nsubband;

		size_t depthrm;

		bool if_alloc_buffer;
		bool if_alloc_bufferT;
		bool if_alloc_dedata;

		std::vector<std::complex<float>> *ptr_dedata;
		std::vector<std::complex<float>> *ptr_buffer;
#ifndef __AVX2__
		std::vector<std::complex<float>> *ptr_bufferT;
#else
		std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> *ptr_bufferT;
#endif

	public:
		std::vector<std::complex<float>> buffer;
		std::vector<std::complex<float>> dedata;
	#ifndef __AVX2__
			std::vector<std::complex<float>> bufferT;
	#else
			std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> bufferT;
	#endif

	public:
		bool ready;
		size_t counter;
		size_t numrm;
		size_t numrmsub;
		size_t nchans;
		size_t ndump;
		double tsamp;
		std::vector<double> frequencies;
		std::vector<double> frequencies_sub;
		size_t maxdepth;
		size_t depthsub;
		double fmax;
		double fmin;
		int maxdelayn;
		size_t offset;
		size_t nsamples;

		std::vector<size_t> mapdm;
		std::vector<size_t> mapsubdm;
		std::vector<int> dmdelayns;
		std::vector<bool> hitdm;
		std::vector<size_t> maprm;
		std::vector<size_t> mapsubrm;
		std::vector<float> rmdelays;
		std::vector<bool> hitrm;
		
#ifndef __AVX2__
		std::vector<std::complex<float>> cache00;
		std::vector<std::complex<float>> cache01;
		std::vector<std::complex<float>> cache10;
		std::vector<std::complex<float>> cache11;
#else
		std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> cache00;
		std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> cache01;
		std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> cache10;
		std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> cache11;
#endif

	public:
		static double dmdelay(double dm, double fh, double fl)
		{
			//return 1/2.41e-4 * dm * (1. / (fl * fl) - 1. / (fh * fh));
			return dm * (fh-fl);

		}

		static double rmdelay(double rm, double fh, double fl)
		{
			//return rm * CONST_C * CONST_C * (1. / (fl * fl) - 1. / (fh * fh));
			return rm * 2 * M_PI / 8 * (fh - fl);
		}
	};

	class DeRMDMX
	{
	public:
		DeRMDMX();
		~DeRMDMX();
		void prepare(DataBuffer<std::complex<float>> &databuffer);
		void prerun(DataBuffer<std::complex<float>> &databuffer);
		void postrun(DataBuffer<std::complex<float>> &databuffer);
		void run(DataBuffer<std::complex<float>> &databuffer);
		void run(size_t i, size_t j);

		void get_subdata(double dm, double rm, DataBuffer<std::complex<float>> &subdata)
		{
			bool succ = false;
			for (size_t k=0; k<dermdms.size(); k++)
			{
				double dms = dermdms[k].dms;
				double dme = dermdms[k].dms + dermdms[k].get_nchans() * dermdms[k].ddm;

				double rms = dermdms[k].rms;
				double rme = dermdms[k].rms + dermdms[k].get_numrm() * dermdms[k].drm;

				if (dms <= dm && dme > dm && rms <= rm && rme > rm)
				{
					dermdms[k].get_subdata(dm, rm, subdata);
					succ = true;
				}
			}

			if (!succ)
				BOOST_LOG_TRIVIAL(error) << "dm or rm is out of range";
		}

		void alloc(size_t k)
		{
			if (cdownsamples[k].td == 1 && cdownsamples[k].fd == 1)
			{
				bufferTs[0].resize(numrmsub * nchans * nsamples0, 0.);
				dedatas[0].resize(numrmsub * nchans * nsamples0, 0.);
			}
			else
			{
				bufferTs[k].resize(dermdms[k].get_numrmsub() * dermdms[k].get_nchans() * dermdms[k].get_nsamples(), 0.);
				dedatas[k].resize(dermdms[k].get_numrmsub() * dermdms[k].get_nchans() * dermdms[k].get_nsamples(), 0.);
			}
		}
		void free(size_t k)
		{
			if (cdownsamples[k].td == 1 && cdownsamples[k].fd == 1)
			{
				bufferTs[0].clear();
				bufferTs[0].shrink_to_fit();

				dedatas[0].clear();
				dedatas[0].shrink_to_fit();
			}
			else
			{
				bufferTs[k].clear();
				bufferTs[k].shrink_to_fit();

				dedatas[k].clear();
				dedatas[k].shrink_to_fit();
			}
		}
		
		void close()
		{
			for (size_t k=0; k<dermdms.size(); k++)
			{
				cdownsamples[k].close();
				dermdms[k].close();

				buffers[k].clear();
				buffers[k].shrink_to_fit();

				bufferTs[k].clear();
				bufferTs[k].shrink_to_fit();

				dedatas[k].clear();
				dedatas[k].shrink_to_fit();
			}
		}

	
	public:
		size_t nsubband;
		double dms;
		double dme;
		double rms;
		double rme;
		double rm_smear; // in phase 0-1
		double drm_init; // optional
		double rm_boost;
		size_t dm_boost; // optional
		double ddm_init; // optional

		size_t depthrm;
	
	private:
		std::vector<std::tuple<double, double, size_t, size_t>> ddplan;
		std::vector<std::tuple<double, double, size_t>> drplan;
		size_t numrm;
		size_t numrmsub;
		size_t nchans;
		size_t nsamples0;
		std::vector<DeRMDM> dermdms;
		std::vector<cDownsample> cdownsamples;
		std::vector<std::vector<std::complex<float>>> buffers;
#ifndef __AVX2__
		std::vector<std::vector<std::complex<float>>> bufferTs;
#else
		std::vector<std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>>> bufferTs;
#endif
		std::vector<std::vector<std::complex<float>>> dedatas;
	};
}

#endif /* DERMDM_H */
