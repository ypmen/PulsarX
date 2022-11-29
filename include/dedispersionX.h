/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-04 08:08:56
 * @modify date 2022-08-04 08:08:56
 * @desc [ultra efficient dedispersion for folding]
 */


#ifndef DEDISPERSIONX_H
#define DEDISPERSIONX_H

#include <vector>
#include <list>
#include "databuffer.h"
#include "downsample.h"
#include "dedisperse.h"

#ifdef __AVX2__
#include <boost/align/aligned_allocator.hpp>
typedef __attribute__(( aligned(32))) float aligned_float;
#endif

namespace Pulsar
{
	class TreeDedispersion
	{
	public:
		TreeDedispersion();
		~TreeDedispersion();
		void alloc_dedata()
		{
			dedata.resize(nchans * nsamples, 0.);
		}
		void free_dedata()
		{
			dedata.clear();
			dedata.shrink_to_fit();
		}
		void alloc_buffer()
		{
			buffer.resize(nsamples * nchans, 0.);
		}
		void free_buffer()
		{
			buffer.clear();
			buffer.shrink_to_fit();
		}
		void alloc_bufferT()
		{
			bufferT.resize(nchans * nsamples, 0.);
		}
		void free_bufferT()
		{
			bufferT.clear();
			bufferT.shrink_to_fit();
		}
		void close();
		void prepare(DataBuffer<float> &databuffer);
		void update_hit(const std::vector<double> &vdm);
		void run(DataBuffer<float> &databuffer);
		void run();
		void get_subdata(double dm, DataBuffer<float> &subdata);
		void get_subdata_tem(double dm, DataBuffer<float> &subdata)
		{
			subdata.resize(ndump, nsubband);
			subdata.tsamp = tsamp;
			subdata.frequencies = frequencies_sub;
		}
		void resize_cache()
		{
			cache0.resize(num_threads * nsamples, 0.);
			cache1.resize(num_threads * nsamples, 0.);
		}

	public:
		bool is_ready(){return ready;}
		size_t get_nchans(){return nchans;}
		size_t get_offset(){return offset;}
		double get_tsamp(){return tsamp;}
		size_t get_ndump(){return ndump;}
		size_t get_nsamples(){return nsamples;}
		size_t get_counter(){return counter;}

		void update_nchans(size_t n){nchans = n;}
		void update_ndump(size_t n){ndump = n;}
		void update_nsamples(size_t n){nsamples = n;offset = nsamples-ndump;}
		void update_tsamp(double t){tsamp = t;}
		void update_frequencies(const std::vector<double> &f){frequencies = f;}

		template<typename Iterator>
		void update_hit(Iterator begin, Iterator end)
		{
			std::vector<bool> hit0(nchans, false);

			for (auto dm=begin; dm!=end; ++dm)
			{
				size_t dmid = (*dm - dms) / ddm;
				hit0[map[dmid]] = true;
			}

			hit.resize((maxdepth + 1) * nchans, false);

			for (size_t i=0; i<hit0.size(); i++)
			{
				if (hit0[i])
				{
					for (size_t k=0; k<maxdepth+1; k++)
					{
						for (size_t j=0; j<nchans; j+=(nchans>>k))
						{
							hit[k * nchans + (i + j) % nchans] = true;
						}
					}
				}
			}
		}

		void hit_all()
		{
			hit.resize((maxdepth + 1) * nchans, false);
			std::fill(hit.begin(), hit.end(), true);
		}

	private:
		void update_delay();
		void update_delay_rec(std::vector<size_t> &vdmid, double &freq, size_t depth, size_t ichan);
		void update_map();
		void transform(size_t depth, size_t ichan);

	public:
		size_t nsubband;
		double dms;
		double ddm;

		bool if_alloc_buffer;
		bool if_alloc_bufferT;
		bool if_alloc_dedata;

		std::vector<float> *ptr_dedata;
		std::vector<float> *ptr_buffer;
#ifndef __AVX2__
		std::vector<float> *ptr_bufferT;
#else
		std::vector<float, boost::alignment::aligned_allocator<float, 32>> *ptr_bufferT;
#endif
	
	private:
		std::vector<float> dedata;
		std::vector<float> buffer;
#ifndef __AVX2__
		std::vector<float> bufferT;
#else
		std::vector<float, boost::alignment::aligned_allocator<float, 32>> bufferT;
#endif

	private:
		bool ready;
		size_t counter;
		size_t nchans;
		size_t nchans_orig;
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
		std::vector<size_t> map;
		std::vector<size_t> mapsub;
		std::vector<int> delayn;
		std::vector<bool> hit;
#ifndef __AVX2__
		std::vector<float> cache0;
		std::vector<float> cache1;
#else
		std::vector<float, boost::alignment::aligned_allocator<float, 32>> cache0;
		std::vector<float, boost::alignment::aligned_allocator<float, 32>> cache1;
#endif

	public:
		static double dmdelay(double dm, double fh, double fl)
		{
			return 4.148741601e3*dm*(1./(fl*fl)-1./(fh*fh));
		}
	};

	class DedispersionX
	{
	public:
		DedispersionX();
		~DedispersionX();
		void prepare(DataBuffer<float> &databuffer);
		void prerun(DataBuffer<float> &databuffer);
		void postrun(DataBuffer<float> &databuffer);
		void run(DataBuffer<float> &databuffer);
		void run(size_t k);
		void alloc(size_t k)
		{
			if (downsamples[k].td == 1 && downsamples[k].fd == 1)
			{
				bufferTs[0].resize(nchans * nsamples0, 0.);
				dedatas[0].resize(nchans * nsamples0, 0.);
			}
			else
			{
				bufferTs[k].resize(treededispersions[k].get_nchans() * treededispersions[k].get_nsamples(), 0.);
				dedatas[k].resize(treededispersions[k].get_nchans() * treededispersions[k].get_nsamples(), 0.);
			}
		}
		void free(size_t k)
		{
			if (downsamples[k].td == 1 && downsamples[k].fd == 1)
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
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				downsamples[k].close();
				treededispersions[k].close();

				buffers[k].clear();
				buffers[k].shrink_to_fit();

				bufferTs[k].clear();
				bufferTs[k].shrink_to_fit();

				dedatas[k].clear();
				dedatas[k].shrink_to_fit();
			}
		}

		void get_subdata(double dm, DataBuffer<float> &subdata)
		{
			bool succ = false;
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					treededispersions[k].get_subdata(dm, subdata);
					succ = true;
				}
			}

			if (!succ)
				BOOST_LOG_TRIVIAL(error) << "dm is out of range";
		}

		void get_subdata_tem(double dm, DataBuffer<float> &subdata)
		{
			bool succ = false;
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					treededispersions[k].get_subdata_tem(dm, subdata);
					succ = true;
				}
			}

			if (!succ)
				BOOST_LOG_TRIVIAL(error) << "dm is out of range";
		}
		
		void get_dmslist(std::vector<double> &dmslist)
		{
			for (auto d=treededispersions.begin(); d!=treededispersions.end(); ++d)
			{
				dmslist.push_back(d->dms);
			}
		}

		size_t get_offset(double dm)
		{
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					if (downsamples[k].td == 1 && downsamples[k].fd == 1)
					{
						size_t nsamples = buffers[0].size()/nchans;
						size_t ndump = treededispersions[k].get_ndump();
						return nsamples - ndump;
					}
					else
					{
						return treededispersions[k].get_offset();
					}
				}
			}

			BOOST_LOG_TRIVIAL(error) << "dm is out of range";

			return -1;
		}

		size_t get_ndump(double dm)
		{
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					return treededispersions[k].get_ndump();
				}
			}

			BOOST_LOG_TRIVIAL(error) << "dm is out of range";

			return -1;
		}

		double get_tsamp(double dm)
		{
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					return treededispersions[k].get_tsamp();
				}
			}

			BOOST_LOG_TRIVIAL(error) << "dm is out of range";

			return 0.;
		}

		size_t get_id(double dm)
		{
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					return k;
				}
			}

			BOOST_LOG_TRIVIAL(error) << "dm is out of range";

			return -1;
		}

		size_t get_counter(double dm)
		{
			for (size_t k=0; k<treededispersions.size(); k++)
			{
				double dms = treededispersions[k].dms;
				double dme = treededispersions[k].dms + treededispersions[k].get_nchans() * treededispersions[k].ddm;
				if (dms <= dm && dme > dm)
				{
					return treededispersions[k].get_counter();
				}
			}

			BOOST_LOG_TRIVIAL(error) << "dm is out of range";

			return -1;
		}

		size_t get_nsegments(){return downsamples.size();}

		double get_dms(size_t k){return treededispersions[k].dms;}
		double get_ddm(size_t k){return treededispersions[k].ddm;}

		bool is_hit(size_t k){return hit[k];}
	
	public:
		size_t nsubband;
		std::list<double> dmlist;
		size_t dm_boost; // optional
		double ddm_init; // optional
	
	private:
		size_t nchans;
		size_t nchans_orig;
		size_t nsamples0;
		std::vector<TreeDedispersion> treededispersions;
		std::vector<Downsample> downsamples;
		std::vector<std::vector<float>> buffers;
#ifndef __AVX2__
		std::vector<std::vector<float>> bufferTs;
#else
		std::vector<std::vector<float, boost::alignment::aligned_allocator<float, 32>>> bufferTs;
#endif
		std::vector<std::vector<float>> dedatas;
		std::vector<bool> hit;

	public:
		static size_t get_maxds(double maxdm, double tsamp, double fmax, double fmin, size_t nchans, double ddm_init=0., double dm_boost=0.)
		{
			double ddm = ddm_init;
			if (ddm == 0.)
				ddm = tsamp / TreeDedispersion::dmdelay(1., fmax, fmin);
			size_t n = 0;
			double dms = 0.;
			do
			{
				dms += nchans * ddm;
				if (dms >= dm_boost)
				{
					ddm *= 2;
					n++;
				}
			} while (dms <= maxdm);

			n = n == 0 ? 1 : n;

			return 1<<(n-1);
		}
	};
}


#endif /* DEDISPERSIONX_H */
