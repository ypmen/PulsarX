/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-08 21:55:15
 * @modify date 2022-08-08 21:55:15
 * @desc [description]
 */

#include "dedispersionX.h"
#include "utils.h"

#ifdef _OPENMP
	#include <omp.h>
#endif

using namespace Pulsar;

TreeDedispersion::TreeDedispersion()
{
	nsubband = 0;
	nchans = 0;
	nchans_orig = 0;
	ndump = 0;
	tsamp = 0.;
	dms = 0.;
	ddm = 0.;

	if_alloc_buffer = true;
	if_alloc_bufferT = true;
	if_alloc_dedata = true;

	ptr_dedata = NULL;
	ptr_buffer = NULL;
	ptr_bufferT = NULL;

	ready = false;
	fmax = 0.;
	fmin = 0.;
	maxdelayn = 0;
	offset = 0;
	nsamples = 0;

	maxdepth = 0;
	depthsub = 0;

	counter = 0;
}

TreeDedispersion::~TreeDedispersion()
{
}

void TreeDedispersion::close()
{
	dedata.clear();
	dedata.shrink_to_fit();

	buffer.clear();
	buffer.shrink_to_fit();

	bufferT.clear();
	bufferT.shrink_to_fit();

	map.clear();
	map.shrink_to_fit();

	mapsub.clear();
	mapsub.shrink_to_fit();

	delayn.clear();
	delayn.shrink_to_fit();

	hit.clear();
	hit.shrink_to_fit();

	cache0.clear();
	cache0.shrink_to_fit();

	cache1.clear();
	cache1.shrink_to_fit();
}

void TreeDedispersion::prepare(DataBuffer<float> &databuffer)
{
	tsamp = databuffer.tsamp;
	ndump = databuffer.nsamples;

	nchans_orig = databuffer.nchans;

	nchans = std::pow(2, std::ceil(std::log2(nchans_orig)));
	frequencies.resize(nchans);
	double fch1 = databuffer.frequencies.front();
	double foff = (databuffer.frequencies.back() - databuffer.frequencies.front()) / (nchans_orig - 1);

	for (size_t j=0; j<nchans; j++)
	{
		frequencies[j] = fch1 + j * foff;
	}

	fmax = *std::max_element(frequencies.begin(), frequencies.end());
	fmin = *std::min_element(frequencies.begin(), frequencies.end());

	maxdelayn = std::ceil(dmdelay(dms + nchans * ddm, fmax, fmin) / tsamp);
	maxdelayn = ceil(1.*maxdelayn/ndump)*ndump;

	nsamples = ndump + maxdelayn;

	if (if_alloc_dedata)
		dedata.resize(nchans * nsamples, 0.);
	if (if_alloc_bufferT)
		bufferT.resize(nchans * nsamples, 0.);
	if (if_alloc_buffer)
		buffer.resize(nsamples * nchans, 0.);
	
	resize_cache();

	size_t tmp = 1;
	while (tmp < nchans)
	{
		maxdepth++;
		if (tmp < nsubband) depthsub++;

		tmp = tmp << 1;
	}

	update_delay();
	update_map();

	offset = nsamples-ndump;

	ready = true;

	std::vector<std::pair<std::string, std::string>> meta = {
			{"dm start", std::to_string(dms)},
			{"dm step", std::to_string(ddm)},
			{"dm end", std::to_string(dms + nchans * ddm)},
			{"dm smearing", std::to_string(dmdelay(0.5*ddm, fmax, fmin))},
			{"maximum dm smearing in channel", std::to_string(dmdelay(dms + nchans * ddm, fmax, fmin) / nchans)},
			{"tsamp", std::to_string(tsamp)},
			{"offset of first sample", std::to_string(offset)},
			{"buffer size", std::to_string(nsamples)},
			{"dump size", std::to_string(ndump)},
			{"maximum depth", std::to_string(maxdepth)},
			{"subband depth", std::to_string(depthsub)}
		};
	format_logging("Dedispersion Info", meta);
}

void TreeDedispersion::update_delay()
{
	delayn.resize((maxdepth + 1) * nchans, 0);
	frequencies_sub.resize(nsubband, 0.);

	double freq;
	update_delay_rec(map, freq, 0, 0);
}

void TreeDedispersion::update_delay_rec(std::vector<size_t> &vdmid, double &freq, size_t depth, size_t ichan)
{
	if (depth == maxdepth)
	{
		vdmid.push_back(0);
		freq = frequencies[ichan];

		if ((1<<depth) == nsubband)
		{
			if (ichan == 0)
			{
				size_t ndm_sub = nchans >> depth;
				size_t ddm_sub = 1 << depth;
				mapsub.resize(ndm_sub, 0);
				for (size_t i=0; i<ndm_sub; i++)
				{
					mapsub[i] = vdmid[i] / ddm_sub;
				}
			}

			frequencies_sub[ichan] = freq;
		}

		return;
	}

	std::vector<size_t> vdmid0, vdmid1;
	double freq0, freq1;
	update_delay_rec(vdmid0, freq0, depth + 1, 2 * ichan);
	update_delay_rec(vdmid1, freq1, depth + 1, 2 * ichan + 1);

	size_t size0 = vdmid0.size();
	size_t size1 = vdmid1.size();
	size_t size = size0 + size1;

	vdmid.resize(size, 0);
	for (size_t i = 0; i < size0; i++)
	{
		vdmid[i] = vdmid0[i];
	}
	for (size_t i = 0; i < size1; i++)
	{
		vdmid[size0 + i] = vdmid1[i] + (1 << depth);
	}

	size_t ndm = (nchans >> (depth + 1));
	double tmpfh = std::max(freq0, freq1);
	double tmpfl = std::min(freq0, freq1);

	freq = tmpfh;

	for (size_t idm=0; idm<ndm; idm++)
	{
		double dm0 = dms + vdmid[idm] * ddm;
		double dm1 = dms + vdmid[ndm + idm] * ddm;

		delayn[depth * nchans + ichan * (ndm * 2) + idm] = std::round(dmdelay(dm0, tmpfh, tmpfl) / tsamp);
		delayn[depth * nchans + ichan * (ndm * 2) + ndm + idm] = std::round(dmdelay(dm1, tmpfh, tmpfl) / tsamp);
	}

	if ((1<<depth) == nsubband)
	{
		if (ichan == 0)
		{
			size_t ndm_sub = nchans >> depth;
			size_t ddm_sub = 1 << depth;
			mapsub.resize(ndm_sub, 0);
			for (size_t i=0; i<ndm_sub; i++)
			{
				mapsub[i] = vdmid[i] / ddm_sub;
			}
		}

		frequencies_sub[ichan] = freq;
	}
}

void TreeDedispersion::update_map()
{
	std::vector<size_t> map_copy(map.begin(), map.end());
	for (size_t i=0; i<map_copy.size(); i++)
	{
		map[map_copy[i]] = i;
	}

	std::vector<size_t> mapsub_copy(mapsub.begin(), mapsub.end());
	for (size_t i=0; i<mapsub_copy.size(); i++)
	{
		mapsub[mapsub_copy[i]] = i;
	}
}

void TreeDedispersion::transform(size_t depth, size_t ichan)
{
#ifndef __AVX2__
	std::vector<float> *temp = bufferT.empty() ? ptr_bufferT : &bufferT;
#else
	std::vector<float, boost::alignment::aligned_allocator<float, 32>> *temp = bufferT.empty() ? ptr_bufferT : &bufferT;
#endif

	if (depth == maxdepth)
	{
		return;
	}

	transform(depth + 1, 2 * ichan);
	transform(depth + 1, 2 * ichan + 1);

	size_t ndm = (nchans >> (depth + 1));
	if (frequencies.front() > frequencies.back())
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (size_t idm=0; idm<ndm; idm++)
		{
			int thread_id = omp_get_thread_num();

			if (hit[depth * nchans + ichan * (ndm * 2) + idm])
			{
				size_t delayn0 = delayn[depth * nchans + ichan * (ndm * 2) + idm];
				std::rotate_copy(temp->begin() + (2 * ichan + 1) * ndm * nsamples + idm * nsamples, temp->begin() + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + delayn0, temp->begin() + (2 * ichan + 1) * ndm * nsamples + (idm + 1) * nsamples, cache0.begin() + thread_id * nsamples);
			}

			if (hit[depth * nchans + ichan * (ndm * 2) + ndm + idm])
			{
				size_t delayn1 = delayn[depth * nchans + ichan * (ndm * 2) + ndm + idm];
				std::rotate_copy(temp->begin() + (2 * ichan + 1) * ndm * nsamples + idm * nsamples, temp->begin() + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + delayn1, temp->begin() + (2 * ichan + 1) * ndm * nsamples + (idm + 1) * nsamples, cache1.begin() + thread_id * nsamples);
			}

			if (hit[depth * nchans + ichan * (ndm * 2) + ndm + idm])
			{
				for (size_t i=0; i<nsamples; i++)
				{
					(*temp)[(2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] = (*temp)[(2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache1[i + thread_id * nsamples];
				}
			}

			if (hit[depth * nchans + ichan * (ndm * 2) + idm])
			{
				for (size_t i=0; i<nsamples; i++)
				{
					(*temp)[(2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] = (*temp)[(2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache0[i + thread_id * nsamples];
				}
			}
		}
	}
	else
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (size_t idm=0; idm<ndm; idm++)
		{
			int thread_id = omp_get_thread_num();

			if (hit[depth * nchans + ichan * (ndm * 2) + idm])
			{
				size_t delayn0 = delayn[depth * nchans + ichan * (ndm * 2) + idm];
				std::rotate_copy(temp->begin() + (2 * ichan + 0) * ndm * nsamples + idm * nsamples, temp->begin() + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + delayn0, temp->begin() + (2 * ichan + 0) * ndm * nsamples + (idm + 1) * nsamples, cache0.begin() + thread_id * nsamples);
			}

			if (hit[depth * nchans + ichan * (ndm * 2) + ndm + idm])
			{
				size_t delayn1 = delayn[depth * nchans + ichan * (ndm * 2) + ndm + idm];
				std::rotate_copy(temp->begin() + (2 * ichan + 0) * ndm * nsamples + idm * nsamples, temp->begin() + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + delayn1, temp->begin() + (2 * ichan + 0) * ndm * nsamples + (idm + 1) * nsamples, cache1.begin() + thread_id * nsamples);
			}

			if (hit[depth * nchans + ichan * (ndm * 2) + idm])
			{
				for (size_t i=0; i<nsamples; i++)
				{
					(*temp)[(2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] = (*temp)[(2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] + cache0[i + thread_id * nsamples];
				}
			}

			if (hit[depth * nchans + ichan * (ndm * 2) + ndm + idm])
			{
				for (size_t i=0; i<nsamples; i++)
				{
					(*temp)[(2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] = (*temp)[(2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] + cache1[i + thread_id * nsamples];
				}
			}
		}
	}
}

void TreeDedispersion::run(DataBuffer<float> &databuffer)
{
	size_t nspace = nsamples - ndump;

	for (size_t i=0; i<ndump; i++)
	{
		for (size_t j=0; j<nchans_orig; j++)
		{
			buffer[(i + nspace) * nchans + j] = databuffer.buffer[i * nchans_orig + j];
		}
	}

	transpose_pad<float>(bufferT.data(), buffer.data(), nsamples, nchans);

	for (size_t j=0; j<nsubband; j++)
	{
		transform(depthsub, j);
	}

	transpose_pad<float>(dedata.data(), bufferT.data(), nsubband, nchans / nsubband * nsamples);

	for (size_t i=0; i<nspace; i++)
	{
		for (size_t j=0; j<nchans_orig; j++)
		{
			buffer[i * nchans + j] = buffer[(i + ndump) * nchans + j];
		}
	}

	counter += ndump;
}

void TreeDedispersion::run()
{
	transpose_pad<float>(ptr_bufferT->data(), ptr_buffer->data(), nsamples, nchans);

	for (size_t j=0; j<nsubband; j++)
	{
		transform(depthsub, j);
	}

	transpose_pad<float>(ptr_dedata->data(), ptr_bufferT->data(), nsubband, nchans / nsubband * nsamples);

	counter += ndump;
}

void TreeDedispersion::get_subdata(double dm, DataBuffer<float> &subdata)
{
	std::vector<float> *temp = dedata.empty() ? ptr_dedata : &dedata;

	subdata.resize(ndump, nsubband);
	subdata.tsamp = tsamp;
	subdata.frequencies = frequencies_sub;

	double ddm_sub = ddm * nsubband;
	size_t dmid = (dm - dms) / ddm_sub;

	size_t k = mapsub[dmid];

	std::copy(temp->begin() + k * nsamples * nsubband, temp->begin() + (k * nsamples + ndump) * nsubband, subdata.buffer.begin());

	subdata.counter += ndump;
}

/* =======================================================================================================================================*/

DedispersionX::DedispersionX()
{
	nsubband = 0;
	dm_boost = 0.;
	ddm_init = 0.;

	nchans = 0;
	nchans_orig = 0;
	nsamples0 = 0;
}

DedispersionX::~DedispersionX()
{
}

void DedispersionX::prepare(DataBuffer<float> &databuffer)
{
	double tsamp = databuffer.tsamp;
	nchans_orig = databuffer.nchans;
	nchans = std::pow(2, std::ceil(std::log2(nchans_orig)));

	std::vector<double> frequencies(nchans, 0.);
	double fch1 = databuffer.frequencies.front();
	double foff = (databuffer.frequencies.back() - databuffer.frequencies.front()) / (nchans_orig - 1);

	for (size_t j=0; j<nchans; j++)
	{
		frequencies[j] = fch1 + j * foff;
	}

	double fmax = *std::max_element(frequencies.begin(), frequencies.end());
	double fmin = *std::min_element(frequencies.begin(), frequencies.end());

	double ddm = ddm_init;
	if (ddm == 0.)
		ddm = tsamp / TreeDedispersion::dmdelay(1., fmax, fmin);

	ddm_init = ddm;

	double dms = 0.;

	std::list<double> dmlist_sort = dmlist;
	dmlist_sort.sort();

	double maxdm = dmlist_sort.back();

	nsamples0 = 0;

	DataBuffer<float> *d = databuffer.get();
	size_t ds = 1;
	do
	{
		double dme = dms + nchans * ddm;

		Downsample downsample(ds, 1);
		downsample.prepare(*d);
		if (ds == 1)
		{
			downsample.closable = true;
			downsample.close();
		}

		TreeDedispersion treededispersion;
		treededispersion.nsubband = nsubband;
		treededispersion.dms = dms;
		treededispersion.ddm = ddm;
		treededispersion.if_alloc_buffer = false;
		treededispersion.if_alloc_bufferT = false;
		treededispersion.if_alloc_dedata = false;

		if (dmlist_sort.front() < dme)
		{
			std::list<double>::iterator begin = dmlist_sort.begin();
			std::list<double>::iterator end = dmlist_sort.begin();
			while ((end != dmlist_sort.end()) && (*end < dme))
			{
				std::advance(end, 1);
			}

			treededispersion.prepare(downsample);
			treededispersion.update_hit(begin, end);

			hit.push_back(true);

			dmlist_sort.erase(begin, end);
		}
		else
		{
			hit.push_back(false);
		}

		if (ds == 1)
		{
			nsamples0 = treededispersion.get_nsamples() > nsamples0 ? treededispersion.get_nsamples() : nsamples0;
		}

		downsamples.push_back(downsample);
		treededispersions.push_back(treededispersion);
		
		d = downsamples.back().get();

		dms += nchans * ddm;

		if (dms >= dm_boost)
		{
			ddm *= 2;
			ds = 2;
		}
	} while (dms <= maxdm);

	// allocate memory for treededispersions
	buffers.resize(treededispersions.size());
	bufferTs.resize(treededispersions.size());
	dedatas.resize(treededispersions.size());
	
	for (size_t k=0; k<treededispersions.size(); k++)
	{
		if (downsamples[k].td == 1 && downsamples[k].fd == 1)
		{
			if (buffers[0].empty() && bufferTs[0].empty() && dedatas[0].empty())
			{
				buffers[0].resize(nsamples0 * nchans, 0.);
			}
			treededispersions[k].update_nsamples(nsamples0);
			treededispersions[k].resize_cache();
			treededispersions[k].ptr_buffer = &buffers[0];
			treededispersions[k].ptr_dedata = &dedatas[0];
			treededispersions[k].ptr_bufferT = &bufferTs[0];
		}
		else
		{
			buffers[k].resize(treededispersions[k].get_nsamples() * treededispersions[k].get_nchans(), 0.);

			treededispersions[k].ptr_buffer = &buffers[k];
			treededispersions[k].ptr_dedata = &dedatas[k];
			treededispersions[k].ptr_bufferT = &bufferTs[k];
		}
	}

	std::string hit_str;
	for (auto h=hit.begin(); h!=hit.end(); ++h) hit_str += *h ? "1 " : "0 ";
	std::vector<std::pair<std::string, std::string>> meta = {
			{"maximum frequency", std::to_string(fmax)},
			{"minimum frequency", std::to_string(fmin)},
			{"number of channels", std::to_string(nchans)},
			{"initial dm step", std::to_string(ddm_init)},
			{"number of dm ranges", std::to_string(downsamples.size())},
			{"dm range hit", hit_str},
		};
	format_logging("Dedispersion Outline", meta);
}

void DedispersionX::prerun(DataBuffer<float> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform data buffering";

	bool closable_bak = databuffer.closable;

	databuffer.closable = false;

	DataBuffer<float> *d = databuffer.get();
	for (size_t k=0; k<downsamples.size(); k++)
	{
		d = downsamples[k].run(*d);
	}

	bool enable = true;
	for (size_t k=0; k<downsamples.size(); k++)
	{
		if (treededispersions[k].is_ready())
		{
			if (downsamples[k].td == 1 && downsamples[k].fd == 1)
			{
				if (enable)
				{
					std::vector<float> &buffer = buffers[0];
					size_t nsamples = buffer.size() / nchans;
					size_t ndump = treededispersions[k].get_ndump();
					size_t nspace = nsamples - ndump;

					for (size_t i=0; i<ndump; i++)
					{
						for (size_t j=0; j<nchans_orig; j++)
						{
							buffer[(i + nspace) * nchans + j] = databuffer.buffer[i * nchans_orig + j];
						}
					}

					enable = false;
				}
			}
			else
			{
				std::vector<float> &buffer = buffers[k];
				size_t nsamples = buffer.size() / nchans;
				size_t ndump = treededispersions[k].get_ndump();
				size_t nspace = nsamples - ndump;

				for (size_t i=0; i<ndump; i++)
				{
					for (size_t j=0; j<nchans_orig; j++)
					{
						buffer[(i + nspace) * nchans + j] = downsamples[k].buffer[i * nchans_orig + j];
					}
				}
			}
		}
	}

	databuffer.closable = closable_bak;

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}

void DedispersionX::postrun(DataBuffer<float> &databuffer)
{
	bool enable = true;
	for (size_t k=0; k<downsamples.size(); k++)
	{
		if (treededispersions[k].is_ready())
		{
			if (downsamples[k].td == 1 && downsamples[k].fd == 1)
			{
				if (enable)
				{
					std::vector<float> &buffer = buffers[0];
					size_t nsamples = buffer.size() / nchans;
					size_t ndump = treededispersions[k].get_ndump();
					size_t nspace = nsamples - ndump;
					for (size_t i=0; i<nspace; i++)
					{
						for (size_t j=0; j<nchans_orig; j++)
						{
							buffer[i * nchans + j] = buffer[(i + ndump) * nchans + j];
						}
					}

					enable = false;
				}
			}
			else
			{
				std::vector<float> &buffer = buffers[k];
				size_t nsamples = buffer.size() / nchans;
				size_t ndump = treededispersions[k].get_ndump();
				size_t nspace = nsamples - ndump;
				for (size_t i=0; i<nspace; i++)
				{
					for (size_t j=0; j<nchans_orig; j++)
					{
						buffer[i * nchans + j] = buffer[(i + ndump) * nchans + j];
					}
				}
			}
		}
	}

	if (databuffer.closable) databuffer.close();
}

void DedispersionX::run(DataBuffer<float> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform dedispersion";

	bool closable_bak = databuffer.closable;

	databuffer.closable = false;

	DataBuffer<float> *d = databuffer.get();
	for (size_t k=0; k<downsamples.size(); k++)
	{
		d = downsamples[k].run(*d);
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
	for (size_t k=0; k<downsamples.size(); k++)
	{
		if (treededispersions[k].is_ready())
		{
			if (downsamples[k].td == 1 && downsamples[k].fd == 1)
			{
				treededispersions[k].run(databuffer);
			}
			else
			{
				treededispersions[k].run(downsamples[k]);
			}
		}
	}

	databuffer.closable = closable_bak;

	if (databuffer.closable) databuffer.close();

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}

void DedispersionX::run(size_t k)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform dedispersion";

	if (treededispersions[k].is_ready())
	{
		treededispersions[k].run();
	}

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}