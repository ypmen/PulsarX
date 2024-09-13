/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-10-14 10:10:24
 * @modify date 2022-10-14 10:10:24
 * @desc [description]
 */

#include <iostream>

#include "logging.h"
#include "dermdm.h"
#include "utils.h"
#include "assert.h"

using namespace Pulsar;

DeRMDM::DeRMDM()
{
	dms = 0.;
	ddm = 0.;
	rms = 0.;
	drm = 0.;

	depthrm = 0;

	if_alloc_buffer = true;
	if_alloc_bufferT = true;
	if_alloc_dedata = true;

	numrm = 0;
	numrmsub = 0;

	ready = false;
	counter = 0;
	nchans = 0;
	ndump = 0;
	tsamp = 0.;
	maxdepth = 0;
	depthsub = 0;
	fmax = 0.;
	fmin = 0.;
	maxdelayn = 0;
	offset = 0;
	nsamples = 0;
}

DeRMDM::~DeRMDM()
{
}

void DeRMDM::close()
{
	dedata.clear();
	dedata.shrink_to_fit();

	buffer.clear();
	buffer.shrink_to_fit();

	bufferT.clear();
	bufferT.shrink_to_fit();

	mapdm.clear();
	mapdm.shrink_to_fit();

	maprm.clear();
	maprm.shrink_to_fit();

	mapsubdm.clear();
	mapsubdm.shrink_to_fit();

	mapsubrm.clear();
	mapsubrm.shrink_to_fit();

	dmdelayns.clear();
	rmdelays.shrink_to_fit();

	hitdm.clear();
	hitdm.shrink_to_fit();

	hitrm.clear();
	hitrm.shrink_to_fit();

	cache00.clear();
	cache00.shrink_to_fit();

	cache01.clear();
	cache01.shrink_to_fit();

	cache10.clear();
	cache10.shrink_to_fit();

	cache11.clear();
	cache11.shrink_to_fit();

}

void DeRMDM::prepare(DataBuffer<std::complex<float>> &databuffer)
{
	tsamp = databuffer.tsamp;
	ndump = databuffer.nsamples;

	nchans = databuffer.nchans;
	frequencies = databuffer.frequencies;

	fmax = *std::max_element(frequencies.begin(), frequencies.end());
	fmin = *std::min_element(frequencies.begin(), frequencies.end());

	maxdelayn = std::ceil(dmdelay(dms + nchans * ddm, fmax, fmin) / tsamp);
	maxdelayn = ceil(1.*maxdelayn/ndump)*ndump;

	nsamples = ndump + maxdelayn;

	numrm = 1 << depthrm;
	assert(numrm <= nchans);

	size_t tmp = 1;
	while (tmp < nchans)
	{
		maxdepth++;
		if (tmp < nsubband) depthsub++;

		tmp = tmp << 1;
	}

	assert(depthrm >= depthsub);
	numrmsub = 1 << (depthrm - depthsub);

	if (if_alloc_dedata)
		dedata.resize(numrmsub * nchans / nsubband * nsamples * nsubband, 0.);
	if (if_alloc_bufferT)
		bufferT.resize(numrmsub * nsubband * nchans / nsubband * nsamples, 0.);
	if (if_alloc_buffer)
		buffer.resize(nsamples * nchans, 0.);

	cache00.resize(nsamples);
	cache01.resize(nsamples);
	cache10.resize(nsamples);
	cache11.resize(nsamples);

	update_dmdelay();
	update_rmdelay();
	update_dmmap();
	update_rmmap();

	offset = nsamples-ndump;

	ready = true;

	std::vector<std::pair<std::string, std::string>> meta = {
			{"dm start", std::to_string(dms)},
			{"dm step", std::to_string(ddm)},
			{"dm end", std::to_string(dms + nchans * ddm)},
			{"dm smearing", std::to_string(dmdelay(0.5*ddm, fmax, fmin))},
			{"maximum dm smearing in channel", std::to_string(dmdelay(dms + nchans * ddm, fmax, fmin) / nchans)},
			{"rm start", std::to_string(rms)},
			{"rm step", std::to_string(drm)},
			{"rm end", std::to_string(rms + numrm * drm)},
			{"rm smearing", std::to_string(rmdelay(0.5*drm, fmax, fmin))},
			{"maximum rm smearing in channel", std::to_string(rmdelay(rms + numrm * drm, fmax, fmin) / nchans)},
			{"tsamp", std::to_string(tsamp)},
			{"offset of first sample", std::to_string(offset)},
			{"buffer size", std::to_string(nsamples)},
			{"dump size", std::to_string(ndump)},
			{"maximum depth", std::to_string(maxdepth)},
			{"subband depth", std::to_string(depthsub)},
			{"RM depth", std::to_string(depthrm)},
		};
	format_logging("DeRMDM Info", meta);
}

void DeRMDM::run(DataBuffer<std::complex<float>> &databuffer)
{
	size_t nspace = nsamples - ndump;

	for (size_t i=0; i<ndump; i++)
	{
		for (size_t j=0; j<nchans; j++)
		{
			buffer[(i + nspace) * nchans + j] = databuffer.buffer[i * nchans + j];
		}
	}

	transpose_pad<std::complex<float>>(bufferT.data(), buffer.data(), nsamples, nchans);

	for (size_t j=0; j<nsubband; j++)
	{
		double freq = 0.;
		transform(freq, depthsub, j);
	}

	for (size_t k=0; k<numrmsub; k++)
	{
		transpose_pad<std::complex<float>>(dedata.data() + k * nchans * nsamples, bufferT.data() + k * nchans * nsamples, nsubband, nchans / nsubband * nsamples);
	}

	for (size_t i=0; i<nspace; i++)
	{
		for (size_t j=0; j<nchans; j++)
		{
			buffer[i * nchans + j] = buffer[(i + ndump) * nchans + j];
		}
	}

	counter += ndump;
}

void DeRMDM::run()
{
	transpose_pad<std::complex<float>>(ptr_bufferT->data(), ptr_buffer->data(), nsamples, nchans);

	for (size_t j=0; j<nsubband; j++)
	{
		double freq = 0.;
		transform(freq, depthsub, j);
	}

	for (size_t k=0; k<numrmsub; k++)
	{
		transpose_pad<std::complex<float>>(ptr_dedata->data() + k * nchans * nsamples, ptr_bufferT->data() + k * nchans * nsamples, nsubband, nchans / nsubband * nsamples);
	}

	counter += ndump;
}

void DeRMDM::get_subdata(double dm, double rm, DataBuffer<std::complex<float>> &subdata)
{
	std::vector<std::complex<float>> *temp = dedata.empty() ? ptr_dedata : &dedata;

	subdata.resize(ndump, nsubband);
	subdata.tsamp = tsamp;
	subdata.frequencies = frequencies_sub;

	double ddm_sub = ddm * nsubband;
	size_t dmid = (dm - dms) / ddm_sub;

	size_t dmk = mapsubdm[dmid];

	double drm_sub = drm * nsubband;
	size_t rmid = (rm - rms) / drm_sub;

	size_t rmk = mapsubrm[rmid];

	std::copy(temp->begin() + rmk * nchans/nsubband * nsamples * nsubband + dmk * nsamples * nsubband, temp->begin() + rmk * nchans/nsubband * nsamples * nsubband + (dmk * nsamples + ndump) * nsubband, subdata.buffer.begin());

	subdata.counter += ndump;
}

void DeRMDM::update_dmdelay()
{
	dmdelayns.resize((maxdepth + 1) * nchans, 0);
	frequencies_sub.resize(nsubband, 0.);

	double freq;
	update_dmdelay_rec(mapdm, freq, 0, 0);
}

void DeRMDM::update_dmdelay_rec(std::vector<size_t> &vdmid, double &freq, size_t depth, size_t ichan)
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
				mapsubdm.resize(ndm_sub, 0);
				for (size_t i=0; i<ndm_sub; i++)
				{
					mapsubdm[i] = vdmid[i] / ddm_sub;
				}
			}

			frequencies_sub[ichan] = freq;
		}

		return;
	}

	std::vector<size_t> vdmid0, vdmid1;
	double freq0, freq1;
	update_dmdelay_rec(vdmid0, freq0, depth + 1, 2 * ichan);
	update_dmdelay_rec(vdmid1, freq1, depth + 1, 2 * ichan + 1);

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

		dmdelayns[depth * nchans + ichan * (ndm * 2) + idm] = std::round(dmdelay(dm0, tmpfh, tmpfl) / tsamp);
		dmdelayns[depth * nchans + ichan * (ndm * 2) + ndm + idm] = std::round(dmdelay(dm1, tmpfh, tmpfl) / tsamp);
	}

	if ((1<<depth) == nsubband)
	{
		if (ichan == 0)
		{
			size_t ndm_sub = nchans >> depth;
			size_t ddm_sub = 1 << depth;
			mapsubdm.resize(ndm_sub, 0);
			for (size_t i=0; i<ndm_sub; i++)
			{
				mapsubdm[i] = vdmid[i] / ddm_sub;
			}
		}

		frequencies_sub[ichan] = freq;
	}
}

void DeRMDM::update_dmmap()
{
	std::vector<size_t> map_copy(mapdm.begin(), mapdm.end());
	for (size_t i=0; i<map_copy.size(); i++)
	{
		mapdm[map_copy[i]] = i;
	}

	std::vector<size_t> mapsub_copy(mapsubdm.begin(), mapsubdm.end());
	for (size_t i=0; i<mapsub_copy.size(); i++)
	{
		mapsubdm[mapsub_copy[i]] = i;
	}
}

void DeRMDM::update_rmdelay()
{
	rmdelays.resize((depthrm + 1) * numrm, 0);

	double freq;
	update_rmdelay_rec(maprm, freq, 0, 0);
}

void DeRMDM::update_freq_rec(double &freq, size_t depth, size_t ichan)
{
	if (depth == maxdepth)
	{
		freq = frequencies[ichan];
		return;
	}

	double freq0, freq1;
	update_freq_rec(freq0, depth + 1, 2 * ichan);
	update_freq_rec(freq1, depth + 1, 2 * ichan + 1);

	freq = std::max(freq0, freq1);
}

void DeRMDM::update_rmdelay_rec(std::vector<size_t> &vrmid, double &freq, size_t depth, size_t ichan)
{
	if (depth == depthrm)
	{
		vrmid.push_back(0);
		update_freq_rec(freq, depth, ichan);

		if ((1<<depth) == nsubband)
		{
			if (ichan == 0)
			{
				size_t nrm_sub = numrm >> depth;
				size_t drm_sub = 1 << depth;
				mapsubrm.resize(nrm_sub, 0);
				for (size_t i=0; i<nrm_sub; i++)
				{
					mapsubrm[i] = vrmid[i] / drm_sub;
				}
			}
		}

		return;
	}

	std::vector<size_t> vrmid0, vrmid1;
	double freq0, freq1;
	update_rmdelay_rec(vrmid0, freq0, depth + 1, 2 * ichan);
	update_rmdelay_rec(vrmid1, freq1, depth + 1, 2 * ichan + 1);

	size_t size0 = vrmid0.size();
	size_t size1 = vrmid1.size();
	size_t size = size0 + size1;

	vrmid.resize(size, 0);
	for (size_t i = 0; i < size0; i++)
	{
		vrmid[i] = vrmid0[i];
	}
	for (size_t i = 0; i < size1; i++)
	{
		vrmid[size0 + i] = vrmid1[i] + (1 << depth);
	}

	size_t nrm = (numrm >> (depth + 1));
	double tmpfh = std::max(freq0, freq1);
	double tmpfl = std::min(freq0, freq1);

	freq = tmpfh;

	for (size_t irm=0; irm<nrm; irm++)
	{
		double rm0 = rms + vrmid[irm] * drm;
		double rm1 = rms + vrmid[nrm + irm] * drm;

		rmdelays[depth * numrm + ichan * (nrm * 2) + irm] = rmdelay(rm0, tmpfh, tmpfl);
		rmdelays[depth * numrm + ichan * (nrm * 2) + nrm + irm] = rmdelay(rm1, tmpfh, tmpfl);
	}

	if ((1<<depth) == nsubband)
	{
		if (ichan == 0)
		{
			size_t nrm_sub = numrm >> depth;
			size_t drm_sub = 1 << depth;
			mapsubrm.resize(nrm_sub, 0);
			for (size_t i=0; i<nrm_sub; i++)
			{
				mapsubrm[i] = vrmid[i] / drm_sub;
			}
		}
	}
}

void DeRMDM::update_rmmap()
{
	std::vector<size_t> map_copy(maprm.begin(), maprm.end());
	for (size_t i=0; i<map_copy.size(); i++)
	{
		maprm[map_copy[i]] = i;
	}

	std::vector<size_t> mapsub_copy(mapsubrm.begin(), mapsubrm.end());
	for (size_t i=0; i<mapsub_copy.size(); i++)
	{
		mapsubrm[mapsub_copy[i]] = i;
	}
}

void DeRMDM::transform(double &freq, size_t depth, size_t ichan)
{
#ifndef __AVX2__
	std::vector<std::complex<float>> *temp = bufferT.empty() ? ptr_bufferT : &bufferT;
#else
	std::vector<std::complex<float>, boost::alignment::aligned_allocator<std::complex<float>, 32>> *temp = bufferT.empty() ? ptr_bufferT : &bufferT;
#endif

	if (depth == maxdepth)
	{
		freq = frequencies[ichan];
		return;
	}

	double freq0, freq1;

	transform(freq0, depth + 1, 2 * ichan);
	transform(freq1, depth + 1, 2 * ichan + 1);

	double tmpfh = std::max(freq0, freq1);
	double tmpfl = std::min(freq0, freq1);

	freq = tmpfh;

	size_t ndm = nchans >> (depth + 1);
	size_t nrm = numrm >> (depth + 1);

	assert(frequencies.front() >= frequencies.back());

	if (nrm > 0)
	{
		for (size_t irm = 0; irm < nrm; irm++)
		{
			float rmdelay0 = rmdelays[depth * numrm + ichan * (2 * nrm) + irm];

			std::complex<float> rmdelayc0;
			rmdelayc0.real(std::cos(-rmdelay0));
			rmdelayc0.imag(std::sin(-rmdelay0));

			float rmdelay1 = rmdelays[depth * numrm + ichan * (2 * nrm) + nrm + irm];

			std::complex<float> rmdelayc1;
			rmdelayc1.real(std::cos(-rmdelay1));
			rmdelayc1.imag(std::sin(-rmdelay1));

			for (size_t idm = 0; idm < ndm; idm++)
			{
				// shift
				if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
				{
					size_t dmdelayn = dmdelayns[depth * nchans + ichan * (2 * ndm) + idm];

					std::rotate_copy(temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + dmdelayn, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + (idm + 1) * nsamples, cache00.begin());
				}

				if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
				{
					size_t dmdelayn = dmdelayns[depth * nchans + ichan * (2 * ndm) + ndm + idm];

					std::rotate_copy(temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + dmdelayn, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + (idm + 1) * nsamples, cache01.begin());
				}

				if (hitrm[depth * numrm + ichan * (2 * nrm) + nrm + irm])
				{
					if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
					{
						std::transform(cache00.begin(), cache00.end(), cache10.begin(), [&rmdelayc1](std::complex<float> &x){return x * rmdelayc1;});
					}

					if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
					{
						std::transform(cache01.begin(), cache01.end(), cache11.begin(), [&rmdelayc1](std::complex<float> &x){return x * rmdelayc1;});
					}
				}

				if (hitrm[depth * numrm + ichan * (2 * nrm) + irm])
				{
					if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
					{
						std::transform(cache00.begin(), cache00.end(), cache00.begin(), [&rmdelayc0](std::complex<float> &x){return x * rmdelayc0;});
					}

					if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
					{
						std::transform(cache01.begin(), cache01.end(), cache01.begin(), [&rmdelayc0](std::complex<float> &x){return x * rmdelayc0;});
					}
				}

				// summation

				if (hitrm[depth * numrm + ichan * (2 * nrm) + nrm + irm])
				{
					if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
					{
						for (size_t i=0; i<nsamples; i++)
						{
							(*temp)[(nrm + irm) * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] = (*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache11[i];
						}
					}

					if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
					{
						for (size_t i=0; i<nsamples; i++)
						{
							(*temp)[(nrm + irm) * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] = (*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache10[i];
						}
					}
				}

				if (hitrm[depth * numrm + ichan * (2 * nrm) + irm])
				{
					if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
					{
						for (size_t i=0; i<nsamples; i++)
						{
							(*temp)[irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] = (*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache01[i];
						}
					}

					if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
					{
						for (size_t i=0; i<nsamples; i++)
						{
							(*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] = (*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache00[i];
						}
					}
				}
			}
		}
	}
	else
	{
		float delayrm = rmdelay(rms, tmpfh, tmpfl);

		std::complex<float> rmdelayc;
		rmdelayc.real(std::cos(-delayrm));
		rmdelayc.imag(std::sin(-delayrm));

		size_t irm = 0;

		for (size_t idm = 0; idm < ndm; idm++)
		{
			if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
			{
				size_t dmdelayn = dmdelayns[depth * nchans + ichan * (2 * ndm) + idm];

				std::rotate_copy(temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + dmdelayn, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + (idm + 1) * nsamples, cache00.begin());
				std::transform(cache00.begin(), cache00.end(), cache00.begin(), [&rmdelayc](std::complex<float> &x){return x * rmdelayc;});
			}

			if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
			{
				size_t dmdelayn = dmdelayns[depth * nchans + ichan * (2 * ndm) + ndm + idm];

				std::rotate_copy(temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + dmdelayn, temp->begin() + irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + (idm + 1) * nsamples, cache01.begin());
				std::transform(cache01.begin(), cache01.end(), cache01.begin(), [&rmdelayc](std::complex<float> &x){return x * rmdelayc;});
			}

			if (hitdm[depth * nchans + ichan * (2 * ndm) + ndm + idm])
			{
				for (size_t i=0; i<nsamples; i++)
				{
					(*temp)[irm * nchans * nsamples + (2 * ichan + 1) * ndm * nsamples + idm * nsamples + i] = (*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache01[i];
				}
			}

			if (hitdm[depth * nchans + ichan * (2 * ndm) + idm])
			{
				for (size_t i=0; i<nsamples; i++)
				{
					(*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] = (*temp)[irm * nchans * nsamples + (2 * ichan + 0) * ndm * nsamples + idm * nsamples + i] + cache00[i];
				}
			}
		}
	}
}

//

DeRMDMX::DeRMDMX()
{
	nsubband = 0;
	dms = 0.;
	dme = 0.;
	rms = 0.;
	rme = 0.;
	rm_smear = 0.1;
	rm_boost = 0.;
	dm_boost = 0.;
	ddm_init = 0.;

	numrm = 0;
	numrmsub = 0;
	nchans = 0;
	nsamples0 = 0;
}

DeRMDMX::~DeRMDMX()
{

}

void DeRMDMX::prepare(DataBuffer<std::complex<float>> &databuffer)
{
	double tsamp = databuffer.tsamp;
	std::vector<double> frequencies = databuffer.frequencies;
	nchans = databuffer.nchans;

	assert((nchans & (nchans - 1)) == 0);

	double fmax = *std::max_element(frequencies.begin(), frequencies.end());
	double fmin = *std::min_element(frequencies.begin(), frequencies.end());

	double ddm = ddm_init;
	if (ddm == 0.)
		ddm = tsamp / DeRMDM::dmdelay(1., fmax, fmin);

	ddm_init = ddm;

	double drm = drm_init;
	if (drm == 0.)
		drm = 2. * M_PI * rm_smear / DeRMDM::rmdelay(1., fmax, fmin);

	numrm = 1 << depthrm;

	size_t depthsub = 0;
	size_t tmp = 1;
	while (tmp < nchans)
	{
		if (tmp < nsubband) depthsub++;

		tmp = tmp << 1;
	}

	assert(depthrm >= depthsub);
	numrmsub = 1 << (depthrm - depthsub);

	// calculate ddplan
	ddplan.clear();
	double dm_start = dms;
	double dm_step = ddm;
	double dm_end = dm_start + nchans * dm_step;
	size_t ds = 1;
	ddplan.push_back(std::make_tuple(dm_start, dm_step, nchans, ds));
	while (dm_end <= dme)
	{
		if (dm_end > dm_boost)
		{
			dm_step *= 2;
			ds *= 2;
		}

		dm_start = dm_end;
		dm_end = dm_start + nchans * dm_step;

		ddplan.push_back(std::make_tuple(dm_start, dm_step, nchans, ds));
	}
	// calculate drplan
	drplan.clear();
	double rm_start = rms;
	double rm_step = drm;
	double rm_end = rm_start + numrm * rm_step;
	drplan.push_back(std::make_tuple(rm_start, rm_step, numrm));
	for (size_t k = 0; k < nchans/numrm-1; k++)
	{
		if (rm_end > rme) break;

		rm_start = rm_end;
		rm_end = rm_start + numrm * rm_step;
		drplan.push_back(std::make_tuple(rm_start, rm_step, numrm));
	}
	while (rm_end <= rme)
	{
		if (rm_end > rme) break;

		if (rm_end > rm_boost) 
		{
			rm_step *= 2;
		}

		for (size_t k = 0; k < nchans/numrm; k++)
		{
			rm_start = rm_end;
			rm_end = rm_start + numrm * rm_step;
			drplan.push_back(std::make_tuple(rm_start, rm_step, numrm));

			if (rm_end > rme) break;
		}
	}

	// initial downsamples dermdms
	for (auto ddp = ddplan.begin(); ddp != ddplan.end(); ++ddp)
	{
		double dm_start = std::get<0>(*ddp);
		double dm_step = std::get<1>(*ddp);
		size_t ndm = std::get<2>(*ddp);
		size_t ds = std::get<3>(*ddp);

		cDownsample cdownsample(ds, 1);
		cdownsample.prepare(databuffer);

		if (ds == 1)
		{
			cdownsample.closable = true;
			cdownsample.close();
		}

		cdownsamples.push_back(cdownsample);

		for (auto drp = drplan.begin(); drp != drplan.end(); ++drp)
		{
			double rm_start = std::get<0>(*drp);
			double rm_step = std::get<1>(*drp);
			size_t nrm = std::get<2>(*drp);

			DeRMDM dermdm;
			dermdm.dms = dm_start;
			dermdm.ddm = dm_step;
			dermdm.rms = rm_start;
			dermdm.drm = rm_step;
			dermdm.nsubband = nsubband;
			dermdm.depthrm = depthrm;
			dermdm.if_alloc_buffer = false;
			dermdm.if_alloc_bufferT = false;
			dermdm.if_alloc_dedata = false;

			dermdm.prepare(cdownsample);
			dermdm.dmhit_all();
			dermdm.rmhit_all();

			if (ds == 1)
			{
				nsamples0 = dermdm.get_nsamples() > nsamples0 ? dermdm.get_nsamples() : nsamples0;
			}

			dermdms.push_back(dermdm);
		}
	}

	// allocate memory for treededispersions
	buffers.resize(cdownsamples.size());
	bufferTs.resize(cdownsamples.size());
	dedatas.resize(cdownsamples.size());

	for (size_t i=0; i<ddplan.size(); i++)
	{
		// alloc buffer memory
		if (cdownsamples[i].td == 1 && cdownsamples[i].fd == 1)
		{
			if (buffers[0].empty() && bufferTs[0].empty() && dedatas[0].empty())
			{
				buffers[0].resize(nsamples0 * nchans, 0.);
			}
		}
		else
		{
			buffers[i].resize(dermdms[i * drplan.size()].get_nsamples() * dermdms[i * drplan.size()].get_nchans(), 0.);
		}

		for (size_t j=0; j<drplan.size(); j++)
		{
			if (cdownsamples[i].td == 1 && cdownsamples[i].fd == 1)
			{
				dermdms[i * drplan.size() + j].update_nsamples(nsamples0);
				dermdms[i * drplan.size() + j].resize_cache();
				dermdms[i * drplan.size() + j].ptr_buffer = &buffers[0];
				dermdms[i * drplan.size() + j].ptr_dedata = &dedatas[0];
				dermdms[i * drplan.size() + j].ptr_bufferT = &bufferTs[0];
			}
			else
			{
				dermdms[i * drplan.size() + j].ptr_buffer = &buffers[i];
				dermdms[i * drplan.size() + j].ptr_dedata = &dedatas[i];
				dermdms[i * drplan.size() + j].ptr_bufferT = &bufferTs[i];
			}
		}
	}
	
	std::vector<std::pair<std::string, std::string>> meta = {
			{"maximum frequency", std::to_string(fmax)},
			{"minimum frequency", std::to_string(fmin)},
			{"number of channels", std::to_string(nchans)},
			{"initial dm step", std::to_string(ddm_init)},
			{"number of dm ranges", std::to_string(ddplan.size())},
			{"initial rm step", std::to_string(drm)},
			{"number of rm ranges", std::to_string(drplan.size())},
		};
	format_logging("DeRMDM Outline", meta);
}

void DeRMDMX::prerun(DataBuffer<std::complex<float>> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform data buffering";

	bool closable_bak = databuffer.closable;

	databuffer.closable = false;

	DataBuffer<std::complex<float>> *d = databuffer.get();
	for (size_t m=0; m<cdownsamples.size(); m++)
	{
		d = cdownsamples[m].run(*d);
	}

	bool enable = true;
	for (size_t m=0; m<cdownsamples.size(); m++)
	{
		bool is_ready = false;
		for (size_t l=0; l<drplan.size(); l++) is_ready = is_ready || dermdms[m * drplan.size() + l].is_ready();
		if (is_ready)
		{
			if (cdownsamples[m].td == 1 && cdownsamples[m].fd == 1)
			{
				if (enable)
				{
					std::vector<std::complex<float>> &buffer = buffers[0];
					size_t nsamples = buffer.size() / nchans;
					size_t ndump = dermdms[m * drplan.size()].get_ndump();
					size_t nspace = nsamples - ndump;
					for (size_t i=0; i<ndump; i++)
					{
						for (size_t j=0; j<nchans; j++)
						{
							buffer[(i + nspace) * nchans + j] = databuffer.buffer[i * nchans + j];
						}
					}

					enable = false;
				}
			}
			else
			{
				std::vector<std::complex<float>> &buffer = buffers[m];
				size_t nsamples = buffer.size() / nchans;
				size_t ndump = dermdms[m * drplan.size()].get_ndump();
				size_t nspace = nsamples - ndump;

				for (size_t i=0; i<ndump; i++)
				{
					for (size_t j=0; j<nchans; j++)
					{
						buffer[(i + nspace) * nchans + j] = cdownsamples[m].buffer[i * nchans + j];
					}
				}
			}
		}
	}

	databuffer.closable = closable_bak;

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}

void DeRMDMX::postrun(DataBuffer<std::complex<float>> &databuffer)
{
	bool enable = true;
	for (size_t m=0; m<cdownsamples.size(); m++)
	{
		bool is_ready = false;
		for (size_t l=0; l<drplan.size(); l++) is_ready = is_ready || dermdms[m * drplan.size() + l].is_ready();
		if (is_ready)
		{
			if (cdownsamples[m].td == 1 && cdownsamples[m].fd == 1)
			{
				if (enable)
				{
					std::vector<std::complex<float>> &buffer = buffers[0];
					size_t nsamples = buffer.size() / nchans;
					size_t ndump = dermdms[m * drplan.size()].get_ndump();
					size_t nspace = nsamples - ndump;

					for (size_t i=0; i<nspace; i++)
					{
						for (size_t j=0; j<nchans; j++)
						{
							buffer[i * nchans + j] = buffer[(i + ndump) * nchans + j];
						}
					}

					enable = false;
				}
			}
			else
			{
				std::vector<std::complex<float>> &buffer = buffers[m];
				size_t nsamples = buffer.size() / nchans;
				size_t ndump = dermdms[m * drplan.size()].get_ndump();
				size_t nspace = nsamples - ndump;

				for (size_t i=0; i<nspace; i++)
				{
					for (size_t j=0; j<nchans; j++)
					{
						buffer[i * nchans + j] = buffer[(i + ndump) * nchans + j];
					}
				}
			}
		}
	}

	if (databuffer.closable) databuffer.close();
}

void DeRMDMX::run(DataBuffer<std::complex<float>> &databuffer)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform defaraday and dedispersion";

	bool closable_bak = databuffer.closable;

	databuffer.closable = false;

	DataBuffer<std::complex<float>> *d = databuffer.get();
	for (size_t i=0; i<cdownsamples.size(); i++)
	{
		d = cdownsamples[i].run(*d);
	}

	for (size_t i=0; i<ddplan.size(); i++)
	{
		for (size_t j=0; j<drplan.size(); j++)
		{
			if (dermdms[i * drplan.size() + j].is_ready())
			{
				if (cdownsamples[i].td == 1 && cdownsamples[i].fd == 1)
				{
					dermdms[i * drplan.size() + j].run(databuffer);
				}
				else
				{
					dermdms[i * drplan.size() + j].run(cdownsamples[i]);
				}
			}
		}
	}

	databuffer.closable = closable_bak;

	if (databuffer.closable) databuffer.close();

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}

void DeRMDMX::run(size_t i, size_t j)
{
	BOOST_LOG_TRIVIAL(debug)<<"perform defaraday and dedispersion";

	if (dermdms[i * drplan.size() + j].is_ready())
	{
		dermdms[i * drplan.size() + j].run();
	}

	BOOST_LOG_TRIVIAL(debug)<<"finished";
}