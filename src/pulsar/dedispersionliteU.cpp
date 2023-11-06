/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 17:37:38
 * @modify date 2020-06-07 17:37:38
 * @desc [description]
 */

#include "dedispersionliteU.h"
#include "dedisperse.h"

using namespace Pulsar;

DedispersionLiteU::DedispersionLiteU()
{
	maxdm = 0.;
	groupsize = 0;

	nsubband = 0;
	counter = 0;
	offset = 0;
	nchans = 0;
	tsamp = 0.;
	ndump = 0;
	nsamples = 0;

	mean_var_ready = false;
}

DedispersionLiteU::~DedispersionLiteU(){}

void DedispersionLiteU::prepare(DataBuffer<float> &databuffer)
{
	assert(nsubband>0);

	ndump = databuffer.nsamples;
	tsamp = databuffer.tsamp;
	nchans = databuffer.nchans;
	frequencies = databuffer.frequencies;
	
	int nch = ceil(nchans*1./nsubband);
	nsubband = ceil(nchans*1./nch);
	frequencies_sub.resize(nsubband, 0.);

	means.resize(nsubband, 0.);
	vars.resize(nsubband, 0.);

	vector<int> fcnt(nsubband, 0);
	for (long int j=0; j<nchans; j++)
	{
		frequencies_sub[j/nch] += frequencies[j];
		fcnt[j/nch] += 1;
	}
	for (long int j=0; j<nsubband; j++)
	{
		frequencies_sub[j] /= fcnt[j];
	}
	
	double fmin = 1e6;
	double fmax = 0.;
	for (long int j=0; j<nchans; j++)
	{
		fmax = frequencies[j]>fmax? frequencies[j]:fmax;
		fmin = frequencies[j]<fmin? frequencies[j]:fmin;
	}
	long int maxdelayn = ceil(dmdelay(maxdm, fmax, fmin)/tsamp);
	nsamples = ceil(1.*maxdelayn/ndump)*ndump + ndump;

	/* memory test */
	try
	{
		float *test = new float [(long int)2*ndump*nchans+(long int)2*nsamples*nchans+(long int)2*groupsize*ndump*nsubband];
		delete [] test;
	}
	catch (std::bad_alloc & exception) 
	{ 
		std::cerr<<"Error: out of memory "<<std::fixed<<std::setprecision(1)<<((long int)2*ndump*nchans+(long int)2*nsamples*nchans+(long int)2*groupsize*ndump*nsubband)*4./1024/1024/1024<<" GB required"<<std::endl; 
		exit(-1);
	}

	buffer.resize(nsamples*nchans, 0.);

	offset = nsamples-ndump;
}

void DedispersionLiteU::updatedm(const std::vector<double> &dms)
{
	if (vdm == dms) return;

	vdm.clear();
	buffersub.clear();
	delayn.clear();

	vdm = dms;

	buffersub.resize(vdm.size()*ndump*nsubband, 0.);

	delayn.resize(vdm.size()*nchans, 0);
	int ndm = vdm.size();
	for (long int j=0; j<nchans; j++)
	{
		for (long int k=0; k<ndm; k++)
		{
			delayn[j*ndm+k] = round(dmdelay(vdm[k], *std::max_element(frequencies.begin(), frequencies.end()), frequencies[j])/tsamp);
		}
	}
}

void DedispersionLiteU::prerun(DataBuffer<float> &databuffer)
{
	assert(databuffer.nsamples == ndump);

	int nspace = nsamples-ndump;

	for (long int i=0; i<ndump; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			buffer[(i+nspace)*nchans+j] = databuffer.buffer[i*nchans+j];
		}
	}

	databuffer.isbusy = false;
	if (databuffer.closable) databuffer.close();

	bufferT.resize(nchans*nsamples, 0.);

	transpose_pad<float>(&bufferT[0], &buffer[0], nsamples, nchans);

	if (databuffer.mean_var_ready)
	{
		int nch = ceil(nchans*1./nsubband);
		std::fill(means.begin(), means.end(), 0.);
		std::fill(vars.begin(), vars.end(), 0.);
		for (long int j=0; j<nchans; j++)
		{
			means[j/nch] += databuffer.weights[j] * databuffer.means[j];
			vars[j/nch] += databuffer.weights[j] * databuffer.vars[j];
		}
		mean_var_ready = databuffer.mean_var_ready;
	}

	counter += ndump;
}

void DedispersionLiteU::run()
{
	int ndm = vdm.size();
	int nch = ceil(nchans*1./nsubband);

	vector<float> buffersubT(nsubband*vdm.size()*ndump, 0.);

	for (long int j=0; j<nchans; j++)
	{
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
		for (long int k=0; k<ndm; k++)
		{
			for (long int i=0; i<ndump; i++)
			{
				buffersubT[(j/nch)*ndm*ndump+k*ndump+i] += bufferT[j*nsamples+i+delayn[j*ndm+k]];
			}
		}
	}
	
	transpose_pad<float>(&buffersub[0], &buffersubT[0], nsubband, ndm*ndump);   
}

void DedispersionLiteU::postrun()
{
	int nspace = nsamples-ndump;

	for (long int i=0; i<nspace; i++)
	{
		for (long int j=0; j<nchans; j++)
		{
			buffer[i*nchans+j] = buffer[(i+ndump)*nchans+j];
		}
	}
}
