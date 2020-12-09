/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 17:37:38
 * @modify date 2020-06-07 17:37:38
 * @desc [description]
 */

#include "dedispersionlite.h"
#include "dedisperse.h"

using namespace Pulsar;

DedispersionLite::DedispersionLite()
{
    nsubband = 0;
    counter = 0;
    offset = 0;
    nchans = 0;
    tsamp = 0.;
    ndump = 0;
    nsamples = 0;
}

DedispersionLite::~DedispersionLite(){}

void DedispersionLite::prepare(DataBuffer<float> &databuffer)
{
    assert(nsubband>0);
    assert(!vdm.empty());

    ndump = databuffer.nsamples;
    tsamp = databuffer.tsamp;
    nchans = databuffer.nchans;
    frequencies = databuffer.frequencies;
    
    int nch = ceil(nchans/nsubband);
    frequencies_sub.resize(nsubband, 0.);

    buffertim.resize(vdm.size()*ndump, 0.);
    buffersub.resize(vdm.size()*ndump*nsubband, 0.);
    buffersubT.resize(nsubband*vdm.size()*ndump, 0.);

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
    long int maxdelayn = ceil(dmdelay(*max_element(vdm.begin(), vdm.end()), fmax, fmin)/tsamp);
    nsamples = ceil(1.*maxdelayn/ndump)*ndump + ndump;
    buffer.resize(nsamples*nchans, 0.);
    bufferT.resize(nchans*nsamples, 0.);

    delayn.resize(vdm.size()*nchans, 0);
    int ndm = vdm.size();
    for (long int j=0; j<nchans; j++)
    {
        for (long int k=0; k<ndm; k++)
        {
            delayn[j*ndm+k] = round(dmdelay(vdm[k], fmax, frequencies[j])/tsamp);
        }
    }

    offset = nsamples-ndump;
}

void DedispersionLite::run(DataBuffer<float> &databuffer)
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

    transpose_pad<float>(&bufferT[0], &buffer[0], nsamples, nchans);

    int ndm = vdm.size();
    int nch = ceil(nchans/nsubband);

    fill(buffersubT.begin(), buffersubT.end(), 0.);
    fill(buffertim.begin(), buffertim.end(), 0.);

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
                buffertim[k*ndump+i] += bufferT[j*nsamples+i+delayn[j*ndm+k]];
            }
        }
    }
    
    transpose_pad<float>(&buffersub[0], &buffersubT[0], nsubband, ndm*ndump);

    for (long int i=0; i<nspace; i++)
    {
        for (long int j=0; j<nchans; j++)
        {
            buffer[i*nchans+j] = buffer[(i+ndump)*nchans+j];
        }
    }

    counter += ndump;
}
