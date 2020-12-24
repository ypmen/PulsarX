/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 18:44:33
 * @modify date 2020-06-07 18:44:33
 * @desc [description]
 */

#include <string.h>

#include "dedisperse.h"
#include "archivelite.h"
#include "constants.h"

extern "C" {
void sgesv_( int* n, int* nrhs, float* a, int* lda, int* ipiv,
                float* b, int* ldb, int* info );
}

using namespace Pulsar;

#define NSBLK 1024

ArchiveLite::ArchiveLite()
{
    start_mjd = 0.;
    sub_mjd = 0.;
    f0 = 0.;
    f1 = 0.;
    acc = 0.;
    dm = 0.;
    snr = 0.;
    nbin = 0;
    nchan = 0;
    npol = 0;
    tbin = 0.;
}

ArchiveLite::ArchiveLite(const ArchiveLite &arch)
{
    start_mjd = arch.start_mjd;
    ref_epoch = arch.ref_epoch;
    f0 = arch.f0;
    f1 = arch.f1;
    acc = arch.acc;
    dm = arch.dm;
    snr = arch.snr;
    nbin = arch.nbin;
    nchan = arch.nchan;
    npol = arch.npol;
    tbin = arch.tbin;
    frequencies = arch.frequencies;    
    profiles = arch.profiles;
}

ArchiveLite & ArchiveLite::operator=(const ArchiveLite &arch)
{
    start_mjd = arch.start_mjd;
    ref_epoch = arch.ref_epoch;
    f0 = arch.f0;
    f1 = arch.f1;
    acc = arch.acc;
    dm = arch.dm;
    snr = arch.snr;
    nbin = arch.nbin;
    nchan = arch.nchan;
    npol = arch.npol;
    tbin = arch.tbin;
    frequencies = arch.frequencies;    
    profiles = arch.profiles;

    return *this;
}

ArchiveLite::~ArchiveLite(){}

void ArchiveLite::resize(int np, int nc, int nb)
{
    npol = np;
    nchan = nc;
    nbin = nb;
    frequencies.resize(nchan, 0.);

    sub_int.resize(npol, nchan, nbin);
}

void ArchiveLite::prepare(DataBuffer<float> &databuffer)
{
    frequencies = databuffer.frequencies;

    sub_mjd = start_mjd + 0.5*databuffer.tsamp;

    if (f1 == 0)
    {
        f1 = acc/CONST_C*f0;
    }

    if (acc == 0)
    {
        acc = f1/f0*CONST_C;
    }
}

bool ArchiveLite::runDspsr(DataBuffer<float> &databuffer)
{
    if (databuffer.counter <= 0)
        return false;

    sub_int.tsubint = databuffer.nsamples*databuffer.tsamp;
    MJD start_time = sub_mjd;
    MJD end_time = sub_mjd + (databuffer.nsamples-1)*databuffer.tsamp;
    MJD epoch = get_epoch(start_time, end_time, ref_epoch);
    sub_int.offs_sub = (epoch-start_mjd).to_second();
    
    double phi = 0.;
    double f = 0.;

    vector<int> hits(nbin, 0);
    vector<float> profilesTPF(nbin*npol*nchan, 0.);
    vector<float> profilesPFT(npol*nchan*nbin, 0.);

    fill(hits.begin(), hits.end(), 0);
    vector<int> binplan(databuffer.nsamples);
    for (long int i=0; i<databuffer.nsamples; i++)
    {
        if (i%NSBLK == 0)
        {
            phi = get_phase(sub_mjd+i*databuffer.tsamp, ref_epoch);
            f = get_ffold(sub_mjd+(i+NSBLK*0.5-0.5)*databuffer.tsamp, ref_epoch);
        }

        phi -= floor(phi);
		int ibin = int(phi*nbin);
        binplan[i] = ibin;
        hits[ibin]++;
        phi += f*databuffer.tsamp;
    }

    fill(profilesTPF.begin(), profilesTPF.end(), 0.);
    for (long int i=0; i<databuffer.nsamples; i++)
    {
        for (long int j=0; j<databuffer.nchans; j++)
        {
            profilesTPF[binplan[i]*databuffer.nchans+j] += databuffer.buffer[i*databuffer.nchans+j];
        }
    }

    transpose_pad<float>(&profilesPFT[0], &profilesTPF[0], nbin, npol*nchan);

    for (long int ipol=0; ipol<npol; ipol++)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
        for (long int ichan=0; ichan<nchan; ichan++)
        {
            double mean = 0;
            long int count = 0;
            for (long int ibin=0; ibin<nbin; ibin++)
            {
                sub_int.data[ipol*nchan*nbin+ichan*nbin+ibin] = profilesPFT[ipol*nchan*nbin+ichan*nbin+ibin];
                if (hits[ibin] != 0)
                {
                    sub_int.data[ipol*nchan*nbin+ichan*nbin+ibin] /= hits[ibin];
                    mean += sub_int.data[ipol*nchan*nbin+ichan*nbin+ibin];
                    count++;
                }
            }
            mean /= count;
            for (long int ibin=0; ibin<nbin; ibin++)
            {
                if (hits[ibin] == 0)
                {
                    sub_int.data[ipol*nchan*nbin+ichan*nbin+ibin] = mean;
                }
            }
        }
    }

    sub_int.ffold = get_ffold(epoch, ref_epoch);
    
    sub_mjd += sub_int.tsubint;

    profiles.push_back(sub_int);

    return true;
}

bool ArchiveLite::runTRLSM(DataBuffer<float> &databuffer)
{
    if (databuffer.counter <= 0)
        return false;

    sub_int.tsubint = databuffer.nsamples*databuffer.tsamp;
    MJD start_time = sub_mjd;
    MJD end_time = sub_mjd + (databuffer.nsamples-1)*databuffer.tsamp;
    MJD epoch = get_epoch(start_time, end_time, ref_epoch);
    sub_int.offs_sub = (epoch-start_mjd).to_second();
    sub_int.ffold = abs(get_ffold(epoch, ref_epoch));

    double phi = 0.;
    double f = 0.;

    vector<float> mxWTW(nbin*nbin, 0.);
    vector<float> vWTd_T(nbin*databuffer.nchans, 0.);

    for (long int i=0; i<databuffer.nsamples; i++)
    {
        if (i%NSBLK == 0)
        {
            phi = get_phase(sub_mjd+(i-0.5)*databuffer.tsamp, ref_epoch);
            f = get_ffold(sub_mjd+(i+NSBLK*0.5-0.5)*databuffer.tsamp, ref_epoch);
        }

        double low_phi = phi;
        double high_phi = phi + f*databuffer.tsamp;
        phi = high_phi;

        if (low_phi > high_phi)
        {
            double tmp = low_phi;
            low_phi = high_phi;
            high_phi = tmp;
        }

        long int low_phin = floor(low_phi*nbin);
        long int high_phin = floor(high_phi*nbin);
        long int nphi = high_phin-low_phin+1;

        assert(nphi<=nbin);

        if (nphi == 1)
        {
            float vWli0 = 1.;

            vWli0 = (high_phi-low_phi)*nbin;

            long int l=low_phin%nbin;
            l = l<0 ? l+nbin:l;
            {
                mxWTW[l*nbin+l] += vWli0*vWli0;
            }

            for (long int j=0; j<databuffer.nchans; j++)
            {       
                vWTd_T[l*databuffer.nchans+j] += vWli0*databuffer.buffer[i*databuffer.nchans+j];
            }
        }
        else if (nphi == 2)
        {
            float vWli0=1., vWli1=1.;

            vWli0 = 1.-(low_phi*nbin-floor(low_phi*nbin));
            vWli1 = high_phi*nbin-floor(high_phi*nbin);

            long int l=low_phin%nbin;
            l = l<0 ? l+nbin:l;
            long int m = high_phin%nbin;
            m = m<0 ? m+nbin:m;

            mxWTW[l*nbin+l] += vWli0*vWli0;
            mxWTW[l*nbin+m] += vWli0*vWli1;
            mxWTW[m*nbin+l] += vWli1*vWli0;
            mxWTW[m*nbin+m] += vWli1*vWli1;

            for (long int j=0; j<databuffer.nchans; j++)
            {
                vWTd_T[l*databuffer.nchans+j] += vWli0*databuffer.buffer[i*databuffer.nchans+j];
                vWTd_T[m*databuffer.nchans+j] += vWli1*databuffer.buffer[i*databuffer.nchans+j];
            }
        }
        else
        {
            vector<float> vWli(nphi, 1.);
            vector<int> binplan(nphi, 0);

            vWli[0] = 1.-(low_phi*nbin-floor(low_phi*nbin));
            vWli[nphi-1] = high_phi*nbin-floor(high_phi*nbin);

            for (long int l=0; l<nphi; l++)
            {
                binplan[l] = (low_phin+l)%nbin;
                binplan[l] = binplan[l]<0 ? binplan[l]+nbin:binplan[l];
            }

            for (long int l=0; l<nphi; l++)
            {
                for (long int m=0; m<nphi; m++)
                {
                    mxWTW[binplan[l]*nbin+binplan[m]] += vWli[l]*vWli[m];
                }
            }

            for (long int l=0; l<nphi; l++)
            {
                for (long int j=0; j<databuffer.nchans; j++)
                {
                    vWTd_T[binplan[l]*databuffer.nchans+j] += vWli[l]*databuffer.buffer[i*databuffer.nchans+j];
                }
            }
        }
    }

    for (long int l=0; l<nbin; l++)
    {
        for (long int m=0; m<nbin; m++)
        {
            mxWTW[l*nbin+m] /= (databuffer.nsamples*databuffer.tsamp*sub_int.ffold);
        }
        mxWTW[l*nbin+l] += 1;
    }

    transpose_pad<float>(&sub_int.data[0], &vWTd_T[0], nbin, npol*nchan);

    int n = nbin;
    int nrhs = databuffer.nchans;
    vector<int> ipiv(n);
    int info;

    sgesv_(&n, &nrhs, &mxWTW[0], &n, &ipiv[0], &sub_int.data[0], &n, &info);
    
    sub_mjd += sub_int.tsubint;

    profiles.push_back(sub_int);

    return true;
}
