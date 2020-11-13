/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 18:44:33
 * @modify date 2020-06-07 18:44:33
 * @desc [description]
 */

#include "dedisperse.h"
#include "archivelite.h"

using namespace Pulsar;

ArchiveLite::ArchiveLite()
{
    start_mjd = 0.;
    sub_mjd = 0.;
    f0 = 0.;
    f1 = 0.;
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
    f0 = arch.f0;
    f1 = arch.f1;
    dm = arch.dm;
    snr = arch.snr;
    nbin = arch.nbin;
    nchan = arch.nchan;
    npol = arch.npol;
    tbin = arch.tbin;
    hits = arch.hits;
    profilesTPF = arch.profilesTPF;
    profilesPFT = arch.profilesPFT;
    frequencies = arch.frequencies;    
    profiles = arch.profiles;
}

ArchiveLite & ArchiveLite::operator=(const ArchiveLite &arch)
{
    start_mjd = arch.start_mjd;
    f0 = arch.f0;
    f1 = arch.f1;
    dm = arch.dm;
    snr = arch.snr;
    nbin = arch.nbin;
    nchan = arch.nchan;
    npol = arch.npol;
    tbin = arch.tbin;
    hits = arch.hits;
    profilesTPF = arch.profilesTPF;
    profilesPFT = arch.profilesPFT;
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
    hits.resize(nbin, 0);
    profilesTPF.resize(nbin*npol*nchan, 0.);
    profilesPFT.resize(npol*nchan*nbin, 0.);

    sub_int.resize(npol, nchan, nbin);
}

void ArchiveLite::prepare(DataBuffer<float> &databuffer)
{
    frequencies = databuffer.frequencies;

    sub_mjd = start_mjd + 0.5*databuffer.tsamp;
}

bool ArchiveLite::run(DataBuffer<float> &databuffer)
{
    if (databuffer.counter <= 0)
        return false;

    sub_int.tsubint = databuffer.nsamples*databuffer.tsamp;
    MJD start_time = sub_mjd;
    MJD end_time = sub_mjd + (databuffer.nsamples-1)*databuffer.tsamp;
    MJD epoch = get_epoch(start_time, end_time);
    sub_int.offs_sub = (epoch-start_mjd).to_second();
    
    double phi = get_phase(sub_mjd);
    double f = get_ffold(epoch, ref_epoch);

    fill(hits.begin(), hits.end(), 0);
    vector<int> binplan(databuffer.nsamples);
    for (long int i=0; i<databuffer.nsamples; i++)
    {
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

    transpose_pad(&profilesPFT[0], &profilesTPF[0], nbin, npol*nchan);

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

    sub_int.ffold = f;
    
    sub_mjd += sub_int.tsubint;

    profiles.push_back(sub_int);

    return true;
}
