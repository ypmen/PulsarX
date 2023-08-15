/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 18:44:33
 * @modify date 2020-06-07 18:44:33
 * @desc [description]
 */

#include <string.h>

#include "dedisperse.h"
#include "psrfits.h"
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
	nblock = 1;
	nbin = 0;
	nchan = 0;
	npol = 0;
	tbin = 0.;
	iblock = 0;

	fref = 1000.;
	use_t2pred = false;
	mean_var_ready = false;
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
	nblock = arch.nblock;
	nbin = arch.nbin;
	nchan = arch.nchan;
	npol = arch.npol;
	tbin = arch.tbin;
	frequencies = arch.frequencies;    
	profiles = arch.profiles;
	sub_mjd = arch.sub_mjd;
	sub_int = arch.sub_int;
	iblock = arch.iblock;

	fref = arch.fref;
	use_t2pred = arch.use_t2pred;
	means = arch.means;
	vars = arch.vars;
	mean_var_ready = arch.mean_var_ready;
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
	nblock = arch.nblock;
	nbin = arch.nbin;
	nchan = arch.nchan;
	npol = arch.npol;
	tbin = arch.tbin;
	frequencies = arch.frequencies;    
	profiles = arch.profiles;
	sub_mjd = arch.sub_mjd;
	sub_int = arch.sub_int;
	iblock = arch.iblock;

	fref = arch.fref;
	use_t2pred = arch.use_t2pred;
	means = arch.means;
	vars = arch.vars;
	mean_var_ready = arch.mean_var_ready;

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

	means.resize(frequencies.size(), 0.);
	vars.resize(frequencies.size(), 0.);

	if (use_t2pred)
	{
		fref = 0.5*(*min_element(frequencies.begin(), frequencies.end())+*max_element(frequencies.begin(), frequencies.end()));
		f0 = 1./pred.get_pfold(ref_epoch.to_day(), fref);
		f1 = pred.get_fdfold(ref_epoch.to_day(), fref);
	}

	sub_mjd = start_mjd + 0.5*databuffer.tsamp;

	if (f1 == 0)
	{
		f1 = -acc/CONST_C*f0;
	}

	if (acc == 0)
	{
		acc = -f1/f0*CONST_C;
	}

	mxWTW.resize(nbin*nbin, 0.);
	vWTd_T.resize(nbin*databuffer.nchans, 0.);

	hits.resize(nbin, 0);
	profilesTPF.resize(nbin*npol*nchan, 0.);
	profilesPFT.resize(npol*nchan*nbin, 0.);
}

bool ArchiveLite::runDspsr(DataBuffer<float> &databuffer)
{
	if (databuffer.counter <= 0)
		return false;

	sub_int.tsubint = nblock*databuffer.nsamples*databuffer.tsamp;
	MJD start_time = sub_mjd;
	MJD end_time = sub_mjd + (nblock*databuffer.nsamples-1)*databuffer.tsamp;
	MJD epoch = get_epoch(start_time, end_time, ref_epoch);
	sub_int.offs_sub = (epoch-start_mjd).to_second();
	
	double phi = 0.;
	double f = 0.;

	vector<int> binplan(databuffer.nsamples);
	for (long int i=iblock*databuffer.nsamples; i<(iblock+1)*databuffer.nsamples; i++)
	{
		if (i%NSBLK == 0)
		{
			phi = get_phase(sub_mjd+i*databuffer.tsamp, ref_epoch);
			f = get_ffold(sub_mjd+(i+NSBLK*0.5-0.5)*databuffer.tsamp, ref_epoch);
		}

		phi -= floor(phi);
		int ibin = int(phi*nbin);
		binplan[i-iblock*databuffer.nsamples] = ibin;
		hits[ibin]++;
		phi += f*databuffer.tsamp;
	}

	for (long int i=0; i<databuffer.nsamples; i++)
	{
		for (long int j=0; j<databuffer.nchans; j++)
		{
			profilesTPF[binplan[i]*databuffer.nchans+j] += databuffer.buffer[i*databuffer.nchans+j];
		}
	}

	if (++iblock == nblock)
	{
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

		iblock = 0;
		fill(hits.begin(), hits.end(), 0);
		fill(profilesTPF.begin(), profilesTPF.end(), 0.);
		fill(profilesPFT.begin(), profilesPFT.end(), 0.);
	}

	return true;
}

bool ArchiveLite::runTRLSM(DataBuffer<float> &databuffer)
{
	if (databuffer.counter <= 0)
		return false;

	sub_int.tsubint = nblock*databuffer.nsamples*databuffer.tsamp;
	MJD start_time = sub_mjd;
	MJD end_time = sub_mjd + (nblock*databuffer.nsamples-1)*databuffer.tsamp;
	MJD epoch = get_epoch(start_time, end_time, ref_epoch);
	sub_int.offs_sub = (epoch-start_mjd).to_second();
	sub_int.ffold = abs(get_ffold(epoch, ref_epoch));

	double phi = 0.;
	double f = 0.;

	for (long int i=iblock*databuffer.nsamples; i<(iblock+1)*databuffer.nsamples; i++)
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
				vWTd_T[l*databuffer.nchans+j] += vWli0*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
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
				vWTd_T[l*databuffer.nchans+j] += vWli0*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
				vWTd_T[m*databuffer.nchans+j] += vWli1*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
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
					vWTd_T[binplan[l]*databuffer.nchans+j] += vWli[l]*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
				}
			}
		}
	}

	if (++iblock == nblock)
	{
		for (long int l=0; l<nbin; l++)
		{
			for (long int m=0; m<nbin; m++)
			{
				mxWTW[l*nbin+m] /= (nblock*databuffer.nsamples*databuffer.tsamp*sub_int.ffold);
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

		iblock = 0;
		fill(mxWTW.begin(), mxWTW.end(), 0.);
		fill(vWTd_T.begin(), vWTd_T.end(), 0.);
	}

	return true;
}

bool ArchiveLite::runPresto(DataBuffer<float> &databuffer)
{
	if (databuffer.counter <= 0)
		return false;

	sub_int.tsubint = nblock*databuffer.nsamples*databuffer.tsamp;
	MJD start_time = sub_mjd;
	MJD end_time = sub_mjd + (nblock*databuffer.nsamples-1)*databuffer.tsamp;
	MJD epoch = get_epoch(start_time, end_time, ref_epoch);
	sub_int.offs_sub = (epoch-start_mjd).to_second();
	sub_int.ffold = abs(get_ffold(epoch, ref_epoch));

	double phi = 0.;
	double f = 0.;

	for (long int i=iblock*databuffer.nsamples; i<(iblock+1)*databuffer.nsamples; i++)
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

			for (long int j=0; j<databuffer.nchans; j++)
			{       
				vWTd_T[l*databuffer.nchans+j] += vWli0*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
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

			for (long int j=0; j<databuffer.nchans; j++)
			{
				vWTd_T[l*databuffer.nchans+j] += vWli0*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
				vWTd_T[m*databuffer.nchans+j] += vWli1*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
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
				for (long int j=0; j<databuffer.nchans; j++)
				{
					vWTd_T[binplan[l]*databuffer.nchans+j] += vWli[l]*databuffer.buffer[(i-iblock*databuffer.nsamples)*databuffer.nchans+j];
				}
			}
		}
	}

	if (databuffer.mean_var_ready)
	{
		double factor = databuffer.tsamp * nbin * f;
		for (long int j=0; j<databuffer.nchans; j++)
		{
			means[j] += databuffer.means[j] * (databuffer.nsamples / nbin) * factor;
			vars[j] += databuffer.vars[j] * (databuffer.nsamples / nbin) * factor * factor;
		}
		mean_var_ready = databuffer.mean_var_ready;
	}

	if (++iblock == nblock)
	{
		transpose_pad<float>(&sub_int.data[0], &vWTd_T[0], nbin, npol*nchan);
		
		sub_mjd += sub_int.tsubint;

		if (mean_var_ready)
		{
			sub_int.means = means;
			sub_int.vars = vars;
			sub_int.mean_var_ready = mean_var_ready;
		}

		profiles.push_back(sub_int);

		iblock = 0;
		fill(vWTd_T.begin(), vWTd_T.end(), 0.);

		std::fill(means.begin(), means.end(), 0.);
		std::fill(vars.begin(), vars.end(), 0.);
	}

	return true;
}

void ArchiveLite::read_archive(const std::string &fname)
{
	Psrfits arch;
	arch.filename = fname;
	arch.open();
	arch.primary.load(arch.fptr);
	arch.load_mode();
	arch.subint.load(arch.fptr);

	assert(arch.subint.mode == Integration::FOLD);
	assert(arch.subint.nsubint > 0);

	start_mjd = arch.primary.start_mjd;
	ref_epoch = arch.primary.start_mjd+arch.subint.integrations[arch.subint.nsubint/2].offs_sub;
	double period = arch.subint.integrations[arch.subint.nsubint/2].folding_period;
	f0 = 1./period;
	if (arch.subint.nsubint)
		f1 = (1./arch.subint.integrations[arch.subint.nsubint-1].folding_period-1./arch.subint.integrations[0].folding_period)/(arch.subint.integrations[arch.subint.nsubint-1].offs_sub-arch.subint.integrations[0].offs_sub);
	else
		f1 = 0.;
	acc = 0.;
	dm = arch.primary.chan_dm;
	snr = 0;

	nbin = arch.subint.nbin;
	nchan = arch.subint.nchan;
	npol = arch.subint.npol;
	tbin = arch.subint.tbin;
	frequencies.resize(nchan);
	std::memcpy(frequencies.data(), arch.subint.integrations[0].frequencies, sizeof(double)*nchan);

	IntegrationLite it;
	for (long int k=0; k<arch.subint.nsubint; k++)
	{
		it.ffold = 1./arch.subint.integrations[k].folding_period;
		it.tsubint = arch.subint.integrations[k].tsubint;
		it.offs_sub = arch.subint.integrations[k].offs_sub;
		int np = arch.subint.integrations[k].npol;
		int nc = arch.subint.integrations[k].nchan;
		int nb = arch.subint.integrations[k].nbin;
		it.resize(np, nc, nb);
		for (long int ip=0; ip<np; ip++)
		{
			for (long int ic=0; ic<nc; ic++)
			{
				for (long int ib=0; ib<nb; ib++)
				{
					it.data[ip*nc*nb+ic*nb+ib] = ((short *)arch.subint.integrations[k].data)[ip*nc*nb+ic*nb+ib]*arch.subint.integrations[k].scales[ip*nc+ic]+arch.subint.integrations[k].offsets[ip*nc+ic];
				}
			}
		}

		profiles.push_back(it);
	}

	/* dedispersion */
	double fmin = 1e6;
	double fmax = 0.;
	for (long int j=0; j<nchan; j++)
	{
		fmax = frequencies[j]>fmax? frequencies[j]:fmax;
		fmin = frequencies[j]<fmin? frequencies[j]:fmin;
	}

	std::vector<float> profile(npol*nchan*nbin, 0.);

	for (auto it=profiles.begin(); it!=profiles.end(); ++it)
	{
		for (long int k=0; k<npol; k++)
		{
			for (long int j=0; j<nchan; j++)
			{
				int delayn = round(Pulsar::DedispersionLite::dmdelay(dm, fmax, frequencies[j])*it->ffold*nbin);
				delayn = (-delayn)%nbin;
				if (delayn<0) delayn += nbin;

				for (long int i=0; i<delayn; i++)
				{
					profile[k*nchan*nbin+j*nbin+i] = it->data[k*nchan*nbin+j*nbin+(i-delayn+nbin)];
				}
				for (long int i=delayn; i<nbin; i++)
				{
					profile[k*nchan*nbin+j*nbin+i] = it->data[k*nchan*nbin+j*nbin+(i-delayn)];
				}
			}
		}
		std::swap(it->data, profile);
	}
}