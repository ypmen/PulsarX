/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-07 15:54:10
 * @modify date 2020-06-07 15:54:10
 * @desc [description]
 */

#ifndef ARCHIVELITE_H
#define ARCHIVELITE_H

#include <fstream>
#include <vector>

#include "dedispersionlite.h"
#include "predictor.h"
#include "kepler.h"
#include "mjd.h"

using namespace std;

namespace Pulsar
{
	class IntegrationLite
	{
	public:
		IntegrationLite()
		{
			ffold = 0.;
			tsubint = 0.;
			offs_sub = 0.;
			npol = 0;
			nchan = 0;
			nbin = 0;

			mean_var_ready = false;
		}
		~IntegrationLite(){}
		void resize(int np, int nc, int nb)
		{
			npol = np;
			nchan = nc;
			nbin = nb;
			data.clear();
			data.resize(npol*nchan*nbin, 0.);
		}
	public:
		double ffold;
		double tsubint;
		double offs_sub;

		int npol;
		int nchan;
		int nbin;

		vector<float> data;
		std::vector<double> means;
		std::vector<double> vars;
		bool mean_var_ready;
	};

	class ArchiveLite
	{
	public:
		ArchiveLite();
		~ArchiveLite();
		void close()
		{
			for (auto p=profiles.begin(); p!=profiles.end(); ++p)
			{
				p->data.clear();
				p->data.shrink_to_fit();
			}
		}
		void prepare(DataBuffer<float> &databuffer);
		bool runDspsr(DataBuffer<float> &databuffer);
		bool runTRLSM(DataBuffer<float> &databuffer);
		bool runPresto(DataBuffer<float> &databuffer);
		void resize(int np, int nc, int nb);
		void dump2bin(const string &rootname)
		{
			ofstream outfile;
			outfile.open(rootname+".spft", ofstream::binary);

			for (auto pro=profiles.begin(); pro!=profiles.end(); ++pro)
			{
				outfile.write((char *)(&((*pro).data[0])), sizeof(float)*npol*nchan*nbin);
			}
			outfile.close();
		}
		void read_archive(const std::string &fname);
	private:
		double get_phase(MJD mjd, MJD &mjdref)
		{
			if (use_t2pred)
			{
				return pred.get_fphase(mjd.to_day(), fref);
			}
			else if (use_kepler)
			{
				return orb.get_fphase(mjd.to_day(), mjdref.to_day());
			}
			else
			{
				long double t = (mjd-mjdref).to_second();
				double phi = f0*t + 0.5*f1*t*t + f2*t*t*t/6.;
				phi -= floor(phi);
				return phi;
			}
		}
		double get_ffold(MJD mjd, MJD &mjdref)
		{
			if (use_t2pred)
			{
				double pfold = pred.get_pfold(mjd.to_day(), fref);
				return 1./pfold;
			}
			else if (use_kepler)
			{
				return orb.get_ffold(mjd.to_day());
			}
			else
			{
				long double t = (mjd-mjdref).to_second();
				double f = f0 + f1*t + 0.5*f2*t;
				return f;
			}
		}
		MJD get_epoch(MJD &start_time, MJD &end_time, MJD &ref_epoch)
		{
			long double t = (start_time.to_second()+end_time.to_second())/2.;
			t -= ref_epoch.to_second();
			long int phi = floor(f0*t+0.5*f1*t*t+f2*t*t*t/6.);
			long double tepoch = 0.;
			if (f1 != 0 || f2 != 0)
			{
				if (f2 == 0)
				{
					double tmp = f0*f0+2*f1*phi;
					tmp = tmp<0? 0:tmp;
					long double tepoch1 = (sqrt(tmp)-f0)/f1;
					long double tepoch2 = (-sqrt(f0*f0+2*f1*phi)-f0)/f1;
					tepoch = abs(t-tepoch1)<abs(t-tepoch2) ? tepoch1:tepoch2;
				}
				else
				{
					double a = phi / f0;
					double b = -0.5 * f1 / f0;
					double c = -f2 / f0 / 6.;

					int n = 8;
					tepoch = t;
					long double dt = 1.;
					while (std::abs(dt) > 1e-9 && --n >= 0)
					{
						dt = a - tepoch  + b * tepoch * tepoch + c * tepoch * tepoch * tepoch;
						tepoch += dt;
					}
				}
			}
			else
				tepoch = phi/f0;

			return MJD((tepoch+ref_epoch.to_second())/86400.);
		}
	public:
		MJD start_mjd;
		MJD ref_epoch;
		double f0;
		double f1;
		double f2;
		double acc;
		double dm;
		double snr;
		long int nblock;
		bool use_t2pred;
		Predictors pred;
		bool use_kepler;
		Kepler orb;
		bool dedispersed;
	public:
		int nbin;
		int nchan;
		int npol;
		double tbin;
		vector<double> frequencies;
		vector<IntegrationLite> profiles;
		MJD sub_mjd;
		IntegrationLite sub_int;
		std::vector<double> means;
		std::vector<double> vars;
		bool mean_var_ready;
	private:
		double fref;
		int iblock;
		vector<float> mxWTW;
		vector<float> vWTd_T;
		vector<int> hits;
		vector<float> profilesTPF;
		vector<float> profilesPFT;
	};
}

#endif /* ARCHIVELITE_H */
