/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-11-02 10:54:12
 * @modify date 2020-11-02 10:54:12
 * @desc [description]
 */

#ifndef ARCHIVEWRITER_H
#define ARCHIVEWRITER_H

#include <string>
#include <vector>

#include "psrfits.h"
#include "mjd.h"
#include "integration.h"

using namespace std;

#include "dedispersionlite.h"
#include "archivelite.h"
#include "gridsearch.h"

class ArchiveWriter
{
public:
	ArchiveWriter();
	~ArchiveWriter();
	void close();
	void prepare();
	void prepare(Pulsar::ArchiveLite &arch, Pulsar::GridSearch &grid);
	void run(Pulsar::ArchiveLite &arch, Pulsar::GridSearch &grid);
	void run(vector<float> &profile, int np, int nc, int nb, double fold_period, double tsubint, double offs_sub);
	void run(vector<float> &profiles, int ns, int np, int nc, int nb, vector<double> &fold_periods, vector<double> &tsubints, vector<double> &offs_subs);
public:
	void get_dedispersion(vector<float> &profiles, Pulsar::ArchiveLite &arch, Pulsar::GridSearch &grid)
	{
		double fmin = 1e6;
		double fmax = 0.;
		for (long int j=0; j<nchan; j++)
		{
			fmax = frequencies[j]>fmax? frequencies[j]:fmax;
			fmin = frequencies[j]<fmin? frequencies[j]:fmin;
		}

		profiles.resize(grid.nsubint*1*nchan*nbin, 0.);

		for (long int l=0; l<grid.nsubint; l++)
		{
			for (long int k=0; k<1; k++)
			{
				for (long int j=0; j<nchan; j++)
				{
					int delayn = round(Pulsar::DedispersionLite::dmdelay(-dm, fmax, frequencies[j])*(grid.f0+grid.f1*grid.tsuboff[l])*nbin);
					delayn = (-delayn)%nbin;
					if (delayn<0) delayn += nbin;

					for (long int i=0; i<delayn; i++)
					{
						profiles[l*1*nchan*nbin+k*nchan*nbin+j*nbin+i] = grid.profiles[l*nchan*nbin+j*nbin+(i-delayn+nbin)];
					}
					for (long int i=delayn; i<nbin; i++)
					{
						profiles[l*1*nchan*nbin+k*nchan*nbin+j*nbin+i] = grid.profiles[l*nchan*nbin+j*nbin+(i-delayn)];
					}
				}
			}
		}
	}

public:
	string rootname;
	string template_file;
	Integration::Mode mode;
	int ibeam;
	string src_name;
	string ra;
	string dec;
public:
	MJD start_mjd;
	int npol;
	int nchan;
	int nbin;
	double tbin;
	double dm;
	vector<double> frequencies;
public:
	int nsubint;
private:
	Integration it;
	Psrfits fits;
};

#endif /* ARCHIVEWRITER_H */
