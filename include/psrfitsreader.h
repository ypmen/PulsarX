/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-20 22:11:40
 * @modify date 2022-08-20 22:11:40
 * @desc [description]
 */

#ifndef PSRFITSREADER_H
#define PSRFITSREADER_H

#include "psrdatareader.h"
#include "psrfits.h"

class PsrfitsReader : public PSRDataReader
{
public:
	PsrfitsReader();
	~PsrfitsReader();
	void check();
	void read_header();
	size_t read_data(DataBuffer<float> &databuffer, size_t ndump, bool virtual_reading = false);
	MJD get_start_mjd_curfile()
	{
		Integration tmpit;
		psf[idmap[ifile_cur]].subint.load_integration(psf[idmap[ifile_cur]].fptr, 0, tmpit);
		return psf[idmap[ifile_cur]].primary.start_mjd + (tmpit.offs_sub - 0.5 * psf[idmap[ifile_cur]].subint.nsblk * psf[idmap[ifile_cur]].subint.tbin);
	}
	double get_tsamp_curfile(){return psf[idmap[ifile_cur]].subint.tbin;}
	size_t get_count_curfile(){return ns_psfn;}
	size_t get_count(){return count;}
	size_t get_ifile(){return ifile_cur;}
	size_t get_ifile_ordered(){return idmap[ifile_cur];}
	void get_filterbank_template(Filterbank &filtem);

private:
	Integration it;
	Integration it8;
	size_t ntot;
	size_t count;
	size_t ns_psfn;
	std::vector<Psrfits> psf;
	long int ifile_cur;
	long int isubint_cur;
	long int isample_cur;
	bool update_file;
	bool update_subint;
};


#endif /* PSRFITSREADER_H */
