/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-05-30 10:44:10
 * @modify date 2021-05-30 10:44:10
 * @desc [description]
 */

#ifndef FILTERBANKWRITER_H
#define FILTERBANKWRITER_H

#include <fstream>

#include "databuffer.h"
#include "filterbank.h"

class FilterbankWriter
{
public:
	FilterbankWriter();
	FilterbankWriter(const FilterbankWriter &filwriter);
	FilterbankWriter & operator=(const FilterbankWriter &filwriter);
	FilterbankWriter(const Filterbank &temp);
	~FilterbankWriter();
	void prepare(DataBuffer<float> &databuffer);
	void run(DataBuffer<float> &databuffer);
public:
	float outstd;
	float outmean;
	Filterbank fil;
	std::ofstream outfile;
};

#endif /* FILTERBANKWRITER_H */
