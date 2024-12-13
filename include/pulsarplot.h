/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-06-14 15:54:14
 * @modify date 2020-06-14 15:54:14
 * @desc [description]
 */

#ifndef PULSARPLOT_H
#define PULSARPLOT_H

#include <map>

#include "utils.h"
#include "gridsearch.h"
#include "dedispersionliteU.h"

namespace Pulsar
{
	class PulsarPlot
	{
	public:
		PulsarPlot();
		~PulsarPlot();
		void plot(const ArchiveLite &archive, GridSearch &gridsearch, std::map<std::string, std::string> &obsinfo, int id, const string &rootname, bool plotx=false, bool save2fits=false);        
	public:
		static void get_dm_chisq_curve(std::vector<float> &vchisq, const std::vector<float> &vddm, const std::vector<double> &frequencies, double f0, double boxphw, int nbin);
		static void get_f0_chisq_curve(std::vector<float> &vchisq, const std::vector<float> &vdf0, const std::vector<double> &tsuboff, double f0, double boxphw, int nbin);
		static void get_f1_chisq_curve(std::vector<float> &vchisq, const std::vector<float> &vdf1, const std::vector<double> &tsuboff, double f0, double boxphw, int nbin);
	};
}

#endif /* PULSARPLOT_H */
