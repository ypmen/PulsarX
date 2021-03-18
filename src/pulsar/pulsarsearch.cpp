/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-05-19 09:44:47
 * @modify date 2020-05-19 09:44:47
 * @desc [description]
 */

#include <fstream>
#include <vector>

#include "pulsarsearch.h"

using namespace std;

PulsarSearch::PulsarSearch()
{
    td = 1;
    fd = 1;
    bswidth = 0.1;
    threMask = 7;
    bandlimit = 10;
    bandlimitKT = 10.;
    threKadaneT = 7;
    threKadaneF = 10;
    widthlimit = 10e-3;
    dms = 0;
    ddm = 1;
    ndm = 1000;
    
    ibeam = 1;
    id = 1;

    outmean = 0.;
    outstd = 0.;
    outnbits = 0;
}

PulsarSearch::~PulsarSearch(){}

void PulsarSearch::prepare(DataBuffer<float> &databuffer)
{
    downsample.td = td;
    downsample.fd = fd;
    downsample.prepare(databuffer);
    downsample.close();

    equalize.prepare(downsample);
    equalize.close();

    baseline.width = bswidth;
    baseline.prepare(equalize);
    baseline.close();

    rfi.prepare(baseline);
    rfi.close();

    if (format == "presto") outnbits = 32;

    dedisp.dms = dms;
    dedisp.ddm = ddm;
    dedisp.ndm = ndm;
    dedisp.ndump = rfi.nsamples;
    dedisp.rootname = rootname;
    dedisp.prepare(rfi);
    dedisp.preparedump(fildedisp, outnbits, format);
}

void PulsarSearch::run(DataBuffer<float> &databuffer)
{
    downsample.open();
    downsample.run(databuffer);
    databuffer.close();
    
    equalize.open();
    equalize.run(downsample);
    downsample.close();

    baseline.open();
    baseline.run(equalize);
    equalize.close();

    rfi.open();
    rfi.zap(baseline, zaplist);
	
    for (auto irfi = rfilist.begin(); irfi!=rfilist.end(); ++irfi)
	{
        if ((*irfi)[0] == "mask")
        {
            rfi.mask(rfi, threMask, stoi((*irfi)[1]), stoi((*irfi)[2]));
        }
        else if ((*irfi)[0] == "kadaneF")
        {
            rfi.kadaneF(rfi, threKadaneF*threKadaneF, widthlimit, stoi((*irfi)[1]), stoi((*irfi)[2]));
        }
        else if ((*irfi)[0] == "kadaneT")
        {
            rfi.kadaneT(rfi, threKadaneT*threKadaneT, bandlimitKT, stoi((*irfi)[1]), stoi((*irfi)[2]));
        }
		else if ((*irfi)[0] == "zdot")
        {
			rfi.zdot(rfi);
        }
		else if ((*irfi)[0] == "zero")
        {
			rfi.zero(rfi);
        }
	}
    baseline.close();

    dedisp.run(rfi, rfi.nsamples);
    rfi.close();

    databuffer.open();
    dedisp.rundump(outmean, outstd, outnbits, format);
}

void plan(variables_map &vm, vector<PulsarSearch> &search)
{
    PulsarSearch sp;

    sp.td = vm["td"].as<int>();
	sp.fd = vm["fd"].as<int>();

    sp.bswidth = vm["baseline"].as<float>();

	vector<string> rfi_opts;
	if (vm.count("rfi"))
	{
        rfi_opts = vm["rfi"].as<vector<string>>();
        for (auto opt=rfi_opts.begin(); opt!=rfi_opts.end(); ++opt)
        {
            if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
            {
                vector<string> temp{*opt, *(opt+1), *(opt+2)};       
                sp.rfilist.push_back(temp);
                advance(opt, 2);
            }
            else if (*opt == "zap")
            {
                sp.zaplist.push_back(pair<double, double>(stod(*(opt+1)), stod(*(opt+2))));
                advance(opt, 2);
            }
            else if (*opt=="zero" or *opt=="zdot")
            {
                vector<string> temp{*opt};
                sp.rfilist.push_back(temp);
            }
        }
	}

    sp.threMask = vm["threMask"].as<float>();
	sp.bandlimit = vm["bandlimit"].as<double>();
    sp.bandlimitKT = vm["bandlimitKT"].as<double>();
    sp.threKadaneF = vm["threKadaneF"].as<float>();
    sp.threKadaneT = vm["threKadaneT"].as<float>();
    sp.widthlimit = vm["widthlimit"].as<double>();

    sp.dms = vm["dms"].as<double>();
	sp.ddm = vm["ddm"].as<double>();
	sp.ndm = vm["ndm"].as<int>();

    sp.outmean = vm["mean"].as<float>();
    sp.outstd = vm["std"].as<float>();
    sp.outnbits = vm["nbits"].as<int>();

    sp.format = vm["format"].as<string>();

    if (vm.count("ddplan"))
    {
        string filename = vm["ddplan"].as<string>();
        string line;
        ifstream ddplan(filename);
        int id = 0;
        while (getline(ddplan, line))
        {
            boost::trim(line);
			if(!isdigit(line[0])) continue;

            vector<string> parameters;
            boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);
            
            sp.td = stoi(parameters[0]);
            sp.fd = stoi(parameters[1]);
            sp.dms = stod(parameters[2]);
            sp.ddm = stod(parameters[3]);
            sp.ndm = stol(parameters[4]);
            
            for (auto opt=parameters.begin()+5; opt!=parameters.end(); ++opt)
            {
                if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
                {
                    vector<string> temp{*opt, *(opt+1), *(opt+2)};       
                    sp.rfilist.push_back(temp);
                    advance(opt, 2);
                }
                else if (*opt == "zap")
                {
                    sp.zaplist.push_back(pair<double, double>(stod(*(opt+1)), stod(*(opt+2))));
                    advance(opt, 2);
                }
                else if (*opt=="zero" or *opt=="zdot")
                {
                    vector<string> temp{*opt};
                    sp.rfilist.push_back(temp);
                }
            }

            sp.id = ++id;
            search.push_back(sp);
        }
        ddplan.close();
    }
    else
    {
        sp.id = 1;
        search.push_back(sp);
    }
}
