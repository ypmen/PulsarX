/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-05-28 11:43:17
 * @modify date 2021-05-28 11:43:17
 * @desc [description]
 */

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <streambuf>

#include "json.hpp"
#include "filmaker.h"

using json = nlohmann::json;

FilMaker::FilMaker()
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
    filltype = "mean";

    telescope_id = -1;
    src_raj = 0.;
    src_dej = 0.;
    ibeam = 1;
    id = 1;

    outmean = 10.;
    outstd = 3.;
    outnbits = 8;
}

FilMaker::~FilMaker(){}

void FilMaker::prepare(DataBuffer<float> &databuffer)
{
    downsample.td = td;
    downsample.fd = fd;
    downsample.prepare(databuffer);
    downsample.close();
    downsample.closable = true;

    equalize.prepare(downsample);
    equalize.close();
    equalize.closable = true;

    baseline.width = bswidth;
    baseline.prepare(equalize);
    baseline.close();
    baseline.closable = true;

    rfi.filltype = filltype;
    rfi.prepare(baseline);
    rfi.close();
    rfi.closable = true;

    stringstream ss_id;
    ss_id<<std::setw(2)<<std::setfill('0')<<id;
    string s_id = ss_id.str();

    filwriter.fil.filename = rootname + "_" + s_id + ".fil";
    std::strcpy(filwriter.fil.source_name, source_name.c_str());
    filwriter.fil.telescope_id = telescope_id;
    filwriter.fil.src_raj = src_raj;
    filwriter.fil.src_dej = src_dej;
    filwriter.fil.ibeam = ibeam;
    filwriter.fil.nbits = outnbits;
    filwriter.outmean = outmean;
    filwriter.outstd = outstd;
    filwriter.prepare(rfi);
}

void FilMaker::run(DataBuffer<float> &databuffer)
{
    DataBuffer<float> *data = downsample.run(databuffer);
    
    data = equalize.run(*data);

    data = baseline.run(*data);

    data = rfi.zap(*data, zaplist);
    if (rfi.isbusy) rfi.closable = false;
	
    for (auto irfi = rfilist.begin(); irfi!=rfilist.end(); ++irfi)
	{
        if ((*irfi)[0] == "mask")
        {
            data = rfi.mask(*data, threMask, stoi((*irfi)[1]), stoi((*irfi)[2]));
            if (rfi.isbusy) rfi.closable = false;
        }
        else if ((*irfi)[0] == "kadaneF")
        {
            data = rfi.kadaneF(*data, threKadaneF*threKadaneF, widthlimit, stoi((*irfi)[1]), stoi((*irfi)[2]));
            if (rfi.isbusy) rfi.closable = false;
        }
        else if ((*irfi)[0] == "kadaneT")
        {
            data = rfi.kadaneT(*data, threKadaneT*threKadaneT, bandlimitKT, stoi((*irfi)[1]), stoi((*irfi)[2]));
            if (rfi.isbusy) rfi.closable = false;
        }
		else if ((*irfi)[0] == "zdot")
        {
			data = rfi.zdot(*data);
            if (rfi.isbusy) rfi.closable = false;
        }
		else if ((*irfi)[0] == "zero")
        {
			data = rfi.zero(*data);
            if (rfi.isbusy) rfi.closable = false;
        }
	}

    if (!databuffer.isbusy) data->closable = true;

    filwriter.run(*data);
}

void plan(variables_map &vm, std::vector<FilMaker> &filmakers)
{
    FilMaker fm;

    fm.bswidth = vm["baseline"].as<vector<float>>().back();

	std::vector<std::string> rfi_opts;
	if (vm.count("rfi"))
	{
        rfi_opts = vm["rfi"].as<std::vector<std::string>>();
        for (auto opt=rfi_opts.begin(); opt!=rfi_opts.end(); ++opt)
        {
            if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
            {
                std::vector<std::string> temp{*opt, *(opt+1), *(opt+2)};       
                fm.rfilist.push_back(temp);
                std::advance(opt, 2);
            }
            else if (*opt == "zap")
            {
                fm.zaplist.push_back(std::pair<double, double>(std::stod(*(opt+1)), std::stod(*(opt+2))));
                std::advance(opt, 2);
            }
            else if (*opt=="zero" or *opt=="zdot")
            {
                std::vector<std::string> temp{*opt};
                fm.rfilist.push_back(temp);
            }
        }
	}

    fm.threMask = vm["threMask"].as<float>();
	fm.bandlimit = vm["bandlimit"].as<double>();
    fm.bandlimitKT = vm["bandlimitKT"].as<double>();
    fm.threKadaneF = vm["threKadaneF"].as<float>();
    fm.threKadaneT = vm["threKadaneT"].as<float>();
    fm.widthlimit = vm["widthlimit"].as<double>();

    fm.outmean = vm["mean"].as<float>();
    fm.outstd = vm["std"].as<float>();
    fm.outnbits = vm["nbits"].as<int>();

    if (vm.count("filplan"))
    {
        std::string filename = vm["filplan"].as<std::string>();
        std::ifstream filplan(filename);
        std::string filplan_str((std::istreambuf_iterator<char>(filplan)), std::istreambuf_iterator<char>());
       
        json filplan_json =json::parse(filplan_str);

        int id = 0;
        for (auto p=filplan_json.begin(); p!=filplan_json.end(); ++p)
        {
            fm.td = (*p)["time_downsample"];
            fm.fd = (*p)["frequency_downsample"];
            fm.bswidth = (*p)["baseline_width"];
            fm.outmean = (*p)["dataout_mean"];
            fm.outstd = (*p)["dataout_std"];
            fm.outnbits = (*p)["dataout_nbits"];
            std::string rfi_flags = (*p)["rfi_flags"];

            std::vector<std::string> parameters;
            boost::split(parameters, rfi_flags, boost::is_any_of("\t "), boost::token_compress_on);

            for (auto opt=parameters.begin(); opt!=parameters.end(); ++opt)
            {
                if (*opt=="mask" or *opt=="kadaneF" or *opt=="kadaneT")
                {
                    std::vector<std::string> temp{*opt, *(opt+1), *(opt+2)};       
                    fm.rfilist.push_back(temp);
                    std::advance(opt, 2);
                }
                else if (*opt == "zap")
                {
                    fm.zaplist.push_back(std::pair<double, double>(std::stod(*(opt+1)), std::stod(*(opt+2))));
                    std::advance(opt, 2);
                }
                else if (*opt=="zero" or *opt=="zdot")
                {
                    std::vector<std::string> temp{*opt};
                    fm.rfilist.push_back(temp);
                }
            }

            fm.id = ++id;
            filmakers.push_back(fm);
        }

        filplan.close();
    }
    else
    {
        fm.id = 1;
        filmakers.push_back(fm);
    }
}