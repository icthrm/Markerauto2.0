#include "opts.h"
#include <iostream>
#include <fstream>
#include "dataf/dataf.h"
#include "dataf/calibration.h"
#include "mrcimg/mrc2img.h"
#include <string>
#include<unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mrcimg/img_util.h>

using namespace std;

int main(int argc, char **argv)
{
	struct options opts;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}
	
	std::cout<<"MRC to PGM."<<std::endl;
	
	util::MrcStack mrcs;
    
	std::string mrcfile = opts.input;
	if(!mrcs.Open(mrcfile.c_str())){
		return 0;
	}
	
	int dos = mrcfile.rfind('/');
    if(dos != std::string::npos){
        mrcfile.erase(0, dos);
    }
    
    dos = mrcfile.rfind('.');
    if(dos != std::string::npos){
        mrcfile.erase(dos, mrcfile.size());
    }
	
	if(access(opts.output.c_str(),0) == -1){
        if (mkdir(opts.output.c_str(),0777)){
            printf("creat folder failed!!!\n");
			return 0;
        }
    }
	
	std::stringstream ss(opts.idxs);
	int idx_s, idx_e; char ch;
	
	ss>>idx_s>>ch>>idx_e;
	for(int i = idx_s; i <= idx_e; i++){
		cv::Mat img = mrcs.GetStackImage(i);
        
		util::ConvertTo1(img, true);

// 		std::cout<<img<<std::endl;
		std::stringstream dest(mrcfile);
		dest<<opts.output<<"/"<<mrcfile<<i<<".pgm";
		util::SaveImage(img, dest.str().c_str());
	}
	
    mrcs.Close();	
}

