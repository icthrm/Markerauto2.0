#include "option.h"
using namespace std;
string io_option::getfilename(string filepath){
    int i,index;
    string name_file;
    for(i=0;i<filepath.length();i++){
        if(filepath[i]=='/') index=i;
    }
    for(i=1;i<filepath.length()-index-4;i++){
        name_file=name_file+filepath[i+index];
    }
    return name_file;
}
void io_option::saveimg(string name,cv::Mat img){
    cv::imwrite(ouputpath+"/"+name+"_"+inputfilename+".jpg",img);
}