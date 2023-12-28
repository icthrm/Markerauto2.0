#include <iostream>
#include <fstream>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
class io_option{
public:
    bool show_img = 0;
    bool save_img = 0;
    bool show_plt = 0;
    string inputfile;
    string ouputpath;
    string inputfilename;

    string getfilename(string filepath);
    void saveimg(string name,cv::Mat img);
};