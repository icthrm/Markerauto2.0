#include"BaseFun.h"

cv::Mat normalize1(cv::Mat img){
/**************
 * @图像归一化函数
 * @输入一个图像，返回归一化后的图像
 ********/
    cv::Mat imginfo;
    cv::Point minIdx,maxIdx;//最值坐标
    double minValue,maxValue;//最值
    cv::minMaxLoc(img,&minValue,&maxValue,&minIdx,&maxIdx);//最值函数
    img.convertTo(imginfo,CV_32FC1,1/(maxValue-minValue),-minValue/(maxValue-minValue));
    return imginfo;
}
double roundness(cv::Mat label_img){
    //计算圆度
    cv::Mat img;
    double a,b;
    std::vector<std::vector<Point>> contours;
    std::vector<Vec4i> hierarchy;

    label_img.convertTo(img,CV_8UC1,1,0);
    findContours(img,contours,hierarchy,cv::RETR_TREE,cv::CHAIN_APPROX_SIMPLE,Point(0,0));//寻找轮廓

    a=cv::contourArea(contours[0])*4*M_PI;
    b = pow(cv::arcLength(contours[0],true), 2);
    if (b == 0) return 0;
    return a/b;
}
void find_local_peak(cv::Mat img, int m,int n,int m_w,int n_w,int &out_m,int &out_n){
    int row,col;
    cv::Mat sub_img;
    cv::Point minloc,maxloc;
    double minval,maxval;

    row=m+m_w;
    col=n+n_w;
    if(row>img.rows) row=img.rows;
    if(col>img.cols) col=img.cols;

    sub_img= img(cv::Range(m,row),cv::Range(n,col));
    cv::minMaxLoc(sub_img,&minval,&maxval,&minloc,&maxloc);

    out_m = m + maxloc.y;
    out_n = n + maxloc.x;
}

bool have255(cv::Mat img){
    int i,j,a;
    for(i=0;i<img.rows;i++){
        for(j=0;j<img.cols;j++){
            if(img.at<uchar>(i,j) == 255) a=1;
        }
    }
    if(a==1)
        return true;
    else
        return false;
}

float sum_pix_float(cv::Mat img){
    int i,j;
    float out;
    out=0;
    for(i=0;i<img.rows;i++){
        for(j=0;j<img.cols;j++){
            out+=img.at<float>(i,j);
        }
    }
    return out;
}
