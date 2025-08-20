#include"templateMatch.h"
cv::Mat expandImage(cv::Mat img,int row,int col){
    /*********
     * @拓展图像，需为float
     * *********/
    int m=img.rows;
    int n=img.cols;
    cv::Mat out(m+row*2,n+col*2,CV_32FC1,cv::Scalar(0));
    img.copyTo(out(cv::Range(row,row+m), cv::Range(col,col+n)));
    return out;
}

double corr_p_p(cv::Mat p1,cv::Mat p2){
    int i,j;
    cv::Mat mean1,std1,pp1,pp2;

    p2.convertTo(pp2,CV_32FC1,1,0);

    cv::meanStdDev(p1,mean1,std1);
    pp1=(p1-mean1.at<double>(0))/std1.at<double>(0);

    cv::meanStdDev(p2,mean1,std1);
    pp2=(pp2-mean1.at<double>(0))/std1.at<double>(0);

    for(i=0;i<p1.rows;i++){
        for(j=0;j<p1.cols;j++){
            pp1.at<float>(i,j)=pp1.at<float>(i,j)*pp2.at<float>(i,j);
        }
    }
    return cv::mean(pp1)[0];
}
