#include"fundmental.h"
float GetMidValue(cv::Mat input){

    int i,j;
    float out;
    cv::Mat img;
    img=input;
    vector<float> vector(img.cols*img.rows);
    for(i=0;i<img.rows;i++){
        for(j=0;j<img.cols;j++){
            vector[i*img.cols+j]=img.at<float>(i,j);
        }
    }
    //std::nth_element(vector[0],vector[(int)(img.cols*img.rows)/2],vector[img.cols*img.rows-1]);
    std::nth_element(vector.begin(),vector.begin()+(int)(img.cols*img.rows)/2,vector.end());
    out=vector[(int)(img.cols*img.rows)/2];
    return out;

}


cv::Mat insertmat_dim2to3(cv::Mat m3,cv::Mat m2,int d){
    /**************
     * @插入函数，将一个二维图像插入三维图像的其中一层
     * @输入一个三维图像m3，一个二维图像m2，将m2插入m3的第d层
     * m2的行数列数要与m3每层的行数列数对应且数据类型均为float
     ********/
    int i,j;
    for(i=0;i<m2.rows;i++){
        for(j=0;j<m2.cols;j++){
            m3.at<float>(d,i,j)=m2.at<float>(i,j);
        }
    }

    return m3;
}

void hardvalsmall(cv::Mat w,double t){
    /**************
     * @硬阈值函数，抹平小值
     * @输入一个图像，与阈值，将所有小于阈值的像素点归零
     * @必须是CV_32FC1
     ********/
    int i,j;

    for(i=0;i<w.rows;i++){
        for(j=0;j<w.cols;j++){
            if(w.at<float>(i,j)<t){
                w.at<float>(i,j)=0;
            }
        }
    }
}

cv::Mat hardval2(cv::Mat w,double t){
    Mat out=w.clone();
    out.convertTo(out,CV_8UC1,1,0);
    for(int i=0;i<out.rows;i++){
        for(int j=0;j<out.cols;j++){
            if(out.at<uchar>(i,j)<0) out.at<uchar>(i,j)=-out.at<uchar>(i,j);
        }
    }
    out=out>=t;
    return out;
}

cv::Mat ToOne(cv::Mat img){
/**************
 * @图像线性归一化函数
 * @输入一个图像，返回归一化后的图像
 ********/
    cv::Mat imginfo;
    cv::Point minIdx,maxIdx;//最值坐标
    double minValue,maxValue;//最值
    cv::minMaxLoc(img,&minValue,&maxValue,&minIdx,&maxIdx);//最值函数
    img.convertTo(img,CV_32FC1,1/(maxValue-minValue),-minValue/(maxValue-minValue));
    return img;

}
void removeSmall(cv::Mat img,double threshold){
    double ae;
    std::vector<std::vector<Point>> contours;
    std::vector<Vec4i> hierarchy;

    img.convertTo(img,CV_8UC1,1,0);
    img=img>200;//大于200的置为255
    findContours(img,contours,hierarchy,cv::RETR_TREE,cv::CHAIN_APPROX_SIMPLE,Point(0,0));//寻找轮廓
    for(int i=0;i<contours.size();i++){
        ae=contourArea(contours[i]);
        if(ae<threshold){
            cv::drawContours(img,contours,i,Scalar(0), -1);
        }
    }


    return;
}

cv::Mat avgImg(cv::Mat img,int k){
    /**************
     * @图像滑动平均函数
     * @输入一个图像与平滑和边长（奇数），返回滑动平均并归一化后的图像
     ********/
    assert(k%2 != 0);//平滑核为偶数时终止程序

    cv::Mat filterImg,filterinfo;

    cv::blur(img,filterImg,cv::Size(k,k));//均值滤波
    ToOne(filterImg);//归一化

    return filterImg;
}
cv::Mat Img_in(cv::Mat img){
    //图像取反函数，必须是CV_8UC1
    cv::Mat outimg;
    img.convertTo(outimg,CV_8UC1,-1,255);
    return outimg;
}
cv::Mat Img_in2(cv::Mat img){
    //图像取反函数，原图像范围0-1
    cv::Mat outimg;
    img.convertTo(outimg,CV_32FC1,-1,1);
    return outimg;
}

void drawPoint(cv::Mat img,cv::Mat res){
    /*******
     * @按照检测结果在图中画出点
     * @res中的顺序为r,x,y,rmse,sigma
    *******/
    int i,x,y,r;
    cv::Point xy;

    for(i=0;i<res.rows;i++){
        r=res.at<float>(i,0)+2;
        xy.y=res.at<float>(i,1);
        xy.x=res.at<float>(i,2);
        cv::circle(img,xy,r,cv::Scalar(255,255,0));
    }

}


void sortMat(Mat &stats, int colId){
    //根据指定列以行为单位排序

    Mat sorted_index;
    cv::sortIdx(stats(cv::Range(0,stats.rows),cv::Range(colId,colId+1)), sorted_index, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);

    Mat sorted_stats = stats.clone();
    int row_num = sorted_index.rows;
    for(int i = 0; i < row_num; i++){
        int _idx = sorted_index.at<int>(i, 0);
        sorted_stats.row(i) = stats.row(_idx) + 0;//必须加0否则会出很难debug的错误
    }
    stats = sorted_stats;
    return;
}
void drawimg(Mat img,string name,int time_wait){
    imshow(name,img);
    waitKey(time_wait);
}
