#include"wavelet.h"
extern const int global_tempout;
extern const int global_timeout;

cv::Mat fullZero2(cv::Mat k,int lk){
     /*************
      * @卷积核补零函数
      * @输入当前卷积核与其长度，将每两个元素之间舔一个0后返回
      *************/
     int i;
     cv::Mat out=cv::Mat::zeros(1,2*lk-1,CV_32FC1);
     for(i=0;i<lk;i++){
          out.at<float>(0,2*i)=k.at<float>(0,i);
     }
     return out;
}

cv::Mat atrous_cv(cv::Mat ai,cv::Mat k){
     /**************
     * @atrous变换处理函数
     * @输入一个图像ai与数组k，返回处理结果，所用卷积核为数组与自身转置的克罗内克积
     ********/
     int i,j;
     cv::Mat ker,out;
     ker=cv::Mat::zeros(k.cols,k.cols,CV_32FC1);
     for(i=0;i<k.cols;i++){
         ker.row(i)=k.at<float>(0,i)*k;
     }//克罗内克积得卷积核
     //cout<<ker<<endl;
     ai.convertTo(ai,CV_64FC1,1,0);
     ker.convertTo(ker,CV_64FC1,1,0);
     cv::filter2D(ai,out,-1,ker,Point(-1, -1), 0,BORDER_REFLECT);
     out.convertTo(out,CV_32FC1,1,0);

     return out;
}

cv::Mat wavelet(cv::Mat a,cv::Mat A){
     /********
      * @小波图像计算函数
      * @a，A为尺度为i，i+1的近似图像，返回尺度i+1经过阈值化的图像
      *******/
     cv::Mat w,temp;
     float temp1;
     double t;

     w=a-A;
     temp1=GetMidValue(w);
     temp=temp1*cv::Mat::ones(a.rows,a.cols,CV_32FC1)-w;
     temp=abs(temp);
     temp1=GetMidValue(temp);
     t=temp1*3/0.67;
     //cout<<"t is "<<t<<endl;
     //if(outpro == 1) cv::imwrite("beforehaidvalsmall.jpg",w*255);
     hardvalsmall(w,t);
     //if(outpro == 1) cv::imwrite("afterhaidvalsmall.jpg",w*255);

     return w;
}

cv::Mat waveletprocess2(cv::Mat img,int J){
     /**************
     * @小波变换处理函数
     * @输入一个图像，与小波变换深度，返回小波相关图像
     * sigma为第一幅图的方差,同时传回
     ********/

     assert(J>1);//小波变换深度小于一时终止程序

     cv::Mat inimg,a,k,temp,temp1,temp2,temp3,temp4,aaa,kk;
     double mean1;
     int m,n,i,j;

     //atrous核初始化
     k=(Mat_<float>(1,5)<<1.0 / 16, 1.0 / 4, 3.0 / 8, 1.0 / 4, 1.0 / 16);


     //Ai初始化，用于记录所有近似图像
     //int size_Ai[3]={J,img.rows,img.cols};
     //Mat Ai(3,size_Ai,CV_32FC1,cv::Scalar(0));

     //W初始化，用以记录所有小波图像
     //int size_Wi[3]={J-1,img.rows,img.cols};
     //Mat Wi(3,size_Wi,CV_32FC1,cv::Scalar(0));

     img.convertTo(img,CV_8UC1,255,0);//转为uchar
     inimg=Img_in(img);//图像取反



     inimg.convertTo(inimg,CV_8UC1,1,0);
     a=inimg.clone();
     a.convertTo(a,CV_32FC1,1,0);
     m=a.rows;
     n=a.cols;
     //Ai=insertmat_dim2to3(Ai,a,0);//记录取反后的原始图像

     //计算近似图像组
     temp=a.clone();
     for(i=1;i<J;i++){
          if(i == J-1)
               temp1=temp.clone();
          temp2=atrous_cv(temp,k);
          temp=temp2;
          kk=fullZero2(k,k.cols);
          k=kk;
     }

     temp2=wavelet(temp1,temp);
     return temp2;
}
