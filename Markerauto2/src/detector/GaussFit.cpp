#include"GaussFit.h"
cv::Mat Gauss_A(cv::Mat loc){
    int i;
    cv::Mat A(4,4,CV_32FC1,cv::Scalar(0));
    for(i=0;i<loc.rows;i++){
        A.at<float>(0,0)+=(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1))*(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));

        A.at<float>(0,1)+=loc.at<float>(i,0)*(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));

        A.at<float>(0,2)+=loc.at<float>(i,1)*(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));

        A.at<float>(0,3)+=(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));


        A.at<float>(1,1)+=loc.at<float>(i,0)*loc.at<float>(i,0);

        A.at<float>(1,2)+=loc.at<float>(i,0)*loc.at<float>(i,1);

        A.at<float>(1,3)+=loc.at<float>(i,0);

        A.at<float>(2,2)+=loc.at<float>(i,1)*loc.at<float>(i,1);

        A.at<float>(2,3)+=loc.at<float>(i,1);

        A.at<float>(3,3)++;
    }

    A.at<float>(1,0)=A.at<float>(0,1);
    A.at<float>(2,0)=A.at<float>(0,2);
    A.at<float>(2,1)=A.at<float>(1,2);
    A.at<float>(3,0)=A.at<float>(0,3);
    A.at<float>(3,1)=A.at<float>(1,3);
    A.at<float>(3,2)=A.at<float>(2,3);

    return A;

}
cv::Mat Gauss_b(cv::Mat info){
    int i;
    cv::Mat b(4,1,CV_32FC1,cv::Scalar(0));

    for(i=0;i<info.rows;i++){
        b.at<float>(0,0)+=(info.at<float>(i,0)*info.at<float>(i,0)+info.at<float>(i,1)*info.at<float>(i,1))*log(info.at<float>(i,2));

        b.at<float>(0,1)+=info.at<float>(i,0)*log(info.at<float>(i,2));

        b.at<float>(0,2)+=info.at<float>(i,1)*log(info.at<float>(i,2));

        b.at<float>(0,3)+=log(info.at<float>(i,2));
    }
    return b;
}




cv::Mat compute_center_Gauss(cv::Mat pic){
    int m,n,i,j;
    float sigma,x0,y0,para_A;
    cv::Mat paraDir,loc,info,X,Y,A,b,x;

    m=pic.rows;
    n=pic.cols;

    X.create(m*n,1,CV_32FC1);
    Y.create(m*n,1,CV_32FC1);
    loc.create(m*n,2,CV_8UC1);
    info.create(m*n,3,CV_32FC1);
    paraDir.create(4,1,CV_32FC1);

    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            info.at<float>(i*n+j,0)=(i+0.0)/(m+0.0);
            info.at<float>(i*n+j,1)=(j+0.0)/(n+0.0);
            info.at<float>(i*n+j,2)=pic.at<float>(i,j);
        }
    }
    A = Gauss_A(info);
    b = Gauss_b(info);
    cv::solve(A,b,x,cv::DECOMP_LU);//因为只有4*4，所以用LU

    sigma=sqrt(0.5/abs(x.at<float>(0,0))*m*n);
    x0=-0.5*x.at<float>(0,1)/x.at<float>(0,0)*m;
    y0=-0.5*x.at<float>(0,2)/x.at<float>(0,0)*n;
    para_A=exp(x.at<float>(0,3)-(x.at<float>(0,2)*x.at<float>(0,2)+x.at<float>(0,1)*x.at<float>(0,1))/(4.0*x.at<float>(0,1)));

    paraDir.at<float>(0,0)=sigma;
    paraDir.at<float>(0,1)=x0;
    paraDir.at<float>(0,2)=y0;
    paraDir.at<float>(0,3)=para_A;

    return paraDir;
}
cv::Mat gaussian(float A,float x0,float y0,float sigma,int size_x,int size_y){
    int i,j;
    cv::Mat f(size_x,size_y,CV_32FC1,cv::Scalar(0));
    for(i=0;i<size_x;i++){
        for(j=0;j<size_y;j++){
            f.at<float>(i,j)=A*exp(-(i-x0+0.0)*(i-x0+0.0)/(2.0*sigma*sigma)-(j-y0+0.0)*(j-y0+0.0)/(2.0*sigma*sigma));
        }
    }

    return f;
}
float compute_gauss_error(cv::Mat pic,cv::Mat dir){
    // 判断gaussian拟合结果的好坏
    // :param picture: 被拟合子图
    // :param x0: 拟合参数x0
    // :param y0: 拟合参数y0
    // :param sigma: 拟合参数sigma
    // :param para_A: 拟合参数A
    // :return:均方根误差
    int i,j;
    cv::Mat f;
    float rmse = 0;
    f = gaussian(dir.at<float>(0,3), dir.at<float>(0,1), dir.at<float>(0,2), dir.at<float>(0,0), pic.rows, pic.cols);
    f=f-pic;
    for(i=0;i<pic.rows;i++){
        for(j=0;j<pic.cols;j++){
                rmse+=f.at<float>(i,j);
        }
    }
    rmse = sqrt(rmse/(pic.rows*pic.cols+0.0));

    return rmse;
}
