#include "hicluster.h"

cv::Mat compute_distance(cv::Mat loc){
    /**********
     * @计算点之间的距离，loc为坐标
     ******/
    int n,i,j;
    cv::Mat distance;

    n=loc.rows;
    distance=cv::Mat(n,n,CV_32FC1,cv::Scalar(0));
    for(i=0;i<n-1;i++){
        for(j=i+1;j<n;j++){
            distance.at<float>(i,j)=sqrt((loc.at<float>(i,0)-loc.at<float>(j,0))*(loc.at<float>(i,0)-loc.at<float>(j,0))+(loc.at<float>(i,1)-loc.at<float>(j,1))*(loc.at<float>(i,1)-loc.at<float>(j,1)));
        }
    }

    return distance;
}

// cv::Mat hicluster(cv::Mat img,double r0){
//     /**************
//      * @聚类函数
//      * @输入一个图像img与初始半径，调用mlpack库的meanshift聚类算法，返回聚类后的结果
//      * @聚类结果存储为cv::Mat,每一维是一个类，为n*2的矩阵，第一个数为这一类的点的个数，第二行为中心点，第三行开始为这个类的坐标，其余为0
//      ********/
//
//
//     int i,j,num,len,most;
//     img.convertTo(img,CV_64FC1,1,0);
//     arma::mat data(reinterpret_cast<double*>(img.data),img.rows,img.cols);
//
//     //大于0的点的坐标
//     arma::uvec  meanShiftData_temp = arma::find(data>0);
//     arma::mat meanShiftData(meanShiftData_temp.n_rows,2);
//
//     for(i=0;i<meanShiftData_temp.n_rows;i++){
//         meanShiftData(i,1)=meanShiftData_temp(i) % data.n_cols;
//         meanShiftData(i,0)=(meanShiftData_temp(i) - (meanShiftData_temp(i) % data.n_cols))/data.n_cols;
//     }
//
//     //MeanShift聚类
//     mlpack::MeanShift<> meanShift(r0);
//     arma::Row<size_t> assignments;
//     arma::mat centroids;
//     meanShift.Cluster((arma::mat) trans(meanShiftData), assignments, centroids);
//     //centroids.print();
//
//     //对每个类进行计数
//     num=assignments.max();
//     arma::mat count_ass;
//     count_ass.zeros(num+1,1);
//     for(i=0;i<(int)assignments.n_elem;i++){
//         count_ass((int)assignments(i),0)++;
//     }
//     most=count_ass.max();
//
//     //将结果按规律输入进clu
//     arma::cube clu;
//     clu.zeros(num+1,most+2,2);
//     arma::uvec  temp1;
//     for(i=0;i<1+num;i++){
//         clu(i,0,0)=count_ass(i);
//         clu(i,1,0)=(int)centroids(0,i);
//         clu(i,1,1)=(int)centroids(1,i);
//         temp1 = arma::find(assignments==i);
//         len=temp1.n_elem;
//         for(j=0;j<len;j++){
//             clu(i,j+2,0)=meanShiftData((int)temp1(j),0);
//             clu(i,j+2,1)=meanShiftData((int)temp1(j),1);
//         }
//     }
//
//     //将cube类型的clu转化为最终的cv::Mat型输出结果
//     int size_out[3]={num+1,most+2,2};
//     cv::Mat out(3,size_out,CV_64FC1,cv::Scalar(0));
//     for(i=0;i<num+1;i++){
//         for(j=0;j<most+2;j++){
//             out.at<double>(i,j,0)=clu(i,j,0);
//             out.at<double>(i,j,1)=clu(i,j,1);
//         }
//     }
//
//     return out;
// }
