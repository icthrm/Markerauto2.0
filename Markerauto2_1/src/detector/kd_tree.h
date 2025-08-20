#include "detector.h"
// #include "mrcimg/mrc2img.h"


// using namespace std;
// using namespace cv;

class Node_kd{

public:
    Node_kd();
    ~Node_kd();

    int dimension;
    cv::Mat x;  //特征坐标，切分的坐标点（坐标集合中中间那个坐标）
    Node_kd *left = NULL;  // 左孩子
    Node_kd *right = NULL;  //右孩子
    Node_kd *parents = NULL;  //父节点
    bool flag_self = false;  //是否被访问过，标记

    void setleft(Node_kd *left);
    void clear_flag(Node_kd *node);  //将标记清零
    void construct(int d,cv::Mat data,Node_kd *node,int layer);  //构建kd树
    float distance(cv::Mat a,cv::Mat b);  //欧氏距离

    cv::Mat search(Node_kd *node,cv::Mat p,cv::Mat L,int K);


};
