#include"kd_tree.h"


inline void sortMat(cv::Mat &stats, int colId){
    //根据指定列以行为单位排序

    cv::Mat sorted_index;
    cv::sortIdx(stats(cv::Range(0,stats.rows),cv::Range(colId,colId+1)), sorted_index, cv::SORT_EVERY_COLUMN + cv::SORT_ASCENDING);

    cv::Mat sorted_stats = stats.clone();
    int row_num = sorted_index.rows;
    for(int i = 0; i < row_num; i++){
        int _idx = sorted_index.at<int>(i, 0);
        sorted_stats.row(i) = stats.row(_idx) + 0;//必须加0否则会出很难debug的错误
    }
    stats = sorted_stats;
    return;
}

Node_kd::Node_kd(){
    Node_kd *left = NULL;  // 左孩子
    Node_kd *right = NULL;  //右孩子
    Node_kd *parents = NULL;  //父节点
    bool flag_self = false;  //是否被访问过，标记
}

void Node_kd::setleft(Node_kd *left1){
    this->left=left1;
}
Node_kd::~Node_kd(){

}

void Node_kd::clear_flag(Node_kd *node){
    node->flag_self=false;
    if(node->left) Node_kd::clear_flag(node->left);
    if(node->right) Node_kd::clear_flag(node->right);
}

void Node_kd::construct(int d,cv::Mat data,Node_kd *node,int layer){
    /***********
    :type d: int
    d是向量的维数
    :type data: cv::Mat
    data是所有向量构成的列表
    :type node: Node_kd
    node是当前进行运算的结点
    :type layer: int
    layer是当前kd树所在层数
    *****************/
    int middle;
    cv::Mat dataleft,dataright;
    Node_kd *left_node,*right_node;

    node->dimension = layer % d;  // 防止维数越界
    // 如果只有一个元素，说明到了叶子结点，该分支结束
    if (data.rows == 1){
        node->x = data;  // 该data中间那个就是唯一的一个坐标
        return;
    }
    if (data.rows == 0){
        //没有代表的数据就作为一个空叶子结点
        return;
    }

    //1,data中的数据按layer%N维进行排序
    sortMat(data,layer % d);

    //2,计算中间点的索引，偶数则取中间两位中较大的一位,记为该结点的特征坐标
    middle = (int)data.rows / 2 ; // 除法取整
    node->x = data(cv::Range(middle,middle+1),cv::Range(0,data.cols));


    //3，划分data
    dataleft = data(cv::Range(0,middle),cv::Range(0,data.cols));
    dataright = data(cv::Range(middle+1,data.rows),cv::Range(0,data.cols));


   // cout<<dataleft<<endl<<dataright<<endl;


    //4,左孩子结点
    left_node=new Node_kd;
    node->left = left_node;
    node->left->parents = node;
    construct(d, dataleft, left_node, layer + 1);
    //cout<<node->left->x<<endl;

    //5，右孩子结点
    right_node=new Node_kd;
    node->right = right_node;
    node->right->parents = node;
    construct(d, dataright, right_node, layer + 1);
}

float Node_kd::distance(cv::Mat a,cv::Mat b){
    int i;
    float dis = 0;

    for(i=0;i<a.cols;i++){
        dis += (a.at<float>(0,i) - b.at<float>(0,i)) * (a.at<float>(0,i) - b.at<float>(0,i));
    }
    return sqrt(dis);
}

void change_L(cv::Mat &L,cv::Mat x,cv::Mat p,int K){
    /**************
    判断并进行是否将该点加入近邻点列表
    :type L: Mat
    L是近邻点列表
    :type x: Mat
    x是判断是否要加入近邻列表的向量
    :type p: Mat
    p是目标向量
    :type K:int
    K是近邻列表的最大元素个数
    ********************/
    //cout<<"input"<<L<<endl;
    int i,j,index;
    Node_kd dis;
    cv::Mat newL,dislist;
    if (L.rows < K){
        newL=cv::Mat::zeros(L.rows+1,x.cols,CV_32FC1);
        for(i=0;i<L.rows;i++){
            for(j=0;j<x.cols;j++){
                newL.at<float>(i,j)=L.at<float>(i,j);
            }
        }
        newL.row(L.rows)=x+0.0;
        L=newL;

        //cout<<"output"<<L<<endl;
        return;
    }

    dislist = cv::Mat::zeros(K,1,CV_32FC1);
    for (i=0;i<K;i++){
        dislist.at<float>(i,0)=dis.distance(p, L(cv::Range(i,i+1),cv::Range(0,L.cols)));
    }
    index=0;
    for(i=1;i<K;i++){
        if(dislist.at<float>(index,0)<dislist.at<float>(i,0))
            index=i;
    }
    if (dis.distance(p, x) < dislist.at<float>(index,0))
        // 若x和p之间的距离小于L到p中最远的点，就用x替换此最远点
        L.row(index) = x+0.0;
    //L=newL;
    //cout<<"output"<<L<<endl;

    return ;
}

cv::Mat Node_kd::search(Node_kd *node,cv::Mat p,cv::Mat L,int K){
    /************************
    :type List: Mat
    :type node: Node_kd
    类Node，整个树的框架，里面包含父子结点信息，以及每个父子结点含有的坐标点
    :type p: Mat
    目标坐标
    :type L: Mat
    L为有k个座位的列表，用于保存已搜寻到的最近点
    :type K: int
    K为近邻个数
    *************************/
    int i;
    float max_dislist,temp;
    cv::Mat dislist,out;
    Node_kd *n,dis;

    //根据p的坐标值和每个点的切分轴向下搜索,先到达底部结点
    n=new Node_kd;
    n = node;  // 用n来记录结点的位置，先从顶部开始,直到叶子结点，循环完的n为叶子节点
    //cout<<"111n "<<n->x<<endl;
    //cout<<"dimension "<<n->dimension<<endl;
    while (true){
        // 若到达了叶子结点则退出循环
        if ((n->left == NULL) && (n->right == NULL))
            break;
        //cout<<n->dimension<<endl;
        if (n->x.at<float>(0,n->dimension) > p.at<float>(0,n->dimension)){
            //cout<<n->x.at<float>(0,n->dimension)<<">"<<p.at<float>(0,n->dimension)<<endl;
            n = n->left;
            //cout<<n->x<<" left"<<endl;
        }else{
            //cout<<n->x.at<float>(0,n->dimension)<<">"<<p.at<float>(0,n->dimension)<<endl;
            n = n->right;
            //cout<<n->x<<" right"<<endl;
        }

    }

    n->flag_self = true;  // 标记为已访问过
    if(!(n->x.rows == 0))  // 若为空叶子结点，则不必记录数值
        change_L(L, n->x, p, K);  // 若符合插入条件，就插入，不符合就不插入

    //cout<<"n.x "<<n->x<<endl;
    while(true){
        // 若当前结点是根结点则输出L算法完成
        if(n->parents == NULL){
            if (L.rows < K)
                std::cout<<"K值超过数据总量"<<std::endl;
            return L;
        }else{
            //cout<<L<<endl;
            //当前结点不是根结点，向上爬一格
            n = n->parents;
            //cout<<n->x<<endl;
            while (n->flag_self == true){
                //若当前结点被访问过，就一直向上爬，到没被访问过的结点为止
                //若向上爬时遇到了已经被访问过的根结点，说明另一边已经搜索过了搜索结束
                if ((n->parents == NULL) && (n->flag_self)){

                    if (L.rows < K)
                        std::cout<<"K值超过数据总量1"<<std::endl;
                    return L;
                }
                n = n->parents;
            }
            //此时n未被访问过,将其标记为访问过
            n->flag_self = true;

            //如果此时 L 里不足 k 个点，则将节点特征加入 L；
            //如果 L 中已满 k 个点，且当前结点与 p 的距离小于与L的最大距离，
            //则用节点特征替换掉 LL 中离最远的点。
            change_L(L, n->x, p, K);
            // 计算p和当前节点切分线的距离。如果该距离小等于于 LL 中最远的距离或者 LL 中不足 kk 个点，则切分线另一边或者 切分线上可能有更近的点，因此在当前节点的另一个枝从 (一) 开始执行。

            dislist = cv::Mat::zeros(L.rows,1,CV_32FC1);
            for (i=0;i<L.rows;i++){
                dislist.at<float>(i,0)=dis.distance(p, L(cv::Range(i,i+1),cv::Range(0,L.cols)));
            }

            max_dislist=0;
            for(i=0;i<dislist.rows;i++){
                if(dislist.at<float>(i,0) > max_dislist)
                    max_dislist=dislist.at<float>(i,0);
            }
            temp=p.at<float>(n->dimension,0) - n->x.at<float>(n->dimension,0);
            temp=temp*temp;
            if ((temp < max_dislist) || (L.rows < K)){
                if (n->left->flag_self == false){
                    out=dis.search(n->left, p, L, K);
                    return out;
                }else{
                    out=dis.search(n->right, p, L, K);
                    return out;
                }
            //如果该距离大于等于 L 中距离 p 最远的距离并且 L 中已有 k 个点，则在切分线另一边不会有更近的点，重新执行(三)
            }
        }
    }
    return L;
}
