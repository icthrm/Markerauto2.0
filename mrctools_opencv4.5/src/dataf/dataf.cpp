#include "dataf.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iostream>
#include "keypoint.h"
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>

void util::FeatureStack::WriteStackByFolder(const char* folderpath) const
{
    if(access(folderpath,0) == -1) {		//create file folder
        mkdir(folderpath,0777);
    }
    std::ostringstream ooo;
    ooo <<folderpath<<"/attributes";
    std::ofstream out(ooo.str().c_str());
    out<<"W:"<<width<<" H:"<<height<<" R:"<<ratio<<"\n";
    out<<"Z:"<<size<<"\n";
    out.close();
    for(int i = 0; i < size; i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        try {
            util::WriteFeaturesByName(vfeatp[i], oss.str().c_str());
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

void util::FeatureStack::ReadStackByFolder(const char* folderpath)
{
    std::ostringstream iii;
    iii <<folderpath<<"/attributes";
    std::ifstream in(iii.str().c_str());
    char ch;
    if(!in.good()) {
        ex::EX_THROW("Can't Open folder");
    }
    in>>ch>>ch>>width>>ch>>ch>>height>>ch>>ch>>ratio;
    in>>ch>>ch>>size;
    in.close();
    if(vfeatp) {
        delete [] vfeatp;
    }
    vfeatp = new std::vector<util::feature>[size];

    for(int i = 0; i < size; i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        try {
            util::ReadFeaturesByName(&vfeatp[i], oss.str().c_str());
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

void util::FiducialStack::WriteFidsByFile(const char* filename) const {
    std::ofstream out(filename);
    if(!out.good()) {
        ex::EX_THROW("Can't Create File");
    }
    std::cout <<std::setprecision(8)<<std::endl;
    out<<width<<"\t"<<height<<"\t"<<size<<"\t"<<ratio<<std::endl;
    for(int i = 0 ; i < size ; i++) {
        out<<"frame "<<i<<std::endl;
        out<<"num = "<<vfidp[i].size() <<std::endl;
        for(int j = 0 ; j < vfidp[i].size() ; j++) {
            float x = vfidp[i][j].x;
            float y = vfidp[i][j].y;
            out<<x<<"\t"<<y<<std::endl;
        }
    }
    return;
}

bool util::FiducialStack::ReadFidsByFile(const char* filename) {
    std::string s;
    std::stringstream ss;
    std::string first_str;
    char ch;

    std::ifstream fin(filename);

    if(!fin.good()) {
        std::cout<<"Unable to open "<<filename<<std::endl;
        return false;
    }

    getline(fin , s);

    ss.str("");
    ss<<s;
    ss>>width>>height>>size>>ratio;

    if(vfidp) {
        delete [] vfidp;
    }
    vfidp = new std::vector<util::point2d>[size];

    while(getline(fin ,s)) {
        ss.clear();
        ss<<s;
        ss>>first_str;
        int feats_num;
        if(first_str == "frame") {
            int frame;
            ss>>frame;
            getline(fin ,s);
            ss.clear();
            ss<<s;
            ss>>first_str;

            if(first_str == "num") {
                ss>>ch>>feats_num;
                for(int i=0; i< feats_num ; i++) {
                    util::point2d pt;
                    fin>>pt.x>>pt.y;
                    vfidp[frame].push_back(pt);
                }
            }
        }
    }
    return true;
}

void util::ImgMatchVector::ReadPairs(img_match* imatch, std::istream& in)
{
    char ch;
    in>>imatch->idx1>>imatch->idx2;
    while(in.peek()!='#' && in.good()) {
        pair p;
        in>>ch>>p.first.x>>ch>>p.first.y>>ch>>ch>>ch>>p.second.x>>ch>>p.second.y>>ch;
        imatch->pairs.push_back(p);
        in.ignore(4096, '\n');
    }
}


void util::ImgMatchVector::PrintPairs(int index, std::ostream& o) const
{
    o<<(*match_vector)[index].idx1<<" "<<(*match_vector)[index].idx2<<std::endl;
    for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
        o<<"("<<((*match_vector)[index].pairs)[i].first.x<<","<<((*match_vector)[index].pairs)[i].first.y<<")&"
         <<"("<<((*match_vector)[index].pairs)[i].second.x<<","<<((*match_vector)[index].pairs)[i].second.y<<")\n";
    }
    o<<"#";
}

void util::ImgMatchVector::CoordinateTransform(int width, int height)
{
    for(int index = 0; index < (*match_vector).size(); index++) {
        for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
            ((*match_vector)[index].pairs)[i].first.x -= width/2;
            ((*match_vector)[index].pairs)[i].first.y -= height/2;
            ((*match_vector)[index].pairs)[i].second.x -= width/2;
            ((*match_vector)[index].pairs)[i].second.y -= height/2;
        }
    }
}

void util::ImgMatchVector::PreRotate(float angle)
{
    for(int index = 0; index < (*match_vector).size(); index++) {
        for(int i = 0; i < (*match_vector)[index].pairs.size(); i++) {
            ((*match_vector)[index].pairs)[i].first.x = cos(angle)*((*match_vector)[index].pairs)[i].first.x-sin(angle)*((*match_vector)[index].pairs)[i].first.y;
            ((*match_vector)[index].pairs)[i].first.y = sin(angle)*((*match_vector)[index].pairs)[i].first.x+cos(angle)*((*match_vector)[index].pairs)[i].first.y;
            ((*match_vector)[index].pairs)[i].second.x = cos(angle)*((*match_vector)[index].pairs)[i].second.x-sin(angle)*((*match_vector)[index].pairs)[i].second.y;
            ((*match_vector)[index].pairs)[i].second.y = sin(angle)*((*match_vector)[index].pairs)[i].second.x+cos(angle)*((*match_vector)[index].pairs)[i].second.y;
        }
    }
}

void util::ImgMatchVector::WriteVectorByFolder(const char* folderpath) const
{
    if(access(folderpath,0) == -1) {		//create file folder
        mkdir(folderpath,0777);
    }
    std::ostringstream ooo;
    ooo <<folderpath<<"/attributes";
    std::ofstream out(ooo.str().c_str());
    out<<"Z:"<<match_vector->size()<<"\n";
    out.close();
    for(int i = 0; i < match_vector->size(); i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ofstream o(oss.str().c_str());
        try {
            PrintPairs(i, o);
        } catch(ex::Exception& e) {
            EX_TRACE("%s\n", e.Msg())
        }
    }
}

void util::ImgMatchVector::ReadVectorByFolder(const char* folderpath)
{
    Clear();
    std::cout <<std::setprecision(8)<<std::endl;
    std::ostringstream iii;
    iii <<folderpath<<"/attributes";
    std::ifstream in(iii.str().c_str());
    char ch;
    int _size;
    in>>ch>>ch>>_size;
    in.close();

    for(int i = 0; i < _size; i++) {
        std::ostringstream oss;
        oss <<folderpath<<"/"<<i<<".txt";
        std::ifstream in(oss.str().c_str());
        if(!in.good()) {
            ex::EX_THROW("Can't Open File");
        }
        img_match& imatch = MallocNewMatch();
        ReadPairs(&imatch, in);
        in.close();
    }
}

util::NodePlane::~NodePlane()
{
    _lnode* ptr = xlist_begin.nxt;
    while(ptr != &xlist_end) {
        if(ptr->ptr) {
            _node* cur = ptr->ptr;
            while(cur) {
                _node* tmp = cur;
                cur = cur->_m2234;
                delete tmp;
            }
        }
        _lnode* tmp = ptr;
        ptr = ptr->nxt;
        delete tmp;
    }
    
    xlist_begin.nxt = &xlist_end;
    xlist_end.pre = &xlist_begin;
}

util::NodePlane::_node* util::NodePlane::BaseAddNode(const util::_point& point)
{
//     assert(point.x >= 0 && point.y >= 0);
    _lnode* ptr = (_m2307->x<=point.x?_m2307:&xlist_begin), * tmpp;
    _node* cur;
    while(ptr != &xlist_end) {
        if(ISEQUAL(ptr->x, point.x)) {
            _m2307 = ptr;
            if(ptr->ptr->p.y > point.y) {
                _node* tmpp = ptr->ptr;
                ptr->ptr = new _node();
                ptr->ptr->p = point;
                ptr->ptr->_m2234 = tmpp;
                ptr->ptr->ref = this->information;
                ptr->ptr->pre = NULL;
                ptr->ptr->nxt = NULL;
                return ptr->ptr;
            }
            else if(ISEQUAL(ptr->ptr->p.y, point.y)) {
                return ptr->ptr;
            }
            else {
                cur = ptr->ptr;
            }
            while(cur->_m2234 != NULL) {
                if(cur->_m2234->p.y > point.y) {
                    _node* tmpp = cur->_m2234;
                    cur->_m2234 = new _node();
                    cur->_m2234->p = point;
                    cur->_m2234->_m2234 = tmpp;
                    cur->_m2234->ref = this->information;
                    cur->_m2234->pre = NULL;
                    cur->_m2234->nxt = NULL;
                    return cur->_m2234;
                }
                else if(cur->_m2234->p.y < point.y) {
                    cur = cur->_m2234;
                }
                else {
                    return cur->_m2234;
                }
            }
            cur->_m2234 = new _node();
            cur = cur->_m2234;
            cur->p = point;
            cur->_m2234 = NULL;
            cur->ref = this->information;
            cur->pre = NULL;
            cur->nxt = NULL;
            return cur;
        }
        else if(ptr->x < point.x &&ptr->nxt->x <= point.x) {
            ptr = ptr->nxt;
        }
        else if(ptr->x < point.x &&point.x < ptr->nxt->x){
			cur = new _node(point);
			
            tmpp = ptr->nxt;
            _m2307 = ptr->nxt = new _lnode(cur);
//             ptr->nxt->x = point.x;
            ptr->nxt->pre = ptr;
            ptr->nxt->nxt = tmpp;
            tmpp->pre = ptr->nxt;
//             cur = ptr->nxt->ptr = new _node();
//             cur->p = point;
            cur->_m2234 = NULL;
            cur->pre = NULL;
            cur->nxt = NULL;
            cur->ref = this->information;
            return cur;
        }
    }
}

util::NodePlane::_node* util::NodePlane::AddNode(const util::_point& point)
{
    if(!pre_visited) {
        pre_visited = BaseAddNode(point);
        return pre_visited;
    }
    if(ISEQUAL(point.x, pre_visited->p.x) && ISEQUAL(point.y, pre_visited->p.y)) {
        return pre_visited;
    }
    else if(ISEQUAL(point.x, pre_visited->p.x) && point.y > pre_visited->p.y) {
        _node* cur = pre_visited;
        while(cur->_m2234 != NULL) {
			if(ISEQUAL(point.y, cur->_m2234->p.y)){
				pre_visited = cur->_m2234;
                return cur->_m2234;
			}
            else if(cur->_m2234->p.y > point.y) {
                _node* tmpp = cur->_m2234;
                cur->_m2234 = new _node();
                cur->_m2234->p = point;
                cur->_m2234->_m2234 = tmpp;
                cur->_m2234->ref = this->information;
                cur->_m2234->pre = NULL;
                cur->_m2234->nxt = NULL;

                pre_visited = cur->_m2234;
                return cur->_m2234;
            }
            else if(cur->_m2234->p.y < point.y){
                cur = cur->_m2234;
            }
//             else{
// 				pre_visited = cur->_m2234;
//                 return cur->_m2234;
// 			}
        }
        cur->_m2234 = new _node();
        cur = cur->_m2234;
        cur->p = point;
        cur->_m2234 = NULL;
        cur->ref = this->information;
        cur->pre = NULL;
        cur->nxt = NULL;
        pre_visited = cur;
        return cur;
    }
    else {
        pre_visited = BaseAddNode(point);
        return pre_visited;
    }
}

int util::NodePlane::Print(std::ostream& o) const
{
    int count = 0;
    _lnode* ptr = xlist_begin.nxt;
    while(ptr != &xlist_end) {
        o<<"#"<<ptr->x<<": ";
        _node* cur = ptr->ptr;
        while(cur != NULL) {
            assert(ISEQUAL(cur->p.x, ptr->x));
            o<<cur->p.y<<" ";
            count++;
            cur = cur->_m2234;
        }
        o<<std::endl;
        ptr = ptr->nxt;
    }
    o<<std::endl;
    return count;
}

int util::TrackSpace::PrintXYPlane(int index, std::ostream& o) const
{
    assert(index >= 0 && index < size);
    return xyplane_array[index].Print(o);
}

void util::TrackSpace::CoordinateTransform(int width, int height)
{
	float w_2 = width*.5, h_2 = height*.5;
    for(int i = 0; i < size; i++) {
        NodePlane::NodeReader reader = xyplane_array[i].GetReader();
        NodePlane::_node* pnode;
        while((pnode = reader.Next()) != NULL) {
			pnode->p.x -= w_2;
			pnode->p.y -= h_2;
        }
    }
}

void util::TrackSpace::PreRotate(float angle)
{
	float cosa = cos(angle), sina = sin(angle);
	for(int i = 0; i < size; i++) {
        NodePlane::NodeReader reader = xyplane_array[i].GetReader();
        NodePlane::_node* pnode;
        while((pnode = reader.Next()) != NULL) {
			_point p = pnode->p;
			pnode->p.x = cosa*p.x-sina*p.y;
			pnode->p.y = sina*p.x+cosa*p.y;
        }
    }
}

//add by wanglianshan ---chongxie
void util::TrackSpace::WriteIMODfidModel(util::TrackSpace& trackspace, int width, int height, int num_frames, const char* filename)
{
    std::vector<std::vector<pt3> > T;
    int zerotilt = num_frames / 2 + num_frames % 2 - 1;
    const unsigned int minTrajectoryLength =(10 > (int)(num_frames * 0.1)) ? 10 : (int)(num_frames * 0.1);

    for(int i = 0; i < trackspace.Size(); i++) {
        util::TrackSpace::Iterator itr = trackspace.Z_Iterator(i);

        while(!itr.IsNULL()) {

            util::TrackSpace::TrackNode node_itr = util::TrackSpace::TrackNode(itr);
            std::vector<pt3> Points ;
            bool isStart = node_itr.IsBegin();
            if(isStart) {
                while(!node_itr.IsNULL()) {
                    pt3 point3d;
                    point3d.x = node_itr.X();
                    point3d.y = node_itr.Y();
                    point3d.z = node_itr.Z();
                    Points.push_back(point3d);
                    node_itr++;
                }
            }

            if(Points.size() > minTrajectoryLength) {
                T.push_back(Points);
            }

            Points.clear();

            itr++;
        }
    }

    std::ofstream out(filename);
    out << "imod 1" << std::endl;
    out << "max " << width << " " << height << " " << num_frames<< std::endl;
    out << "offsets 0 0 0" << std::endl;
    out << "angles 0 0 0" << std::endl;
    out << "scale 1 1 1" << std::endl;
    out << "mousemode 1" << std::endl;
    out << "drawmode 1" << std::endl;
    out << "b&w_level 0,200" << std::endl;
    out << "resolution 3" << std::endl;
    out << "threshold 128" << std::endl;
    out << "pixsize 1" << std::endl;
    out << "units pixels" << std::endl;
    out << "symbol circle"<< std::endl;//so imod displays all the markers at once
    out << "size 7"<<std::endl;

    out << "object 0 " << T.size() << " 0" << std::endl;
    out << "name" << std::endl;
    out << "color 0 1 0 0" << std::endl;
    out << "open" << std::endl;
    out << "linewidth 1" << std::endl;
    out << "surfsize  0" << std::endl;
    out << "pointsize 0" << std::endl;
    out << "axis      0" << std::endl;
    out << "drawmode  1" << std::endl;

    for(int i=0 ; i< T.size() ; i++) {
        out << "contour " << i << " 0 " << T[i].size() << std::endl;
        for(int j =0; j< T[i].size() ; j++) {
            out << T[i][j].x<< " " <<T[i][j].y<< " " << T[i][j].z<< std::endl;
        }
    }
}
