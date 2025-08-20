#ifndef DATAFORMAT_H__
#define DATAFORMAT_H__

#include <vector>
#include <list>
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <string>
#include <climits>
#include <cmath>
#include "keypoint.h"
#include "util/exception.h" 


namespace util {
    
typedef _point point2d;

inline float L2(const point2d& p1, const point2d& p2)
{
	float difx = p1.x-p2.x;
	float dify = p1.y-p2.y;
	
	return difx*difx+dify*dify;
}

inline float CalTriangleArea(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3)
{
    float s = ((p3.x-p1.x)*(p2.y-p1.y) - (p2.x-p1.x)*(p3.y-p1.y))*.5;

    return s>0 ? s : -s;
}


inline bool IsPolygonConvex(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3, const util::point2d& p4)
{
	float z1, z2, z3, z4;
 
	z1 = ((p2.x - p1.x) * (p4.y - p1.y) - (p4.x - p1.x) * (p2.y - p1.y));
	z2 = ((p4.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p4.y - p1.y));
	z3 = ((p4.x - p2.x) * (p3.y - p2.y) - (p3.x - p2.x) * (p4.y - p2.y));
	z4 = ((p3.x - p2.x) * (p1.y - p2.y) - (p1.x - p2.x) * (p3.y - p2.y));
 
	return (z1 * z2 > 0) && (z3 * z4 > 0);
}

inline float CalQuadrangleArea(const util::point2d& p1, const util::point2d& p2, const util::point2d& p3, const util::point2d& p4)
{
    return (CalTriangleArea(p1, p2, p3)+CalTriangleArea(p2, p3, p4)+CalTriangleArea(p1, p2, p4)+CalTriangleArea(p1, p3, p4))*.5;
}
// static void Add(FiducialStack fidstack1, FiducialStack fidstack2, FiducialStack* new_fidstack);
    
/** @brief Stack of fiducial markers*/
class FiducialStack {
private:
    int size;
    std::vector<util::point2d>* vfidp;
    int width, height;
    float ratio;

public:
    FiducialStack():size(0), vfidp(NULL), width(-1), height(-1), ratio(1) {}

    FiducialStack(int _size):size(_size) {
        vfidp = new std::vector<util::point2d>[size];
    }

    ~FiducialStack() {
        if(vfidp) {
            delete [] vfidp;
        }
    }

    void ReSize(int  _size) {
        if(size != _size) {
            delete [] vfidp;
            size = _size;
            vfidp = new std::vector<util::point2d>[size];
        }
    }
    void SetRatio(float r) {
        ratio = r;
    }

    void Clear() {
        delete [] vfidp;
        vfidp = new std::vector<util::point2d>[size];
    }

    void SetWxH(int _width, int _height) {
        width = _width;
        height = _height;
    }

    void Release() {
        if(vfidp) {
            delete [] vfidp;
        }
        vfidp = NULL;
    }

    std::vector<util::point2d>& V(int index) {
        assert(index>=0&&index<size);
        return vfidp[index];
    }
    const std::vector<util::point2d>& V(int index) const {
        assert(index>=0&&index<size);
        return vfidp[index];
    }

    int Size() const {
        return size;
    }
    int Width() const {
        return width;
    }
    int Height() const {
        return height;
    }
    float Ratio() const {
        return ratio;
    }
    
//     void Add(FiducialStack fidstack1, FiducialStack fidstack2){
//         std::cout<<"ok"<<std::endl;
// 
//     }
    void WriteFidsByFile(const char* filename) const;
    bool ReadFidsByFile(const char* filename);
    static void div(const std::vector< int >& m_index_high, const std::vector< int >& m_index_low, FiducialStack* fidstack_low, FiducialStack* fidstack_high, FiducialStack* fidstack);
    static void comb(const std::vector< int >& m_index_high, const std::vector< int >& m_index_low, FiducialStack* fidstack_low, FiducialStack* fidstack_high, FiducialStack* fidstack);
    static void Add(FiducialStack* fidstack1, FiducialStack* fidstack2, FiducialStack* new_fidstack);
};

// void Add(FiducialStack fidstack1, FiducialStack fidstack2, FiducialStack* new_fidstack);

typedef std::pair<util::_point, util::_point> pair;

struct img_match {
    int idx1, idx2;		//the index of first image and sencond image; which is matched
    std::vector<pair> pairs;
    size_t size() const {
        return pairs.size();
    }
};


/** @brief MatchPairStack; the output of keymatch process*/
class ImgMatchVector {
private:
    std::vector<img_match>* match_vector;
public:
    ImgMatchVector()
    {
        match_vector = new std::vector<img_match>();
    }

    ~ImgMatchVector()
    {
        if(match_vector) {
            delete match_vector;
        }
    }

    void Clear()
    {
        match_vector->clear();
    }

    void Release()
    {
        if(match_vector) {
            delete match_vector;
        }
        match_vector = NULL;
    }

    int Size() const {
        return match_vector->size();
    }

    img_match& operator[](int index)
    {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }

    const img_match& operator[](int index) const
    {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }

    img_match& V(int index) {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }
    const img_match& V(int index) const {
        assert(index>=0&&index<match_vector->size());
        return (*match_vector)[index];
    }

    img_match& MallocNewMatch()
    {
        img_match data;
        match_vector->push_back(data);
        return (*match_vector)[match_vector->size()-1];
    }

    void PushBack(const img_match& data)
    {
        match_vector->push_back(data);
    }
    
    img_match* GetMatchSetWithIdx(int idx1, int idx2, bool& noex){							//!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for(int i = 0; i < match_vector->size(); i++){
			if((*match_vector)[i].idx1 == idx1 && (*match_vector)[i].idx2 == idx2){
				noex = true;
				return &((*match_vector)[i]);
			}
			if((*match_vector)[i].idx1 == idx2 && (*match_vector)[i].idx2 == idx1){
				noex = false;
				return &((*match_vector)[i]);
			}
		}
		return NULL;
	}

    static void ReadPairs(img_match* pairs, std::istream& in);
    void CoordinateTransform(int width, int height);
    void PreRotate(float angle);		//angle in degree
    void PrintPairs(int index, std::ostream& o) const;
    void WriteVectorByFolder(const char* folderpath) const;
    void ReadVectorByFolder(const char* folderpath);
};

/** @brief the class node stores the point in an progect image and the relations between the points*/
struct _info {
    float angle;
    int z;
};

#define ISEQUAL(x, y)		fabs((x)-(y)) < 0.01

class NodePlane {
    friend class TrackSpace;
private:
	struct _node{
		_point p;
		_node* pre;
		_node* nxt;
		void* ref;
		long extra_idx;
		_node():pre(NULL),nxt(NULL),extra_idx(-1){
		}
		_node(const _point& pp):p(pp),pre(NULL),nxt(NULL),extra_idx(-1){}
		
		bool operator > (const _node& n) const{
			return p.x+p.y > n.p.x+n.p.y;
		}
		
		bool operator < (const _node& n) const{
			return p.x+p.y < n.p.x+n.p.y;
		}
		
		bool operator == (const _node& n) const{
			return ISEQUAL(p.x, n.p.x) && ISEQUAL(p.y, n.p.y);
		}
		
		bool operator != (const _node& n) const{
			return !(ISEQUAL(p.x, n.p.x) && ISEQUAL(p.y, n.p.y));
		}
	};
	
private:
	std::list<_node> node_list;
    void* information;
    std::list<_node>::iterator pre_visited;

public:
    class NodeReader {
        friend class NodePlane;
    private:
		std::list<_node>& node_list;
		std::list<_node>::iterator itr;
		std::list<_node>::iterator end;
		
    private:
        NodeReader(std::list<_node>& nlist):node_list(nlist){itr = node_list.begin(); end = node_list.end();}
        
    public:
        _node* Next()
        {
			_node* ptr;
			if(itr != end){
				ptr = &(*itr);
				itr++;
			}
			else{
				ptr = NULL;
			}
            
            return ptr;
        }
    };
	
public:
    NodePlane(){
        pre_visited = node_list.begin();
    }
    
    ~NodePlane(){}
	
    void SetInfo(void* info){
        information = info;
    }

    _node* AddNode(int x, int y){
        return AddNode(_point(x,y));
    }

    _node* AddNode(const _point& point);

    static void Link(_node* pnode1, _node* pnode2){
        _node* tmp;
        if(!pnode1->nxt) {
            tmp = pnode1;
        }
        else {
            tmp = pnode1;
            while(tmp->nxt && ((_info*)(tmp->nxt->ref))->z < ((_info*)(pnode2->ref))->z) {
                tmp = tmp->nxt;
            }
        }

        _node* tmp2;

        if(!pnode2->pre) {
            tmp2 = pnode2;
        }
        else {
            tmp2 = pnode2;
            while(tmp2->pre){
                tmp2 = tmp2->pre;
				if(*tmp2 == *pnode1){
					return;			//the path has already covered pnode1;
				}
				if(((_info*)(tmp2->ref))->z+3 < ((_info*)(pnode1->ref))->z){		//it is not a safe path
					return;
				}
            }
        }
//         tmp->nxt = tmp2;
//         tmp2->pre = tmp;
        
        if(!tmp->nxt && ((_info*)(tmp->ref))->z < ((_info*)(tmp2->ref))->z){
            tmp->nxt = tmp2;
            tmp2->pre = tmp;
        }
        else if(!tmp->nxt && ((_info*)(tmp->ref))->z > ((_info*)(tmp2->ref))->z
                && tmp->pre && ((_info*)(tmp->pre->ref))->z < ((_info*)(tmp2->ref))->z
                && tmp2->nxt && ((_info*)(tmp->ref))->z < ((_info*)(tmp2->nxt->ref))->z) {
            tmp->nxt = tmp2->nxt;
            tmp2->nxt->pre = tmp;
            tmp->pre->nxt = tmp2;
            tmp2->pre = tmp->pre;
            tmp2->nxt = tmp;
            tmp->pre = tmp2;
        }
        else if((((_info*)(tmp->ref))->z < ((_info*)(pnode2->ref))->z)
                && tmp->nxt && (((_info*)(tmp->nxt->ref))->z > ((_info*)(pnode2->ref))->z) && pnode2->pre == NULL) {
            pnode2->nxt = tmp->nxt;
            tmp->nxt->pre = pnode2;
            tmp->nxt = pnode2;
            pnode2->pre = tmp;
        }
        else{
// 			if(!(tmp->p.x == tmp2->p.x && tmp->p.y == tmp2->p.y)){
// 				std::cout<<((_info*)(tmp->ref))->z<<" "<<((_info*)(tmp2->ref))->z<<std::endl;
// 			}
		}
    }

    NodeReader GetReader()
    {
        NodeReader nreader(node_list);
        return nreader;
    }

public:
    int Print(std::ostream& o) const;
};

union pt2 {
    struct {
        double x;
        double y;
    };
    double v[2];
    const pt2& operator = (const pt2& pt_) {
        x = pt_.x;
        y = pt_.y;
        return *this;
    }
};

union pt3 {
    struct {
        double x;
        double y;
        double z;
    };
    double v[3];

    const pt3& operator = (const pt3& pt_) {
        x = pt_.x;
        y = pt_.y;
        z = pt_.z;
        return *this;
    }
};

class TrackSpace {
private:
    util::NodePlane* xyplane_array;	//xlist stores the x coordinate data;and, xyplane_array sotres the z coordinate data
    std::vector<pt3> pt3_v;
    _info* info_array;			//info_array sotres the z coordinate data
    int size;

public:
    class Iterator {
        friend class TrackSpace;
        friend class TrackNode;
    private:
        std::list<util::NodePlane::_node>* pnodes;
        std::vector<pt3>& pt3_v;
        std::list<util::NodePlane::_node>::iterator itr;
		
    private:
        Iterator(TrackSpace& trspce):pt3_v(trspce.pt3_v) {}

    public:
        const float& X()
        {
            return (*itr).p.x;
        }

        const float& Y()
        {
            return (*itr).p.y;
        }

        const int& Z()
        {
            return ((_info*)(*itr).ref)->z;
        }

        const float& Angle()
        {
            return ((_info*)(*itr).ref)->angle;
        }

        /** @brief ++:like i++*/
        Iterator operator++(int)
        {
            Iterator tmp = *this;
            itr++;
            return tmp;
        }

        /** @brief ++:like ++i*/
        Iterator& operator ++()
        {
            itr++;
            return *this;
        }

        /** @brief --:like i--*/
        Iterator operator--(int)
        {
            Iterator tmp = *this;
            itr--;
            return tmp;
        }

        /** @brief --:like --i*/
        Iterator& operator --()
        {
            itr--;
            return *this;
        }

        bool operator ==(const Iterator& it) const
        {
            return pnodes == it.pnodes && itr == it.itr;
        }

        bool operator !=(const Iterator& it) const
        {
            return pnodes != it.pnodes || itr != it.itr;
        }

        const Iterator& operator =(const Iterator& it)
        {
            pnodes = it.pnodes;
            itr = it.itr;
            return *this;
        }

        bool IsNULL() {
            return itr == pnodes->end();
        }

        void SetExtraAsZVectorIndex(long idx) {
            (*itr).extra_idx = idx;
        }
    };

    class TrackNode {
    private:
        util::NodePlane::_node* pnode;
        std::vector<pt3>& pt3_v;
    public:
        TrackNode(const Iterator& it):pt3_v(it.pt3_v)
        {
            pnode = &(*(it.itr));
        }

        TrackNode(const TrackNode& tracknode):pnode(tracknode.pnode), pt3_v(tracknode.pt3_v) {}

        bool IsBegin() const
        {
            return pnode->pre == NULL;
        }

        bool IsRBegin() const
        {
            return pnode->nxt == NULL;
        }

        const float& X() const
        {
            return pnode->p.x;
        }

        const float& Y() const
        {
            return pnode->p.y;
        }

        const int& Z() const
        {
            return ((_info*)(pnode->ref))->z;
        }

        const float& Angle() const
        {
            return ((_info*)(pnode->ref))->angle;
        }

        TrackNode operator++(int)
        {
            TrackNode tmp = *this;
            pnode = pnode->nxt;
            return tmp;
        }

        /** @brief ++:like ++i*/
        TrackNode& operator ++()
        {
            pnode = pnode->nxt;
            return *this;
        }

        /** @brief --:like i--*/
        TrackNode operator--(int)
        {
            TrackNode tmp = *this;
            pnode = pnode->pre;
            return tmp;
        }

        /** @brief --:like --i*/
        TrackNode& operator --()
        {
            pnode = pnode->pre;
            return *this;
        }

        const TrackNode& operator = (const TrackNode& tracknode) {
            pnode = tracknode.pnode;
            pt3_v = tracknode.pt3_v;
            return *this;
        }

        bool IsNULL() const {
            return !pnode;
        }

        void AddExtra(const double pt3_[3]) {
            pt3 tmp;
            tmp.x = pt3_[0];
            tmp.y = pt3_[1];
            tmp.z = pt3_[2];
            pt3_v.push_back(tmp);
            pnode->extra_idx = pt3_v.size()-1;
            NodePlane::_node* tmpp = pnode;
            while((tmpp = tmpp->nxt) != NULL) {
                tmpp->extra_idx = pnode->extra_idx;
            }
            tmpp = pnode;
            while((tmpp = tmpp->pre) != NULL) {
                tmpp->extra_idx = pnode->extra_idx;
            }
        }

        void AddExtraToNxt(const double pt3_[3]) {
            pt3 tmp;
            tmp.x = pt3_[0];
            tmp.y = pt3_[1];
            tmp.z = pt3_[2];
            pt3_v.push_back(tmp);
            pnode->extra_idx = pt3_v.size()-1;
            NodePlane::_node* tmpp = pnode;
            while((tmpp = tmpp->nxt) != NULL) {
                tmpp->extra_idx = pnode->extra_idx;
            }
        }

        void AddExtraToPre(const double pt3_[3]) {
            pt3 tmp;
            tmp.x = pt3_[0];
            tmp.y = pt3_[1];
            tmp.z = pt3_[2];
            pt3_v.push_back(tmp);
            pnode->extra_idx = pt3_v.size()-1;
            NodePlane::_node* tmpp = pnode;
            while((tmpp = tmpp->pre) != NULL) {
                tmpp->extra_idx = pnode->extra_idx;
            }
        }

        bool HasExtra() const {
            return pnode->extra_idx != -1;
        }

        void GetExtraAsZVectorIndex(long* idx) {
            *idx = pnode->extra_idx;
        }

        long VecIndex() {
            return pnode->extra_idx;
        }

        void GetExtra(double pt3_[3]) const {
            memcpy(pt3_, pt3_v[pnode->extra_idx].v, sizeof(double)*3);
        }
    };
public:
    TrackSpace(const ImgMatchVector& imv, std::vector<float> angles) {}
    TrackSpace():xyplane_array(NULL), info_array(NULL), size(0) {}
    ~TrackSpace() {
        if(xyplane_array) {
            Release();
        }
    }

    void Release() {
        assert(xyplane_array);
        delete [] xyplane_array;
        delete [] info_array;
        xyplane_array = NULL;
        info_array = NULL;
        pt3_v.clear();
    }

    int Size() const {
        return size;
    }

    float Z_Angle(int z) const {
        return info_array[z].angle;
    }

    void InsertMatchVector(const ImgMatchVector& imv)
	{
		 for(int i = 0; i < imv.Size(); i++){
            NodePlane::_node* pnode1 = NULL, * pnode2 = NULL;
            int idx1 = imv[i].idx1;
            int idx2 = imv[i].idx2;
			
            for(int j = 0; j < imv[i].size(); j++){
                pnode1 = xyplane_array[idx1].AddNode((imv[i].pairs)[j].first);
                pnode2 = xyplane_array[idx2].AddNode((imv[i].pairs)[j].second);
                if(idx1 < idx2){						//WARNING
                    NodePlane::Link(pnode1, pnode2);
                }
                else if(idx1 > idx2){
                    NodePlane::Link(pnode2, pnode1);
                }
            }
        }
	}
    
    void Create(const ImgMatchVector& imv, const std::vector<float>& angles)
    {
        assert(!xyplane_array);
        EX_TIME_BEGIN("Creating TrackSpace")

        size = angles.size();
        xyplane_array = new util::NodePlane[size];		//create z coordinate
        info_array = new _info[size];

        for(int i = 0; i < size; i++) {
            info_array[i].angle = angles[i];
            info_array[i].z = i;
            xyplane_array[i].SetInfo(&(info_array[i]));
        }
        
        InsertMatchVector(imv);

        EX_TIME_END("Creating TrackSpace")
    }

    Iterator Z_Iterator(int z)
    {
        assert(z >= 0&&z < size);
        Iterator citr(*this);
        citr.pnodes =  &(xyplane_array[z].node_list);
        citr.itr = xyplane_array[z].node_list.begin();
        citr.pt3_v = pt3_v;
        return citr;
    }

    int Z_Size(int z)
    {
        assert(z >= 0&&z < size);
        return xyplane_array[z].node_list.size();
    }

    int PrintXYPlane(int index, std::ostream& o) const;
	void CoordinateTransform(int width, int height);
	void PreRotate(float angle);		//angle in degree
	static void WriteIMODfidModel(util::TrackSpace& trackspace, int width, int height, int num_frames, const char* filename);
};

static bool ReadAnglesByName(const char* name, std::vector<float>* angles)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }

    while(in.good()) {
        float angle;
        in>>angle;
        if(in.fail()) {
            break;
        }
        angles->push_back(angle);
    }
    in.close();
    return true;
}

}
#endif
