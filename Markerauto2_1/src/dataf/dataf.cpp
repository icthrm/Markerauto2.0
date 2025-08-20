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

// void util::FiducialStack::Add(util::FiducialStack fidstack1, util::FiducialStack fidstack2)
// {
//     int size1, size2;
//     size1=fidstack1.Size();
//     size2=fidstack2.Size();
//     size=size1+size2;
//     for(int i=0;i<size1;i++)
//     {
//         std::vector<util::point2d>& fids = vfidp[i];
//         for(int j=0;j<fidstack1.V(i).size();j++)
//         {
//             util::point2d fid;
//             fid.x=fidstack1.V(i)[j].x;
//             fid.y=fidstack1.V(i)[j].y;
//             fids.push_back(fid);
//         }
//         
//     }
//     for(int i=0;i<size2;i++)
//     {
//         std::vector<util::point2d>& fids = vfidp[i];
//         for(int j=0;j<fidstack2.V(i).size();j++)
//         {
//             util::point2d fid;
//             fid.x=fidstack1.V(i)[j].x;
//             fid.y=fidstack1.V(i)[j].y;
//             fids.push_back(fid);
//         }
//         
//     }
// //     std::vector<util::point2d>& fids = vfidp[1];
// }
// void util::FiducialStack::A(int a, int b)
// {
//     std::cout<<"a+b"<<a+b<<std::endl;
// }
void util::FiducialStack::div(const std::vector<int>& m_index_high, const std::vector<int>& m_index_low, util::FiducialStack* fidstack_low, util::FiducialStack* fidstack_high, util::FiducialStack* fidstack)
{
    fidstack_low->ReSize(m_index_low.size());
    int num=fidstack->Size()-m_index_low.size();
    fidstack_high->ReSize(num);
    int num_low=0, num_high=0;
    for(int i=0;i<fidstack->Size();i++){
//         std::vector<util::point2d>& fids = fidstack_low->V(num_low);
//         std::vector<util::point2d>& fids2 = fidstack_high->V(num_high);
//         std::vector<int>::iterator itr = std::find(m_index.begin(), m_index.end(), i);
        if(std::count(m_index_low.begin(), m_index_low.end(), i)){
            std::vector<util::point2d>& fids = fidstack_low->V(num_low);
            for(int j=0;j<fidstack->V(i).size();j++)
            {
                util::point2d fid;
                fid.x=fidstack->V(i)[j].x;
                fid.y=fidstack->V(i)[j].y;
                fids.push_back(fid);
            }
//             std::cout<<"num_low:"<<num_low<<std::endl;
            num_low++;
        }
        else{
            std::vector<util::point2d>& fids2 = fidstack_high->V(num_high);
            for(int j=0;j<fidstack->V(i).size();j++)
            {
                util::point2d fid;
                fid.x=fidstack->V(i)[j].x;
                fid.y=fidstack->V(i)[j].y;
                fids2.push_back(fid);
            }
//             std::cout<<"num_high:"<<num_high<<std::endl;
            num_high++;
        }
    }
    
}

void util::FiducialStack::comb(const std::vector<int>& m_index_high, const std::vector<int>& m_index_low, util::FiducialStack* fidstack_low, util::FiducialStack* fidstack_high, util::FiducialStack* fidstack)
{
    int num;
    num=fidstack_high->Size()+fidstack_low->Size();
    fidstack->ReSize(num);
    int num_low=0, num_high=0;
    for(int i=0;i<fidstack->Size();i++){
        if(std::count(m_index_low.begin(), m_index_low.end(), i)){
            std::vector<util::point2d>& fids = fidstack_low->V(num_low);
            for(int j=0;j<fids.size();j++)
            {
                util::point2d fid;
                fid.x=fids[j].x;
                fid.y=fids[j].y;
                fidstack->V(i).push_back(fid);
            }
//             std::cout<<"num_low:"<<num_low<<std::endl;
            num_low++;
        }
        else{
            std::vector<util::point2d>& fids2 = fidstack_high->V(num_high);
            for(int j=0;j<fids2.size();j++)
            {
                util::point2d fid;
                fid.x=fids2[j].x;
                fid.y=fids2[j].y;
                fidstack->V(i).push_back(fid);
            }
//             std::cout<<"num_high:"<<num_high<<std::endl;
            num_high++;
        }
    }
    
}




void util::FiducialStack::Add(util::FiducialStack* fidstack1, util::FiducialStack* fidstack2, util::FiducialStack *new_fidstack)
{
    int size1, size2;
//     std::cout<<"ok"<<std::endl;
    size1=fidstack1->Size();
    size2=fidstack2->Size();
    new_fidstack->ReSize(size1+size2);
//     std::cout<<"fid_size1:"<<fidstack1.Size()<<std::endl;
//     std::cout<<"fid_size2:"<<fidstack2.Size()<<std::endl;
//     std::cout<<"fid_size:"<<new_fidstack->Size()<<std::endl;
    int num=0;
    int center=0.5*fidstack2->Size();
    for(int i=0;i<new_fidstack->Size();i++)
    {
        std::vector<util::point2d>& fids = new_fidstack->V(i);
//         std::vector<util::point2d>::iterator itr;
//         for(itr=fidstack1.begi;)
        if(i<0.5*fidstack2->Size())
        {
            for(int j=0;j<fidstack2->V(i).size();j++)
            {
                util::point2d fid;
                fid.x=fidstack2->V(i)[j].x;
                fid.y=fidstack2->V(i)[j].y;
                fids.push_back(fid);
            }
        }
        else{
            if(num<fidstack1->Size()){
                for(int j=0;j<fidstack1->V(num).size();j++)
                {
                    util::point2d fid;
                    fid.x=fidstack1->V(num)[j].x;
                    fid.y=fidstack1->V(num)[j].y;
                    fids.push_back(fid);
                }
                num++;
            }
            else{
                for(int j=0;j<fidstack2->V(center).size();j++)
                {
                    util::point2d fid;
                    fid.x=fidstack2->V(center)[j].x;
                    fid.y=fidstack2->V(center)[j].y;
                    fids.push_back(fid);
                }
                center++;
            }
                
        }
//         if(i<fidstack1->Size()){
//             for(int j=0;j<fidstack1->V(i).size();j++)
//             {
//                 util::point2d fid;
//                 fid.x=fidstack1->V(i)[j].x;
//                 fid.y=fidstack1->V(i)[j].y;
//                 fids.push_back(fid);
//             }
//         }
//         else {
//             for(int j=0;j<fidstack2->V(num).size();j++)
//             {
//                 util::point2d fid;
//                 fid.x=fidstack2->V(num)[j].x;
//                 fid.y=fidstack2->V(num)[j].y;
//                 fids.push_back(fid);
//             }
//             num++;
//         }
    }
//     std::cout<<"ok"<<std::endl;
//     std::vector<util::point2d>& fids = vfidp[1];
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

util::NodePlane::_node* util::NodePlane::AddNode(const util::_point& point)
{
	_node newnode(point);
	newnode.ref = information;
	bool found = false;
	
	std::list<_node>::iterator itr;
	
	if(pre_visited != node_list.begin() && *pre_visited > newnode){
		itr = node_list.begin();
	}
	else{
		itr = pre_visited;
	}
	
	for(; itr != node_list.end(); itr++){
		if(*itr == newnode){
			pre_visited = itr;
			found = true;
			break;
		}
		if(*itr > newnode){
			break;
		}
	}
	
	if(!found){
		pre_visited = node_list.insert(itr, newnode);
	}
	return &(*pre_visited);
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
