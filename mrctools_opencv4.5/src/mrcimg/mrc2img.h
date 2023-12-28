#ifndef MRC2IMG_H__
#define MRC2IMG_H__

#include "mrcheader.h"
#include <iostream>
#include <string>
#include <fstream>
#include "opencv2/core.hpp"

namespace util{
    
/** @brief MrcStack; only for MRC of float format. */
class MrcStack{
private:
    MRCheader* header;
    std::string name;
    std::ifstream in;
	std::ofstream out;
	cv::Mat* slices;
	bool is_cashed;

private:
	cv::Mat ReadStackImage(size_t index);
	
private:
	void FixTitlePadding(char *label);
	
    bool AnalysLabel(const std::string& label , std::string& key , std::string&  value );
	
public:
	enum Mode{MODE_BYTE = 0, MODE_SHORT = 1, MODE_FLOAT = 2};
	
public:
    float inplane_rotation;//////////
	
public:
    MrcStack():header(nullptr), slices(nullptr), is_cashed(false){}
    
    ~MrcStack(){if(header){delete header;}}
    
    bool Open(const char* filename);
    
    void Close(){FreeCache(); if(header){in.close(); out.close(); delete header; header = nullptr;}}
    
    cv::Mat GetStackImage(int index);
	
	void DoCaching();
	
	void FreeCache();
	
	cv::Mat operator [](int idx);
	
	cv::Mat const& operator [] (int idx) const;
    
    int Size() const{return header->nz;}
    
    int Z() const{return header->nz;}
    
    void SetZ(int nz){header->nz = nz;}
    
    void SetSize(int nz){header->nz = nz;}
    
    int Width() const{return header->nx;}
    
    int X() const{return header->nx;}
    
    void SetX(int nx){header->nx = nx;}
    
    void SetWidth(int nx){header->nx = nx;}
    
    int Height() const{return header->ny;}
    
    int Y() const{return header->ny;}
    
    void SetY(int ny){header->ny = ny;}
    
    void SetHeight(int ny){header->ny = ny;}
    
    void SetVolumeSize(int nx, int ny, int nz){header->nx = nx; header->ny = ny; header->nz = nz;}
    
//     void FillWithValue(float pxval=0);
    
    const MRCheader& Header() const{return *header;}
    
    MRCheader& AllocHeader() {if(!header) header = new MrcHeader(); return *header;}
    
    MRCheader& GetHeaderReference(){return *header;}
    
    const char* Name() const{return name.c_str();}
    
    void PrintHeader(std::ostream& o = std::cout) const;
	
	void CopyToNewStack(MrcStack& nmrcs) const;
	
	const MrcStack& operator = (const MrcStack& _mrcs);
	
	void SetHeader(Mode mode, float amin, float amean, float amax);
	
	void SetName(const char* __name);
	
	void WriteHeaderToFile(const char* __name = NULL);
	
	void AppendStackImageToFile(const cv::Mat* img);
};

}

#endif
