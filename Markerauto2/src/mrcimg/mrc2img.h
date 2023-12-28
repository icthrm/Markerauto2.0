#ifndef MRC2IMG_H__
#define MRC2IMG_H__

#include "mrcheader.h"
#include <iostream>
#include <string>
#include <fstream>
// #include "opencv2/core.hpp"
#include <opencv2/core/core.hpp>
// #include <cxcore.hpp>
// #include <cv.hpp>


namespace util{
    
/** @brief MrcStack; only for MRC of float format. */
class MrcStack{
private:
    MRCheader* header;
    std::string name;
    std::ifstream in;
    cv::Mat* slices;
	bool is_cashed;

private:
	cv::Mat ReadStackImage(size_t index);
	
public:
    MrcStack():header(nullptr), slices(nullptr), is_cashed(false){}
    
    ~MrcStack(){if(header){delete header;}}
    
    bool Open(const char* filename);
    
    void Close(){if(header){in.close(); delete header; header = nullptr;} FreeCache();}
    
    cv::Mat GetStackImage(int index);
	
	void DoCaching();
// 	
	void FreeCache();
// 	
	cv::Mat const& operator [](int idx);
	
	cv::Mat  const& operator [] (int idx) const;
    
    int Size() const{return header->nz;}
    
    int Width() const{return header->nx;}
    
    int Height() const{return header->ny;}
    
    const char* Name() const{return name.c_str();}
    
    void PrintHeader(std::ostream& o = std::cout) const;
};

}

#endif
