#include "mrcheader.h"
#include <fstream>
#include <climits>
#include <cassert>

bool MrcReadHead(MrcHeader* header, const char* filename)
{
    std::ifstream in;
    in.open(filename);
    if(!in.good()){
	return false;
    }
    
    in.read((char*)header, sizeof(MrcHeader));
    if(!(header->cmap[0]=='M'&&header->cmap[1]=='A'&&header->cmap[2]=='P')){
	in.close();
	return false;       
    }
    
    in.close();
    return true;
}

bool MrcWriteHead(const MrcHeader& header, const char* filename)
{
    std::ofstream out;
    out.open(filename, std::ios::out|std::ios::ate);
    out.seekp(std::ios::beg);
    out.write(reinterpret_cast<const char*>(&header),sizeof(MrcHeader));
    out.close();
}

void MrcInitHead(MrcHeader *head)
{
    head->nx=1;
    head->ny=1;
    head->nz=1;

    head->mode=MRC_MODE_BYTE;

    head->nxstart=0;
    head->nystart=0;
    head->nzstart=0;

    head->mx=1;
    head->my=1;
    head->mz=1;

    head->xlen=1;
    head->ylen=1;
    head->zlen=1;

    head->alpha=90;
    head->beta=90;
    head->gamma=90;

    head->mapc=1;
    head->mapr=2;
    head->maps=3;

    head->amin=0;
    head->amax=255;
    head->amean=128;

    head->ispg=1;
    head->nsymbt=0;

    head->next=0;

    head->creatid=1000;
    head->cmap[0]='M';
    head->cmap[1]='A';
    head->cmap[2]='P';
    head->stamp[0]='D';
}

/**slcN couts from 0 to N-1, so if you want to read the first slice slcN shoud be 0*/
bool MrcReadSlice(const char* filename, int slcN, char axis, float *slcdata)
{
    
    std::ifstream in(filename);
    if(!in.good()){
	return false;
    }
    MRCheader head;    
    in.read((char*)&head, sizeof(MrcHeader));
    
    int psize;
    
    switch(head.mode){
    case MRC_MODE_BYTE :
        psize=sizeof(unsigned char);
        break;
    case MRC_MODE_SHORT :
        psize=sizeof(short);
        break;
    case MRC_MODE_FLOAT :
        psize=sizeof(float);
        break;
    }
    
    short buf_short;
    unsigned char buf_byte;
    float buf_float;
    char* vbuf;

    switch(head.mode){
    case MRC_MODE_BYTE:
	vbuf = (char*)&buf_byte;
	break;

    case MRC_MODE_SHORT:
	vbuf = (char*)&buf_short;
	break;

    case MRC_MODE_FLOAT:
	vbuf = (char*)&buf_float;
	break;
    }
    
    switch(axis){
    case 'x':
    case 'X':
        in.seekg(sizeof(MrcHeader)+head.next+slcN*psize,std::ios::beg);
	
	for(int i = 0; i < head.ny*head.nz; i++){
	    in.read(vbuf,psize);
	    switch(head.mode){
	    case MRC_MODE_BYTE:
		slcdata[i] = *((unsigned char*)vbuf);
		break;

	    case MRC_MODE_SHORT:
		slcdata[i] = *((short*)vbuf);
		break;

	    case MRC_MODE_FLOAT:
		slcdata[i] = *((float*)vbuf);
		break;
	    }
	    in.seekg((head.nx-1)*psize,std::ios::cur);
	}
        break;

    case 'y':
    case 'Y':
	for(int k = 0; k < head.nz; k++){
	    in.seekg(sizeof(MrcHeader)+head.next+(k*head.nx*head.ny+head.nx*slcN)*psize,std::ios::beg);

	    for(int i = 0; i < head.nx; i++){
		in.read(vbuf,psize);
		switch(head.mode){
		case MRC_MODE_BYTE:
		    slcdata[k*head.nx+i] = *((unsigned char*)vbuf);
		    break;

		case MRC_MODE_SHORT:
		    slcdata[k*head.nx+i] = *((short*)vbuf);
		    break;

		case MRC_MODE_FLOAT:
		    slcdata[k*head.nx+i] = *((float*)vbuf);
		    break;
		}
	    }
	}
        break;

    case 'z':
    case 'Z':
	in.seekg(sizeof(MrcHeader)+head.next+slcN*head.nx*head.ny*psize,std::ios::beg);

	for(int i = 0; i < head.nx*head.ny; i++){
	    in.read(vbuf,psize);
	    switch(head.mode){
	    case MRC_MODE_BYTE:
		slcdata[i] = *((unsigned char*)vbuf);
		break;

	    case MRC_MODE_SHORT:
		slcdata[i] = *((short*)vbuf);
		break;

	    case MRC_MODE_FLOAT:
		slcdata[i] = *((float*)vbuf);
		break;
	    }
	}
        break;
    }
    return true;
}

bool MrcAppData(int modetype, int size, const void* data, const char* filename)
{
    int psize;
    
    switch(modetype){
    case MRC_MODE_BYTE :
        psize=sizeof(unsigned char);
        break;
    case MRC_MODE_SHORT :
        psize=sizeof(short);
        break;
    case MRC_MODE_FLOAT :
        psize=sizeof(float);
        break;
    }
	
    std::ofstream out;
    out.open(filename, std::ios::out|std::ios::app);
    out.write((char*)data, size*psize);
    out.close();
}

void MrcUpdateHead(MrcHeader* header, void* data)
{
    unsigned char* vchar;
    short* vshort;
    float* vfloat;
    
    float amin = INT_MAX, amax = INT_MIN, amean = 0; 
    
    switch(header->mode){
    case MRC_MODE_BYTE :
        vchar = (unsigned char*)data;
	for(int i = 0; i < header->nx*header->ny*header->nz; i++){
	    amean += vchar[i];
	    if(vchar[i] < amin){
		amin = vchar[i];
	    }
	    if(vchar[i] > amax){
		amax = vchar[i];
	    }
	}
        break;
    case MRC_MODE_SHORT :
        vshort = (short*)data;
	for(int i = 0; i < header->nx*header->ny*header->nz; i++){
	    amean += vshort[i];
	    if(vshort[i] < amin){
		amin = vshort[i];
	    }
	    if(vshort[i] > amax){
		amax = vshort[i];
	    }
	}
        break;
    case MRC_MODE_FLOAT :
        vfloat = (float*)data;
	for(int i = 0; i < header->nx*header->ny*header->nz; i++){
	    amean += vfloat[i];
	    if(vfloat[i] < amin){
		amin = vfloat[i];
	    }
	    if(vfloat[i] > amax){
		amax = vfloat[i];
	    }
	}	
        break;
    }
    header->amin=amin;
    header->amax=amax;
    header->amean=amean/header->nx*header->ny*header->nz;    
}

void MrcSave(const MrcHeader& header, const void* data, const char* filename)
{
    MrcWriteHead(header, filename);
    MrcAppData(header.mode, header.nx*header.ny*header.nz, data, filename);
}

void MrcConvertTo1(MrcHeader* header, void* data)
{
    assert(header->amax-header->amin != 0);
    for(int i = 0; i < header->nx*header->ny*header->nz; i++){
	switch(header->mode){
	case MRC_MODE_BYTE :
	    ((unsigned char*)data)[i] = (((unsigned char*)data)[i]-header->amin)/(header->amax-header->amin);
	    break;
	case MRC_MODE_SHORT :
	    ((short*)data)[i] = (((short*)data)[i]-header->amin)/(header->amax-header->amin);
	    break;
	case MRC_MODE_FLOAT :
	    ((float*)data)[i] = (((float*)data)[i]-header->amin)/(header->amax-header->amin);
	    break;
	}
    }
    header->amean = (header->amean-header->amin)/(header->amax-header->amin);
    header->amin = 0;
    header->amax = 1;
}

void MrcDataConTo1(const MrcHeader& header, void* data)
{
    assert(header.amax-header.amin != 0);
    for(int i = 0; i < header.nx*header.ny*header.nz; i++){
	switch(header.mode){
	case MRC_MODE_BYTE :
	    ((unsigned char*)data)[i] = (((unsigned char*)data)[i]-header.amin)/(header.amax-header.amin);
	    break;
	case MRC_MODE_SHORT :
	    ((short*)data)[i] = (((short*)data)[i]-header.amin)/(header.amax-header.amin);
	    break;
	case MRC_MODE_FLOAT :
	    ((float*)data)[i] = (((float*)data)[i]-header.amin)/(header.amax-header.amin);
	    break;
	}
    }
}