#include "mrc2img.h"


bool util::MrcStack::Open(const char* filename) {
    in.open(filename);
    if(!in.good()) {
        return false;
    }
    
    if(header){
		delete header;
    }
    
    header = new MRCheader();

    in.read((char*) header, sizeof(MrcHeader));
//     if ( ! ( header->cmap[0]=='M'&&header->cmap[1]=='A'&&header->cmap[2]=='P' ) ) {
//         in.close();
//         delete header;
//         header = NULL;
//         return false;
//     }

    std::string key , value;
    inplane_rotation = 0;
    for(int i = 0 ; i< MRC_NLABELS ; i++){
        header->labels[i][MRC_LABEL_SIZE] = 0;
		
        if(i < header->nlabl) {
            FixTitlePadding(header->labels[i]);
            if(AnalysLabel(header->labels[i] , key , value))
            {
                if(key == "Tilt axis rotation angle") {
                    std::istringstream ss(value);
                    float rot;
                    ss>>rot;
                    inplane_rotation = rot;
                }
            }
        }
    }

    name = std::string(filename);
    return true;
}

cv::Mat util::MrcStack::ReadStackImage(size_t index)
{
    assert(header != NULL && index >= 0 && index < header->nz);
    assert(in.good());

    const int& nx = header->nx;			//the coordinate is different(have not been tested)
    const int& ny = header->ny;

    cv::Size size(nx, ny);
    cv::Mat img = cv::Mat::zeros(size, CV_32FC1);

    int bufsize = nx*ny;

    switch(header->mode) {
    case MRC_MODE_BYTE: {
        unsigned char tmpbuf[bufsize];
        for(size_t y = 0; y < ny; y++) {
            in.seekg(sizeof(MrcHeader) +header->next+ (index*nx*ny+y*nx) *sizeof(unsigned char), std::ios::beg);
            in.read((char*) & (tmpbuf[y*nx]), nx*sizeof(unsigned char));
        }

        unsigned char* src = tmpbuf;
        for(size_t y = 0; y < img.size().height; y++) {
            float* start = (float*)(img.data+y*img.step);
            for(int x = 0; x < img.size().width; x++) {
                *start++ = *src++/255.0f;
            }
        }

        break;
    }
    case MRC_MODE_SHORT: {
        short* tmpbuf = new short[bufsize];
        for(size_t y = 0; y < ny; y++) {
            in.seekg(sizeof(MrcHeader) +header->next+ (index*nx*ny+y*nx) *sizeof(short), std::ios::beg);
            in.read((char*) & (tmpbuf[y*nx]),nx*sizeof(short));
        }

        short* src = tmpbuf;
        for(size_t y = 0; y < img.size().height; y++) {
            float* start = (float*)(img.data+y*img.step);
            for(int x = 0; x < img.size().width; x++) {
                *start++ = *src++;
            }
        }

        delete [] tmpbuf;
        break;
    }
//     case MRC_MODE_SHORT: {
//         short tmpbuf[bufsize];
//         for(size_t y = 0; y < ny; y++) {
//             in.seekg(sizeof(MrcHeader) +header->next+ (index*nx*ny+y*nx) *sizeof(short), std::ios::beg);
//             in.read((char*) & (tmpbuf[y*nx]),nx*sizeof(short));
//         }
//
//         short* src = tmpbuf;
//         for(size_t y = 0; y < img.size().height; y++) {
//             float* start = (float*)(img.data+y*img.step);
//             for(int x = 0; x < img.size().width; x++) {
//                 *start++ = *src++;
//             }
//         }
//
//         break;
//     }
    case MRC_MODE_FLOAT: {
        for(size_t y = 0; y < ny; y++) {
            in.seekg(sizeof(MrcHeader)+header->next+(index*nx*ny+y*nx)*sizeof(float), std::ios::beg);
            in.read((char*)(img.data+y*img.step),nx*sizeof(float));
        }

        break;
    }
    default:
        break;
//         delete img;
//         return NULL;
    }

    return img;
}

void util::MrcStack::DoCaching()
{
    if(!slices) {
        slices = new cv::Mat[header->nz];
    }

    for(int i = 0; i < header->nz; i++) {
        slices[i] = ReadStackImage(i);
    }

    is_cashed = true;
}

void util::MrcStack::FreeCache()
{
    if(!slices) {
        return;
    }

    delete [] slices;
    slices = NULL;
    is_cashed = false;
}

cv::Mat util::MrcStack::GetStackImage(int index) {
    if(!is_cashed) {
        return ReadStackImage(index);
    }
    else{
        return slices[index].clone();
    }
}

cv::Mat util::MrcStack::operator[](int idx)
{
    return slices[idx];
}

cv::Mat const& util::MrcStack::operator[](int idx) const
{
    return slices[idx];
}

void util::MrcStack::CopyToNewStack(util::MrcStack& nmrcs) const
{
    nmrcs.FreeCache();

    if(header) {
        if(!nmrcs.header) {
            nmrcs.header = new MrcHeader();
        }

        memcpy(nmrcs.header, header, sizeof(MrcHeader));
    }

    nmrcs.name = name;
    nmrcs.in.open(nmrcs.name.c_str());
    nmrcs.is_cashed = is_cashed;

    if(is_cashed) {
        if(nmrcs.slices){
            delete [] nmrcs.slices;
        }
        nmrcs.slices = new cv::Mat[header->nz];
        for(int i = 0; i < header->nz; i++) {
            slices[i].copyTo(nmrcs.slices[i]);
        }
    }
}

const util::MrcStack& util::MrcStack::operator=(const util::MrcStack& _mrcs)
{
    _mrcs.CopyToNewStack(*this);
    return *this;
}

void util::MrcStack::SetHeader(util::MrcStack::Mode mode, float amin, float amean, float amax)
{
    header->amin = amin;
    header->amean = amean;
    header->amax = amax;
    header->mode = mode;
}

void util::MrcStack::SetName(const char* __name)
{
    name = std::string(__name);
}

void util::MrcStack::WriteHeaderToFile(const char* __name)
{
    if(__name) {
        name = std::string(__name);
    }
    
    out.open(name.c_str(), std::ios::binary);
    out.write((char*)header, sizeof(MRCheader));

    char nullbuf[header->next];

    out.write(nullbuf, header->next);
}

void util::MrcStack::AppendStackImageToFile(const cv::Mat* img)
{
    for(size_t y = 0; y < header->ny; y++) {
        out.write((char*)(img->ptr()+y*img->step), header->nx*sizeof(float));
    }
}

static inline bool IsSpace(char c)
{
    if (' ' == c || '\t' == c)
        return true;
    return false;
}

static inline void Trim(std::string & str)
{
    if (str.empty()) {
        return;
    }
    int i, start_pos, end_pos;
    for (i = 0; i < str.size(); ++i) {
        if (!IsSpace(str[i])) {
            break;
        }
    }
    if (i == str.size()) { // 全部是空白字符串
        str = "";
        return;
    }

    start_pos = i;

    for (i = str.size() - 1; i >= 0; --i) {
        if (!IsSpace(str[i])) {
            break;
        }
    }
    end_pos = i;

    str = str.substr(start_pos, end_pos - start_pos + 1);
}


bool util::MrcStack::AnalysLabel(const std::string& label, std::string& key, std::string& value)
{
    int pos;
    if((pos = label.find('=')) == -1)
        return false;
    key = label.substr( 0 , pos);
    value = label.substr(pos+1 ,label.size() - (pos + 1));
    Trim(key);
    if (key.empty()) {
        return false;
    }
    Trim(value);
    return true;
}

void util::MrcStack::FixTitlePadding(char *label)
{
    int len;
    label[MRC_LABEL_SIZE] = 0x00;
    len = strlen(label);
    if (len < MRC_LABEL_SIZE)
        memset(&label[len], ' ', MRC_LABEL_SIZE - len);
}


void util::MrcStack::PrintHeader(std::ostream& o) const {
    o<<"\n nx: "<<header->nx;         /*  # of Columns                  */
    o<<"\n ny: "<<header->ny;         /*  # of Rows                     */
    o<<"\n nz: "<<header->nz;         /*  # of Sections.                */
    o<<"\n mode: "<<header->mode;       /*  given by #define MRC_MODE...  */

    o<<"\n nxstart: "<<header->nxstart;    /*  Starting point of sub image.  */
    o<<"\n nystart: "<<header->nystart;
    o<<"\n nzstart: "<<header->nzstart;

    o<<"\n mx: "<<header->mx;         /* Grid size in x, y, and z       */
    o<<"\n my: "<<header->my;
    o<<"\n mz: "<<header->mz;

    o<<"\n xlen: "<<header->xlen;       /* length of x element in um.     */
    o<<"\n ylen: "<<header->ylen;       /* get scale = xlen/nx ...        */
    o<<"\n zlen: "<<header->zlen;

    o<<"\n alpha: "<<header->alpha;      /* cell angles, ignore */
    o<<"\n beta: "<<header->beta;
    o<<"\n gamma: "<<header->gamma;

    o<<"\n mapc: "<<header->mapc;       /* map coloumn 1=x,2=y,3=z.       */
    o<<"\n mapr: "<<header->mapr;       /* map row     1=x,2=y,3=z.       */
    o<<"\n maps: "<<header->maps;       /* map section 1=x,2=y,3=z.       */

    o<<"\n amin: "<<header->amin;
    o<<"\n amax: "<<header->amax;
    o<<"\n amean: "<<header->amean;

    o<<"\n ispg: "<<header->ispg;       /* image type */
    o<<"\n nsymbt: "<<header->nsymbt;     /* space group number */


    /* 64 bytes */

    o<<"\n next: "<<header->next;
    o<<"\n sizeof header: "<<sizeof(MRCheader);
    o<<"\n creatid: "<<header->creatid;  /* Creator id, hvem = 1000, DeltaVision = -16224 */


    o<<"\n blank: "<<header->blank;

    o<<"\n nint: "<<header->nint;
    o<<"\n nreal: "<<header->nreal;
    o<<"\n sub: "<<header->sub;
    o<<"\n zfac: "<<header->zfac;

    o<<"\n min2: "<<header->min2;
    o<<"\n max2: "<<header->max2;
    o<<"\n min3: "<<header->min3;
    o<<"\n max3: "<<header->max3;
    o<<"\n min4: "<<header->min4;
    o<<"\n max4: "<<header->max4;


    o<<"\n idtype: "<<header->idtype;
    o<<"\n lens: "<<header->lens;
    o<<"\n nd1: "<<header->nd1;     /* Devide by 100 to get o<<header-> value. */
    o<<"\n nd2: "<<header->nd2;
    o<<"\n vd1: "<<header->vd1;
    o<<"\n vd2: "<<header->vd2;
    o<<"\n tiltangles: "<<header->tiltangles[0]<<" "<<header->tiltangles[1]<<" "
     <<header->tiltangles[2]<<" "<<header->tiltangles[3]<<" "<<header->tiltangles[4]<<" "
     <<header->tiltangles[5]<<" ";  /* 0,1,2 = original:  3,4,5 = current */


    o<<"\n xorg: "<<header->xorg;
    o<<"\n yorg: "<<header->yorg;
    o<<"\n zorg: "<<header->zorg;
    o<<"\n cmap: "<<header->cmap[0]<<header->cmap[1]<<header->cmap[2]<<header->cmap[3];
    o<<"\n stamp: "<<header->stamp[0]<<header->stamp[1]<<header->stamp[2]<<header->stamp[3];
    o<<"\n rms: "<<header->rms;

    for(int i = 0; i < header->nlabl; i++) {
        o<<"\n labels["<<i<<"]: "<<header->labels[i];
    }
    o.flush();
}


