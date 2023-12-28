#include "opts.h"
#include <iostream>
#include <fstream>
#include "dataf/dataf.h"
#include "modelmatch/match_core.h"
#include "detector/detector.h"
#include "nbundle/bundle_core.h"

// #define TESTDEVELOP

using namespace std;
// using namespace ann_1_1_char;


int main(int argc, char **argv)
{

    int tmp_dia=24;

    cout<<CV_VERSION<<endl;
    EX_TRACE("\nMARKERAUTO --version 1.6.3\n")

    struct options opts;
    opts.diameter = -1;
    opts.verbose = 0;
    opts.rotation_angle = 0;
	opts.testmode = false;
    
    if(GetOpts(argc, argv, &opts) <= 0) {
        EX_TRACE("***WRONG INPUT.\n");
        return -1;
    }
    
    vector<float> angles;
    vector<float> angles_heigh;
    vector<float> angles_new;
    vector<float> angles_po;
    vector<float> angles_neg;
    vector<int> m_index2;
    vector<int> m_index_high;
    vector<int> m_index_low;
    if(!util::ReadAnglesByName(opts.inputangle, &angles)) {
        std::cout<<"Can't open tilt angle file."<<endl;
        return -1;
    }
    
    int angele=45;
    for(int i=0;i<angles.size();i++)
    {
        m_index2.push_back(i);
        if(abs(angles[i])<angele) 
        {
            angles_new.push_back(angles[i]);
            m_index_low.push_back(i);
        }
        else{ 
            angles_heigh.push_back(angles[i]);
            m_index_high.push_back(i);
        }
        if(angles[i]>angele) angles_po.push_back(angles[i]);
        if(angles[i]<-angele) angles_neg.push_back(angles[i]);
        
    }
//     std::cout<<m_index_high.size()<<std::endl;
//     std::cout<<m_index_low.size()<<std::endl;
//     for(int i=0;i<angles_new.size();i++)
//     {
//         std::cout<<"new:"<<angles_new[i]<<std::endl;
//         std::cout<<"index:"<<m_index[i]<<std::endl;
//     }
    
    util::MrcStack mrcs;
    mrcs.Open(opts.input);

    util::FiducialStack fidstack;
//     fidstack.ReSize(angles_new.size());
//     fidstack.ReadFidsByFile("fids.txt");

#ifndef TESTDEVELOP
    EX_TIME_BEGIN("\n%sDo DetectorMain", _WAVE)
//     Detector::DetectorMain2(mrcs, &fidstack, opts.diameter, m_index_low, 1);
    Detector::DetectorMain(mrcs, &fidstack, opts.diameter, m_index_low, 1);
    fidstack.WriteFidsByFile("fids.txt");
//  if(1) {
//         Detector::Test(mrcs, fidstack, opts.diameter);
//     }
    if(opts.verbose >= 1) {
        Detector::Test(mrcs, fidstack, opts.diameter);
    }

    EX_TIME_END("Do DetectorMain")
#else
    fidstack.ReadFidsByFile("fids.txt");
#endif

    int sum=0;
    for(int i=0;i<fidstack.Size();i++)
    {
        std::vector<util::point2d>& fids = fidstack.V(i);
        sum=sum+fids.size();
    }
    std::cout<<sum<<std::endl;

//     std::cout<<m_index_low.size()<<std::endl;
    util::FiducialStack fidstack_low, fidstack_high;
    fidstack_low.SetWxH(fidstack.Width(), fidstack.Height());;
    fidstack_high.SetWxH(fidstack.Width(), fidstack.Height());;
    util::FiducialStack::div(m_index_high, m_index_low, &fidstack_low, &fidstack_high, &fidstack);
//     std::cout<<fidstack.Size()<<std::endl;
//     fidstack_low.WriteFidsByFile("fids_low.txt");
//     fidstack_high.WriteFidsByFile("fids_high.txt");
    
    util::ImgMatchVector imvector;
    HSetVector hset;
// #undef TESTDEVELOP
#ifndef TESTDEVELOP
    cv::Mat tmplt;
    util::SeriesReadFromFile(&tmplt, "avgtmplt");
    EX_TIME_BEGIN("\n%sDo MatchMain", _WAVE)
//     std::cout<<tmplt.size().width<<std::endl;
//     ModelMatch::MatchMain(fidstack_low, angles_new, &imvector, &hset, 0.85*tmp_dia, true, opts.testmode);	//0.85
    ModelMatch::MatchMain(fidstack_low, angles_new, &imvector, &hset, 0.85*tmplt.size().width, true, opts.testmode);	//0.85   11.34 23  19
    imvector.WriteVectorByFolder("matches");
    hset.WriteVectorByFolder("transmx");
//     ModelMatch::MatchMain(fidstack, angles, &imvector_all, &hset_all, 0.5*tmplt.size().width, true, opts.testmode);	//0.85    
//     imvector_all.WriteVectorByFolder("matches_all");
//     hset_all.WriteVectorByFolder("transmx_all");

    if(opts.verbose >= 1) {
        ModelMatch::Test(mrcs, imvector, 0.5f, "matches_ill");
    }

    EX_TIME_END("Do MatchMain")
#else
    imvector.ReadVectorByFolder("matches");
    hset.ReadVectorByFolder("transmx");
#endif
    
    util::TrackSpace trackspace;
    trackspace.Create(imvector, angles_new);

    util::FiducialStack addedfsk;
    util::ImgMatchVector addedimv;
    
#ifndef TESTDEVELOP
    Detector::LocalDetectorMain(mrcs, trackspace, hset, -1, &addedfsk, &addedimv, m_index_low);
    addedfsk.WriteFidsByFile("addfids.txt");
    addedimv.WriteVectorByFolder("addimvec");
#else
    addedfsk.ReadFidsByFile("addfids.txt");
    addedimv.ReadVectorByFolder("addimvec");
#endif
// 	Detector::Test(mrcs, addedfsk, opts.diameter, "addedfsk");
// 	ModelMatch::Test(mrcs, addedimv, 0.5f, "addedimv_ill");
    

    trackspace.InsertMatchVector(addedimv);


    trackspace.CoordinateTransform(fidstack.Width(), fidstack.Height());

    if(opts.rotation_angle < -0.01 || opts.rotation_angle > 0.01) {
        EX_TRACE("Do pre-rotation of series...\n")
        trackspace.PreRotate(-DEG2RAD(opts.rotation_angle));
    }
    
    std::vector<mx::pproj_params> cameras;
//     std::vector<mx::pproj_params> cameras_tmp;
    std::vector<mx::pproj_params> cameras_tmp;
    std::vector<v3_t> points;
    EX_TIME_BEGIN("\n%sDo BundleMain", _WAVE)
//     fidstack.Width()
    PPBundleApp::BundleMain(trackspace, fidstack.Width()/*1024*//*featsk.Width()*/, fidstack.Height() /*1024*//*featsk.Height()*/, &cameras, &cameras_tmp, &points, true);
    PPBundleApp::PrintCamerasAsIMOD(cameras, -DEG2RAD(opts.rotation_angle), 1, opts.outputxf, "xtiltangle.txt", opts.outputangle, "invalid.txt");
    EX_TIME_END("Do BundleMain")
    
    std::vector<mx::pproj_params> cameras_new;
//     std::vector<util::point2d> fid;
    util::FiducialStack addfids;
    util::FiducialStack proaddfids;
    addfids.SetWxH(mrcs.Width(), mrcs.Height());
    PPBundleApp::CreateMain(angles_heigh, &cameras_new, points, &addfids);
    
    
//     addfids.WriteFidsByFile("new_fids.txt");
    
//     int sum=0;
//     for(int i=0;i<1;i++)
//     {
//         std::vector<util::point2d>& fids = fidstack.V(i);
//         if(fids.size()<800)
//         {
//
//         }
//     }
    util::ImgMatchVector imvector2;
    HSetVector hset2;
    ModelMatch::MatchSameAngleMain(fidstack_high, addfids, angles_heigh, &imvector2, &hset2, 0.85*tmplt.size().width, true, opts.testmode);	//0.85
//     ModelMatch::MatchSameAngleMain(fidstack_high, addfids, angles_heigh, &imvector2, &hset2, 0.85*tmp_dia, true, opts.testmode);	//0.85
//     imvector2.WriteVectorByFolder("matches2");
//     hset2.WriteVectorByFolder("transmx2");
    util::FiducialStack pro_fids;
    util::ImgMatchVector imvector3;
    HSetVector hset3;

    util::FiducialStack allfids;
    util::FiducialStack fidstack_new;
    Detector::TransMain(mrcs, hset2, &fidstack_high, &addfids, &proaddfids, &pro_fids, &allfids, -1, m_index_high);
//     util::FiducialStack::comb(m_index_high, m_index_low, &fidstack_low, &fidstack_high, &fidstack_new);      //test
    util::FiducialStack::comb(m_index_high, m_index_low, &fidstack_low, &allfids, &fidstack_new);
    proaddfids.WriteFidsByFile("proadd.txt");
//     Detector::TransMain2(imvector3, hset2, &addfids, &pro_fids);
//     fidstack_new.WriteFidsByFile("fids_new");
    fidstack_new.SetWxH(mrcs.Width(), mrcs.Height());
    vector<float> angles_final;
    util::ImgMatchVector imvector_new;
    HSetVector hset_new;
    EX_TIME_BEGIN("\n%sDo MatchMain", _WAVE)
    ModelMatch::MatchMain(fidstack_new, angles, &imvector_new, &hset_new, 0.85*tmplt.size().width, true, opts.testmode);	//0.85
//     ModelMatch::MatchMain(fidstack_new, angles, &imvector_new, &hset_new, 0.85*tmp_dia, true, opts.testmode);	//0.85  23
    imvector_new.WriteVectorByFolder("matches_after");
    hset_new.WriteVectorByFolder("transmx_after");
    EX_TIME_END("Do MatchMain")
//  if(1) {
//         ModelMatch::Test(mrcs, imvector_new, 0.5f, "matches_ill");
//     }

    util::TrackSpace trackspace2;
    trackspace2.Create(imvector_new, angles);

    util::FiducialStack addedfsk2;
    util::ImgMatchVector addedimv2;
#ifndef TESTDEVELOP
    Detector::LocalDetectorMain(mrcs, trackspace2, hset_new, -1, &addedfsk2, &addedimv2, m_index2);
    addedfsk2.WriteFidsByFile("addfids2.txt");
    addedimv2.WriteVectorByFolder("addimvec2");
#else
    addedfsk.ReadFidsByFile("addfids.txt");
    addedimv.ReadVectorByFolder("addimvec");
#endif

    trackspace2.InsertMatchVector(addedimv2);



    trackspace2.CoordinateTransform(fidstack_new.Width(), fidstack_new.Height());

    if(opts.rotation_angle < -0.01 || opts.rotation_angle > 0.01) {
        EX_TRACE("Do pre-rotation of series...\n")
        trackspace2.PreRotate(-DEG2RAD(opts.rotation_angle));
    }

    std::vector<mx::pproj_params> cameras2;
    std::vector<v3_t> points2;
    ModelMatch::InitCamera(angles, hset_new, cameras, &cameras_tmp, m_index_low);

    EX_TIME_BEGIN("\n%sDo BundleMain", _WAVE)
    PPBundleApp::BundleMain(trackspace2, 0/*1024*//*featsk.Width()*/, 0/*1024*//*featsk.Height()*/, &cameras2, &cameras_tmp, &points2, false);
    PPBundleApp::PrintCamerasAsIMOD(cameras2, -DEG2RAD(opts.rotation_angle), 1, opts.outputxf, "xtiltangle2.txt", opts.outputangle, "invalid2.txt");
    EX_TIME_END("Do BundleMain")
    mrcs.Close();
}
