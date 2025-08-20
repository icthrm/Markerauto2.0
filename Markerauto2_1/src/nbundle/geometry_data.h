#ifndef GEOMETRY_H__
#define GEOMETRY_H__

#include <vector>
#include <cassert>
#include <algorithm>
#include <list>
#include <utility>
#include <string.h>
#include "micros.h"


namespace bundle{

struct Keypoint{
public:
    float m_x, m_y;        	///> Subpixel location of keypoint. 
    int m_extra;  		///> 4 bytes of extra storage */
    int m_track;  		///> Track index this point corresponds to 
    
public:
    Keypoint(){
        m_x = m_y = 0.0;
        m_extra = -1;
        m_track = -1;
    }

    Keypoint(float x, float y) : m_x(x), m_y(y), m_extra(-1), m_track(-1){}

    virtual ~Keypoint(){}
};

/** @brief Data struct for matches */
struct KeypointMatch{
public:
    int m_idx1, m_idx2;
public:
    KeypointMatch(){}

    KeypointMatch(int idx1, int idx2) :m_idx1(idx1), m_idx2(idx2){}
};

typedef std::pair<unsigned long, unsigned long> MatchIndex;

struct AdjListElem{
public:
    bool operator< (const AdjListElem &other) const{
        return m_index < other.m_index;
    }

    unsigned int m_index;
    std::vector<KeypointMatch> m_match_list;
};

typedef std::vector<AdjListElem> MatchAdjList;

class MatchTable
{
private:
    std::vector<MatchAdjList> m_match_lists;

public:
    MatchTable(){ }

    MatchTable(int num_images){
        m_match_lists.resize(num_images);
    }

    void SetMatch(MatchIndex idx){
        if (Contains(idx))
            return;  // already set

        /* Using vector */
        AdjListElem e;
        e.m_index = idx.second;
        MatchAdjList &l = m_match_lists[idx.first];
        MatchAdjList::iterator p = std::lower_bound(l.begin(), l.end(), e);
        l.insert(p, e);
    }

    void AddMatch(MatchIndex idx, KeypointMatch m){
        assert(Contains(idx));
        GetMatchList(idx).push_back(m);
    }

    void ClearMatch(MatchIndex idx){ // but don't erase!
        if (Contains(idx)){
            GetMatchList(idx).clear();
        }
    }

    void RemoveMatch(MatchIndex idx){
        if (Contains(idx)){
            std::vector<KeypointMatch> &match_list = GetMatchList(idx);
            match_list.clear();
            AdjListElem e;
            e.m_index = idx.second;
            MatchAdjList &l = m_match_lists[idx.first];
            std::pair<MatchAdjList::iterator, MatchAdjList::iterator> p = std::equal_range(l.begin(), l.end(), e);

            assert(p.first != p.second); // l.end());

            l.erase(p.first, p.second);
        }
    }

    unsigned int GetNumMatches(MatchIndex idx){
        if (!Contains(idx))
            return 0;
	return GetMatchList(idx).size();
    }

    std::vector<KeypointMatch> &GetMatchList(MatchIndex idx){
        AdjListElem e;
        e.m_index = idx.second;
        MatchAdjList &l = m_match_lists[idx.first];
        std::pair<MatchAdjList::iterator, MatchAdjList::iterator> p = std::equal_range(l.begin(), l.end(), e);

        assert(p.first != p.second); // l.end());

        return (p.first)->m_match_list;
    }
    
    const std::vector<KeypointMatch> &GetMatchList(MatchIndex idx) const{
        AdjListElem e;
        e.m_index = idx.second;
        const MatchAdjList &l = m_match_lists[idx.first];
        std::pair<MatchAdjList::const_iterator, MatchAdjList::const_iterator> p = std::equal_range(l.begin(), l.end(), e);

        assert(p.first != p.second); // l.end());

        return (p.first)->m_match_list;
    }

    bool Contains(MatchIndex idx) const{
        AdjListElem e;
        e.m_index = idx.second;
        const MatchAdjList &l = m_match_lists[idx.first];
        std::pair<MatchAdjList::const_iterator,
        MatchAdjList::const_iterator> p = std::equal_range(l.begin(), l.end(), e);

        return (p.first != p.second); // l.end());
    }

    void RemoveAll(){
        int num_lists = m_match_lists.size();

        for (int i = 0; i < num_lists; i++){
            m_match_lists[i].clear();
        }
    }

    unsigned int GetNumNeighbors(unsigned int i){
        return m_match_lists[i].size();
    }

    MatchAdjList &GetNeighbors(unsigned int i){
        return m_match_lists[i];
    }

    MatchAdjList::iterator Begin(unsigned int i){
        return m_match_lists[i].begin();
    }

    MatchAdjList::iterator End(unsigned int i){
        return m_match_lists[i].end();
    }
};

/* Return the match index of a pair of images */
inline MatchIndex GetMatchIndex(int i1, int i2){
    return MatchIndex((unsigned long) i1, (unsigned long) i2);
}

typedef std::pair<int,int> ImageKey;
typedef std::vector<ImageKey> ImageKeyVector;

/* Data for tracks */
struct TrackData{
    ImageKeyVector m_views;
    int m_extra;
    
public:
    TrackData() : m_extra(-1){}
    TrackData(ImageKeyVector views) : m_views(views), m_extra(-1){}
};

/* Data for 3D points */
struct PointData{
    double m_pos[3];  /* 3D position of the point */

    ImageKeyVector m_views;  /* View / keys corresponding to this point */
    bool m_fixed;      /* Should this point be fixed during bundle adjustment? */
    
public:
    PointData(){ m_fixed = false; }
};

class ImageData{
public:
    std::vector<Keypoint> m_keys;               ///> Keypoints in this image    
    std::vector<int> m_visible_points;  	///> Indices of points(track) visible in this system 
    std::vector<int> m_visible_keys;
    
    bool m_ignore_in_bundle;  			///> if ignore this image during bundle adjustment 
    float m_angle_in_mrc;
      
public:
    ImageData(){
	m_ignore_in_bundle = false;
    }
    
    /** @brief add a keypoint to image(may cause some question, axis system)*/
    void AddKey(float x, float y)				//should be careful
    {
	Keypoint kps(x, y);
	
	m_keys.push_back(kps);
    }
    
    /** @brief return the size of m_keys*/
    int GetNumKeys(){ return m_keys.size();}
    
    /** @brief set the extra ref in m_keys */
    void SetTracks(){
	int num_tracks = (int)m_visible_points.size();
	
	for (int i = 0; i < num_tracks; i++){
	    int tr = m_visible_points[i];
	    int key = m_visible_keys[i];
	    
	    assert(key < (int)m_keys.size());

	    m_keys[key].m_track = tr;
	}
    }
};

class GeometryData{
public:
    /* Geometry data */
    std::vector<TrackData> m_track_data;   		///> Information about the detected 3D tracks
    std::vector<ImageData> m_image_data;   		///> Image data     
    MatchTable m_matches;
    
public:
    int NumImages() const{ return (int)m_image_data.size();}
    
    /** @brief Get keys */
    Keypoint& GetKey(int img, int key){ return m_image_data[img].m_keys[key];}    
    const Keypoint& GetKey(int img, int key) const{ return m_image_data[img].m_keys[key];}  
    int GetNumKeys(int img) const{ return (int)m_image_data[img].m_keys.size();}
    
    /** @brief Operation wrap of MatchTable */  
    bool HasImagesMatch(int i1, int i2){ return m_matches.Contains(MatchIndex(i1, i2));}
    int GetNumMatches(int i1, int i2){ return m_matches.GetNumMatches(MatchIndex(i1, i2));}
    float OriAngles(int img_idx) const{ return m_image_data[img_idx].m_angle_in_mrc;}
    
};

}

#endif
