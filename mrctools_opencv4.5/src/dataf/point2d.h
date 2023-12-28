#ifndef POINT2D_H
#define POINT2D_H

#include <iostream>
#include <vector>

using namespace std;

class Point2D
{
public:
    float x, y;
    int markerType;//needed to identify to which marker size belongs to when there are different markers sizes in the same sample
    int frameID;
    Point2D()
    {
    }
    Point2D(float a, float b, int ID)
    {
        x = a;
        y = b;
        frameID= ID;
    }
    Point2D(float a, float b, int ID, int markerType_)
    {
        x = a;
        y = b;
        markerType=markerType_;
        frameID= ID;
    }
    //famatnote: whenever you have vectors or pointers you need a copy constructor
    Point2D & operator=(const Point2D &p)
    {
        if (this!=&p)
        {
            x=p.x;
            y=p.y;
            markerType=p.markerType;
            frameID=p.frameID;
        }
        return *this;

    }
    Point2D(const Point2D &p)//copy constructor
    {
        //cout<<"I am in point2D constructor copy"<<endl;

        x=p.x;
        y=p.y;
        markerType=p.markerType;
        frameID=p.frameID;
    }
    //famatnote: whenever we have pointers or vectors we need a destructor
    ~Point2D()
    {
        //cout<<"I am in Point2D destructor and index0 size is="<<index0.size()<<endl;

    }

    friend istream &operator>>(istream &is, Point2D &p);
    friend ostream &operator<<(ostream &os,
                               const Point2D &p);
};

#endif