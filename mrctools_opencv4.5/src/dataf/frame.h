#ifndef FRAME_H
#define FRAME_H
#include "point2d.h"
#include <vector>
using namespace std;
class Point2D;
class frame
{
 public:
  frame();
  frame(int frameID_ );
//   frame(int frameID, int width, int height);
//   vector<Point2D*> p; // the markers found in this frame
  int frameID;
//   int width, height;
  bool discard;//if true indicates we should not use this frame
};

typedef struct{
    double A11;
    double A12;
    double A21;
    double A22;
    double x_shift;
    double y_shift;
}XfParam;

#endif