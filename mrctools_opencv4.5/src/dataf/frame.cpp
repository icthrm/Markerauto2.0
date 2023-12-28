#include "frame.h"

// frame::frame(int new_frameID, int new_width, int new_height)
// {
//   frameID = new_frameID;
//   width = new_width;
//   height = new_height;
//   discard = false;
// }

frame::frame()
{
    discard = false;
}

frame::frame(int frameID_)
{
    frameID = frameID_;
    discard = false;
}
