#ifndef MICRO_H__
#define MICRO_H__

/** @brief predefined micros*/
#define CLAMP(x,mn,mx) (((x) < mn) ? mn : (((x) > mx) ? mx : (x)))

#ifndef MAX
#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#define MAX4(a,b,c,d) (MAX(MAX(a,b), MAX(c,d)))
#define MIN4(a,b,c,d) (MIN(MIN(a,b), MIN(c,d)))

#define DEG2RAD(d) ((d) * (M_PI / 180.0))
#define RAD2DEG(r) ((r) * (180.0 / M_PI))

#define SQ(x) ((x) * (x))
#define CB(x) ((x) * (x) * (x))

/* Return the index of entry (x, y) in a 2D array */
#define INDEX(x, y, w) ((y) * (w) + (x))

#define SGN(x) ((x) < 0 ? (-1) : (1))

#endif