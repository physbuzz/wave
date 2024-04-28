#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

namespace Utils{

//Some string manipulation functions for saving files. pad_int(1234,5) returns "01234".
inline std::string pad_int(int arg, int padcount) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(padcount) << arg;
    return ss.str();
}

//Returns a file name in the form of "prefix00###suffix". For example "image0032.bmp"
inline std::string getFilename(std::string prefix, int num, int padcount, std::string suffix) {
    return prefix + pad_int(num, padcount) + suffix;
}

inline float randf(float min, float max) {
    return min + float(rand())*(max-min)/RAND_MAX;
}



class BresenhamLine {
    int x1,y1,x2,y2, dx,dy, step;
    float m,b;
    //int dx,dy,sx,sy,err,e2;
public:
    BresenhamLine(int x1,int y1,int x2,int y2) : x1(x1),y1(y1),x2(x2),y2(y2) {
        dx=x2-x1;
        dy=y2-y1;
        if(std::abs(dx)>std::abs(dy)){
            step=dx<0?-1:1; //calculate the direction from x1 to x2
            m=dy/float(dx);//calculate the slope as rise/run
            b=y1-m*x1; //calculate the y-intercept
        } else {
            if(dy==0){
                step=0;
                m=0;
                b=0;
            } else { //if dy=0 and abs(dx)<=abs(dy), we've drawn a dot and so need to do nothing.
                step=dy<0?-1:1; //calculate the direction from y1 to y2
                m=dx/float(dy); //calculate the slope as run/rise
                b=x1-m*y1; //calculate the x-intercept
            }
        }
    }
    class BLIterator {
        BresenhamLine *line;
        int x,y;
        bool end;
    public:
        BLIterator(BresenhamLine *line, bool end=false) : line(line), x(line->x1),y(line->y1), end(end)  { }

        std::pair<int,int> operator*() { return std::make_pair(x,y); }
        bool operator!=(const BLIterator &other) {
            return !(end && other.end); //|| (other.x!=x && other.y!=y), if we wanted to check for that.
        }
        BLIterator &operator++() {
            if(line->step==0){
                end=true;
            } else if(std::abs(line->dx)>std::abs(line->dy)){
                if(x==line->x2)
                    end=true;
                x+=line->step;
                y=std::round(line->m*x+line->b);
            } else {
                if(y==line->y2)
                    end=true;
                y+=line->step;
                x=std::round(line->m*y+line->b);
            }
            return *this;
        }
    };
    BLIterator begin() { return BLIterator(this);}
    BLIterator end() { return BLIterator(this,true);};
};



}
#endif
