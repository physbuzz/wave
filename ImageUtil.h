#ifndef IMAGEUTIL_H
#define IMAGEUTIL_H

#include <vector>
#include <fstream>
#include <string>
#include <iostream>

//MISC UTILS
inline 
int floatToRGB(double r,double g,double b,double a=0){
	int a1=(int)(255*a);
	int r1=(int)(255*r);
	int g1=(int)(255*g);
	int b1=(int)(255*b);
	a1=a1>0?(a1>255?255:a1):0;
	r1=r1>0?(r1>255?255:r1):0;
	g1=g1>0?(g1>255?255:g1):0;
	b1=b1>0?(b1>255?255:b1):0;
	return (a1<<24)|(r1<<16)|(g1<<8)|b1;
}

inline 
int intToRGB(int r,int g,int b,int a=0){
	a=a>0?(a>255?255:a):0;
	r=r>0?(r>255?255:r):0;
	g=g>0?(g>255?255:g):0;
	b=b>0?(b>255?255:b):0;
	return (a<<24)|(r<<16)|(g<<8)|b;
}

//writes all the binary data given to a file.
inline
void writeCharFile(std::string name,std::vector<char> data){
	std::ofstream output(name.c_str(),std::ofstream::out|std::ofstream::binary);
	output.write(&data[0],sizeof(char)*data.size());
	output.close();
}

/*Example usage of Image:

Image x(1280,800)
x.put(100,100,floatToRGB(1,0,0));
x.save("myBitmap.bmp");

Saves a bitmap with a single red dot on it. For some weird reason 
the vertical axis is flipped, but that's not important right now. */
class Image{
	int width;
	int height;

	std::vector<char> data;
public:
	Image(int w,int h);
	Image(std::vector<double> arr,int w,int h);
	void createFromArray(std::vector<double> arr,int width,int height);
	void put(int x,int y,int c);
	int get(int x,int y);
	void save(std::string filename);
    void line(int x, int y, int x2, int y2,int c);
    int getWidth();
    int getHeight();
};

/*
It's convenient sometimes to just keep a 2D array of doubles
while treating it as an "image".
 */
class DoubleImage{
	int width;
	int height;

	std::vector<double> data;
public:
	int getWidth();
	int getHeight();
	DoubleImage();
	DoubleImage(int w,int h);
	DoubleImage(std::vector<double> arr,int w,int h);
	void createFromArray(std::vector<double> arr,int width,int height);
	const std::vector<double> getData();
	void put(int x,int y,double c);
	void increase(int x,int y,double c);
    void zeroImage();
	double get(int x,int y);
	double &at(int x,int y);
	double getMax();
	double getMin();
	void unitStretch(); //scale and translate so that the lowest value is 0 and the highest is 1.
};


inline double DoubleImage::getMax(){
	double max=data[0];
	for(int x=0;x<width;x++) 
		for(int y=0;y<height;y++) 
			if(data[y*width+x]>max)
				max=data[y*width+x];
	return max;
}
inline double DoubleImage::getMin(){
	double min=data[0];
	for(int x=0;x<width;x++) 
		for(int y=0;y<height;y++) 
			if(data[y*width+x]<min)
				min=data[y*width+x];
	return min;
}
inline void DoubleImage::unitStretch(){
	double min=getMin();
	double max=getMax();
	if(max-min>0)
		for(int x=0;x<width;x++) 
			for(int y=0;y<height;y++) 
				data[y*width+x]=(data[y*width+x]-min)/(max-min);
}

inline 
Image::Image(int w,int h): width(w),height(h),data(std::vector<char>(w*h*4,0)){
}

inline
Image::Image(std::vector<double> arr,int w,int h) {
	createFromArray(arr,w,h);
}

inline
void Image::createFromArray(std::vector<double> arr,int width,int height) {
	data=std::vector<char>(width*height*4,0);
	this->width=width;
	this->height=height;
	for(int x=0;x<width;x++) {
		for(int y=0;y<height;y++) {
			double a=arr[y*width+x];
			put(x,y,floatToRGB(a,a,a));
		}
	}
}

inline
int Image::getWidth(){return width;}

inline
int Image::getHeight(){return height;}

inline
void Image::put(int x,int y,int c){
	if(x>=0&&x<width && y>=0&&y<height)
		*((int *)&data[(x+y*width)*4])=c;

}

inline
void Image::line(int x, int y, int x2, int y2, int c){
    //Bresenham line algorithm.
    int w=x2-x;
    int h=y2-y;
    int dx1=0, dy1=0, dx2=0, dy2=0;

    if(w<0)
        dx1=-1;
    else if(w>0)
        dx1=1;

    if(h<0)
        dy1=-1;
    else if(h>0)
        dy1=1;

    if(w<0)
        dx2=-1;
    else if(w>0)
        dx2=1;

    int longest=std::abs(w);
    int shortest=std::abs(h);
    if (!(longest>shortest))
    {
        longest=std::abs(h);
        shortest=std::abs(w);
        if (h<0)
            dy2=-1;
        else if(h>0)
            dy2=1;
        dx2=0;
    }
    int numerator=longest>>1;
    for(int i=0;i<=longest;i++)
    {
        put(x,y,c);
        numerator+=shortest;
        if(!(numerator<longest))
        {
            numerator-=longest;
            x+=dx1;
            y+=dy1;
        }
        else
        {
            x+=dx2;
            y+=dy2;
        }
    }
}

inline
int Image::get(int x,int y){
	if(x>=0&&x<width && y>=0&&y<height)
		return *((int *)&data[(x+y*width)*4]);
	else
		return 0;
}

inline
void Image::save(std::string filename){
	int k=width*height;
	int s=4*k;
	int filesize=54+s;

	//disgusting code assuming little endian
	std::vector<char> data2(filesize);
	data2[0]='B';
	data2[1]='M';
	*((int *)&data2[2])=filesize;
	*((int *)&data2[10])=54;

	*((int *)&data2[14])=40; //header size
	*((int *)&data2[18])=width; //bitmap width in pixels
	*((int *)&data2[22])=height; //bitmap height in pixels
	*((int *)&data2[26])=1; //"Must be set to 1."
	*((int *)&data2[28])=32; //bits per pixel
	*((int *)&data2[30])=0; //compression method (0=none)
	*((int *)&data2[34])=s; //size of raw bitmap data in bytes
	*((int *)&data2[38])=100; //horizontal rez in pixels/meter
	*((int *)&data2[42])=100; //vertical rez in pixels/meter
	*((int *)&data2[46])=0; //"the number of colors in the color palette, or 0 to default to 2n."
	*((int *)&data2[50])=0; //"the number of important colors used, or 0 when every color is important; generally ignored."
	for(int x=0;x<data.size();x++){
		data2[x+54]=data[x];
	}

	writeCharFile(filename,data2);
}

inline
DoubleImage::DoubleImage() : width(0),height(0),data(){

}

inline
DoubleImage::DoubleImage(int w,int h) : width(w),height(h),data(w*h){

}

inline
DoubleImage::DoubleImage(std::vector<double> arr,int w,int h) {
	createFromArray(arr,w,h);
}

inline
void DoubleImage::createFromArray(std::vector<double> arr,int width,int height) {
	 if(arr.size()<width*height) {
		std::cout<<"DoubleImage given invalid array size"<<std::endl;
		return;
	}
	data=arr;
	this->width=width;
	this->height=height;
}

inline
const std::vector<double> DoubleImage::getData() {
	return data;
}

inline
void DoubleImage::increase(int x,int y,double c) {
	data[y*width+x]+=c;
}

inline 
void DoubleImage::zeroImage(){
    for(int x=0;x<width;x++) 
        for(int y=0;y<height;y++) 
            data[y*width+x]=0;
}

inline
void DoubleImage::put(int x,int y,double c) {
	data[y*width+x]=c;
}

inline
double DoubleImage::get(int x,int y) {
	return data[y*width+x];
}

inline
double &DoubleImage::at(int x,int y) {
	return data[y*width+x];
}

inline
int DoubleImage::getWidth() {
	return width;
}

inline
int DoubleImage::getHeight() {
	return height;
}
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

typedef struct rgb {
  float r, g, b;
} RGB;

typedef struct hsl {
  float h, s, l;
} HSL;

/*
 * Converts an RGB color value to HSL. Conversion formula
 * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
 * Assumes r, g, and b are contained in the set [0, 255] and
 * returns HSL in the set [0, 1].
 */
HSL rgb2hsl(float r, float g, float b) {

  HSL result;

  r /= 255;
  g /= 255;
  b /= 255;

  float max = MAX(MAX(r,g),b);
  float min = MIN(MIN(r,g),b);

  result.h = result.s = result.l = (max + min) / 2;

  if (max == min) {
    result.h = result.s = 0; // achromatic
  }
  else {
    float d = max - min;
    result.s = (result.l > 0.5) ? d / (2 - max - min) : d / (max + min);

    if (max == r) {
      result.h = (g - b) / d + (g < b ? 6 : 0);
    }
    else if (max == g) {
      result.h = (b - r) / d + 2;
    }
    else if (max == b) {
      result.h = (r - g) / d + 4;
    }

    result.h /= 6;
  }

  return result;

}

////////////////////////////////////////////////////////////////////////

/*
 * Converts an HUE to r, g or b.
 * returns float in the set [0, 1].
 */
float hue2rgb(float p, float q, float t) {

  if (t < 0)
    t += 1;
  if (t > 1)
    t -= 1;
  if (t < 1./6)
    return p + (q - p) * 6 * t;
  if (t < 1./2)
    return q;
  if (t < 2./3)
    return p + (q - p) * (2./3 - t) * 6;

  return p;

}

////////////////////////////////////////////////////////////////////////

/*
 * Converts an HSL color value to RGB. Conversion formula
 * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
 * Assumes h, s, and l are contained in the set [0, 1] and
 * returns RGB in the set [0, 255].
 */
RGB hsl2rgb(float h, float s, float l) {

  RGB result;

  if(0 == s) {
    result.r = result.g = result.b = l; // achromatic
  }
  else {
    float q = l < 0.5 ? l * (1 + s) : l + s - l * s;
    float p = 2 * l - q;
    result.r = hue2rgb(p, q, h + 1./3) * 255;
    result.g = hue2rgb(p, q, h) * 255;
    result.b = hue2rgb(p, q, h - 1./3) * 255;
  }

  return result;

}

#endif
