#include <vector>
#include <cassert>
#include <stdexcept>
#include <climits>
#include "utils.h"
#include "ImageUtil.h"
#include "math.h"

#define BOUNDARY_EPS 1.0e-2

template<typename Float>
class Wave {
public:
    int n,m;
    std::vector<Float> phicur, philast, cellMass, fieldMass2, txx, tyy, txy, damping, dirich;

    Wave(int n, int m) : n(n),m(m),
    phicur(n*m,Float(0)), philast(n*m,Float(0)), cellMass(n*m,Float(1)), fieldMass2(n*m,Float(0)),
    txx((n-1)*m,Float(1)), tyy(n*(m-1),Float(1)), txy((n-1)*(m-1),Float(0)),
    damping(n*m,Float(0)), dirich(n*m,Float(1)) {
        //assert(n>1 && m>1);
        //assert(size_t(n)*size_t(m)<size_t(INT_MAX));
        if(n<=1||m<=1){
            throw std::invalid_argument("Grid dimensions must be greater than 1.");
        }
        if (size_t(n) * size_t(m) >= size_t(INT_MAX)) {
            throw std::overflow_error("Grid size is too large and may cause overflow.");
        }
    }

    void timestep(Float dt){
        std::vector<Float> phinext(n*m,Float(0));
        bool isDirichlet=true;

        for(int j=0;j<m;j++){
            for(int i=0;i<n;i++){
                int index=j*n+i;
                Float force=0;
                // If we are inside a boundary, do nothing.
                if(dirich[index]<=BOUNDARY_EPS){
                    phinext[index]=0;
                    continue;
                }
                //First do it poorly!
                /* 4pt stencil for tension forces */
                //i+1
                if(i==n-1){
                    //exception 1, we're on the boundary of the simulation domain.
                    if(isDirichlet){
                        //force+=-txx[index-1]*phicur[index];
                        force+=-txx[j*(n-1)+i-1]*phicur[index];
                    }
                } else if(dirich[index+1]<=0){
                    //Exception 2, we have a Dirichlet boundary inside the simulation domain.
                    force+=txx[j*(n-1)+i]*phicur[index]*(dirich[index+1]/dirich[index]-Float(1));
                } else {
                    //Generic case: simple force from tension.
                    force+=txx[j*(n-1)+i]*(phicur[index+1]-phicur[index]);
                }

                //i-1
                if(i==0){
                    if(isDirichlet){
                        force+=-txx[j*(n-1)+i]*phicur[index];
                    }
                } else if(dirich[index-1]<=0){
                    force+=txx[j*(n-1)+i]*phicur[index]*(dirich[index-1]/dirich[index]-Float(1));
                } else {
                    force+=txx[j*(n-1)+i]*(phicur[index-1]-phicur[index]);
                }

                //j+1
                if(j==m-1){
                    if(isDirichlet){
                        force+=-tyy[index-n]*phicur[index];
                        //Index sanity check: Maximum(index-n)=Maximum(j*n+i-n)=(m-1)*n+(n-1)-n=(m-1)*n-1,
                        //which is the size we initialized tyy to.  
                    }
                } else if(dirich[index+n]<=0){
                    force+=tyy[index]*phicur[index]*(dirich[index+n]/dirich[index]-Float(1));
                    //Sanity check: Maximum(index+n)=Maximum(j*n+i+n)=(m-2)*n+(n-1)+n=m*n-1,
                    //which is the size of phicur.
                } else {
                    force+=tyy[index]*(phicur[index+n]-phicur[index]);
                }

                //j-1
                if(j==0){
                    if(isDirichlet){
                        force+=-tyy[index]*phicur[index];
                    }
                } else if(dirich[index-n]<=0){
                    force+=tyy[index]*phicur[index]*(dirich[index-n]/dirich[index]-Float(1));
                } else {
                    force+=tyy[index]*(phicur[index-n]-phicur[index]);
                }
                
                /* 4pt stencil for the tension diagonals
                 * */

                //i+1, j+1
                
                if(i==n-1 || j==m-1){
                    //exception 1, we're on the boundary of the simulation domain.
                    //Let's ignore this case.
                } else if(dirich[index+n+1]<=0){
                    //Exception 2, we have a Dirichlet boundary inside the simulation domain.
                    force+=0.5*txy[j*(n-1)+i]*phicur[index]*(dirich[index+n+1]/dirich[index]-Float(1));
                } else {
                    //Generic case: cross-force from tension.
                    force+=0.5*txy[j*(n-1)+i]*(phicur[index+n+1]-phicur[index]);
                }

                //i-1,j-1
                
                if(i==0 || j==0){
                    //ignore this case
                } else if(dirich[index-n-1]<=0){
                    force+=0.5*txy[(j-1)*(n-1)+i-1]*phicur[index]*(dirich[index-n-1]/dirich[index]-Float(1));
                } else {
                    force+=0.5*txy[(j-1)*(n-1)+i-1]*(phicur[index-n-1]-phicur[index]);
                }

                //i+1,j-1
                if(i==n-1 || j==0){
                    //ignore this case
                } else if(dirich[index-n+1]<=0){
                    force+=-0.5*txy[(j-1)*(n-1)+i]*phicur[index]*(dirich[index-n+1]/dirich[index]-Float(1));
                } else {
                    //j*n+i-(n-1)
                    //min: 
                    force+=-0.5*txy.at((j-1)*(n-1)+i)*(phicur[index-n+1]-phicur[index]);
                }

                //i-1,j+1
                if(i==0 || j==m-1){
                    //ignore this case
                } else if(dirich[index+n-1]<=0){
                    force+=-0.5*txy[j*(n-1)+i-1]*phicur[index]*(dirich[index+n-1]/dirich[index]-Float(1));
                } else {
                    force+=-0.5*txy[j*(n-1)+i-1]*(phicur[index+n-1]-phicur[index]);
                }


                //Old way of adding friction and timestepping:
                //force+=-damping[index]*(phicur[index]-philast[index])/dt - phicur[index]*fieldMass2[index];
                //phinext[index]=2*phicur[index]-philast[index]+dt*dt*force/cellMass[index];

                //New way:
                phinext[index]=(2*phicur[index]-philast[index]+dt*dt*force/cellMass[index]-phicur[index]*fieldMass2[index]
                    +Float(0.5)*damping[index]*dt*philast[index])/(1+Float(0.5)*damping[index]*dt);
            }
        }
        philast=phicur;
        phicur=phinext;
    }
    void saveImage(std::string fname,int cropping=0){
        assert(n-2*cropping>=0 && m-2*cropping>=0);
        Image outimg(n-2*cropping,m-2*cropping);
        for(int j=cropping;j<m-cropping;j++){
            for(int i=cropping;i<n-cropping;i++){
                int index=j*n+i;
                /*
                Float amp=phicur[index];
                int red=amp>0?int(amp*255):0;
                int blue=amp<0?int(-amp*255):0;
                outimg.put(i-cropping,j-cropping,intToRGB(red,0,blue));*/
                Float amp=phicur[index];
                int c=std::abs(amp)*255;
                outimg.put(i-cropping,j-cropping,intToRGB(c,c,c));
            }
        }
        outimg.save(fname);
    }

    void setPhiSoft(int i, int j, Float to, Float weight){
        if(weight>1)
            phicur.at(j*n+i)=to;
        else if(weight>0)
            phicur.at(j*n+i)=phicur.at(j*n+i)*(Float(1)-weight)+weight*to;
    }
};


float waveInitialConditions(float x, float y, float s, float w, float fade, float t){
    float a=x/(s*s*w);
    float b=y/s;
    return std::exp(-0.5*b*b/(1+a*a))*std::cos(0.5*a*b*b/(1+a*a)-t*w+x*w)*(t<fade?(t/fade):1); 
}

float waveWeight(float x, float y, float s, float w, float fade, float t){
    float a=x/(s*s*w);
    float b=y/s;
    return std::exp(-0.5*b*b/(1+a*a)); 
}


float g00func(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    return (r>rs+cutoff)?(1.0f/(1-rs/r)):(1.0f/(1-rs/(rs+cutoff)));
}
float Gxxfunc(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    if(r<=rs+cutoff)
        r=rs+cutoff;
    return 1.0f-rs*x*x/(r*r*r);
}
float Gyyfunc(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    if(r<=rs+cutoff)
        r=rs+cutoff;
    return 1.0f-rs*y*y/(r*r*r);
}
float Gxyfunc(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    if(r<=rs+cutoff)
        r=rs+cutoff;
    return -rs*x*y/(r*r*r);
}

int main(){

    int imgx=1280, imgy=800;
    int dampingSize=90;
    float LSim=18;
    float dx=LSim/imgx;
    float maxDamping=5.0;




    int nx=imgx+2*dampingSize,ny=imgy+2*dampingSize;
    float L=dx*nx;

    Wave<float> wave(nx,ny);

    float myCellMass=1.0f*(dx*dx);
    float rs=2.0f;
    float cutoff=0.01;
    float dampingcutoff=0.25;
    for(int j=0;j<ny;j++){
        for(int i=0;i<nx;i++){


            float x=(float(i)-nx/2.0f)/nx*L, y=(float(j)-ny/2.0f)/nx*L;
            float r=std::sqrt(x*x+y*y);

            float g00=g00func(x,y,rs,cutoff);
            float Gxx=Gxxfunc(x+0.5*dx,y,rs,cutoff);
            float Gxy=Gxyfunc(x+0.5*dx,y+0.5*dx,rs,cutoff);
            float Gyy=Gyyfunc(x,y+0.5*dx,rs,cutoff);

            wave.cellMass.at(nx*j+i)=g00*dx*dx; 
            if(i<nx-1)
                wave.txx.at((nx-1)*j+i)=Gxx;
            if(j<ny-1)
                wave.tyy.at(nx*j+i)=Gyy;
            if(i<nx-1&&j<ny-1)
                wave.txy.at((nx-1)*j+i)=Gxy;
            if(r<rs+dampingcutoff){
                float damping1=maxDamping*3;
                float damping2=maxDamping*20;
                wave.damping[nx*j+i]=
                    (r>rs)? (1.0f-(r-rs)/(dampingcutoff))*damping1:
                    (1-r/rs)*(damping2-damping1)+damping1;
            }

            /*
            float reps=std::sqrt(x*x+y*y+0.2*0.2);
            wave.cellMass[nx*j+i]=myCellMass*(1.0f+2.0/reps);
            wave.txx[(nx-1)*j+i]=myCellMass*(1.0f+2.0/reps);*/

            //if(std::sqrt(x*x+y*y)<3.0f)
                //wave.cellMass[nx*j+i]*=3;
            //wave.dirich[nx*j+i]=std::sqrt(x*x+y*y)-0.707f;
            //float r=std::sqrt(x*x+y*y);
            //float r0=2.0f;
            //wave.dirich[nx*j+i]=r-r0/3;
            //wave.damping[nx*j+i]=r>r0?0:(r0-r)/r0*1.0;
        }
    }
    assert(dampingSize<nx && dampingSize<ny);
    for(int k=0;k<dampingSize;k++){
        for(int j=0;j<ny;j++){
            wave.damping[j*nx+k]+=(float(1)-float(k)/dampingSize)*maxDamping;
            wave.damping[j*nx+nx-k-1]+=(float(1)-float(k)/dampingSize)*maxDamping;
        }

        for(int i=0;i<nx;i++){
            wave.damping[k*nx+i]+=(float(1)-float(k)/dampingSize)*maxDamping;
            wave.damping[(ny-k-1)*nx+i]+=(float(1)-float(k)/dampingSize)*maxDamping;
        }
    }

    float dt=0.0065;
    float time=0;
    int nsteps=36000;
    for(int k=0;k<nsteps;k++){
        for(int j=0;j<ny;j++){
            int i=dampingSize;
            float x=(float(i)/nx-0.5)*L, y=(float(j)-ny/2.0f)/nx*L;
            float s=1.0;
            float w=12.0;

            wave.setPhiSoft(dampingSize,j,
                    waveInitialConditions(x,y-3.5,s,w,0.2,time),
                    waveWeight(x,y-3.5,s,w,0.2,time)*10.0f);
            //wave.phicur[j*nx+dampingSize]=waveInitialConditions(x,y,1,5,0.2,time);
        }
        wave.timestep(dt);
        time+=dt;
        if(k%4==0){
            wave.saveImage(getFilename("out",k/4,5,".bmp"),dampingSize);
        }
    }
    return 0;
}













