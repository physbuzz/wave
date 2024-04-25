#include <vector>
#include <cassert>
#include <stdexcept>
#include <climits>
#include "utils.h"
#include "ImageUtil.h"
#include "math.h"
#include "omp.h"

#define BOUNDARY_EPS 1.0e-2





namespace Wave {


namespace Math {

/* Returns the amplitude of a Gaussian beam pointing in the +x direction, 
with width s at its most focused (at the origin), frequency w, and wavenumber k.
It's up to the user to specify w and k correctly. Dispersion (fieldMass>0) is not accounted for.
The maximum amplitude at (x,y)=(0,0) is normalized to 1.
*/

inline float gaussianBeamAmplitude(float x, float y, float t, float s, float w, float k, float phase) {
    float a=x/(s*s*k);
    float b=y/s;
    return std::exp(-0.5*b*b/(1+a*a))*std::cos(0.5*a*b*b/(1+a*a)-t*w+x*k+phase); 

}


inline float gaussianBeamAmplitudeWithDispersion(float x, float y, float t, float m, float s, float w, float k, float phase) {
    std::cerr<<"gaussianBeamAmplitudeWithDispersion not implemented yet."<<std::endl;
    return 0;
}

/* Schwarzschild black hole metric (upper index) components.
Specifically, these are the components of $\sqrt{-g}g^{\mu\nu}$ in the + - - - convention.
So, GXX, GYY, and GXY are the *negatives* of Txx, Tyy, and Txy in my implementation.

GTX and GTY are zero.
 */

float SBHGTT(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    return (r>rs+cutoff)?(1.0f/(1-rs/r)):(1.0f/(1-rs/(rs+cutoff)));
}

float SBHGXX(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    if(r<=rs+cutoff)
        r=rs+cutoff;
    return 1.0f-rs*x*x/(r*r*r);
}

float SBHGYY(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    if(r<=rs+cutoff)
        r=rs+cutoff;
    return 1.0f-rs*y*y/(r*r*r);
}

float SBHGXY(float x, float y, float rs, float cutoff){
    float r=std::sqrt(x*x+y*y);
    if(r<=rs+cutoff)
        r=rs+cutoff;
    return -rs*x*y/(r*r*r);
}


/* Kerr black hole metric (upper index) components.
These are the components of $\sqrt{-g}g^{\mu\nu}$ in the + - - - convention.
GXX, GYY, and GXY are the *negatives* of txx, tyy, and txy in my implementation.
GTX and GTY are the +ttx and tty. I think this is correct because the wave equation is
0=G^{tt} \ddot{\phi}
+\partial_i(G^{ij}\partial_j \phi)
+2 G^{0i}\partial_i \dot{\phi}+(\partial_i G^{0i})\dot{\phi}
Ignore the divergence of G^{0i}, and ignore the Laplacian term. Then this looks like a transport equation
$0=\ddot{\phi}+(2 \vec{G})/G^{00} \dot{\phi} = D \dot{\phi}/Dt$.

This convention works for the relativistic wave equation case, but I'll have to see if the convention is commensurate
with more ordinary cases (waves on the surface of a river).

These functions can be optimized much better! They are output from a mathematica script.
 */

float KERRGTT(float x, float y, float rs, float a) {
    float r=std::sqrt(x*x+y*y);
    return ((std::pow(r,3) + std::pow(a,2)*(r + rs))/(r*(std::pow(a,2) + r*(r - rs))));
}

float KERRGTX(float x, float y, float rs, float a) {
    float r=std::sqrt(x*x+y*y);
    return (-((a*rs*y)/(std::pow(a,2)*r + std::pow(r,3) - std::pow(r,2)*rs)));
}

float KERRGTY(float x, float y, float rs, float a) {
    float r=std::sqrt(x*x+y*y);
    return ((a*rs*x)/(std::pow(a,2)*r + std::pow(r,3) - std::pow(r,2)*rs));
}

float KERRGXX(float x, float y, float rs, float a) {
    float r=std::sqrt(x*x+y*y);
    return (-1 + (-std::pow(a,2) + r*rs)/std::pow(r,2) + ((std::pow(a,4) + 2*std::pow(a,2)*r*(r - rs) + std::pow(r,2)*rs*(-r + rs))*std::pow(y,2))/(std::pow(r,4)*(std::pow(a,2) + r*(r - rs))));
}

float KERRGXY(float x, float y, float rs, float a) {
    float r=std::sqrt(x*x+y*y);
    return (-(((std::pow(a,4) + 2*std::pow(a,2)*r*(r - rs) + std::pow(r,2)*rs*(-r + rs))*x*y)/(std::pow(r,4)*(std::pow(a,2) + r*(r - rs)))));
}

float KERRGYY(float x, float y, float rs, float a) {
    float r=std::sqrt(x*x+y*y);
    return (-((std::pow(a,4)*std::pow(y,2) + r*(r - rs)*(std::pow(r,4) + (2*std::pow(a,2) - r*rs)*std::pow(y,2)))/(std::pow(r,4)*(std::pow(a,2) + r*(r - rs)))));
}

} // end namespace Wave::Math


enum SourceType{
    POINT,
    PLANE,
    GAUSSIAN_BEAM
};
enum SourceStyle{
    SOFT_SIMPLE,
    SOFT_ADVANCED,
    HARD
};

class WaveKernel;

struct KernelSource {
public:
    SourceType type;
    SourceStyle style;
    //KernelSource position
    float x,y;
    //KernelSource direction and wave number
    //float kx,ky;
    //KernelSource frequency (it is up to the user to specify omega^2/c^2=kx^2+ky^2)
    float omega;
    float phase;
    float amplitude;
    
    //slowly apply the source over the first few seconds of sim time.
    float fadeTime;
    //If SourceStyle is SOFT, then we allow the wave to be generated over some distance fadeDist
    //float fadeDist;
    void apply(WaveKernel* wave);
private:
    float fadeFunc(float time) const;
    //Algo: look up the cell and set phi to cos(omega*time)*(time<fadeTime?time/fadeTime:1)
    void applyPointHard(WaveKernel* wave);
    //Algo: Do the same as applyPointHard, but lerp between the old value and the new value.
    void applyPointSoftSimple(WaveKernel* wave);
    //Algo: Treat things as if we added a source term to the wave equation, eg. wave->phi[...]+=(calculated value)
    void applyPointSoftAdvanced(WaveKernel* wave);
    //Algo: Bresenham line along the source plane setting it equal to some cos(omega*t) value.
    void applyPlaneHard(WaveKernel* wave);
    //Algo: same but lerp a bit.
    void applyPlaneSoftSimple(WaveKernel* wave);
    //Algo: apply a source term over fadeDist / in a big rotated rectangle. This might be hard to do fast, so do it slow first and then check.
    void applyPlaneSoftAdvanced(WaveKernel* wave);
    //Algo: Bresenham line, set it equal to the Gaussian beam formula I found.
    void applyGaussianBeamHard(WaveKernel* wave);
    //Algo: Bresenham line, lerp between phi[...] and the hard value.
    void applyGaussianBeamSoftSimple(WaveKernel* wave);
    //Algo: apply a source term to create a Gaussian beam using formulas I derived.
    void applyGaussianBeamSoftAdvanced(WaveKernel* wave);
};
KernelSource createPointKernelSource(float x, float y, float omega, float amplitude=1, float phase=0, float fadeTime=0.1);



enum ColorScheme {
    RED_BLUE,
    GRAYSCALE
};





class WaveKernel {
public:
    int grid_n,grid_m;
    float dt;
    float time_elapsed;
    std::vector<float> phicur, philast, cellMass, fieldMass2, ttx, tty, txx, tyy, txy, damping, dirich;
    std::vector<KernelSource> sources;

    WaveKernel(int n=2, int m=2,float defaultMass=1) : grid_n(n),grid_m(m), dt(0.01), time_elapsed(0),
    phicur(n*m,0), philast(n*m,0), cellMass(n*m,defaultMass), fieldMass2(n*m,0),
    ttx(n*m,0),tty(n*m,0),
    txx(n*m,1), tyy(n*m,1), txy(n*m,0),
    damping(n*m,0), dirich(n*m,1), sources() {
        if(n<=1||m<=1){
            throw std::invalid_argument("Grid dimensions must be greater than 1.");
        }
        if(size_t(n)*size_t(m)>=size_t(INT_MAX)) {
            throw std::overflow_error("Grid size is too large and may cause overflow.");
        }
    }

    //Getters for all the array variables w/ bounds checking through vector.at
    //Can do a const version if we need it later.
    float &getPhi(int i, int j) { return phicur.at(j*grid_n+i); }
    //float getPhi(int i, int j) const { return phicur.at(j*grid_n+i); }
    float &getPhiLast(int i, int j) { return philast.at(j*grid_n+i); }
    float &getCellMass(int i, int j) { return cellMass.at(j*grid_n+i); }
    float &getFieldMass2(int i, int j) { return fieldMass2.at(j*grid_n+i); }
    float &getTxx(int i, int j) { return txx.at(j*grid_n+i);}
    float &getTyy(int i, int j) { return tyy.at(j*grid_n+i); }
    float &getTxy(int i, int j) { return txy.at(j*grid_n+i); }
    float &getDamping(int i, int j) { return damping.at(j*grid_n+i); }
    float &getDirich(int i, int j) { return dirich.at(j*grid_n+i); }
    float &getTtx(int i, int j) { return ttx.at(j*grid_n+i); }
    float &getTty(int i, int j) { return tty.at(j*grid_n+i); }

    void addPointSource(float i, float j, float omega, float amplitude=1, float phase=0, float fadeTime=0.1) {
        if(i<=0 || j<=0 || i>=grid_n || j>=grid_m) {
            std::cerr<<"Point source must be inside the simulation domain."<<std::endl;
            return;
        }
        sources.push_back(createPointKernelSource(i,j,omega,amplitude,phase,fadeTime));
    }

    void timestep(){
        std::vector<float> phinext(grid_n*grid_m,float(0));

        //Dirichlet boundary conditions enforced by always updating the boundary cells to 0.
        #if defined(_OPENMP)
        #pragma omp parallel for
        #endif
        for(int j=1;j<grid_m-1;j++){
            for(int i=1;i<grid_n-1;i++){
                int index=j*grid_n+i;
                float force=0;
                // If we are inside a boundary, do nothing.
                if(dirich[index]<=BOUNDARY_EPS){
                    phinext[index]=0;
                    continue;
                }

                //First do it poorly!
                /* 4pt stencil for tension forces */
                //i+1
                //Exceptional case: we have a Dirichlet boundary inside the simulation domain.
                //  we get the value for phicur[index+1] by requiring (phicur[index+1]-phicur[index])*t+phicur[index]=0 
                //  where (dirich[index+1]-dirich[index])*t+dirich[index]==0
                //Generic case: simple force from tension.
                if(dirich[index+1]<=0)
                    force+=txx[index]*phicur[index]*(dirich[index+1]/dirich[index]-float(1));
                else 
                    force+=txx[index]*(phicur[index+1]-phicur[index]);

                //i-1
                if(dirich[index-1]<=0)
                    force+=txx[index-1]*phicur[index]*(dirich[index-1]/dirich[index]-float(1));
                else 
                    force+=txx[index-1]*(phicur[index-1]-phicur[index]);

                //j+1 
                if(dirich[index+grid_n]<=0)
                    force+=tyy[index]*phicur[index]*(dirich[index+grid_n]/dirich[index]-float(1));
                else
                    force+=tyy[index]*(phicur[index+grid_n]-phicur[index]);

                //j-1
                if(dirich[index-grid_n]<=0)
                    force+=tyy[index-grid_n]*phicur[index]*(dirich[index-grid_n]/dirich[index]-float(1));
                else 
                    force+=tyy[index-grid_n]*(phicur[index-grid_n]-phicur[index]);
                
                //4pt stencil for the tension diagonals
                //i+1,j+1
                if(dirich[index+grid_n+1]<=0)
                    force+=0.5*txy[index]*phicur[index]*(dirich[index+grid_n+1]/dirich[index]-float(1));
                else 
                    force+=0.5*txy[index]*(phicur[index+grid_n+1]-phicur[index]);

                //i-1,j-1
                if(dirich[index-grid_n-1]<=0)
                    force+=0.5*txy[index-grid_n-1]*phicur[index]*(dirich[index-grid_n-1]/dirich[index]-float(1));
                else
                    force+=0.5*txy[index-grid_n-1]*(phicur[index-grid_n-1]-phicur[index]);

                //i+1,j-1
                if(dirich[index-grid_n+1]<=0)
                    force+=-0.5*txy[index-grid_n+1]*phicur[index]*(dirich[index-grid_n+1]/dirich[index]-float(1));
                else
                    force+=-0.5*txy[index-grid_n+1]*(phicur[index-grid_n+1]-phicur[index]);

                //i-1,j+1
                if(dirich[index+grid_n-1]<=0)
                    force+=-0.5*txy[index+grid_n-1]*phicur[index]*(dirich[index+grid_n-1]/dirich[index]-float(1));
                else
                    force+=-0.5*txy[index+grid_n-1]*(phicur[index+grid_n-1]-phicur[index]);
                
                //Note: this ttx/tty convection turn this integrator from truly CTCS to 
                //something not time reversible! Uh, I hope this isn't unconditionally 
                //unstable or anything. Centered difference and downwind schemes are famously unstable
                //for pure convection problems, that's like numerical methods 101.
                //ttx forces due to convection
                if(dirich[index+1]>0)
                    force+=ttx[index]*(phicur[index+1]-philast[index+1])/dt;                
                if(dirich[index-1]>0)
                    force+=-ttx[index-1]*(phicur[index-1]-philast[index-1])/dt;                
                //tty convection
                if(dirich[index+grid_n]>0)
                    force+=tty[index]*(phicur[index+grid_n]-philast[index+grid_n])/dt;                
                if(dirich[index-grid_n]>0)
                    force+=-tty[index-grid_n]*(phicur[index-grid_n]-philast[index-grid_n])/dt;                
                


                //Old way of adding friction and timestepping:
                //force+=-damping[index]*(phicur[index]-philast[index])/dt - phicur[index]*fieldMass2[index];
                //phinext[index]=2*phicur[index]-philast[index]+dt*dt*force/cellMass[index];

                //Add the mass term, if there is one.
                force-=phicur[index]*fieldMass2[index];
                //New way:
                //The equation is cellMass*(phinext+philast-2*phicur)/(dt*dt)==force-damping*cellMass*(phinext-philast)/(2*dt)
                phinext[index]=(2*phicur[index]-philast[index]+dt*dt*force/cellMass[index]
                    +float(0.5)*damping[index]*dt*philast[index])/(1+float(0.5)*damping[index]*dt);
            }
        }
        philast=phicur;
        phicur=phinext;

        time_elapsed+=dt;
        for(auto source : sources){
            source.apply(this);
        }
    }

    void setPhiSoft(int i, int j, float to, float weight){
        if(weight>1)
            phicur.at(j*grid_n+i)=to;
        else if(weight>0)
            phicur.at(j*grid_n+i)=phicur.at(j*grid_n+i)*(float(1)-weight)+weight*to;
    }
};


inline void KernelSource::apply(WaveKernel* wave){
    switch(type) {
        case POINT:
            switch(style){
                case SOFT_SIMPLE:
                    applyPointSoftSimple(wave);
                    break;
                case SOFT_ADVANCED:
                    applyPointSoftAdvanced(wave);
                    break;
                case HARD:
                    applyPointHard(wave);
                    break;
            }
            break;
        case PLANE:
            switch(style){
                case SOFT_SIMPLE:
                    applyPlaneSoftSimple(wave);
                    break;
                case SOFT_ADVANCED:
                    applyPlaneSoftAdvanced(wave);
                    break;
                case HARD:
                    applyPlaneHard(wave);
                    break;
            }
            break;
        case GAUSSIAN_BEAM:
            switch(style){
                case SOFT_SIMPLE:
                    applyGaussianBeamSoftSimple(wave);
                    break;
                case SOFT_ADVANCED:
                    applyGaussianBeamSoftAdvanced(wave);
                    break;
                case HARD:
                    applyGaussianBeamHard(wave);
                    break;
            }
            break;
    }
}



inline float KernelSource::fadeFunc(float time) const {
    return time<=0?0:(time<fadeTime?time/fadeTime:1);
}
//Algo: look up the cell and set phi to cos(omega*time)*(time<fadeTime?time/fadeTime:1)
inline void KernelSource::applyPointHard(WaveKernel* wave){
    if(int(x)<=0 || int(y)<=0 || int(x)>=wave->grid_n || int(y)>=wave->grid_m) {
        std::cerr<<"Point source must be inside the simulation domain."<<std::endl;
        return;
    }
    wave->getPhi(int(x),int(y))=amplitude*std::cos(omega*wave->time_elapsed+phase)*fadeFunc(wave->time_elapsed);
}
//Algo: Do the same as applyPointHard, but lerp between the old value and the new value.
inline void KernelSource::applyPointSoftSimple(WaveKernel* wave){
    std::cerr<<"applyPointSoftSimple not implemented yet."<<std::endl;
}
//Algo: Treat things as if we added a source term to the wave equation, eg. wave->phi[...]+=(calculated value)
inline void KernelSource::applyPointSoftAdvanced(WaveKernel* wave){
    std::cerr<<"applyPointSoftAdvanced not implemented yet."<<std::endl;
}
//Algo: Bresenham line along the source plane setting it equal to some cos(omega*t) value.
inline void KernelSource::applyPlaneHard(WaveKernel* wave){
    std::cerr<<"applyPlaneHard not implemented yet."<<std::endl;
}
//Algo: same but lerp a bit.
inline void KernelSource::applyPlaneSoftSimple(WaveKernel* wave){
    std::cerr<<"applyPlaneSoftSimple not implemented yet."<<std::endl;
}
//Algo: apply a source term over fadeDist / in a big rotated rectangle. This might be hard to do fast, so do it slow first and then check.
inline void KernelSource::applyPlaneSoftAdvanced(WaveKernel* wave){
    std::cerr<<"applyPlaneSoftAdvanced not implemented yet."<<std::endl;
}
//Algo: Bresenham line, set it equal to the Gaussian beam formula I found.
inline void KernelSource::applyGaussianBeamHard(WaveKernel* wave){
    std::cerr<<"applyGaussianBeamHard not implemented yet."<<std::endl;
}
//Algo: Bresenham line, lerp between phi[...] and the hard value.
inline void KernelSource::applyGaussianBeamSoftSimple(WaveKernel* wave){
    std::cerr<<"applyGaussianBeamSoftSimple not implemented yet."<<std::endl;
}
//Algo: apply a source term to create a Gaussian beam using formulas I derived.
inline void KernelSource::applyGaussianBeamSoftAdvanced(WaveKernel* wave){
    std::cerr<<"applyGaussianBeamSoftAdvanced not implemented yet."<<std::endl;
}

KernelSource createPointKernelSource(float x, float y, float omega, float amplitude, float phase, float fadeTime) {
    //Checks of validity have to be done inside Wave::Wave and Wave::WaveKernel, not here.
    KernelSource src;
    src.type=POINT;
    src.style=HARD;
    src.x=x;
    src.y=y;
    src.omega=omega;
    src.fadeTime=fadeTime;
    src.phase=phase;
    src.amplitude=amplitude;
    return src;
}

Image renderToImage(const WaveKernel &ker,int xCropping=0,int yCropping=0, ColorScheme scheme=RED_BLUE){
    if(ker.grid_n-2*xCropping<0 || ker.grid_m-2*yCropping<0)
        throw std::invalid_argument("Cropping is too large for the grid size.");

    Image outimg(ker.grid_n-2*xCropping,ker.grid_m-2*yCropping);
    switch(scheme){
        case RED_BLUE:
            for(int j=yCropping;j<ker.grid_m-yCropping;j++){
                for(int i=xCropping;i<ker.grid_n-xCropping;i++){
                    int index=j*ker.grid_n+i;
                    float amp=ker.phicur[index];
                    int red=amp>0?int(amp*255):0;
                    int blue=amp<0?int(-amp*255):0;
                    outimg.put(i-xCropping,j-yCropping,intToRGB(red,0,blue));
                }
            }
            break;
        case GRAYSCALE:
            for(int j=yCropping;j<ker.grid_m-yCropping;j++){
                for(int i=xCropping;i<ker.grid_n-xCropping;i++){
                    int index=j*ker.grid_n+i;
                    float amp=ker.phicur[index];
                    int c=std::abs(amp)*255;
                    outimg.put(i-xCropping,j-yCropping,intToRGB(c,c,c));
                }
            }
            break;
        default:
            throw std::invalid_argument("Invalid color scheme.");
    }
    //outimg.save(fname);
    return outimg;
}


/* A Configuration class for the different simulation domains we might want to set up.
The WaveKernel has no notion of "delta x" (although it does have a "delta t"), so Wave 
is responsible for abstracting things so that we don't have to deal with the WaveKernel directly.

As a general rule I will keep specific mathematical details out of this class. Instead, we'll have
functions like 
void setupBlackHoleConditions(Wave &sim, ...parameters...);
void setupOceanWaveConditions(Wave &sim, ...parameters...);
void setupDirichletBox(Wave &sim, ...parameters...);
void setupNeumannBox(Wave &sim, ...parameters...);

Once we set the width/height of the WaveKernel and instantiate it, we can't go back and call
setWidth,setHeight, and so on. So here's some state machine structure.


    // Constructor
    Wave(int width, int height);

    // Simulation functions
    Float timestep(Float dt);
    void saveImage(const std::string& filename, int cropping = 0);
    void setPhiSoft(int i, int j, Float to, Float weight);

    // Configuration and setup
    void configureGrid(int width, int height);
    void configureAbsorbingBoundaries(Float size, Float maxDamping);
    void addPlaneWave(Float posX, Float posY, Float kX, Float kY, const std::string& mode);

private:
    int n, m;
    std::vector<Float> phicur, philast, cellMass, fieldMass2, txx, tyy, txy, damping, dirich;

    // Helper functions
    void initializeGrids();
    void checkGridSize();

 */



class Wave {
private:

    //State machine. We need to go from DOMAIN_CONFIG to SIM_CONFIG to RUNNING.
    //DOMAIN_CONFIG allows us to change the simulation size, domain size, etc.
    //SIM_CONFIG allows us to modify simple settings like dt, nSubsteps
    //SIM_CONFIG also lets us add sources and sinks and modify tensions/densities inside the WaveKernel
    //RUNNING will warn us if we try to modify dt,nSubsteps. But it's not as critical a separation as DOMAIN_CONFIG.
    enum WaveState { 
        DOMAIN_CONFIG,
        SIM_CONFIG, 
        RUNNING,
    } state;

    struct DomainConfig {
        //domain size (left,bottom), (top,right). 
        float l,b,r,t;

        //the kernel will have grid_n=nx+2*npadx, grid_m=ny+2*npady;
        //ny is not stored as a variable, it's calculated to enforce a square cell size.
        int nx;
        int npadx,npady;
        //damping scales up linearly on the boundary padding region to maxDamping.
        float maxDamping;
        //speed of light, which affects the mass of each cell (not the tension)
        float speedOfLight;

        float getDx() const {
            return (r-l)/nx;
        }
        int getNy() const {
            return std::ceil((t-b)*nx/(r-l));
        }
        bool isValid() const {
            return t>b && r>l && nx>1;
        }

        //Coordinate transformations from domain <-> grid. 
        //returns/arguments are floats to allow for sub-pixel placement of sources
        //and evaluating eg. iToX(i+0.5).
        float xToI(float x) const { return std::floor((x-l)/getDx())+npadx; }
        float iToX(float i) const { return (i-npadx)*getDx()+l; }
        float yToJ(float y) const { return std::floor((y-b)/getDx())+npady; }
        float jToY(float j) const { return (j-npady)*getDx()+b; }

        DomainConfig() : l(0), b(0), r(1), t(1), nx(256), npadx(20), npady(20), maxDamping(5), speedOfLight(1) {}
    } domain;

    struct TimestepConfig {
        //the actual timestep value is stored in ker.dt, but we want to warn if it hasn't been set.
        bool dtSet;
        int nSubsteps;
        TimestepConfig() : dtSet(false), nSubsteps(1) {}
    } timestep;

    struct ImageConfig {
        std::string folder;
        std::string prefix;
        int digits;
        std::string format;
        ColorScheme colorScheme;
        ImageConfig() : folder("."), prefix("image"), digits(5), format("bmp") {}
        std::string getFilename(int index) const {
            return Utils::getFilename(folder+"/"+prefix,index,digits,"."+format);
        }
    } image;


    

    //if verbosity==0, print nothing except fatal errors.
    int verbosity;

    WaveKernel ker;


    //Only call this after the kernel has been instantiated
    void kernelSetDamping(){
        for(int k=0;k<domain.npadx;k++){
            for(int j=0;j<ker.grid_m;j++){
                ker.getDamping(k,j)+=(float(1)-float(k)/domain.npadx)*domain.maxDamping;
                ker.getDamping(ker.grid_n-k-1,j)+=(float(1)-float(k)/domain.npadx)*domain.maxDamping;
            }
        }
        for(int k=0;k<domain.npady;k++){
            for(int i=0;i<ker.grid_n;i++){
                ker.getDamping(i,k)+=(float(1)-float(k)/domain.npady)*domain.maxDamping;
                ker.getDamping(i,ker.grid_m-k-1)+=(float(1)-float(k)/domain.npady)*domain.maxDamping;
            }
        }
    }

    void instantiateKernel() {
        if(!domain.isValid()) 
            throw std::invalid_argument("Cannot instantiate kernel without a valid domain and nx.");

        // If the padding is negative, add 10% padding (5% on each side)
        if(domain.npadx<0)
            domain.npadx=int(std::ceil(0.05*domain.nx));

        if(domain.npady<0)
            domain.npady=int(std::ceil(0.05*domain.getNy()));

        int grid_n=domain.nx+2*domain.npadx, grid_m=domain.getNy()+2*domain.npady;
        if(verbosity>=2){
            std::cout<<"Instantiating simulation grid.\n";
            std::cout<<"domain=[["<<domain.l<<","<<domain.b<<"],["<<domain.r<<","<<domain.t<<"]]\n";
            std::cout<<"grid_size=["<<grid_n<<","<<grid_m<<"]\n";
            std::cout<<"max_damping="<<domain.maxDamping<<std::endl;
        }

        float dx=domain.getDx();

        //Instantiate kernel. The "mass" of each cell is related to the speed of light by c=sqrt(mu/T),
        //T defaults to 1, so we pass in mu dx^2 = c^2 T dx^2. 
        ker=WaveKernel(grid_n,grid_m,dx*dx*domain.speedOfLight*domain.speedOfLight);
        //Fill the appropriate damping values inside the Kernel
        kernelSetDamping();
        state=SIM_CONFIG;
    }

    void kernelCheck() {
        if(state==DOMAIN_CONFIG){
            instantiateKernel();
        }
    }


    //Check that the kernel is instantiated and lazily instantiate it if it isn't.
public:

    Wave() : state(DOMAIN_CONFIG), domain(), verbosity(2) {}
    void setVerbosity(int v) { verbosity=v; }

///////////////////////////////////////////////////////////////////////////////    
////////////////////////////////// DOMAIN_CONFIG /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////    

    void setDomain(float left, float bottom,float right, float top) {
        if(state!=DOMAIN_CONFIG)
            throw std::runtime_error("Cannot set domain after kernel is instantiated.");
        if(top<=bottom || right<=left)
            throw std::invalid_argument("Invalid domain.");
        domain.b=bottom; domain.l=left; domain.t=top; domain.r=right;
    }

    void setNx(int nx) {
        if(state!=DOMAIN_CONFIG)
            throw std::runtime_error("Cannot set nx after kernel is instantiated.");
        if(nx<=1)
            throw std::invalid_argument("Invalid nx.");
        domain.nx=nx;
    }

    // Set the padding area which linearly scales up to maxDamping.
    // Setting a padding to a negative number will add 5% padding on each side.
    void setDamping(int padx=20, int pady=20, float maxDamping=5) {
        if(state!=DOMAIN_CONFIG)
            throw std::runtime_error("Cannot set damping after kernel is instantiated.");
        domain.npadx=padx;
        domain.npady=pady;
        if(maxDamping>=0)
            domain.maxDamping=maxDamping;
        else 
            throw std::invalid_argument("Invalid maxDamping.");
    }

    //Set speed of light. By default this affects the mass of each cell, not the tensions.
    void setSpeedOfLight(float c) {
        if(state!=DOMAIN_CONFIG)
            throw std::runtime_error("Cannot set speed of light after kernel is instantiated.");
        if(c<=0)
            throw std::invalid_argument("Invalid speed of light.");
        domain.speedOfLight=c;
    }

///////////////////////////////////////////////////////////////////////////////    
////////////////////////////////// SIM_CONFIG /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////    
    void setDt(float dt) {
        kernelCheck();
        if(state==RUNNING) {
            if(verbosity>=1)
                std::cerr<<"Warning: Changing dt while simulation is running."<<std::endl;
        }

        if(dt<=0)
            throw std::invalid_argument("Invalid dt. Needs to be >0.");

        ker.dt=dt;
        timestep.dtSet=true;
    }
    void setNSubsteps(int nsub) {
        kernelCheck();
        if(nsub<=0)
            throw std::invalid_argument("Invalid nSubsteps. Needs to be >=1.");
        timestep.nSubsteps=nsub;
    }
    void addPointSource(float x, float y, float omega, float amplitude=1, float phase=0, float fadeTime=0.1) {
        kernelCheck();
        float i=domain.xToI(x), j=domain.yToJ(y);
        ker.addPointSource(i,j,omega,amplitude,phase,fadeTime);
    }
    //Original plan was to have KernelSource be public facing, but because of the coordinate transformation
    //this is a bad idea. In the future I could implement a Source object, then coordinate transform to a KernelSource.
    /*
    void addSource(const KernelSource &src) {
        if(state==DOMAIN_CONFIG)
            throw std::runtime_error("Cannot add sources before domain is configured");
        
        ker.sources.push_back(src);
    }*/


    void setColorScheme(ColorScheme scheme) {
        image.colorScheme=scheme;
    }
    void setOutputFolder(const std::string &folder) {
        image.folder=folder;
    }
    void setImageFormat(const std::string &prefix, int digits, const std::string &format) {
        image.prefix=prefix;
        image.digits=digits;
        image.format=format;
    }
    void generateImageSeries(float totalTime) {
        kernelCheck();

        if(verbosity>=2)
            std::cout<<"Generating image series.\n";
        int nsteps=int(totalTime/(ker.dt*timestep.nSubsteps));
        for(int k=0;k<nsteps;k++){
            for(int j=0;j<timestep.nSubsteps;j++){
                ker.timestep();
            }
            renderToImage(ker,domain.npadx,domain.npady,image.colorScheme).save(image.getFilename(k));
        }
    }

};

////////////////////////////////////////////////////////////////////////////////    
/////////////////////////////  FUNCTIONS /////////////////////////////
////////////////////////////////////////////////////////////////////////////////    

}





int main() {
    int imgx=640, imgy=480;
    int dampingSize=90;
    float LSim=18;
    float maxDamping=15.0;
    float aspect=float(imgy)/imgx;
    int fps=10;
    int substeps=4;


#if defined(_OPENMP)
    omp_set_num_threads(8);
#endif

    Wave::Wave wave;
    wave.setVerbosity(2);
    wave.setDomain(-LSim/2, -LSim*aspect/2, LSim/2, LSim*aspect/2);
    wave.setNx(imgx);
    wave.setDamping(dampingSize,dampingSize,maxDamping);

    wave.setSpeedOfLight(1.0);

    wave.setDt(0.01);
    wave.setNSubsteps(substeps);

    //Add five random sources with random frequencies.
    srand(time(NULL));
    float sz=1.0;
    for(int k=0;k<60;k++){
        float x=Utils::randf(-sz,sz)-7;
        float y=Utils::randf(-sz,sz);
        float omega=Utils::randf(-4,4);
        float phase=Utils::randf(0,2*M_PI);
        float amp=Utils::randf(0.1,2.0);
        wave.addPointSource(x,y,omega, amp,phase);
    }
    //wave.addPointSource(0,0,10,10);

    wave.setColorScheme(Wave::ColorScheme::RED_BLUE);
    wave.setOutputFolder("out");
    wave.setImageFormat("pointsource_gray_",5,"bmp");
    wave.generateImageSeries(30);
    return 0;
}

    



/* old timestepping code.
    for(int k=0;k<1000;k++){
        wave.timestep(dt);
        time+=dt;
        if(k%10==0){
            wave.saveImage(getFilename("out",k/10,5,".bmp"),dampingSize);
        }
    } */



/*
void setSchwarzschildConditions(Wave::Wave &sim, float rs, float cutoff) {
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
                //wave.damping[nx*j+i]=
                    //(r>rs)? (1.0f-(r-rs)/(dampingcutoff))*damping1:
                    //(1-r/rs)*(damping2-damping1)+damping1;
            }

            float reps=std::sqrt(x*x+y*y+0.2*0.2);
            wave.cellMass[nx*j+i]=myCellMass*(1.0f+2.0/reps);
            wave.txx[(nx-1)*j+i]=myCellMass*(1.0f+2.0/reps);

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
}
    */
   














