Examples of using the API:

"Light Waves" config:
Parabolic mirror

``` python

import wave

imgx=1280
imgy=800
# Create a new simulation configuration with the bottom left at (-4,-4) and the top right at (4,4).
# Speed of light set to 1. Tell the API that we want our configuration options to be printed out.
sim=wave.domain(bl=[-4,-4],tr=[4,4],c=1, verbose=True)

# Set absorbing boundary conditions: a padding region of 0.5 will be added all around the simulation,
# with a damping value slowly ramping up towards the edges of the simulation. 
sim.absorb(size=0.5, max_damping=5)

# set the number of x cells to use in the simulation domain. The number of y cells will be decided 
# automatically. In this case, the number of y cells is also 1280, the delta_x of the simulation is
# 8/1280, and because we add a padding of 1 the whole simulation domain actually consists of 
# ceil(9 *1280/8 )= 1440 cells. (or, a grid of 1440*1440 cells).
sim.setNx(imgx)

# Specify that we want a plane wave source with wavevector 5 starting from the left of the screen
# and going towards the right.
sim.planewave(position=[-4,0], k=[1,0], sh="soft")

# construct a parabolic mirror on the right side of the simulation domain with Dirichlet boundary conditions.
parabola_pos=4
parabola_height=2
parabola_const=0.5
parabola_nseg=10
parabola_dy=2*parabola_height/parabola_nseg

sim.strokeSize(0.1)
mirror=sim.begin_line(mode="dirichlet")
for y_ind in range(parabola_nseg+1):
    y=-parabola_height+y_ind*parabola_dy
    x=parabola_pos - parabola_const*y**2
    mirror.vertex(parabola_pos-parabola_height,y)
mirror.end_line()


# duration=16 seconds (so, since c=1 and our simulation domain has size 8, we simulate for two crossing times.)
# nframes=16*60, so 16 seconds of 60fps video
# 3 substeps per frame. So the total dt will be dt=1/(60*3). 
sim.config_parameters(duration=16, nframes=16*60, nsubsteps=3)

# output images of the format "img/run1_000.png", delete images of format "img/run1_*.png", and overwrite images if they exist.
sim.config_img(outputdirectory="img",fname="run1_",format="png",overwrite=True,clearimg=True)

# Run the simulation.
sim.simulate()
```


Rain drop 42 degree reflection
Waves in a maze
Beam splitter
Total internal reflection
Wifi signal
Fresnel lens
Scattering in a field of water droplets
Importing the example files from falstad's simulation
Incoherent source
Incoherent beam
Reflection-transmission drawing

"Water Waves" config:
Simulate 






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

    void timestep(Float dt)
    void saveImage(std::string fname,int cropping=0)

    void setPhiSoft(int i, int j, Float to, Float weight)
};



# More Thoughts
How would I keep everything in real units?

One benefit is that it would be easy to generate preview videos (at the lowest resolution that is attainable). Then you could just multiply by a global quality factor and get the high res version. 

Image output size needs to be specified in pixels.

For the time config, the user can specify:
    delta_t, number of frames, substeps per frame (my current setup)
    Total time, frames per second of simulation time (float). dt would be chosen automatically.
    Total time, nframes. 

It would be good to be able to 



# List of Parameters
I don't think I'm going to go the route of having a super smart configuration manager. Though it would be 
an interesting thing to try (think: CSS-like). I think I'm just going to say that you 
have to specify my favorite set of variables (imgw,imgh, domain left/right and bottom, numSuperSample=-1,
tension=1 and speedOfLight, totalTime, numFrames, substepsPerFrame=-1). Everything else
is calculated.


- Highest wavenumber k
- imgw, imgh, numSuperSample=1
- domainWidth, domainHeight
- deltat, numFrames, substepsPerFrame
- paddingHoriz, paddingVertical, maxPadding
- deltax, gridn, gridm,
- speedOfLight, density, tension
- cfl<1, kSafetyFactor < 1

``` mathematica
hardConstraints = {deltax==(domainWidth)/(imgw numSuperSample),
gridn==imgw nSuperSample+2 paddingHoriz/deltax,
gridm==imgh nSuperSample+2 paddingVertical/deltax,
domainHeight/domainWidth==imgh/imgw,

speedOfLight deltat/deltax == cfl (* we want CFL<1 *),
speedOfLight==Sqrt[tension/density],
totalTime==deltat numFrames substepsPerFrame,
k == kSafetyFactor pi/deltax (* we want all of our wavenumbers to be <<maxK *)
};

softConstraints = { kSafetyFactor<1, cfl<1, domainLeft<domainRight, domainTop<domainBottom,
}

```

