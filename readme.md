
# Basic Wave Simulator

Example of a Gaussian beam:

![](gitimg/example-gauss.gif)

Example of an obstacle with Dirichlet boundary conditions:

![](gitimg/example-obstacle.gif)


Currently adding nice configuration options. See [a previous commit](https://github.com/physbuzz/wave/commit/8b09db706fece08278c0f010503d1de3669316d0) for
the functional Schwarzschild black hole run config.

# TODO

 - [x] Store some static math functions (like the Schwarzschild G, Kerr G, Newtonian G) somewhere
 - [x] Store the renderer as a static function (with options for color scheme).
 - [x] Write a Wave class for configuring initial conditions
 - [x] decide how I want to handle the transformation between Source position and kernel position
 - [x] Write up at least the hard point source
 - [ ] Code up the previous configs and make sure they compile and run
 - [ ] Let image deal with LodePNG.
 - [ ] Test the Kerr black hole scenario
 - [ ] Find my previous method for CLI and implement some options (changing dt, nx,ny on the fly).
 - [ ] Implement the DirichletLine thing
 - [ ] Litmus test: if I can run all of my previous test scenarios and make videos easily without recompiling, even at different nx and ny, then I have something highly configurable and awesome! From there I could choose to work on more features (probably yes) or optimization (meh?). Loading all of the math.falstad examples would be a cool start!

Niceties:
 - [ ] Allow dt to be changed before domain is initialized
 - [ ] Ensure that Source::fadeTime is always positive / checks are made.

Start writing a Wave class to handle all of the configuration.