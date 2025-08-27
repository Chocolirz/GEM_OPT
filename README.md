# Simulating GEM detector using Garfield++

This is part of the code used for simulating electrical, optical readouts, etc. in thick glass-GEMs for the MIGDAL experiment. Different folders contain different things we want to measure. 

## Projects

- [GemWithRim](https://github.com/Chocolirz/GEM_OPT/tree/main/GemWithRim) is the main program for single-GEM simulations, which supports and generates data for other projects. 
- [Drift](https://github.com/Chocolirz/GEM_OPT/tree/main/Drift) contains results for the collection ratio. 
- [field](https://github.com/Chocolirz/GEM_OPT/tree/main/field) is there to visualise electric field at different positions, field data should be generated from the main program in correct format (.csv or .txt, I mainly used .txt). 
- [Map](https://github.com/Chocolirz/GEM_OPT/tree/main/Map) is there to generate distribution of photons at the focused plane (bottom of the third GEM) and photon collection map for three GEMs.
- [ParallelPlate](https://github.com/Chocolirz/GEM_OPT/tree/main/ParallelPlate) is an independent project built upon simple geometry, it serves as a starting point to get familiar with ```Garfield++```
- [Photon](https://github.com/Chocolirz/GEM_OPT/tree/main/Photon) generates distribution of photons at the ITO plane for a single GEM. 
- [Stack](https://github.com/Chocolirz/GEM_OPT/tree/main/Stack) investigates electrical readout for three GEMs. 


Slides for some results can be found [here](https://indico.stfc.ac.uk/event/1629/) (does not include optical readout for three GEMs). 

## How to run a program

For ```.C``` files, we execute them using ```cmake```. First, check that everything has been set properly in the ```CMakeLists.txt``` file. Next, in its folder, run

```mkdir build```, ```cd build```

and then 

```cmake ..```, ```make -j(number of cores)```

if everything works fine, you will see something like

```built excutable (name)```

that means you can run the project by 

```./(name)```.

For python it's simpler, just use

```python ./(name).py```

or ```python3``` for some environments. 