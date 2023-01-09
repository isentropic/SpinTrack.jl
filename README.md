# SpinTrack.jl
This software is specifically designed to perform for storage ring proton EDM
experiment. 

This software is developed by Zhanibek Omarov under supervision from Yannis K.
Semertzidis and Selcuk Haciomeroglu during Ph.D. program at KAIST and IBS-CAPP. 

## Installation instructions
Install Julia - [link](https://julialang.org/downloads/platform/)

Next, paste this in the Julia prompt (REPL):
```
import Pkg;
Pkg.add(url="https://github.com/isentropic/SpinTrack.jl")
Pkg.add("Plots")
```

If you want Jupyter notebooks, refer to [this link](https://julialang.github.io/IJulia.jl/stable/manual/installation/)

## Getting started
Once the Installation instructions are performed, paste the following to get the
most basic tracking example:

```
using SpinTrack
using Plots

p = symmetric_hybrid_ring()
sol = get_solution(u1_long(1e-6), p);
plot(sol)
```
In case `Plots` library errors, you might want to try a different plotting
[backend](https://docs.juliaplots.org/latest/backends).

## More examples
More additional examples could be found in [examples/](examples/) in jupyter
notebook format. 

## Useful links
1. https://julialang.org/
2. https://docs.julialang.org/en/v1/manual/getting-started/


## TODO items
1. Add comprehensive documentation
2. Add unit tests 
3. Add CI via github actions
4. Add links to the thesis, paper and notes
5. Add more tutorials

## Attribution
Please provide the appropriate attribution if used in your works or find helpful.

1. Bibtex:
``` bibtex
@software{Omarov_SpinTrack_jl_2022,
author = {Omarov, Zhanibek},
month = {8},
title = {{SpinTrack.jl}},
version = {0.0.1},
year = {2022}
}
```
2. Or use the github sidebar "citation" interface
3. Please give a github star for visibility.
