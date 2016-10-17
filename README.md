# polynomial-path-smoothing

This repository contains julia code for implementing polynomial path smoothing based on the work by Richter and Bry.

"poly.jl" Contains old code which is very likely buggy
"notebook.ipynb" is a jupyter notebook which contains an interactive example of the polynomial path smoothing approach.

# Getting started with the notebook
Install Julia and IJulia: https://github.com/JuliaLang/IJulia.jl
Start a Julia session, and run

```
using IJulia;
notebook();
```

Then open the appropriate .ipynb file and you're golden.

# Getting started with other code
The file `poly_helper.jl` contains the functions used to generate smooth paths. To see an example, open the julia CLI and try 
``` include("poly.jl") ```
``` test_multiseg(10) ```
