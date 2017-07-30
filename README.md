# ADOLC and IPOPT wrapper for finite difference or Automatic Differentiation

ADOLC and IPOPT wrapper for finite difference or Automatic Differentiation. 

A converting tool from regular double to adouble is provided.

Based on the example from [Ipopt distribution](http://web.mit.edu/ipopt_v3.11.6/doc/documentation.pdf)

This is just for my own future reference.

[The optimization test case](http://ab-initio.mit.edu/wiki/index.php/NLopt_Tutorial "NLP Example")

## Tested on:

- Windows 10 - mingw64

## How to use:

- Create a class inhering from `(ad_)OptInterface` and define some optimization parameters and functions as x lower and upper bounds; gradient lower and upper bounds; objective function (override) and constraint function (override). In this case it is the `(ad_)optexample1.h/(ad_)myExample1`.

- Create the `main` function, which defines the case example, the Ipopt C++ wrapper class `TNLP`. To run the optimization call `createAppAndRun(nlp);`

- The converting tool from a regular double file to ADOLC capable adouble types is provided. For this is required to mark the doubles that will be activated as `double/*a*/` and also to define objects (classes and functions) that you be updated with a new name `ad_OLDNAME`. This tool requires a file and save a copy with the prefix `ad_`. It can be used as: `python(3) convertAdolc.py -f filename.h(.cpp)`

## TODO:

- Test on debian
- Check how the Ipopt was installed

*(Note: was used the [newcppproject.sh](https://gist.github.com/caiofcm/83d4d3d2370546d846454ff74dea7348) to generate the Makefile)*


