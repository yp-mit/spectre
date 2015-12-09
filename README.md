<img src="https://sites.google.com/site/durisi/software/spectre%20logo%20high-res.jpg?attredirects=0)" alt="spectre logo" align="middle" height="142" width="142"> 
short packet communication toolbox
===========================================



Shannon theory describes fundamental limits of communication and compression systems. Classic closed-form results (such as the well known log(1 + SNR) formula) apply only to the regime of infinite blocklength (infinite packet size/ infinite delay). 

For finite blocklengths, no closed-form results are usually obtainable, but there exist tight upper and lower bounds on fundamental limits, as well as approximations. This repository provides  numerical routines to compute these bounds and these approximations for some popular channel and source models.


Content
--------------------

Achievability bounds, converse bounds and approximation for the following source and channel models:

* AWGN single-antenna channel
* Quasi-static fading multi-antenna channel
* Rayleigh block-fading multi-antenna channel (no CSI)
* Rayleigh block-fading multi-antenna channel (full CSI at receiver)
* Binary memoryless source and Hamming distortion
* Gaussian memoryless source and quadratic distortion

Getting started
------------------

Consult the toolbox [manual](https://sites.google.com/site/durisi/manual.pdf?attredirects=0&d=1). For the latest version,
```
  cd documentation
  make manual.pdf
```



Want to contribute?
-------------------

The toolbox is under development and the participation of additional members of the information and communication theory communities to this endeavor  is warmly welcomed! 

For questions and access permissions: email <fblcode-list@mit.edu>



Contributors (alphabetic order)
---------------------------------------
* [Austin Collins](http://www.mit.edu/~austinc/) (MIT)
* [Giuseppe Durisi](https://sites.google.com/site/durisi/) (Chalmers)
* [Victoria Kostina](http://vkostina.caltech.edu/) (Caltech)
* [Johan Östman](https://www.chalmers.se/en/staff/Pages/johanos.aspx) (Chalmers)
* [Yury Polyanskiy](http://people.lids.mit.edu/yp/homepage/) (MIT)
* [Ido Tal](http://webee.technion.ac.il/people/idotal/) (Technion)
* [Wei Yang](https://sites.google.com/site/weiyangcth/) (Princeton) 
