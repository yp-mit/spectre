This folder contains numerical evaluations of the finite-blocklength lossless source coding bounds for a binary symmetric source model. More details are given in the comments in each `.m` file.

* Point-to-point (almost-lossless) bounds:
    * `p2p_rcu`: (achievability) RCU bound [1, Th. 4]
    * `p2p_han_converse`: Han's converse [2, Lemma 1.3.2]
    * `p2p_ht_converse`: binary hypothesis testing converse [3, Appendix A]
    * `p2p_optimum`: optimum rate
* Slepian-Wolf bounds:
    * `sw_rcu`: (achievability) RCU bound [1, Th. 18]
    * `sw_han_converse`: Han's converse [2, Lemma 7.2.2]
    * `sw_ht_converse`: composite hypothesis testing converse [1, Th. 19]
* Gaussian approximation:
    * `gaussian_approx`: third-order Gaussian approximation
* Helper functions:
    * `binomial_coeff`: compute logarithmic values of binomial coefficients, avoiding overflow
    * `my_logsumexp`: compute logarithmic sum of the exponentials of the given vector, avoiding overflow
    * `my_logdiffexp`: compute logarithmic difference of the exponentials of two given number, avoiding overflow
    * `my_bisect`: called by `sw_ht_converse` to compute the parameters in composite hypothesis tests
    * `p2p_ht_converse`: called by `sw_ht_converse` to get initial values for the parameters in composite hypothesis tests
    
    
[1] S. Chen, M. Effros, and V. Kostina, “Lossless source coding in the point-to-point, multiple access, and random access scenarios,” ArXiv preprint, 2019. Available at [https://arxiv.org/abs/1902.03366](https://arxiv.org/abs/1902.03366). <br />
[2] T. S. Han, Information-Spectrum Methods in Information Theory. Springer-Verlag Berlin Heidelberg, 2003. <br />
[3] V. Kostina and S. Verdu, “Fixed-length lossy compression in the finite blocklength regime,” IEEE Trans. Inf. Theory, vol. 58, no. 6, pp. 3309–3338, Jun. 2012.

