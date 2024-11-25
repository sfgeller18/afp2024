# AFP 2024

**Status:**

- Started working on convertible debt tree from goldman paper. Relevant class is NaiveConvertibleTree in tree.hpp.

- Added code from past project for cross platofrm pseudo-random number generation using for now just MKL but willl add the option for CUDA as well. This will be useful if we try to implement multi-factor Stochastic Process simulators for S, r and credit spreads using Sobol generator. Code currently includes methods for generating distributions with desired 3rd and 4th moments (standard = Gaussian, general = generalized gaussian (arbitrary 4th moment), skewGenNormal = generalized gaussian w/ skew).

