---
title: "Home"
date: 2023-01-22T00:08:19Z
draft: false
mathjax: true
---

# Introduction

This blog has been written as supplementary material to support the [SageMath implementation 
of SQISign](https://github.com/LearningToSQI/SQISign-SageMath). The aim of this implementation is to create an educational resource for the 
interested isogenist or cryptographer and focuses on being clear and verbose. We do not claim
any novelty in the blog posts (see references below and given in the code descriptions
for the papers we followed).
Additionally, the implementation aims (within reason) to allow user-supplied parameters, which
we hope will aid in research and experimentation to enable progress on optimising the protocol.  

The implementation follows the initial 2020 description of [SQISign](https://eprint.iacr.org/2020/1240), 
and in particular does not
include the most recent improvements from
[New algorithms for the Deuring correspondence: Towards practical and secure SQISign signatures](https://eprint.iacr.org/2022/234).
Including these optimisations is planned for future work.

Although the code has been written with plenty of comments, there are certain aspects of SQISign
which require careful attention or special tricks to deal with edge cases. This is the primary focus of 
this blog, and we hope that these pages will act as a guide for anyone else interested in implementing
SQISign themselves.

Implementing SQISign has been a lot of fun. More than anything, we hope that with the SageMath code 
and these pages, more people will be able to learn about and experiment with the protocol.

### Acknowledgements 

[SQISign-SageMath](https://github.com/LearningToSQI/SQISign-SageMath) was implemented by Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer 
and Giacomo Pope. 
We thank Antonin Leroux and Lorenz Panny for their advice and helpful 
discussions while writing the SQISign implementation.

The blog *Learning to SQI* was written by Maria Corte-Real Santos and Giacomo Pope. We thank our 
code co-authors, Jonathan Komada Eriksen and Michael Meyer for their helpful feedback while writing 
this material. 

## Existing Implementations 

When published in 2020, the paper was accompanied by a C/Pari implementation, as well as a
Magma reference implementation:

- [SQISign](https://github.com/SQISign/sqisign), C/Pari (Official),
- [SQISign](https://github.com/SQISign/sqisign-magma), Magma (Official Reference).

Similarly, when the faster description was published in 2022, a new implementation
with various improvements was also released:

- [SQISign](https://github.com/SQISign/sqisign2), C/Pari (Official).

The most recent implementation follows the updated 2023 paper 
*"New algorithms for the Deuring correspondence - Towards practical and secure SQISign signatures"* by 
Luca De Feo, Antonin Leroux, Benjamin Wesolowski and Patrick Longa and is approximately 3-4 times faster
than the original implementation:

- [SQISign (EC 2023)](https://github.com/SQISign/sqisign-ec23), C/Pari (Official).


## References

This blog does not aim to be self-contained and much of the context for the objects and
algorithms discussed come from the SQISign paper and various references.

- [SQISign: compact post-quantum signatures from quaternions and isogenies](https://eprint.iacr.org/2020/1240), Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski (2020).
- [New algorithms for the Deuring correspondence: toward practical and secure SQISign signatures](https://eprint.iacr.org/2022/234), Luca De Feo, Antonin Leroux, Patrick Longa and Benjamin Wesolowski (2022).
- [Quaternion algebras and isogeny-based cryptography](https://www.lix.polytechnique.fr/Labo/Antonin.LEROUX/manuscrit_these.pdf), Antonin Leroux (PhD Thesis) (2022).
- [On the quaternion $\ell$-isogeny path problem](https://arxiv.org/abs/1406.0981), David Kohel, Kristin Lauter, Christophe Petit, Jean-Pierre Tignol (2014).
- [Supersingular isogeny graphs and endomorphism rings: reductions and solutions](https://eprint.iacr.org/2018/371), Kirsten Eisenträger, Sean Hallgren, Kristin Lauter, Travis Morrison Christophe Petit (2018).
- [An improvement to the quaternion analogue of the $\ell$-isogeny path problem](https://crypto.iacr.org/2018/affevents/mathcrypt/page.html), Christophe Petit and Spike Smith (2018).
- [Cryptographic Smooth Neighbors](https://eprint.iacr.org/2022/1439) by Giacomo Bruno, Maria Corte-Real Santos, Craig Costello, Jonathan Komada Eriksen,
Michael Meyer, Michael Naehrig, and Bruno Sterner (2022).
- [Deuring for the People: Supersingular Elliptic Curves with Prescribed Endomorphism Ring in General Characteristic](https://ia.cr/2023/106) Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni (2023).

[Back to Top](#top)