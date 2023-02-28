---
title: "Future Work"
date: 2023-01-01T22:07:49Z
draft: false
mathjax: true
---


## New SQISign

Currently, our implementation follows the original 2020 paper: [SQISign: compact post-quantum signatures from quaternions and isogenies](https://eprint.iacr.org/2020/1240), by Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski.

However, in 2022, an eprint was put online called *New algorithms for the Deuring correspondence*, which has new descriptions of ideal to kernel algorithms
which estimate an approximate two-fold speed up for both keygen and signing. 

At the time of writing this blog, the paper has just been updated ready for EUROCRYPT 2023: [New algorithms for the Deuring correspondence: toward practical and secure SQISign signatures](https://eprint.iacr.org/2022/234), by Luca De Feo, Antonin Leroux, Patrick Longa and Benjamin Wesolowski.

This is a really beautiful paper with some exciting results. Aside from anything else, the new implementation is 3-4x faster
across all of the SQISign protocol!

With the publication of these new results, this SageMath implementation can really only be considered as a prelude 
to implementing modern SQISign, but everyone has to start somewhere.

The main focus of future work will be implementing all of the additional improvements from this new paper so that we can
bring the SageMath implementation of SQISign up to the current standard. However, before we do this, there are a few 
low hanging fruit left in our current implementation that should be addressed. We discuss a few of these below.

## Faster Isogenies

The current bottleneck of our implementation of SQISign is the computation of the large $\ell$-isogenies
when computing the various isogenies needed within `IdealToIsogenySmallFromKLPT()`. In fact, when running
SQISign, approximately 75% of the total computation time is spent within `EllipticCurveIsogenyFactored()`,
which is a wrapper function around SageMath's own `EllipticCurveIsogeny` which has a few performance improvements,
such as using the `velusqrt` algorithm when $\ell > 400$.

When we started working on SQISign, we took a shortcut in following the specification. Rather than write our own
$x$-only isogeny computations, we extended the base field to $\mathbb{F}\_{p^4}$ so that the available torsion 
for our curves $E / \mathbb{F}\_{p^4}$ is $(p^2 - 1)$. This meant sacrificing performance 
(one would expect arithmetic for $\mathbb{F}\_{p^4}$ to be about four times slower than over $\mathbb{F}\_{p^2}$), 
but gaining the ability to use all of SageMath's built-in functions for isogenies and to focus our attention on all the other 
algorithms we would need for SQISign.

Now, with a working implementation of SQISign, our focus can shift to refactoring the code and working with
$x$-only operations, so we can compute isogenies on $E / \mathbb{F}_{p^2}$. This should end up with an approximate 
four times speed up. However, once all the additional code surrounding the multiplication is taken into account, 
we find that the speed-up is not nearly so dramatic.

As a first example, let us take some large prime $\ell$ which divides $(p+1)$. We can compute an isogeny of this degree
over either $\mathbb{F}\_{p^2}$ or $\mathbb{F}\_{p^4}$, and this should guide us to how much time we could save:

```py
sage: p = 73743043621499797449074820543863456997944695372324032511999999999999999999999
sage: l = 6983
sage: assert l.divides(p+1)
sage: # Computation over Fp4, as is currently performed
sage: Fp4.<z4> = GF(p^4)
sage: E4 = EllipticCurve(Fp4, [1,0])
sage: P4 = ((p^2-1) // l) * E4.random_point()
sage: assert P4.order() == l
sage: 
sage: # Compute isogeny using Vélu's formula
sage: %timeit E4.isogeny(P4)
2 s ± 31.1 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
sage: # vélusqrt offers a significant speed up
sage: %timeit E4.isogeny(P4, algorithm="velusqrt")
276 ms ± 8.85 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
sage:
sage:
sage: # Now, perform the "same" computation over Fp2 instead
sage: Fp2.<z2> = GF(p^2)
sage: E2 = EllipticCurve(Fp2, [1,0])
sage: P2 = ((p+1) // l) * E2.random_point()
sage: assert P2.order() == l
sage:
sage: # Compute isogeny using Vélu's formula
sage: %timeit E2.isogeny(P2)
1.65 s ± 38.6 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
sage: # vélusqrt offers a significant speed up
sage: %timeit E2.isogeny(P2, algorithm="velusqrt")
174 ms ± 3.37 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)
```

From this tiny experiment, it seems that the saving from multiplications in $\mathbb{F}_{p^2}$ are not dominant as one may expect due to the other aspects of the SageMath functions.
Performing over $\mathbb{F}\_{p^2}$ for Vélu's formulas seems to save about 20%, while for the $\sqrt{\text{élu}}$
algorithm, the saving is closer to 40%.
We see that there is still time
to be saved, so it is important for us to work on, but how this work should be done is not so obvious.

One option would be to ditch the in-built SageMath isogeny computations and work on our own $x$-only isogeny implementation.
This would require more development and would move away from "SQISign in SageMath" to something closer to 
"SQISign in Python", but it may result in a multiplicative speed up like we are hoping for. We can also specialise
to Montgomery and Edwards curves, for which there are particularly efficient implementations of computing the
action and codomain of isogenies. These are not suitable for general purpose isogenies in SageMath (not every
elliptic curve can even be written in Montgomery form), but it would be suitable for us and align more closely with the C implementation. 


## Friends of Friends of Quaternions

While working on implementing SQISign, a complementary piece of research 
[Deuring for the People: Supersingular Elliptic Curves with Prescribed Endomorphism Ring in General Characteristic](https://ia.cr/2023/106)
by Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni, was published. 

In the paper, there are several nice implementation tricks to improve the general performance of sub-algorithms while 
working with the Deuring correspondence, most of which carry over directly for SQISign. One in particular is noticing 
a simplification to the `ideal_to_kernel()` algorithm: by evaluating the endomorphism action of the *conjugate* of the ideal 
generator, the ideal's kernel is directly computed, avoiding expensive discrete logarithm computations. This is discussed in more
detail on the page [Between Kernels and Ideals](/posts/ideal-kernel/).

Following their description, we have implemented this in our own code. There is still more to gain from implementing 
various other improvements shown in the paper.

Their paper was accompanied with SageMath code for the [Constructive Deuring Correspondence](https://github.com/friends-of-quaternions/deuring).
As well as some clean and easy to read code for the KLPT algorithm and ideal to isogeny computations for generic extension fields, 
this code included $x$-only isogeny computations. Adapting this code for our own uses seems like a perfect starting point 
to improve the performance of the code, but there is an interesting problem to overcome.

In their code, an $x$-only isogeny is computed by first deriving the kernel polynomial using $x$-only computations.
The isogeny itself is computed using Kohel's algorithm from the kernel polynomial using in-built SageMath functions. 
However, when the degree of the isogeny grows (and hence the degree of the kernel polynomial), Kohel's algorithm becomes slow: 

```py
sage: p = 73743043621499797449074820543863456997944695372324032511999999999999999999999
sage: l = 6983
sage: assert l.divides(p+1)
sage: Fp2.<z2> = GF(p^2)
sage: E = EllipticCurve(Fp2, [1,0])
sage: P = (p+1) // l * E.random_point()
sage: # Compute the isogeny using Velu's algorithm
sage: %timeit E.isogeny(P)
1.74 s ± 35.3 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
sage: # Compute the isogeny use Kohel's algorithm 
sage: h = E.isogeny(P).kernel_polynomial()
sage: %timeit E.isogeny(h, check=False)
5.54 s ± 500 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)
```

We see that even though the $x$-only computations would allow us to remain over $\mathbb{F}_{p^2}$, when we have to
work with large $\ell$-degree polynomials for Kohel's, SageMath gets really slow. It is an interesting problem, and 
something which might be able to be addressed within SageMath.

Ultimately, it means that simply adapting the $x$-only code from *Deuring for the People* will slow down, rather 
than speed up our own implementation for SageMath.

## Community Feedback

Most of all, we hope that by sharing this code, we encourage others to work on SQISign and other protocols that use the Deuring 
correspondence. We have done our best to write code which is easy to read and adapt, but all code can be improved, either by 
generalising, making functions more robust or including additional documentation and examples.

If you are interested in helping with this implementation, either through improvements of existing code or the addition of new code,
please get in contact. 
The author's contact details can be found on their personal pages: [Maria](https://www.mariascrs.com), [Giacomo](https://giacomopope.com).

We also warmly welcome issues being filed on the [GitHub repo](https://github.com/LearningToSQI/SQISign-SageMath), 
or pull requests to include new features.

[Back to Top](#top)