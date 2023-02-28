---
title: "Estimating Bounds for KLPT"
date: 2023-01-28T00:08:19Z
draft: false
mathjax: true
---

The KLPT algorithm takes an input left-$\OO_0$ ideal $I$ and computes
an equivalent ideal $J$ with norm dividing some smooth integer $\TT$.
At a high level, this is accomplished with the following steps:

1. Compute an ideal $L \sim I$ with prime norm: $n(L) = N$.
2. Compute $\gamma \in \BB$ with reduced norm: $\Nrd(\gamma) = L_1 N$ with $L_1 \mathrel{|} \TT $.
3. Compute $\mu \in \BB$ such that $\gamma \cdot \mu \in L$.
4. Compute $\nu \in \BB$ with $\Nrd(\nu) = L_2$ with $L_1 L_2 \mathrel{|} \TT $, using $\mu$.
5. Compute $\beta = \gamma \cdot \nu \in L$ with reduced norm $\Nrd(\beta) = L_1 L_2 N$.
6. Compute $J = \chi_L(\beta)$ with $n(J) = L_1 L_2$ such that $J \sim L \sim I$.

In this discussion, we will not detail *how* all these steps work, but rather how
we choose the bounds we place on deriving the values $L_1$ and $L_2$ such that we can heuristically 
rely on the algorithm to return $J$ with high probability. 

Applying KLPT to SQISign, we want $L_1 L_2$ to be as small as possible to ensure that each
step is efficient. Most importantly, we need this to divide the available torsion on the curve
so we do not have to work with extension fields. 
For `SigningKLPT()`, $L_1 L_2 = \ell^e$ is fixed from the parameter set, for security considerations,
but keeping the output norm small is important, as this controls the degree of the response isogeny
$\sigma$, and this computation is by far the most expensive in the SQISign protocol.

## Searching for small primes

Although it is maybe not obvious from the steps outlined above, making $n(J)$ small
is all wrapped up in finding the smallest possible value of $n(L) = N$. 

The trick to solving Step 1 is computing the Minkowski reduced basis of the input
ideal: $I = \langle  \beta_1, \beta_2, \beta_3, \beta_4 \rangle$ and taking random 
linear combinations of these basis elements until we find
some $\alpha \in I$ with norm $\Nrd(\alpha) = n(I) \cdot N$. We can then compute the equivalent 
ideal $L \sim I = \chi_I(\alpha)$ with norm $n(L) = N$.

**Note**: for random input, we would expect all elements of the reduced basis to have
reduced norm $\Nrd(\beta_i) \simeq \sqrt{p}$, and so we expect the norm of $L$ to be
$n(L) = N \simeq \sqrt{p}$. However, in SQISign, we occasionally have that ideals
have particularly small norm. When this happens, the reduced basis is orthogonal with
two elements of exceptionally small reduced norm. In these cases, it can be impossible
to find prime norm ideals with $n(L) < p$. This is discussed in more detail on the page
[Equivalent Prime Norm Ideals](/posts/prime-ideal/).

For this discussion, let us write $p \simeq 2^{2\lambda}$ and assume we have obtained some ideal $L$ 
with prime norm $N \simeq 2^{\lambda + \delta}$, where $\delta$ is some small positive 
value. The question is, how large can $\delta$ be before we are likely to find a solution?
We wish to pick $\delta$ small enough that we have a small output for $n(J)$, but large enough
that we are overwhelmingly likely to find solutions for the equivalent prime norm ideal and
hence for the KLPT algorithm.

## Finding a quaternion algebra element with given norm

The next step is to find $\gamma$ with norm $M = L_1 \cdot N$. This is accomplished with 
`RepresentIntegerHeuristic()` and requires that the norm $M > p$. To find $\gamma$,
we have to be able to find a solution $M' = x^2 + y^2$ using Cornacchia's algorithm. 
On random input, we expect this to take roughly $\log(p)$ attempts, and so a good rule of thumb is to ensure
that $M > p\log(p) \simeq 2^{2\lambda + \epsilon}$. Assuming that $N \simeq 2^{\lambda + \delta}$
we then have a bound for $L_1 \simeq 2^{\lambda + \epsilon - \delta}$ where we expect 
$\epsilon \simeq \log\log (p)$.

Practically, we do this by taking the input target norm $\TT$ and factoring it. We then
randomly take factors from $\TT$ for use in $L_1$ until $M$ is large enough such that we expect
to find a solution for $\gamma$. We then update $\TT^\prime = \TT / L_1$ for the following steps.

## Finding a solution to the strong approximation 

Step four is where $N$ comes back into focus. Our goal is to find some $\nu$ with norm 
$L_2$ such that $L_1 L_2$ divides the target norm $\TT$. From the estimates above, we
then have a maximal size for $L_2$ given $L_1$ computed earlier:

$$
L_{2}^{\text{Max}} = \TT / L_1 \simeq \TT \cdot 2^{\delta - \lambda - \epsilon}.
$$ 

Using the lattice method described in [Strong Approximation Lattice Trick](/posts/strong-approximation-lattice/), 
we know that an approximate minimum bound for
the reduced norm is $\Nrd(\nu) \simeq pN^3$. This means practically we require
$\Nrd(\nu) = L_2 \geq pN^3$ for there to be a reasonable chance at finding a close vector.

Combining these two bounds, we must have that:

$$
L_{2}^{\text{Max}}  \simeq \TT \cdot 2^{\delta - \lambda - \epsilon} \geq pN^3 \simeq 2^{5\lambda + 3\delta},
$$

which allow us to represent our heuristic bounds in terms of the target norm $\TT$:

$$
\TT > 2^{6\lambda + 2\delta + \epsilon} \simeq p^3 \cdot 2^{2\delta + \epsilon}.
$$

## Concrete bounds for SQISign

When running the KLPT algorithm in SQISign, we use it to find ideals with norm that divides the
available odd torsion. Further tricks described on the page [Computing Isogenies from Ideals](/posts/ideal-to-isogeny/) 
allow us to square this such that $\TT = T^2$. We can then estimate our bounds by considering:

$$
\frac{T^2}{p^3} > 2^{2\delta + \epsilon},
$$

where both $p$ and $T$ are fixed constants for the protocol. 
For the SQISign prime $p_{6983}$, we have that $p \simeq 2^{256}$ and $T \simeq 2^{428}$.
Using the bounds above, we have that $T^2 p^{-3} \simeq 2^{89}$. Allowing $\epsilon = \log\log(p) \simeq 8$ 
for solving $\textsf{RepresentInteger}_{\OO_0}$ we find the maximum possible value which allows solutions 
would be $\delta = 40$.

Practically, we find that allowing $\delta \simeq \epsilon \simeq \log\log(p)$ seems to work with good
probability, however fine-tuning these constants will require more testing. Allowing $\delta$ to be smaller
than its maximum makes finding the prime norm a little harder, but makes finding a solution to the strong approximation
a little easier. As we randomly sample $L_i$, the bounds input into the strong approximation are more fuzzy, so
it appears helpful to save some bits by shrinking $\delta$. 

As described on the page [Equivalent Prime Norm Ideals](/posts/prime-ideal/), there are a few additional tricks which help keep the size of $N$ small, so setting
$\delta \simeq \log\log(p)$ seems fairly robust.


## Setting the signing length of SQISign

For the SQISign identification protocol to be zero-knowledge, the output of `response()` must have fixed length $\ell^e$.

The output of `response()` is a degree $\ell^e$ isogeny, which is derived from a cyclic ideal $J$ of norm $n(J) = \ell^e$. 
The ideal $J$ is found by running `SigningKLPT()`, a generalisation of the KLPT algorithm which allows as input an
ideal $I$ with generic left order $\OO$ providing that a connecting ideal $I_\tau$ with left order $\OO_0$ and right order
$\OO$ is known.

We decompose $e$ as $e_1 + e_2$, where $L_1 = \ell^{e_1}$ and $L_2 = \ell^{e_2}$. The goal is to minimise $e$ for 
efficiency of the protocol, while allowing it to be large enough that heuristically we always find an equivalent 
ideal from `SigningKLPT()`.

**Setting** $e_1$: In SQISign, we have a connecting ideal $I_\tau$ of prime norm, so we need not find an equivalent prime norm ideal. 
Following the efficient `keygen()` from appendix D of the SQISign paper, the norm is also particularly 
small: $n(I_\tau) \simeq p^{1/4}$. For the input ideal $I$, we will need an equivalent prime norm ideal $L \sim I$ 
which we expect
to have norm $n(L) \simeq \sqrt{p}$. We can set $e_1$ such that $L_1 \cdot n(L) > p$, and it very likely that we 
find $\gamma \in \BB$ with norm $L_1 \cdot n(L)$. 

**Setting** $e_2$: As we want $e$ to be fixed, we then just allow $e_2 = e - e_1$, but how small can we make $L_2$
while still finding solutions? The size of $L_2$ becomes important when we attempt to solve a generalisation of the 
strong approximation where the modulus is a composite integer: $N = n(L) \cdot n(I_\tau)$. 

Using the lattice trick described in [Strong Approximation Lattice Trick](/posts/strong-approximation-lattice/), 
we expect to find solutions of size $L_2 > pN^3$. Using the estimates
above, we have $L_1 > p^{1/2}$ and $L_2 > p^{13/4}$, so we should have $e > \frac{15}{4} \log_\ell(p)$. For the 256-bit
prime used in SQISign this would restrict us to $e > 960$ for a theoretical bound and $e = 1000$ is a comfortable 
heuristic bound.

In practice, we set this bound for $e$ in [`setup.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/setup.py) as a global variable. Then, when running `SigningKLPT()`, we set a sensible
value for $e_1$ once we find $n(L)$ for $L \sim I$.
Finally, we set $e_2 = e - e_1$ and attempt to solve the strong approximation for fixed $L_2$. If successful, we find
an ideal $J$ of norm $\ell^e$ equivalent to the input ideal $I$.

### Fixing the signing length

In fact, there is a little additional complication in setting $e_2$. When we compute an ideal using the KLPT 
algorithm, the output ideal is not necessarily cyclic. So although we may have some ideal $J$ such that $n(J) = \ell^e$,
once we scale the ideal to ensure $J$ is cyclic (this is necessary such that the corresponding isogeny is cyclic), 
we may find $J'$ has norm $n(J') = \ell^{e-\epsilon_0}$, for some small $\epsilon_0$.

To deal with this, `SigningKLPT()` attempts to "overshoot" $L_2$ and pick some $\epsilon_0$ so we have $e_2 = e - e_1 + \epsilon_0$.
The output from the KLPT algorithm with then be some ideal of norm $n(J) = \ell^{e + \epsilon_0}$ and we hope that the cyclic
ideal $n(J')$ has norm exactly $\ell^e$.

**Note**: we were unable to derive a way to ensure that the output had the correct norm. So ultimately we guess $\epsilon_0$ as best
as we can, and if $J'$ has the wrong norm, we just run `SigningKLPT()` again. This does not add too much overhead to the computation
time, mainly because of how expensive the isogenies are currently.

Our current technique is to take 

$$
\epsilon_0 = 2 \left\lfloor \frac{\log\log(p)}{4} \right\rfloor + \epsilon(\gamma). 
$$

Where $\epsilon(\gamma)$ is a dynamic correction we make depending on the derived $\gamma \in \BB$ with
$\Nrd(\gamma) = L_1 n(L)$. 

First, we write $\gamma$ in the basis of $\OO_0$ and compute $g$, the $\gcd$ of its coefficients.
We then compute the $\epsilon(\gamma) = \gcd(g, L_1)$. This might be easier to follow in code,
so we include the estimate below:

```py
def quaternion_change_basis(γ, O):
    """
    Computes the coefficients of a quaternion γ
    in the basis of a given order O
    """
    O_matrix = Matrix([b.coefficient_tuple() for b in O.basis()])
    γ_vector = vector(γ.coefficient_tuple())
    γO_coeffs = γ_vector * O_matrix.inverse()

    assert γ == sum([a * b for a, b in zip(γO_coeffs, O.basis())])
    return γO_coeffs


def quarternion_basis_gcd(γ, O):
    """
    Computes the gcd of the coefficients of a
    quaternion γ in the basis of a given order O
    """
    γO_coeffs = quaternion_change_basis(γ, O)
    return gcd(γO_coeffs)


def derive_L2_SigningKLPT(γ, L1, e1):
    """
    Given L1 = l^e1 and γ try and compute L2
    so that the output of SigningKLPT has norm
    exactly 2^e
    """
    g = quarternion_basis_gcd(γ, O0)
    extra = 2 * (floor(loglogp / 4) + ZZ(gcd(g, L1).valuation(l)))
    e2 = e - e1 + extra
    return l**e2
```

[Back to Top](#top)