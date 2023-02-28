---
title: "Equivalent Prime Norm Ideals"
date: 2023-01-27T22:05:58Z
draft: false
mathjax: true
---

The first step of the KLPT algorithm is to find an ideal $L$ equivalent to the input ideal
$I$ with prime norm. The naive algorithm to compute equivalent ideals with prime norm is
fairly simple, but getting SQISign to work in practice means being very careful with this step.

Practically, we want an ideal with prime norm, as we will need to perform computations
modulo its norm. When it is prime, computations are easy and efficient.

On this page, we outline how we can naively implement this algorithm and then continue to cover 
the edge cases which appear for SQISign, how they can cause issues and how we can stop these complications
from becoming a problem. 

## Computing equivalent ideals with prime norm

The trick to finding an ideal $L \sim I$ with prime norm is to first find a *good* element 
of the quaternion algebra $\BB$ contained in the ideal, namely an element $\alpha \in I$. Every such element has reduced norm 
$\Nrd(\alpha) = N \cdot n(I)$ for some integer $N$. We are interested in the *scaled norm*: 

$$
q_{I}(\alpha) = \Nrd(\alpha) / n(I),
$$

which corresponds to the norm of the ideal 

$$L = \chi_I(\alpha) = I \frac{\bar{\alpha}}{n(I)}, \quad n(L) = q_I(\alpha),$$ 

equivalent to $I$. 
We see that if we can find an element $\alpha \in I$ of prime scaled norm, we can efficiently compute 
an equivalent ideal with prime norm as $\chi_I(\alpha)$. 

We can enumerate elements $\alpha \in I$ by computing a basis $\langle b_0, b_1, b_2, b_3 \rangle$
of the ideal and looking at linear combinations: $\alpha = \sum c_i b_i$ for $c_i \in \mathbb{Z}$. 
If all we need is an ideal with prime norm,
then heuristically we can compute random linear combinations and take the scaled norms
as random integers. A prime norm will be stumbled upon fairly
quickly.

For SQISign we impose a further restriction. Not only do we want $L$ to have prime norm, but we need the
norm to be small as possible. The size of $n(L)$ determines (roughly) the size of the output of KLPT
itself. The smaller the output of KLPT, the smaller norm ideals we have to work with and ultimately 
the smaller the degree of the isogenies of which we have to compute via the Deuring correspondence. 
For more detailed discussions of
the norms and KLPT bounds see the page [Estimating Bounds for KLPT](/posts/klpt-bounds/).


## Finding small algebra elements

The best we can hope for when looking for algebra elements $\alpha \in I$ with small prime scaled norm is to 
find many elements with small reduced norm and hopefully find one which is prime.

The first thing we do is compute the Minkowski reduced basis of the given ideal. For an ideal $I$ with maximal
left-order $\OO$,
the Minkowski reduced basis $I = \langle \beta_0, \beta_1, \beta_2, \beta_3 \rangle$ satisfies the bound

$$
p^2 \leq 16 q_I(\beta_0)q_I(\beta_1)q_I(\beta_2)q_I(\beta_3) \leq 4p^2.
$$

We can compute this basis in SageMath with the following snippet:

```py
def reduced_basis(I):
    """
    Computes the Minkowski reduced basis of the
    input ideal. 
    """

    def _matrix_to_gens(M, B):
        """
        Converts from a matrix to generators in the quat.
        algebra
        """
        return [sum(c * g for c, g in zip(row, B)) for row in M]

    B = I.basis()
    G = I.gram_matrix()
    U = G.LLL_gram().transpose()

    return _matrix_to_gens(U, B)
```

Note that the reduced basis of an ideal is the basis of the smallest equivalent ideal to $I$. This means
that given an ideal $J \sim I$, we would expect to get the same reduced basis elements inputting either $J$ or $I$.
This will become important later, when we consider what happens when the smallest equivalent ideal has particularly 
small norm.

Given the basis $\beta_0, \dots, \beta_3$, we can now find small norm elements as above by looking at 
random linear combinations of the basis. If we keep the coefficents $c_i \in \mathbb{Z}$ small, the 
output should be small too. However, we can do a little better than this by enumerating short
vectors of the lattice corresponding to the reduced basis. 

We implemented this short vector enumeration in [`lattices.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/lattices.py) in the function `generate_small_norm_quat()`.
Essentially, we take the four basis elements and construct a four by four matrix from them. We then use `fpylll`
to enumerate short vectors in the corresponding lattice, which we return as elements $\alpha \in \BB$.

We then iterate through these elements, compute the scaled norm, and simply return the first element whose 
scaled norm is prime.

For a random input $I$, we expect the basis elements to all be approximately of the size $q_I(\beta_i) \simeq \sqrt{p}$,
and so when enumerating short vectors of the lattice, we expect to find small solutions approximately this size.

Most of the time, this is indeed what happens and the algorithm finishes quickly. However, there are certain steps
in SQISign where we find ourselves in a special case and the norms of the basis elements $\beta_i$ are unbalanced.
These cases can mean it is very unlikely to small prime norm elements $\alpha \in I$. 

**Fool me once...** Before continuing, we mention one other (easy to solve) complication.
Sometimes when running KLPT we find a perfectly good $L$ with prime $n(L)$ but something later in the algorithm
fails and we need to start again by picking a new $L$ with new norm.
In the repeated run, if we were to deterministically return the first smallest prime norm element from our lattice
enumeration, we are likely to repeat the remaining steps and encounter the same error over and over again. 
To stop this, we store the prime norms we have discovered from previous runs
within a set `previous`. When a prime norm is found, if it is contained within `previous` we continue enumerating short vectors
until the next prime is found.  


## Sometimes primes are hard to find

When the norm of the smallest equivalent ideal to $I$ is small, then our LLL reduction inside 
`reduced_basis()` does *too good* of a job at finding short vectors and we find two basis elements with very small 
norm which are also usually orthogonal. 

The first problem with finding very small norm basis elements is that the determinant restriction of the Minkowski reduction
means the remaining two basis elements
have very large norm. As an example of this, take the unit ideal of $\OO_0$ which has a basis $\langle 1, i, (i + j)/2, (1 + k)/2 \rangle$. 
The basis is already Minkowski reduced, and the elements have norm: $(1, 1, (p+1)/4, (p+1)/4)$. 
Taking random linear combinations would give you elements of norm $\Nrd(\alpha) \simeq p$, far too big for SQISign to run as intended.

So, why can we not just take random combinations of the first two elements in these cases. The problem here comes from the orthogonality. If we have 
$\beta_1 = i \beta_0$ then we have $\Nrd(\beta_0) = \Nrd(\beta_1)$. We can then compute the reduced norm of a linear combination of these as:
$\Nrd(c_0\beta_0 + c_1\beta_1) = (c_0^2 + c_1^2)\Nrd(\beta_0)$. Unless we have $\Nrd(\beta_0) = 1$ or one of $c_i = 0$ then the result
will never be prime. 

There is no magic trick to sorting this problem out directly. There will be ideals where we will be very unlikely to be able to find small prime norm elements of the reduction produces bases like this.
The tricks we use in SQISign is simply to avoid having to deal with these
bad ideals by using alternative methods.


### Bad ideals and where to find them

For random ideals, this problem is very unlikely to happen. However, for SQISign, this problem appears when the input ideal $I$ with left order
$\OO_0$ and right order $\OO$ connects two elliptic curves which are close on the isogeny graph. Essentially, the closeness of the curves
means there is a small degree isogeny linking them and via the Deuring correspondence there is some small norm ideal $J \sim I$. This ideal
has the "bad" Minkowski reduced basis and causes an issue in our algorithm. This is discussed in much more detail in the page
[Small Steps from Curves with Small Endomorphisms](/posts/small-steps/) 

## Prime powers

One fix we can apply in a few places within SQISign is that there is an ideal $M$ equivalent to the input ideal $I$
which has norm $\ell^\times$. For our implementation this is always $\ell=2$, a power
of two. 
As mentioned above, we require the norm to be prime so we can perform arithmetic during later stages of the KLPT algorithm.
Despite this, although slightly more complicated, we can take $n(M) = \ell^\times$ instead and most steps are unchanged.
As long as we have $n(M)$ small enough (about $\sqrt{p}$), then we can try and run everything using $N = 2^\times$ and skip
finding an equivalent prime norm ideal all together. 
This slightly complicates the `IdealModConstraint()` and `StrongApproximationHeuristic()` algorithms, but not enough to
stop this being a good fix. 
Essentially, care has to be taken that certain integers can be inverted mod $2^k$, which essentially just means
making sure values are odd and finding new values when we know something is going to break.

## Here's one I made earlier

Sometimes it is the case that we are working with some input ideal $I$ which is equivalent to the ideal $I_\tau$ from `keygen()`,
which has prime norm $n(I_\tau) \simeq p^{1/4}$.
When this is the case, we can skip this finding a prime norm ideal and simply set $L = I_\tau$. 

In our implementation, this is used in two places:

- In `keygen()` the last iteration of `IsogenyFromIdealFromKLPT()` has as input an ideal equivalent to $J_\tau \sim I_\tau$. By passing 
  $I_\tau$ through, we can skip this algorithm for the last call to KLPT, saving us from a bad Ideal.
- In response, we make a call to `IsogenyFromIdealCoprime()` which in part needs us to run KLPT on the connecting ideal $J_\tau \sim I_\tau$.
  Again, we simply supply $I_\tau$ as an optional argument to the KLPT algorithm and skip computing $L$ for this case.

## "Prime" norms

This last modification not only helps find good algebra elements when the basis is particularly reduced, but also improves the 
efficiency of the algorithm for all steps.

The ultimate goal of KLPT is to find an ideal $J \sim I$ which has norm dividing some target norm $\mathcal{T}$. 
Following the [C/Pari implementation](https://github.com/SQISign/sqisign) of SQISign, we can relax the condition that $N = n(L)$ is prime and instead look for 
$N = L_0 N'$ for some prime $N'$ and $L_0 \mathrel{|} \mathcal{T}$. Ultimately, we can absorb this new factor $L_0$ in the 
algorithm, perform computations modulo $N'$, and return an ideal $J \sim I$ such that $n(J) = L_0 L_1 L_2 \mathrel{|} \mathcal{T}$.

Although the resulting norm of $J$ will be increased by $L_0$, we dramatically improve the chance of finding elements
$\alpha \in I$ which are suitable for the algorithm. Concretely, as more $\alpha \in I$ are *good*, we
are more likely to be able to pick among the smallest elements generated by the lattice reduction, ensuring $N'$ as small as possible.
Additionally, as composite values can now be accepted, if the norm of the ideal is small enough that we find an orthogonal
reduced basis, we are much more likely to find a suitable $\alpha \in I$.

This small section gives an overview of how this works and the implementation is contained within all the functions 
needed for the KLPT algorithm in [`KLPT.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/KLPT.py).

### Allowing composite equivalent ideals

Firstly, given a generator $\alpha$ of $L = \OO \langle\alpha, N \rangle$ we can recover a prime norm ideal $L' = \OO \alpha + \OO N' $.
Note that although $L'$ has prime norm, it will not be equivalent to the input ideal $I$. 
Disregarding this problem for now, we continue through most of the KLPT algorithm as normal to find 
$\beta' \in L'$ of reduced norm $\Nrd(\beta') = N' L_1 L_2$, simply by using $L', N'$ 
instead of $L, N$ throughout all the steps algorithm.

However, when $L_0 \neq 1$, our computation will end with some element $\beta' \in L'$ rather than $\beta \in L \sim I$.
So, despite having been able to find and work with the prime norm ideal $L'$, as it is not equivalent to our input ideal $I$ 
we cannot expect the value $J = \chi_L(\beta')$ to be an integral ideal equivalent to the input $I$.

We correct for this problem at the end of the `strong_approximation_wrapper_KLPT()`. The goal is to, given $L$, $L'$ and $\beta' \in L'$, 
compute a new element $\beta \in L$ for which $L_0 L_1 L_2 \mathrel{|} \Nrd(\beta)$. Once we have this, we can construct
an ideal $J$ from $L$ equivalent to $I$ with the correct norm.

First, compute a generator $\delta$ such that $L = \langle \delta, N \rangle$ with the additional constraint that

$$
\gcd(\Nrd(\delta), L_0 L_1 L_2 n(L)^2) = n(L).
$$

Given $\delta$, we can then compute $\beta = \delta \bar{\beta}'$, which will have reduced norm 

$$\Nrd(\beta) = f N L_1 L_2 = f N' L_0 L_1 L_2,$$ 

where $f$ is some unimportant cofactor with $\gcd(f, n(L)) = 1$. Now, as $\beta \in L$ by construction, and 
also has $L_0 L_1 L_2$ as a factor of its reduced norm, we can finish the algorithm. 
Given $\OO$, the left-order of $L$, we can compute $J \sim L \sim I$ from 

$$
J = \OO \beta + \OO L_0 L_1 L_2,
$$

and have a good ideal to return from the KLPT algorithm with a relaxed and faster implementation of finding 
equivalent "prime" norm ideals $L$.

In our implementation, we use this trick for every call of the KLPT algorithm, with the exception of the `SigningKLPT()`.
This is because when we compute our response isogeny, it is vital to control the size of the output ideal $J$ 
as it is fixed to ensure security of SQISign. 


[Back to Top](#top)