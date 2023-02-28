---
title: "Strong Approximation Lattice Trick"
date: 2023-01-26T09:02:23Z
draft: false
mathjax: true
---

## Overview

The KLPT algorithm, which we call `EquivalentSmoothIdealHeuristic()`, takes as input a left $\OO\_0$-ideal $I$ and a smooth integer $\TT$, and computes an equivalent ideal $J$ with norm $\TT$. There are four main steps:
1. Compute an ideal $L \sim I$ of norm $N$ using `EquivalentPrimeIdealHeuristic()`.
2. Next, we run `RepresentIntegerHeuristic()` to obtain $\gamma \in \OO\_0$ of norm $N\TT$.
3. Solving the `IdealModConstraint()` on input $L, \gamma$ derives $(C_0 : D_0) \in \mathbb{P}^1(\mathbb{Z}/N\mathbb{Z})$ such that $\mu_0 = j(C_0 + \omega D_0)$ where $\gamma\mu_0 \in L$. 
4. Run `StrongApproximationHeuristic()` on input $N$, $C_0, D_0$ to obtain $\mu = \lambda\mu_0 + N\mu_1 \in \OO\_0$ of norm dividing $N\TT/n(\gamma)$. 
5. Set $\beta = \gamma\mu$, then output $J = \chi_L(\beta)$.

In practice, we run through the factors of $\TT$ (that are large enough for the algorithms to run successfully), so that $J$ has norm *dividing* $\TT$. 

In this post, we focus on the fourth step: the *strong approximation step*. 


## Strong approximation algorithm: random solutions

In the following, we identify the distinguished quadratic subring $R[\omega]$ of $\BB$.
Define a function $f_{\omega}(x,y)$ to return the norm of the element $z = x + \omega y \in R[\omega]$. 
In SQISign, we have $p \equiv 3 \bmod 4$, so $\omega = i$ and $f_\omega(x,y) = x^2 + y^2$. 
For the rest of this page we assume that $\omega = i$. 

The generic idea of `StrongApproximationHeuristic()` is
to sample $y, z$ in some randomized way and try to see if we can find a solution $t, x$ to the norm equation $f(t, x) = M − pf(y, z)$ 
using Cornacchia's algorithm.

Let $\mathcal{D}(\TT)$ be the set of divisors of $\TT$.  The $\textsf{StrongApproximation}_{\mathcal{D}(\TT)}$ algorithm works as follows.

Setting $\mu_0 = j(C_0 +  iD_0)$ as in `IdealModConstraint()`, we want to find $\mu_1$ such that $\mu = \lambda\mu_0 + N\mu_1$ has norm in $\mathcal{D}(\TT)$. Our $\mu_1$ will be of the form: $$\mu_1 = x + iy + j(z + it),$$ and so we want to find a solution $(\lambda,x,y,z,t)$ to the norm equation 
$$n(\mu) = N^2f(x,y) + pf(\lambda  C_0 + Nz, \lambda D_0 + Nt) = M,$$ where $M \in \mathcal{D}(\TT)$.
Working modulo $N$, we get $$M = pf(\lambda C_0, \lambda D_0) = p\lambda^2 f(C_0, D_0) \bmod N.$$

So, we first compute $\lambda$ such that 

$$
\lambda^2 = \frac{M}{pf(C_0, D_0)}\bmod N.
$$

However, such a $\lambda$ exists if and only if the RHS is a quadratic residue modulo $N$. Therefore, in the first step, we run through $M \in \mathcal{D}(\TT)$ with $M \geq N^4p$ (see [Estimating Bounds for KLPT](/posts/klpt-bounds/)) until $M / (pf(C_0, D_0))$ is a quadratic residue, and then compute the corresponding $\lambda$. 

Next, we generate a random pair $(z, t)$ such that:
$$
M - pf(\lambda C_0 + Nz, \lambda D_0 + Nt) = 0 \bmod {N^2}.
$$ We do this by picking $z$ randomly and then solving a linear equation modulo $N$ to recover $t$. 

Given $\lambda$, we now compute $$M^\prime = \frac{M - pf(\lambda C_0 + Nz, \lambda D_0 + Nt)}{N^2},$$ and determine if the equation $f(x,y) = M^\prime$ has any solutions. If no solution exists we sample another random pair $(z,t)$. 

Once we have found a solution $(x,y)$ for sampled $(z,t)$, we output $$\mu = \lambda j(C_0 + D_0i ) + N(x + i t + j(z + i t)).$$

The solution $\mu$ will have norm $\Nrd(\mu) = M \geq pN^4$. However, we can obtain smaller solutions creating a special lattice, applying lattice reduction
and using carefully chosen integers rather than random ones.

## Efficiently generating short vectors

The aim of this section is to explain how we can reduce the size of the output. We do this by obtaining *good* solutions $(z,t)$ rather sampling a random pair. We define good solutions as the ones corresponding to a small value of $pf(\lambda C_0 + Nz, \lambda D_0 + Nt)$. The key to obtaining such solutions is to observe that the space of solutions modulo $N^2$ is a translated lattice.

To find good solutions, we do the following:
1. Construct the appropriate lattice containing these solutions.
2. Construct a *target vector* of the form $$\begin{bmatrix}
z^\prime \\\ t^\prime
\end{bmatrix} + N\begin{bmatrix}
z_0 \\\ t_0
\end{bmatrix}$$ and find vectors in the lattice close to this vector. Enumerate these close vectors.
3. For each close vector $(z_c,t_c)$, construct $z = \frac{z_c}{N} + z_0$ and $t = \frac{t_c}{N} + t_0$ and if the congruence $$M - pf(\lambda C_0 + Nz, \lambda D_0 + Nt) = 0 \bmod N^2, $$ holds, compute the integer $$M^\prime = \frac{M - pf(\lambda C_0 + Nz, \lambda D_0 + Nt)}{N^2}.$$ If $M^\prime > 0$ continue to the next step, otherwise, try the next close vector. 
4. Using the Cornacchia algorithm, determine whether there exists a solution $(x,y)$ the norm equation $f(x,y) = M^\prime$. If no such solution exists, return to Step 3 and try another close vector. Otherwise, output $\mu = \lambda\mu_0 + N\mu_1$, where $\mu_0 = j(C_0 +  iD_0)$ and $\mu_1 = t + \omega x + j(y + iz)$.

This is implemented in `strong_approximation_construct_lattice()` and `strong_approximation_lattice_heuristic()` 
in [`KLPT.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/KLPT.py).

We now describe steps 1, 2, and 3 in more detail.

### Step 1: constructing the lattice

We first describe how to construct the appropriate lattice. We want 
$$
\begin{aligned}
&M - pf(\lambda C_0 + Nz, \lambda D_0 + Nt) = 0 \bmod N^2\\\
\implies &M - p\lambda^2f(C_0, D_0) - 2p\lambda N \left(  C_0  z +  D_0 t \right) = 0\bmod N^2 \\\
\implies &c_0 - N \left(  c_z z +  c_t t \right) = 0\bmod N^2.
\end{aligned}
$$ where $c_z = 2p\lambda C_0$, $c_t = 2p\lambda D_0$ and $c_0 = M - p\lambda^2f(C_0,D_0)$.
We translate and scale the lattice by $$(z,t) \mapsto \big(z, t - c_0c_t^{-1}\big).$$ Then, dividing through by $N$, we obtain the relation
$$
\begin{aligned}
c_z z + c_t t  = 0\bmod N,
\end{aligned}
$$ for which we get a basis of solutions $(1, - c_zc_t^{-1} \bmod N), (0, N)$. 
We consider the lattice defined by matrix
$$
N\begin{bmatrix}
1 & -c_y\cdot c_z^{-1} \bmod N\\\ 0 & N
\end{bmatrix}.
$$ Here, the factor of $N$ multiplying the lattice matrix comes from the fact that we divided through by $N$. If $(z,t)$ is contained in this lattice, then $(\frac{z}{N}, \frac{t}{N} + c_0c_t^{-1})$ is a good candidate for a solution to $$M - pf(\lambda C_0 + Ny, \lambda D_0 + Nz) = 0 \bmod N^2.$$ 

By considering the discriminant of this lattice and using the Gaussian heuristic, we can show the smallest size of solution we can expect with $M^\prime > 0$ positive is $pN^3$, rather than the $pN^4$ expected by taking random solutions.
Therefore, using this lattice we can then hope to find smaller solutions and correspondingly pick a smaller $M \in \mathcal{D}(\TT)$ of size approximately $pN^3$. 

Though it may be unclear by reading this page in isolation, as $\mu$ has smaller reduced norm, this will decrease the norm of the output ideal and hence the degree of any isogenies computed from ideals generated via KLPT. Ultimately leading to a speed up for `keygen()` and `response()` and `verify_response()` of the SQISign identification protocol (respectively, `signing()` and `verify()` in the signature scheme). 

### Step 2: constructing the target vector

To obtain the target vector, we first try to find a solution to
$$
M - p\lambda^2f(C_0, D_0) - 2\lambda C_0pNz - 2\lambda D_0p Nt = 0\bmod N^2.
$$ The $(z,t)$, however, must lie in $\mathbb{Z}$, and so we set our target vector to be an (integer) approximation of this solution, namely:
$$
\mathbf{v}_t = \begin{bmatrix}
y \\\ z
\end{bmatrix} = - \lambda \cdot \begin{bmatrix}
C_0 \\\ D_0
\end{bmatrix} - N\cdot \begin{bmatrix}
0 \\\ c_0\cdot c_z^{-1} \bmod N
\end{bmatrix},
$$ where $c_z$ as above and $c_0 = M - p\lambda^2f(C_0,D_0)$.
This may not lie on the lattice, but we can hope to find a solution close enough so that the norm is kept close to $pN^3$. By enumerating vectors close to this target vector, we can find good solutions. 

### Step 3: finding vectors close to the target 

We use `generate_close_vectors()` implemented in [`lattices.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/lattices.py) to find vectors close to the target (i.e., within a certain bounded distance from the target). It works as follows.

1. Use `solve_closest_vector_problem()` to compute the vector $\mathbf{v}_{\text{closest}}$ in the lattice closest to the target vector $\mathbf{v}_t$. Set $d$ to be the distance between $\mathbf{v}_t$ and $\mathbf{v}$. 
2. Let $b_0 = \lfloor \frac{M}{p} \rfloor$ and set the norm bound to be $B =  \lfloor b_0 + d \rfloor + 2\sqrt{b_0 + d}$. 
3. Find short vectors of size less than $B$ in the lattice using `generate_short_vectors()`.
4. Using the short vectors $\mathbf{v}\_{\text{short}}$ computed in the step above, output close vectors by computing the sum $\mathbf{v}\_{\text{closest}} + \mathbf{v}\_{\text{short}}$.

### Enumerating short vectors using SageMath

The method described above can really only be efficient if it is easy to enumerate many short vectors of a lattice.

SageMath has a few ways to do this, one of which being associated to the `QuadraticForm` class with `short_vector_list_up_to_length(n)`,
but it seemed to be slow and error prone. A first iteration of our implementation called to the Pari function `qfminim()` which is fairly efficient,
but loads all short vectors into memory so you have to choose how many to look through and compute them all before testing if you have a 
good solution.

We instead used the `fpylll` library which comes with SageMath, and which allows the enumeration of short vectors. This means
we can built an iterator and test each short solution as we find it. We include the code snippet below, as we hope this could
be useful in many other contexts.


```py
from fpylll import IntegerMatrix, Enumeration, EvaluatorStrategy
from fpylll.fplll.gso import MatGSO

def generate_short_vectors(lattice_basis, bound, count=2000):
    """
    Given a lattice `lattice_basis`, an upper bound for target 
    norms and the total count of short vectors wanted, creates 
    a generator of short vectors of the lattice with 
    norm < bound.
    """
    # LLL reduce the lattice_basis
    L = lattice_basis.LLL()

    # Move from Sage world to Fypll world
    A = IntegerMatrix.from_matrix(L)

    # Gram-Schmidt Othogonalization
    G = MatGSO(A)
    _ = G.update_gso()

    # Enumeration class on G with `count` solutions
    E = Enumeration(
        G, nr_solutions=count, strategy=EvaluatorStrategy.BEST_N_SOLUTIONS
    )

    # We need the row count when we call enumerate
    r = L.nrows()

    # If enumerate finds no solutions it raises an error, so we
    # wrap it in a try block
    try:
        # The arguments of enumerate are:
        # E.enumerate(first_row, last_row, max_dist, max_dist_expo)
        short_vectors = E.enumerate(0, r, bound, 0)
    except:
        short_vectors = []
        
    for _, xis in short_vectors:
        # Returns values x1,x2,...xr such that
        # x0*row[0] + ... + xr*row[r] = short vector
        v3 = vector([ZZ(xi) for xi in xis])
        v = v3 * L
        yield v
```

## Edge cases 

### Powers of two

As described in [Equivalent Prime Norm Ideals](/posts/prime-ideal/), occasionally it is required that we input an ideal $L$ with norm $N = 2^k$, rather than some large prime.
When this is the case, we will not be able to invert $c_z$ as $2 \mid \gcd(c_z, N)$. Therefore, we must modify `StrongApproximationHeuristic()` to handle this. 

When constructing the lattice, we double the modulus $N$, and half $c_y$ and $c_z$. More precisely, set the modulus to be $2N$ and $c_y$ and $c_z$ to be $p\lambda C_0$ and $p\lambda D_0$, respectively. We then run Step 3 in the same way. 

Note that care also has to be taken in the previous step with `IdealModConstraint()`, to see
a little more information, see the comments in `check_ideal_mod_constraint()`, which is called to make sure $(C_0 : D_0)$ are suitable values. 

### Composite modulus

Another edge case happens in `SigningKLPT()`, where the input modulus $N$ is the product of two large primes. The Strong Approximation step is actually
no different for this case, but instead extra care is taken before computing the short lattice solutions by making sure $\gamma$ is picked in such a way 
that $\lambda^2$ is indeed a quadratic residue. Note for efficiency, we will compute $\lambda$ as a square root for each of the factors of $N$ and then
recombine them using the Chinese remainder theorem.

For more information on running the KLPT for `SigningKLPT()`, the original [SQISign Paper](https://eprint.iacr.org/2020/1240) discusses this edge case very carefully.


[Back to Top](#top)