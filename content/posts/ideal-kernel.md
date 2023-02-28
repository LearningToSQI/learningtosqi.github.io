---
title: "Between Kernels and Ideals"
date: 2023-01-24T15:20:23Z
draft: false
mathjax: true
---

Central to constructive applications of the Deuring correspondence is 
efficiently converting between equivalent ideals and isogenies. The core of this work 
is converting between a cyclic ideal $I$ of norm $n(I) = D$ 
and a kernel generator $K$ of order $D$. This page is dedicated to discussing how to perform
these conversions and follows the description of $\textsf{IdealToKernel}$ and $\textsf{KernelToIdeal}$, 
described by Algorithms 19 and 20 of 
[Leroux's Thesis](https://www.lix.polytechnique.fr/Labo/Antonin.LEROUX/manuscrit_these.pdf).

Although we have attempted to make this page self-contained, we point an interested
reader to both the first explicit 
formulation in [[GPS17]](https://eprint.iacr.org/2016/1154), 
and the clear and well motivated discussions from 
[Leroux's Thesis: Section 4.2](https://www.lix.polytechnique.fr/Labo/Antonin.LEROUX/manuscrit_these.pdf), as well as the recent publication 
[[EPSV23]](https://eprint.iacr.org/2023/106).

## Evaluating endomorphisms

Central to converting between kernel generators and ideals 
is being able to compute the isomorphism $\xi : \End(E) \to \OO$ and
subsequently using this to evaluate the action of some algebra element 
$\alpha \in \BB$ on points $P \in E$ via the mapping $\xi^{-1}(\alpha)(P)$.

The construction of the isomorphism comes from our knowledge of the basis
of the maximal order $\OO_0$ and the endomorphism ring of the special 
curve $E_0$:

$$
\begin{aligned}
\OO_0 &= \left\langle 1, i, \frac{i + j}{2}, \frac{1 + k}{2} \right\rangle, \\\
\End(E_0) &= \left\langle 1, \iota, \frac{\iota + \pi}{2}, \frac{1 + \iota \pi}{2} \right\rangle, \\\
\end{aligned}
$$

where $\BB = \langle 1, i, j, k \rangle$, and $\iota$ and $\pi$ are the twisting 
and Frobenius endomorphisms respectively:

$$
\iota : (x,y) \mapsto (-x, \sqrt{-1} y), \qquad \pi : (x,y) \mapsto (x^p,y^p).
$$

Concretely, these endomorphisms can be evaluated on any point $P \in E_0$ with the
following functions:

```python
# Global values
p = E0.base().characteristic()
sqrt_minus_one = E0.base()(-1).sqrt(all=True)[0] # Deterministic

def E01(P):
    """
    Identity map, does nothing
    """
    return P


def E0ι(P):
    """
    Returns ι(P) = (-x, √-1 y)
    """
    return E0(-P[0], sqrt_minus_one * P[1])


def E0π(P):
    """
    Returns π(P) = (X^p, Y^p, Z^p)
    """
    return E0(P[0] ** p, P[1] ** p, P[2] ** p)


def E0ιπ(P):
    """
    Returns ιπ(P) = (-X^p, √-1 Y^p, Z^p)
    """
    return E0(-P[0] ** p, sqrt_minus_one * P[1] ** p, P[2] ** p)
```

With the endomorphisms implemented, we now need to look at how we can map algebra elements
to endomorphisms. Generically, any element $\alpha \in \BB$ can be written in the following
way:

$$
\alpha = a_0 + a_1 i + a_2 j + a_3 k, \quad a_i \in \QQ.
$$

However, we will be working with integral, cyclic ideals $I$ with left order $\OO_0$, and we 
will be evaluating either basis elements of $\OO_0$ or the action of generators of the ideal 
$I = \OO_0\langle \alpha, D \rangle$. For these cases, our elements can always be written as:

$$
\alpha = \frac{a_0 + a_1 i + a_2 j + a_3 k}{d}, \quad a_i \in \ZZ, \quad d \in \\{1, 2\\}.
$$

We can then compute the isomorphism $\xi : (1, i, j, k) \mapsto ([1], \iota, \pi, \iota \circ \pi)$,
allowing us to express the endomorphism $\theta_\alpha = \xi^{-1}(\alpha) \in \End(E_0)$:

$$
[d] \circ \theta_\alpha  = [a_0] + [a_1] \circ \iota + [a_2] \circ \pi + [a_3] \circ \iota \circ \pi ,
$$

where $[m]$ denotes the endomorphism of scalar multiplication by $m$. 

When $d = 1$, we can simply evaluate the right-hand side to obtain the action of $\theta_\alpha$
on some point $P \in E_0$.
For the case of $d = 2$, we must first compute some $Q$ such that $P= [2]Q$ such that 
$\theta_\alpha(P) = [2] \circ \theta_\alpha (Q)$.

SageMath allows us to easily compute all $Q$ such that $P = [2]Q$ with the function `P.division_points(2)`.
We only need one such point, so we can simply compute `Q = P.division_points(2)[0]`. 

**Note**: for SQISign, it is important that we can compute such a $Q$ without needing to
perform a field extension to allow the division. This requires that $P$ cannot have maximal 
even torsion. 
Understanding this is particularly important when we are taking 
step sizes in `IdealToIsogenyFromKLPT()` by computing isogenies from ideals of norm $2^\text{step}$.
The fact we may need to divide by two during this conversion
means that if the maximal available torsion is $2^f$, then we are restricted to pick
$\text{step} \leq f - 1$.

### Putting this all together

With all the above discussion, we have a fairly friendly looking implementation of evaluating 
the action of $\xi^{-1}(\alpha)$ on a point $P$ of order $D$:

```python
def eval_endomorphism(α, P, D):
    """
    Evaluates the action of an endomorphism
    f ∈ End(E0) on a point P ∈ E.
    """
    # Basis for the endomorphisms, functions as above
    EndE0 = [E01, E0ι, E0π, E0ιπ]

    # Unpack the coefficients of the generator α such that 
    # α = (w + x*i + y*j + z*k) / d for d,w,x,y,z in ZZ
    d, *α_coeffs = α.denominator_and_integer_coefficient_tuple()

    # For generators of integral ideals, we expect the denominator 
    # to be at most 2
    assert d in (1, 2), "Something is wrong with the input ideal"
    if d == 2:
        # Divide out by two before evaluation if needed
        P = P.division_points(d)[0]

    # Compute the image of α(P)
    # α(P) = [α0] * P + [α1] * ι(P) + [α2] * π(P) + [α3] * ι(π(P))
    αP = sum(ai * θ(P) for ai, θ in zip(α_coeffs, EndE0))

    return αP
```


### Pick a curve

Technically, we can work with any curve $E$ providing that we know
the connecting isogeny $\phi : E_0 \to E$, and in our implementation 
[`deuring.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/deuring.py), we have allowed for this with an optional parameter
`connecting_isogenies` which is the tuple $(\phi, \widehat\phi)$. However,
in SQISign, everything is mapped back to $E_0$ or $\OO_0$ before converting,
so to allow for a clearer discussion we have simply assumed that our ideals 
$I$ have left order $\OO_0$ and our kernel generators are points $P \in E_0$. 

Essentially, the trick for working with generic curves is to take $P \in E$ 
and compute the image $\widehat\phi(P) \in E_0$ to use in the above function. 
Then, after we find $P_\alpha = \theta_\alpha \circ \widehat\phi(P) \in E_0$, we
return $\phi(P_\alpha) \in E$. 

**Note**: The composition $\phi \circ \widehat\phi = [D]$ where 
$D = \deg(\phi)$. So in the above, we have computed:

$$
[D] \circ \xi^{-1}(\alpha) = \phi \circ \theta_\alpha \circ \widehat\phi,
$$

So, to recover the true action of $\xi^{-1}(\alpha)(P)$, we must divide by $[D]$.


## Computing the kernels from ideals

With the ability to evaluate $\xi^{-1}(\alpha)(P)$ on a point $P \in E_0$, obtaining 
a kernel generator from an ideal is now straightforward. For the ease of notation, we
will write $\alpha(P)$ as the action of some $\alpha \in \BB$ on a point $P \in E_0$,
dropping the explicit mention of the isomorphism $\xi$.


For an ideal $I$, let us denote the $I$-torsion subgroup 

$$
E[I] = \\{ P \in E \mathrel{|} \alpha(P) = 0 \\;\\; \forall \alpha \in I \\} ,
$$

which should be understood as the kernel $E_0[I] = \ker(\phi_I)$ for the isogeny
 $\phi_I : E_0 \to E_0 / \langle I \rangle$. For an ideal of norm $D$, the group $E_0[I]$ will have
 order $D$.

We will only be considering cyclic ideals as these correspond to cyclic isogenies.
For clarification, see the page [Working with Cyclic Ideals](/posts/cyclic-ideals).
When $I$ is a cyclic ideal, we have that
$E_0[I]$ is a cyclic group and there is some generator point $K \in E_0[D]$
such that $E_0[I] = \langle K \rangle $.

Additionally, for any cyclic ideal with left order $\OO_0$, we can always find an 
element $\alpha \in I$ with reduced norm $\gcd(\Nrd(\alpha), D^2) = D$. 
We call this element the *generator* of the ideal, and we write:

$$
I = \OO_0\langle \alpha, D \rangle = \OO_0 \alpha + \OO_0 D.
$$

Given the generator $\alpha$, we can then recover $E_0[I] = E_0[\alpha] \cap E_0[D]$ 
by first computing the torsion basis $E_0[D] = \langle P, Q \rangle$ and then computing the kernel generator $K = [a]P + [b]Q \in E_0[D]$.

This is achieved by solving for $a,b$ such that the action of $\alpha(K)$ 
gives the identity point on $E_0$.
We can efficiently solve this by exploiting the fact that $\xi^{-1} (\alpha)$ is an endomorphism, and
hence a group homomorphism, to reduce the problem to finding $a,b$ such that $[a]\alpha(P) = -[b]\alpha(Q)$. We do this by
solving a discrete logarithm problem, which is efficient as our points all have smooth order.

### Prime power norms

When $D$ is a prime power, we have that at least one of $\alpha(P)$ or $\alpha(Q)$
will have order $D$. Without loss of generality, assume $\alpha(Q)$ has order $D$.
We can then use $\alpha(Q)$ as the base point for solving the discrete log and find
the integer $a$ such that $\alpha(P) = [a] \alpha(Q)$. The kernel generator 
$K = P - [a]Q$ can then be returned.

We can solve this in SageMath with the following snippet:

```python
αP, αQ = [eval_endomorphism(α, X, D) for X in (P, Q)]
if not has_order_D(αQ, D):
    αP, αQ = αQ, αP
    P, Q = Q, P
a = αQ.discrete_log(αP)
assert αP == a*αQ
return P - a*Q
```

### Allowing $D$ to be composite

Although working with prime powers allows an easier discussion, by allowing composite $D$
we can perform this algorithm only once, rather than for all primes factoring $D$.
However, to do this, we must first slightly modify
the discrete logarithm steps. 

For generic $D$, we might have that neither 
$\alpha(P)$ nor $\alpha(Q)$ have order $D$, so we cannot solve the discrete logarithm
as above to recover $a$.
Instead, we compute $\alpha(P), \alpha(Q)$ as before and write them in linear combinations
of the torsion basis $E_0[D]$:
$$
\alpha(P) = [a_1] P + [b_1] Q, \qquad \alpha(Q) = [a_2] P + [b_2] Q.
$$
The integers $a_i, b_i$ can be efficiently recovered by solving two bi-dimensional 
discrete logarithm problem. In out implementation, this is the function `BiDLP(R, P, Q, D)`
in the [`utilities.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/utilities.py) file.

Given the integers $a_i, b_i$, we can recover $E_0[I]$ by first computing the kernel basis 
of the matrix:

$$
\mathbf{A} = 
\begin{bmatrix}
a_1 & b_1 \\\
a_2 & b_2 
\end{bmatrix} \in \mathbb{M}_2(\mathbb{Z} / D \mathbb{Z}), \qquad 
\ker(\mathbf{A}) = 
\begin{bmatrix}
u_1 & v_1 \\\
u_2 & v_2 
\end{bmatrix}
$$

and from this, we recover 

$$
E_0[I] = \langle [u_1] P + [v_1] Q,  [u_2] P + [v_2] Q\rangle.
$$

We can perform this computation in SageMath with the following snippet:

```python
# Evaluate α(P) and α(Q)
R, S = [EvalEndomorphism(α, X) for X in (P,Q)]

# Find x,y such that R = xP + yQ
Rx, Ry = BiDLP(D, R, P, Q)

# Find x,y such that R = xP + yQ
Sx, Sy = BiDLP(D, S, P, Q)

# matrix of α on D-torsion w.r.t. basis (P,Q)
mat = matrix([(Rx, Ry), (Sx, Sy)])

# Thanks to Lorenz Panny for the following Sage
# Hack to find kernel of a matrix over ZZ/DZZ
# See: https://trac.sagemath.org/ticket/34862
ker = mat.stack(D*identity_matrix(2)).left_kernel().basis_matrix()[:,:2].change_ring(Zmod(D))

return [u*P + v*Q for u,v in ker if u or v]
```

For our implementation, we want to not just have the kernel 
$E_0[I] = \langle K_1, K_2\rangle$ but rather a single generating point $K$. 
Although there are deterministic ways to solve for this, it seems to require computing 
the order of both $K_1$ and $K_2$, which for large characteristic can be slow.

Instead, we find that a fairly sensible method is to just try linear combinations of $K_i$ until
a point of order $D$ is found.

If we have that either $K_1$ or $K_2$ has order exactly $D$, we can return this as our kernel generator
immediately. In the case where neither has full order, we can construct $K = K_1 + [x] K_2$ and see whether
this has full order. Checking if an element has order exactly $D$ is about ten times faster than computing
its order, so this check is almost always faster.

We implemented this as a small function, which we copy below:

```python
def derive_cyclic_generator(P, Q, D):
    """
    Given generators <P,Q> of a cyclic group
    of order D, find K such that G = <K>
    
    Heuristically, it seems easy to randomly
    find a K this way, and is about 10x faster
    than the deterministic method as we do not
    need to compute the order of P or Q.
    """
    K = P + Q
    for _ in range(1000):
        if has_order_D(K, D):
            return K
        K += Q
    raise ValueError(f"Never found a cyclic generator!")
```

We acknowledge there may be a significantly better way to do this, but this works for now.

### Avoiding discrete logs

The description above is how our code was implemented while developing this project. However, there was a recent publication:
[Deuring for the People: Supersingular Elliptic Curves with Prescribed Endomorphism Ring in General Characteristic](https://ia.cr/2023/106) 
by Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni, which had some tricks to remove the discrete logarithms.

The authors noticed that instead of computing the action of $\alpha$ on the generators of $E_0[D] = \langle P, Q \rangle$
and then solving a discrete logarithm to recover the kernel, one could directly recover 
$E_0[I] = \langle \bar\alpha(P), \bar\alpha(Q) \rangle$
by computing the action of the *conjugate* of the generator.
For more details, please see Section 4.1 of the above reference.

We now have a very simply algorithm for composite $D$. Given an ideal $I$ of norm $D$:

- Compute the generator of the ideal $\alpha$ and take the conjugate $\bar\alpha$
- Compute the torsion basis $E_0[D] = \langle P, Q \rangle$.
- Evaluate the action of $\bar\alpha$ on this basis to obtain $\langle \bar\alpha(P), \bar\alpha(Q) \rangle$.
- (Optionally) take a linear combination of $\langle \bar\alpha(P), \bar\alpha(Q) \rangle$ to find a point $K$ of order $D$.

### Putting it all together

Putting everything together, we get the following clean implementation

```py
def ideal_to_kernel(E, I):
    """
    Given a supersingular elliptic curve E0
    and a cyclic ideal I with left order O0 
    produces a generator K_I ∈ E[n(I)] such that
    ker(ϕ) = ⟨K_I⟩ = ⟨α_bar(P), α_bar(Q)⟩
    """
    assert is_cyclic(I), "Input ideal is not cyclic"

    # Degree of the isogeny we will to compute
    D = ZZ(I.norm())

    # Compute a generator such that I = O<α, D>
    α = ideal_generator(I)

    # Compute the torsion basis of E[D]
    P, Q = torsion_basis(E, D)

    # Evaluate R = α_bar(P)
    α_bar = α.conjugate()
    R = eval_endomorphism(α_bar, P, D)
    
    # If this has full order, we can stop here
    # as R = α_bar(P) generates the kernel
    if has_order_D(R, D):
        return R

    # Same again for S = α_bar(Q)
    S = eval_endomorphism(α_bar, Q, D)
    if has_order_D(S, D):
        return S

    # Neither R or S had full order, so we to find a
    # linear combination of R, S which has order D
    return derive_cyclic_generator(R, S, D)
```



## Computing ideals from kernels

Unlike `ideal_to_kernel()`, which is used throughout SQISign, converting from
a kernel generator to an ideal is only performed during the commitment stage.
In this final section, we
overview how this algorithm is implemented. We follow 
[Leroux's Thesis: Algorithm 20](https://www.lix.polytechnique.fr/Labo/Antonin.LEROUX/manuscrit_these.pdf) 
without modification, so the summary is shorter than the section above.

We know some point $K \in E[D]$ which generates the kernel of a 
cyclic isogeny $\phi$ of degree $D$, and we wish to compute a generator $\alpha$ 
such that $I_\phi = \OO \alpha + \OO D$ for $\OO \cong \End(E)$.

To do this, we start with a basis for $\OO$ which we denote 
$\langle \beta_0, \beta_1, \beta_2, \beta_3 \rangle$,
and we know that our generator can be written as some integral linear 
combination of the basis:
$\alpha = a_0 \beta_0 +  a_1 \beta_1 + a_2 \beta_2 + a_3 \beta_3$ for 
$a_i \in \ZZ$.

As $\ker(\phi) = \langle K \rangle$, if $\alpha$ the generator of the ideal $I_\phi$, 
we must have that $\alpha(K)$ is the identity point on $E$. This is what allows us to
determine the coefficients $a_i$.
Using our isomorphism $\xi$, we map the basis for $\OO$ to a basis of $\End(E)$ such that 
$\xi(\beta_i) = \theta_i$. We then map the above $\alpha$ to some $\theta \in \End(E)$:

$$
\theta = [a_0] \theta_0 + [a_1] \theta_1 + [a_2] \theta_2 + [a_3] \theta_3.
$$
Then, we look for $a_i$ such that $\theta(K)$ is the identity point on $E$.

The method described in Algorithm 20 is to find two basis elements
$\theta_i$ and $\theta_j$ which through their action on $K$ generate the torsion basis, i.e.,
$E[D] = \langle \theta_i(K), \theta_j(K) \rangle$.

Given this basis, we can compute any point $P \in E[D]$ as a linear combination of
$\theta_i(K)$ and $\theta_j(K)$. In particular, for some $\theta_k$ not equal to $\theta_i$ 
or $\theta_j$, we can represent the point $\theta_k(K) = [a] \theta_i(K) + [b] \theta_j(K)$ 
by solving the bi-dimensional discrete logarithm problem.

As we have $\theta_k(K) = [a] \theta_i(K) + [b] \theta_j(K)$, we have derived an endomorphism

$$
\theta = \theta_k - [a] \theta_i - [b] \theta_j,
$$

which maps $K$ to the identity point on $E$. This is precisely what we
were looking for, and we can finish our algorithm by computing

$$
\alpha = \beta_k - a \beta_i - b \beta_j, \quad I_\phi = \OO\langle \alpha, D \rangle.
$$

### Putting it all together

Putting this together into one function, we have the following implementation:

```py
def kernel_to_ideal(P, D):
    """
    Given a point P ∈ E[D] compute the
    ideal I(<P>))
    """
    # Compute a basis β1,β2,β3,β4 of O0
    # with norm coprime to D
    βs = compute_coprime_basis(D)

    # Compute the image of all the points
    # β(P) by acting with θ ≅ β
    θs = [eval_endomorphism(β, P, D) for β in βs]

    # Find θi, θj which generates E[D]
    i, j = find_torsion_basis_EndE(θs, D)
    θi, θj = θs[i], θs[j]

    # Pick k ≠ i,j such that
    k = set([0, 1, 2, 3]).difference([i, j]).pop()
    θk = θs[k]

    # Solve the discrete log
    a, b = BiDLP(θk, θi, θj, D)
    assert a * θi + b * θj == θk

    # Create the Quaternion Algebra element
    α = βs[k] - a * βs[i] - b * βs[j]
    return O0 * α + O0 * D
```

### Finding the torsion basis

To find $E[D] = \langle \theta_i(K), \theta_j(K) \rangle$ it is enough to just
check all pairs $\theta_i(K)$, $\theta_j(K)$ by computing the Weil pairing
$e(\theta_i(K), \theta_j(K))$. If the multiplicative order of this element 
is $D$, then $\theta_i(K)$, $\theta_j(K)$ must be linearly independent.

```py
def find_torsion_basis_EndE(θPs, D):
    """
    Looks for θi, θj such that θi(P), θj(P) generates E[D]
    """
    for i in range(4):
        for j in range(i + 1, 4):
            eθiθj = θPs[i].weil_pairing(θPs[j], D, algorithm="pari")
            if has_order_D(eθiθj, D, multiplicative=True):
                return i, j
    raise ValueError(f"No basis for E[D] found with given point")
```

### Computing a good basis

Another thing we must be careful of is making sure we compute a basis for $\OO$ 
such that the reduced norm of all the basis elements are coprime to $D$.
This is done to ensure that we find the correct generator $\alpha$ and not something
which only corresponds to some subgroup of the kernel.

Currently, we find a *good* basis with the following (fairly hacky) function:

```py
def compute_coprime_basis(D):
    """
    Start with basis <1, i, (i + j) / 2, (1 + k) / 2>
    and find a new basis such that the norm of each basis
    element is coprime to `D`.
    """
    O0_basis = O0.basis()
    θs = []
    for f in O0_basis:
        while True:
            if gcd(f.reduced_norm(), D) == 1:
                θs.append(f)
                break
            f += O0_basis[0] + O0_basis[1]
    return θs
```

Potentially there is much nicer solution than this. If you have a suggestion, 
please let us know! For now, this is working, so hopefully that is good enough.


[Back to Top](#top)





