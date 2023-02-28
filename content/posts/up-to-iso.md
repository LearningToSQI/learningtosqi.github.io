---
title: "Correct up to Isomorphism"
date: 2023-01-29T08:00:00Z
draft: false
mathjax: true
---

Informally speaking, the Deuring correspondence gives a bijection between maximal orders in the quaternion algebra $\BB$ and the ideals connecting them, to endomorphism rings of supersingular curves over $\mathbb{F}\_{p^2}$ and the isogenies between them. 
However, when computing the isogeny from a given ideal under this correspondence, it will only be correct *up to isomorphism*. Conversely, the correspondence yields the correct ideal *up to equivalence*. In the context of isogeny-based cryptography, this does not pose a problem as we are only interested in curves up to isomorphism. However, when working with these objects in practice, it can become an issue. 

To demonstrate this, consider the following situation. Say we have ideals $I$ and $J$ obtained using `kernel_to_ideal()`, corresponding to isogenies $\varphi_I$ and $\varphi_J$, and we would like to compute $IJ$. For the multiplication to be well-defined, we must have that the orders are compatible in the sense that $O_R(I) = O_L(J)$. In particular, these orders must be the *same*, not just isomorphic. However, by the Deuring correspondence, `kernel_to_ideal()` will only return the correct $I$ and $J$ up to equivalence, and as a result we can only be certain that $\OO_R(I) \cong \OO_L(J)$. 

The aim of this blogpost is to highlight generic fixes that we used throughout our code to obtain the correct ideals and isogenies.


### Obtaining the correct ideal

Suppose that we are given two equivalent ideals $I \sim J$. This means that there exists an $\alpha \in \BB^{\times}$ such that $J = I\alpha$. On input $I, J$, our algorithm `left_isomorphism()` outputs such an $\alpha$.

```py
def left_isomorphism(I, J):
    """
    Given two isomorphic left ideals I, J computes
    α such that J = I*α
    """
    B = I.quaternion_algebra()

    if B != J.quaternion_algebra():
        raise ValueError("Arguments must be ideals in the same algebra.")

    if I.left_order() != J.left_order():
        raise ValueError("Arguments must have the same left order.")

    IJ = I.conjugate() * J
    L = reduced_basis(IJ)
    for t in L:
        α = t / I.norm()
        if J == I * α:
            return α

    raise ValueError("Could not find a left isomorphism...")
```

The correctness of `left_isomorphism()` follows from the proof of Lemma 1 in the [original SQISign paper](https://eprint.iacr.org/2020/1240). Indeed, 
if $I$ and $J$ have the same left order, we can identify the ideal $\bar{I}J$ to the principal ideal $\OO_R(I)\beta$. Then, $J = \chi_I(\beta) = I\frac{\beta}{n(I)}$. 
Our code considers $\beta$ in the LLL reduced basis of $\bar{I}J$, so that $\beta$ has reduced norm $n(I)n(J)$, and $I\frac{\beta}{n(I)}$ has norm $n(J)$. 

Now let us turn our attention to multiplying ideals $I$ and $J$ using our algorithm `multiply_ideals()`, which will work if the orders of $I$ and $J$ are compatible. If we only have $\OO_R(I) \cong \OO_L(J)$ but we know a $K \sim I$ with $\OO_R(K) = \OO_L(J)$, we are still able to multiply the ideals $I$, $J$ by supplying extra information to `multiply_ideals()`. 

We first compute $\beta$ such that $I = K \beta$ by running `left_isomorphism(K, I)`, and then set $J' = \beta^{-1}J \beta$. We then compute $IJ$ as $IJ'$. The latter is computable as $I\beta^{-1}J \beta = KJ\beta$ and $K$ and $J$ have compatible orders. 

We note that, by definition, $n(J) = n(J')$. This is important as we usually require the norm of $J$ to have certain properties. For example, $J$ will usually be the output of `EquivalentSmoothIdealHeuristic()`, and we therefore want it to have smooth norm. 

```py
def multiply_ideals(I, J, beta=None):
    """
    Computes I*J when O_R(I) ≃ O_L(J)

    If these orders do not match, we compute
    as isomorphism which takes O_L(J) to O_R(I)
    """
    if I.right_order() != J.left_order():
        if beta is None:
            raise ValueError(
                "Right and left orders, do not match. Must supply an automorphism, beta"
            )

        J = beta ** (-1) * J * beta
        assert I.right_order() == J.left_order(), "Orders does still not match"
    return I * J
```

To see this in action, look at `IdealToIsogenyFromKLPT()` available in the file [`deuring.py`](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/deuring.py). If $I_m$ is the ideal at step $m$ of the filtration, we want to compute $JI_m$. However, to do this we must supply `multiply_ideals()` with the automorphism `left_isomorphism(K, J)`. 
 

### Obtaining the correct isogeny

Suppose we have computed an isogeny $\varphi: E_1 \rightarrow E$, however the codomain is only certain to be correct up to isomorphism. Namely, we *actually* need an isogeny $\varphi': E_1 \rightarrow E'$, where $E \cong E'$. Using SageMath's in-built functions, this is easy to fix: we first compute the isomorphism $\psi: E \rightarrow E'$ using `E.isomorphism_to(E')` and then set $\varphi' = \psi \circ \varphi$. 

If instead we are in a situation where the domain is incorrect, instead of post-composing, we simply *pre-compose* with the appropriate isomorphism.

[Back to Top](#top)


