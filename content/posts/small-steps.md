---
title: "Small Steps from Curves with Small Endomorphisms"
date: 2023-01-24T14:10:44Z
draft: false
mathjax: true
---

This page is used to explain the optional argument `end_close_to_E0` available to `IdealToIsogenyFromKLPT()` and why it is necessary to include such that the algorithm can run successfully in SQISign for both `keygen()` and `response()`. For a more detailed description of the algorithm `IdealToIsogenyFromKLPT()`, read the page [Computing Isogenies from Ideals](/posts/ideal-to-isogeny).

## Motivation

We first compute the ideal filtration

$$
I = \tilde{I}_v \subset \dots \subset \tilde{I}_1 \subset \tilde{I}_0,
$$ 

where each quotient $I\_m = \tilde{I}\_{m}/\tilde{I}\_{m-1}$ has norm $\ell^{\text{step}}$. 

We have that the left 
$\OO\_0$-ideal $JI_{m}$ computed in Step 4 of `IdealToIsogenyFromKLPT()` is equivalent to 
$K_{m-1}I_m = K I_1 \cdots I_m$ of norm $n(K)\ell^{m\cdot \text{step}}$. 
When computing the ideal filtration, if $\ell^{\text{step}}$ does not divide $n(I)$, then one step will be smaller: 
$\ell^{\text{small}} = n(I) \bmod \ell^{\text{step}}$. The value of `end_close_to_E0` determines whether 
this small step should be taken first or last.

The problem of closeness all boils down to the fact that when we find the isogeny corresponding to $JI_m$, 
we use `IdealToIsogenySmallFromKLPT()`, within which we 
run `EquivalentPrimeIdealHeuristic()` via `EquivalentSmoothIdealHeuristic()` with $JI_m$ as input.
An important thing to note about this subroutine is that it only cares about the *equivalence* class of the 
input ideal. This means that even if $JI_m$ has large norm, if there is an equivalent ideal with small norm, 
the algorithm could still fail.

When running the equivalent prime norm ideal sub-algorithm of the KLPT algorithm, part of the heuristic success 
is that the input ideal is random. However, the curve $E_0$ and corresponding order $\OO_0$ are special in that 
the endomorphisms of $E_0$ are particularly small. For some of the inputs we need, the ideals are connected
to this special order have particularly small-norm and our heuristic assumptions are incorrect.

To illustrate the problem generally, imagine an ideal $I$ with left and right orders $\OO_0$ and $\OO$.
Then consider the isogeny $\phi : E_0 \to E$, such that $\OO \cong \End(E)$. When the curves $E_0$ and $E$ 
are connected by a small degree isogeny (i.e., they are close in the isogeny graph) then there will be an ideal $J \sim I$
with particularly small norm $n(J) = \deg(\phi)$. 
When this is the case, it is very unlikely for some $L \sim I \sim J$ with small prime norm to be found.
For more information on why this is a problem, see [Equivalent Prime Norm Ideals](/posts/prime-ideal/).

**Example**: During `keygen()` we compute $\tau' : E_0 \to E_A$ of $\ell$-power norm approximately of degree 
$\deg(\tau') \simeq p$. The ideal equivalent via the Deuring correspondence is $J_\tau$ with 
$n(J_\tau) = \deg(\tau')$. However, we know that there is also an equivalent ideal $I_\tau \sim J_\tau$ of norm
$n(I_\tau) \simeq p^{1/4}$, and so $E_0$ and $E_A$ are actually much closer than the norm of $J_\tau$ suggests.
Running `EquivalentPrimeIdealHeuristic()` with $J_\tau$ as input will compute the reduced basis corresponding 
to the much smaller $I_\tau$ and likely find no prime norm ideal $L$.

For SQISign to run successfully, we need to be able to work around the steps of the algorithm for which the input ideals
have particularly small norm by using clever tricks. The function `IdealToIsogenyFromKLPT()` is used in both key generation 
and in the response. To handle both of these use cases, we need slightly different tricks, and we control this with the 
optional `bool`, `end_close_to_E0` which determines the order of the step sizes for the filtration. We give more detail 
about this for the rest of the section.

## Working around key generation

During `keygen()`, we compute an isogeny $\tau' : E_0 \to E_A$ of $\ell$-power norm. Here we have that *both* ends of the
isogeny are close to $E_0$. The domain of $\tau'$ *is* $E_0$, which is as close as you could ever be and, as explained
above, $E_A$ is connected to $E_0$ by some unknown prime-degree isogeny $\ell \simeq p^{1/4}$. For the intermediate steps
of this algorithm, we will be walking away from $E_0$ and towards $E_A$ with steps of size $n(K) \ell^{m \cdot \text{step}}$.
Essentially, almost every step of this algorithm requires working with ideals which connect curves close to 
$E_0$. Our hope is to make the distance just big enough to make everything work.

The trick for `keygen()` is to counter-intuitively take the smallest possible step from $E_0$ as the first step. 
The idea is to avoid calling `EquivalentPrimeIdealHeuristic()` all together and instead work with input ideals 
of $\ell$ power norm. We can do this so long as the input ideal has norm smaller than $\sqrt{p}$.
For more information, see the page [Equivalent Prime Norm Ideals](/posts/prime-ideal/).

Let us discuss precisely how this is achieved in practice, which boils down to considering the size of the
current step and the norm of the connecting ideal $n(K)$, both of which are powers of $\ell$. 
On input $JI_m$, we can always make an ideal $M$ with norm a power of $\ell$ with the following
steps:

1. Find $\alpha$ such that $K = J\alpha$. Then, $\alpha$ will have reduced norm $\frac{n(K)}{n(J)}$. 
2. Set $M = JI_m \alpha$. 

By construction, $M$ is equivalent to $JI_m$ and will have norm $n(M) = \frac{n(JI_m)n(K)}{n(J)} = n(I_m)n(K) = \ell^\times$. 
Providing that $\ell^\times$ is small enough (about $\sqrt{p}$), we can use $M$ as input for KLPT without needing to
compute an equivalent prime norm ideal. If the norm of $M$ is too large, then the lattice solutions from the Strong Approximation
would become too large. Luckily the crossing-over point of $M$ being too large is exactly when we can rely comfortably on 
`EquivalentPrimeIdealHeuristic()` to work as expected for random ideals.

Concretely, in key generation, we input to our algorithm the ideal $J_{\tau}$ to compute $\tau'$. 
As $J_{\tau}$ is already a left $\OO\_0$-order, the first ideal in the filtration will already be a left $\OO\_0$-ideal 
and the connecting ideal $K$ is set as unit ideal in $\OO\_0$ (so $n(K) = 1$) with corresponding trivial isogeny 
$\varphi_K : E_0 \to E_0$.

The fact that $n(K) = 1$ for the first step is what allows us to use the shortest step first. On the first call, we have 
$n(M) = \ell^{\text{small}}$ which will be very small (for SQISign $\ell^{\text{step}} \simeq p^{1/4}$) and so we use this as
our prime power norm ideal $L$ in the KLPT algorithm. 

For the second step, we will have $K$ such that $n(K) = \ell^{\text{small}}$ and the step of $JI_1$ will be a full step of size
$\ell^{\text{step}}$.
Applying the same trick above, we will have that $n(M) = \ell^{\text{step} + \text{small}}$. Again, this will be small enough
we can use this in place of deriving an equivalent prime norm ideal.

For the third step we are $\ell^{2\cdot \text{step} + \text{small}}$ from $E_0$ and $\ell^\text{step}$ from $E_A$. Remember that as $E_A$ 
is itself close to $E_0$, if we had saved the smallest step for now, we be $\ell^{3\cdot \text{step}}$ from $E_0$ but only 
$\ell^{\text{small}}$ from $E_A$. The closeness to $E_A$ (and hence $E_0$) would mean we are likely to fail. As the norm of 
$K$ grows with each step, we cannot use the trick described above to compute a good $M$ this far into the chain.
Having used the small step early and having this step as large as possible is what saves us during `keygen()`.

For the final step, we input an ideal equivalent to $J_\tau$ itself and it is hopeless to try and use the above trick. Now, $n(K)$
is too big and the curve $E_A$ is too close to $E_0$ to use `EquivalentPrimeIdealHeuristic()`. 
Luckily, we already know an equivalent prime norm ideal: $I_\tau$. We pass this ideal
as an optional parameter and use it in the last step of the function and do not call `EquivalentPrimeIdealHeuristic()` at all.

We see that ultimately three of the four steps for `keygen()` fall apart without proper handling of edge cases. Luckily, 
there are ways to work around each of them. 

## Responding with another trick

For the response isogeny $\sigma : E_A \to E_2$, we only have $E_A$ as being close to $E_0$. However, unlike `keygen()`, we will not have
$n(K) = 1$ as input for the first step, so we cannot hope to use the above trick to do the first steps using $\ell$-power tricks.
Here, the additional input given in this case is $K = J_{\tau} \sim I_{\tau}$, an ideal of norm $\approx p^{1/4}$. At each step of the filtration, 
our ideals will be equivalent to a left $\OO_0$-ideal of norm $N_{\tau}\ell^{m\cdot \text{step}}$

The goal is to make sure that the initial step is large enough for `EquivalentPrimeIdealHeuristic()` run successfully, 
so we set `end_close_to_E0` to be `False` and have the first step be as big as possible. This puts the potentially very 
small last step at the end, which is no problem, as $E_2$ is far enough from $E_0$ that we do not expect for the corresponding 
ideal to have small norm.

There is however a slightly different problem. Before running the ideal chain for `response()`, we call `IdealToIsogenyCoprime()` with $J_\tau$
as input and pass this through to the KLPT algorithm. We know that this is a problematic ideal and so we cannot run this as intended, but 
luckily again we have $I_\tau$. By passing this in as the equivalent prime norm ideal for $J_\tau$, we can skip running 
`EquivalentPrimeIdealHeuristic()` during the call to the KLPT algorithm within `IdealToIsogenyCoprime()`.

[Back to Top](#top)
