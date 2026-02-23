# Chapter 13A: Quotient and Dual Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Linear Transformations (Ch05) · Inner Product Spaces (Ch08)

**Chapter Outline**: From Internal Structure to External Mappings → Definition of Quotient Space $V/W$ → Cosets and Linear Operations → Dimension Formula for Quotient Spaces → Induced Mappings → Dual Space $V^*$ → Linear Functionals → Dual Bases → Annihilators ($W^0$) → Transpose Mappings and Dual Operators → Double Dual Space and Natural Isomorphism → Applications: Introduction to Tensor Analysis, Covariant vs. Contravariant Vectors in Physics

**Extension**: Quotient spaces are tools for "ignoring" information within subspaces, key to reducing dimensionality via invariant subspaces; dual spaces are the essential language for modern geometry, physics (General Relativity), and Functional Analysis.

</div>

Linear algebra studies not only the internal construction of spaces but also the "projections" and "mappings" between them. Quotient spaces allow us to simplify our vision by ignoring specific subspaces (like ignoring phase when studying the amplitude of a wave); dual spaces treat mappings themselves as vectors, providing the mathematical vehicle for measurement, observation, and tensor calculus.

---

## 13A.1 Quotient Space $V/W$

!!! definition "Definition 13A.1 (Cosets and Quotient Spaces)"
    Let $W$ be a subspace of $V$. For any $v \in V$, the **coset** of $v$ relative to $W$ is defined as:
    $$v + W = \{v + w : w \in W\}$$
    The set of all cosets forms the **Quotient Space**, denoted $V/W$.
    - **Addition**: $(u+W) + (v+W) = (u+v) + W$
    - **Scalar Multiplication**: $c(v+W) = (cv) + W$

!!! theorem "Theorem 13A.1 (Dimension of Quotient Spaces)"
    If $V$ is finite-dimensional, then:
    $$\dim(V/W) = \dim V - \dim W$$
    This reflects that the quotient space is essentially the co-dimensional space obtained by "subtracting" the degrees of freedom of the subspace.

---

## 13A.2 Dual Space $V^*$

!!! definition "Definition 13A.2 (Linear Functionals and Dual Spaces)"
    A linear mapping from $V$ to its scalar field $F$ is called a **linear functional**.
    The vector space of all linear functionals on $V$ is the **Dual Space**, denoted $V^*$.

!!! theorem "Theorem 13A.2 (Dual Basis)"
    Let $\{e_1, \ldots, e_n\}$ be a basis for $V$. The corresponding **dual basis** $\{e^1, \ldots, e^n\} \subset V^*$ is defined by:
    $$e^i(e_j) = \delta^i_j = \begin{cases} 1 & i=j \\ 0 & i \neq j \end{cases}$$
    In the finite-dimensional case, $\dim V^* = \dim V$.

---

## 13A.3 Annihilators ($W^0$)

!!! definition "Definition 13A.3 (Annihilator)"
    The **annihilator** of a subspace $W \subseteq V$ is the set of functionals in $V^*$ that vanish on $W$:
    $$W^0 = \{ f \in V^* : f(w) = 0, \forall w \in W \}$$
    **Property**: $\dim W^0 = \dim(V/W) = \dim V - \dim W$.

---

## Exercises

**1. [Quotient Space] In $\mathbb{R}^3$, if $W = \operatorname{span}\{e_1, e_2\}$, describe $V/W$ geometrically.**

??? success "Solution"
    **Geometric Analysis:**
    1. $W$ is the $xy$-plane.
    2. A coset $v+W$ is the translation of the $xy$-plane along the vector $v$.
    3. Note: All vectors falling on the same line perpendicular to the $xy$-plane (i.e., parallel to the $z$-axis) have a difference that belongs to $W$, so they correspond to the same coset.
    **Conclusion**: $V/W$ can be viewed geometrically as the $z$-axis, or more accurately, the set of all planes parallel to the $xy$-plane. Its dimension is $3-2=1$.

**2. [Coset Logic] Prove that $u+W = v+W$ if and only if $u-v \in W$.**

??? success "Solution"
    **Proof:**
    1. If $u+W = v+W$, then $u \in v+W$ (since $u = u+0 \in u+W$).
    2. By definition of a coset, there exists some $w \in W$ such that $u = v + w$.
    3. Subtracting $v$ gives $u - v = w$.
    4. Since $w \in W$, it follows that $u - v \in W$.
    **Intuition**: Two vectors belong to the same coset iff their "difference" lies in the subspace being ignored.

**3. [Dual Basis] Find the dual basis $\{f_1, f_2\}$ for $B = \{(1, 1), (1, 0)\}$ in $\mathbb{R}^2$.**

??? success "Solution"
    **Calculation:**
    1. Let $f_1(x, y) = ax+by$.
       - $f_1(1, 1) = a+b = 1$
       - $f_1(1, 0) = a = 0$
       Solving gives $a=0, b=1$. Thus $f_1(x, y) = y$.
    2. Let $f_2(x, y) = cx+dy$.
       - $f_2(1, 1) = c+d = 0$
       - $f_2(1, 0) = c = 1$
       Solving gives $c=1, d=-1$. Thus $f_2(x, y) = x-y$.
    **Conclusion**: The dual basis is $\{y, x-y\}$.

**4. [Annihilator] Prove that $W^0$ is a subspace of $V^*$.**

??? success "Solution"
    **Axiom Verification:**
    1. **Zero Element**: The zero functional $\mathbf{0}(v)=0$ clearly vanishes on $W$, so $\mathbf{0} \in W^0$.
    2. **Additive Closure**: If $f, g \in W^0$, then for any $w \in W$, $(f+g)(w) = f(w)+g(w) = 0+0=0$, so $f+g \in W^0$.
    3. **Scalar Closure**: If $f \in W^0$, then $(cf)(w) = c(f(w)) = c(0)=0$, so $cf \in W^0$.
    **Conclusion**: The annihilator satisfies the subspace axioms.

**5. [Transpose] How is the transpose $T^*: U^* \to V^*$ of a linear map $T: V \to U$ defined?**

??? success "Solution"
    **Definition:**
    For any functional $g \in U^*$, its image under $T^*$ is a functional in $V^*$ defined by:
    $$T^*(g) = g \circ T$$
    That is, $(T^*g)(v) = g(T(v))$ for $v \in V$.
    **Significance**: This shows how the dual operator "pulls back" observations from the output space to the input space.

**6. [Dimension] Prove $\dim(V/W) = \dim W^0$.**

??? success "Solution"
    **Proof:**
    1. We know $\dim(V/W) = \dim V - \dim W$.
    2. Consider the restriction map from $V^*$ to $W^*$. Its kernel is exactly $W^0$.
    3. By the rank-nullity theorem, $\dim V^* = \dim \operatorname{im} + \dim W^0$. Since the restriction is surjective, the dimension of its image is $\dim W^* = \dim W$.
    4. Therefore, $\dim W^0 = \dim V - \dim W$.
    **Conclusion**: The dimensions are equal.

**7. [Naturality] Why is $V \cong V^{**}$ natural, while $V \cong V^*$ is not?**

??? success "Solution"
    **Analysis:**
    1. **$V \to V^*$**: Requires an explicit choice of basis to define the mapping $e_i \to e^i$. Changing the basis changes the mapping.
    2. **$V \to V^{**}$**: There exists an evaluation map $\phi(v)(f) = f(v)$.
    3. The definition of this map **does not use a basis**. The relationship "vector acts on functional" is invariant under coordinate changes.
    **Conclusion**: Only isomorphisms independent of basis choice are termed "natural isomorphisms."

**8. [Functional] In $P_n$, is the map $f(p) = \int_0^1 p(x) dx$ a linear functional?**

??? success "Solution"
    **Verification:**
    1. The result is a scalar (real number).
    2. Integration is linear: $\int (ap+bq) = a\int p + b\int q$.
    **Conclusion**: Yes, it is a classic linear functional on the space of polynomials.

**9. [Natural Projection] Define the natural projection $\pi: V \to V/W$. What is its kernel?**

??? success "Solution"
    **Definition:**
    $\pi(v) = v + W$.
    **Finding the Kernel:**
    $\ker(\pi) = \{ v \in V : \pi(v) = 0_{V/W} \}$.
    The zero element of the quotient space is $0 + W = W$.
    By coset equality, $v + W = W \iff v \in W$.
    **Conclusion**: The kernel of the natural projection is exactly the subspace $W$.

**10. [Physics] In physics, to which space do covariant vectors usually belong?**

??? success "Solution"
    **Conclusion:**
    Covariant vectors (also called covectors) belong to the **dual space $V^*$**.
    **Reasoning**: Covariant components transform the same way as the basis vectors (accompanying the change of basis), which is the algebraic behavior of linear functionals under a change of coordinates. Standard position vectors (contravariant vectors) belong to the primal space $V$.

## Chapter Summary

Quotient and dual spaces elevate the abstract dimension of linear algebra:

1.  **Information Compression**: Quotient spaces modularize complex structures through equivalence classes, serving as the key to induced transformations and fundamental homomorphism theorems.
2.  **Algebra of Observation**: Dual spaces formalize the act of "measurement," establishing functionals as independent objects of study and laying the groundwork for index juggling and tensor calculus in physics.
3.  **Return to Naturality**: The natural isomorphism of the double dual reveals the deepest symmetry between mathematical objects and their observers, marking the transition from pure computation to categorical thinking.
