# Chapter 13A: Quotient and Dual Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Linear Transformations (Ch05)

**Chapter Outline**: Definition of Quotient Spaces → Cosets and Linear Operations → Dimension Formulas → Induced Mappings → Dual Space $V^*$ → Linear Functionals → Dual Bases → Annihilators ($W^0$) → Transpose of a Linear Map → Double Dual Space and Natural Isomorphism → Category-theoretic Perspective

**Extension**: Quotient spaces are tools for "ignoring" information within subspaces; Dual spaces are the essential language for establishing modern geometry, physics (Tensor Analysis), and Functional Analysis.

</div>

Linear algebra studies not only the internal construction of spaces but also the "projections" and "mappings" between them. Quotient spaces allow us to simplify our vision by ignoring specific subspaces, while dual spaces treat mappings themselves as vectors, providing the mathematical vehicle for measurement, observation, and tensor calculus.

---

## 13A.1 Quotient Space $V/W$

!!! definition "Definition 13A.1 (Cosets and Quotient Spaces)"
    Let $W$ be a subspace of $V$. For any $v \in V$, the **coset** of $v$ relative to $W$ is defined as $v + W = \{v + w : w \in W\}$.
    The set of all cosets forms the **Quotient Space**, denoted $V/W$.
    - **Addition**: $(u+W) + (v+W) = (u+v) + W$
    - **Scalar Multiplication**: $c(v+W) = cv + W$

!!! theorem "Theorem 13A.1 (Dimension of Quotient Spaces)"
    If $V$ is finite-dimensional, then $\dim(V/W) = \dim V - \dim W$.

---

## 13A.2 Dual Space $V^*$

!!! definition "Definition 13A.2 (Linear Functionals and Dual Spaces)"
    A linear mapping from $V$ to its scalar field $F$ is a **linear functional**. The vector space of all linear functionals on $V$ is the **Dual Space**, denoted $V^*$.

!!! theorem "Theorem 13A.2 (Dual Basis)"
    Let $\{e_1, \ldots, e_n\}$ be a basis for $V$. The corresponding **dual basis** $\{e^1, \ldots, e^n\}$ is defined by:
    $$e^i(e_j) = \delta_{ij}$$
    Consequently, for finite-dimensional spaces, $\dim V^* = \dim V$.

---

## 13A.3 Annihilators and Transpose Maps

!!! definition "Definition 13A.3 (Annihilator $W^0$)"
    The **annihilator** of a subspace $W$ is the set of all functionals in $V^*$ that vanish on $W$:
    $$W^0 = \{f \in V^* : f(w) = 0, \forall w \in W\}$$
    **Property**: $\dim W^0 = \dim V - \dim W$.

---

## 13A.4 Natural Isomorphism and the Double Dual

!!! theorem "Theorem 13A.3 (Natural Isomorphism to Double Dual)"
    For a finite-dimensional space $V$, there exists a **natural isomorphism** $\Phi: V \to V^{**}$ that is independent of the choice of basis, defined by:
    $$\Phi(v)(f) = f(v)$$
    This implies that vectors can be viewed as operators acting on functionals.

---

## Exercises

1. **[Quotient] In $\mathbb{R}^3$, if $W = \operatorname{span}\{e_1, e_2\}$, describe $V/W$ geometrically.**

   ??? success "Solution"
       $V/W$ can be viewed as the set of lines perpendicular to the $xy$-plane (each coset is a vertical line). Its dimension is $3-2=1$.

2. **[Coset] Prove that $u+W = v+W$ if and only if $u-v \in W$.**

   ??? success "Solution"
       If $u+W = v+W$, then $u \in v+W$, so $u = v+w$ for some $w \in W$, meaning $u-v = w \in W$.

3. **[Dual Basis] Find the dual basis for $B = \{(1, 1), (1, 0)\}$ in $\mathbb{R}^2$.**

   ??? success "Solution"
       Let the dual basis be $\{f_1, f_2\}$.
       $f_1(1, 1)=1, f_1(1, 0)=0 \implies f_1(x, y) = y$.
       $f_2(1, 1)=0, f_2(1, 0)=1 \implies f_2(x, y) = x-y$.

4. **[Annihilator] Prove that $W^0$ is a subspace of $V^*$.**

   ??? success "Solution"
       If $f, g \in W^0$, then $(f+g)(w) = f(w)+g(w) = 0+0=0$, showing closure under addition. Scalar multiplication follows similarly.

5. **[Transpose] How is the transpose $T^*: U^* \to V^*$ of a linear map $T: V \to U$ defined?**

   ??? success "Solution"
       $T^*(f) = f \circ T$.

6. **[Dimension] Prove $\dim(V/W) = \dim W^0$.**

   ??? success "Solution"
       Both sides are equal to $\dim V - \dim W$.

7. **[Naturality] Why is $V \cong V^{**}$ natural, while $V \cong V^*$ is not?**

   ??? success "Solution"
       Constructing $V \to V^*$ requires an explicit choice of basis to define the mapping; the definition of $V \to V^{**}$ (evaluating a functional at a vector) is entirely basis-independent.

8. **[Functional] In $P_n$, is the map $f(p) = \int_0^1 p(x) dx$ a linear functional?**

   ??? success "Solution"
       Yes. Integration is a linear operation and the result is a scalar.

9. **[Projection] Define the natural projection $\pi: V \to V/W$. What is its kernel?**

   ??? success "Solution"
       $\pi(v) = v + W$. Its kernel $\ker(\pi) = W$.

10. **[Physics] In physics, which space typically corresponds to covariant vectors?**

   ??? success "Solution"
        The dual space $V^*$.

## Chapter Summary

Quotient and dual spaces elevate the abstract dimension of linear algebra:

1.  **Information Compression**: Quotient spaces modularize complex structures through equivalence classes, serving as the key to induced transformations and fundamental homomorphism theorems.
2.  **Algebra of Observation**: Dual spaces formalize the act of "measurement," establishing functionals as independent objects of study and laying the groundwork for index juggling and tensor calculus in physics.
3.  **Return to Naturality**: The natural isomorphism of the double dual reveals the deepest symmetry between mathematical objects and their observers, marking the transition from pure computation to categorical thinking.
