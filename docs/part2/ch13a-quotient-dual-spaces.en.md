# Chapter 13A: Quotient Spaces and Dual Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Linear Transformations (Ch5)

**Chapter Outline**: Equivalence Relations and Partitioning → Definition of Quotient Space $V/W$ → Natural Projection Mapping → Fundamental Theorem of Homomorphism → Dual Space $V^*$ → Dual Basis → Double Dual and Natural Isomorphism → Abstract Definition of Transpose Transformations

**Extension**: Quotient spaces are tools for "dimensional reduction" by ignoring differences within a subspace; Dual spaces are core to functional analysis, revealing symmetry between vectors and observers (linear functionals).

</div>

This chapter pushes the study of vector spaces to a higher level of abstraction. **Quotient spaces** simplify structure by ignoring differences in certain dimensions, while **dual spaces** establish a parallel world composed of "measurement tools" (linear functionals).

---

## 13A.1 Quotient Spaces and Duality

!!! definition "Definition 13A.1 (Quotient Space)"
    Let $W$ be a subspace of $V$. The quotient space $V/W$ is the set of all cosets $v+W$, with dimension $\dim V - \dim W$.

!!! theorem "Theorem 13A.3 (Dual Basis)"
    If $V$ has a basis $\{e_1, \dots, e_n\}$, then $V^*$ has a unique basis $\{\phi_1, \dots, \phi_n\}$ satisfying $\phi_i(e_j) = \delta_{ij}$.

---

## Exercises

1. **[Coset Ops] In $V = \mathbb{R}^2$, let $W$ be the $x$-axis. Describe the geometric meaning of the quotient space $V/W$.**
   ??? success "Solution"
       Each coset $v+W$ represents a line parallel to the $x$-axis. The quotient space $V/W$ essentially compresses all these horizontal lines into points, creating a structure isomorphic to the $y$-axis.

2. **[Dimension] If $\dim V = 5$ and $\dim W = 2$, find $\dim(V/W)$.**
   ??? success "Solution"
       $\dim(V/W) = 5 - 2 = 3$.

3. **[Linear Functional] Verify that $\phi(x, y) = 2x + 3y$ is a linear functional on $\mathbb{R}^2$.**
   ??? success "Solution"
       $\phi(u+v) = 2(x_1+x_2) + 3(y_1+y_2) = \phi(u) + \phi(v)$.
       $\phi(cu) = 2(cx) + 3(cy) = c\phi(u)$.
       Linearity holds, so it is a linear functional.

4. **[Dual Basis Calculation] In $\mathbb{R}^2$, take the basis $e_1 = (1, 0), e_2 = (1, 1)$. Find its dual basis $\phi_1, \phi_2$.**
   ??? success "Solution"
       Let $\phi_1(x, y) = ax+by$. Since $\phi_1(1, 0)=1 \implies a=1$ and $\phi_1(1, 1)=0 \implies 1+b=0 \implies b=-1$. Thus $\phi_1(x, y) = x-y$.
       Similarly, $\phi_2(1, 0)=0 \implies a=0$ and $\phi_2(1, 1)=1 \implies b=1$. Thus $\phi_2(x, y) = y$.

5. **[Annihilator] Define the annihilator $W^0$ of a subspace $W$.**
   ??? success "Solution"
       $W^0 = \{ \phi \in V^* : \phi(w) = 0, \forall w \in W \}$. It is the subspace of all functionals in the dual space that "don't see" $W$.

6. **[Isomorphism Theorem] Prove $V/ \ker(T) \cong \operatorname{Im}(T)$.**
   ??? success "Solution"
       Define $\bar{T}(v+W) = T(v)$. It's easy to show this mapping is well-defined and is a bijective linear map. This shows that the image structure is determined by the remaining structure after quotienting out the kernel.

7. **[Double Dual] Why is a finite-dimensional space $V$ naturally isomorphic to $V^{**}$?**
   ??? success "Solution"
       Define $\psi: V \to V^{**}$ such that $\psi(v)(\phi) = \phi(v)$. Since this map is independent of the choice of basis and the dimensions match, it is a natural isomorphism.

8. **[Transpose] Let $T: V \to W$. How is its transpose $T^*: W^* \to V^*$ defined?**
   ??? success "Solution"
       For $\phi \in W^*$, $T^*(\phi) = \phi \circ T$. This reflects the pullback of observers from the output space to the input space.

9. **[Projection Mapping] Prove that the natural projection $\pi: V \to V/W$ is surjective and $\ker \pi = W$.**
   ??? success "Solution"
       For any $v+W \in V/W$, clearly $\pi(v) = v+W$, so it's surjective.
       $\pi(v) = 0+W \iff v \in W$, so the kernel is $W$.

10. **[Calculation] In $P_1$ (linear polynomials), the basis is $\{1, x\}$. Find the coordinates of the evaluation map $\phi(p) = p(1)$.**
    ??? success "Solution"
        $\phi(1) = 1, \phi(x) = 1$. Thus $\phi = 1\phi_1 + 1\phi_2$. The coordinates are $[1, 1]^T$.

## Chapter Summary

This chapter deepens the algebraic understanding of space structures:

1. **Structural Simplification**: Quotient spaces reveal residual features outside a subspace by "grouping similar terms."
2. **Observational Symmetry**: Dual spaces place "state" and "measurement" on equal footing.
3. **Categorical Links**: Isomorphism theorems and transposes unify spaces, subspaces, and mappings seamlessly.
