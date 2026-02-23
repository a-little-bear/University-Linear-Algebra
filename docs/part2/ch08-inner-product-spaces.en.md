# Chapter 08: Inner Product Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Orthogonality (Ch7) · Fields (Ch00)

**Chapter Outline**: Definition of Inner Products → Axioms → Norms Induced by Inner Products → Cauchy-Schwarz Inequality → Triangle Inequality → Angles between Vectors → Orthogonal Complements → Adjoint of a Linear Operator → Self-Adjoint and Unitary Operators → Projection Theorem

**Extension**: Inner product spaces add "geometry" (lengths and angles) to the "topology" of vector spaces; they are the starting point for Hilbert space theory in analysis.

</div>

While any vector space allows for addition and scaling, an **inner product space** adds the ability to measure lengths and angles. By defining a scalar-valued function $\langle u, v \rangle$, we can determine if vectors are "close" or "perpendicular." This geometric layer is essential for optimization and functional analysis. The **Cauchy-Schwarz Inequality** is the most fundamental result in this area, providing the bound that ensures the consistency of our definitions of length and angle across all dimensions.

---

## 08.1 Definitions and Inequalities

!!! definition "Definition 08.1 (Inner Product)"
    An inner product on $V$ is a map $\langle \cdot, \cdot \rangle: V \times V \to F$ such that:
    1. **Positivity**: $\langle v, v \rangle \ge 0$ and $\langle v, v \rangle = 0 \iff v = 0$
    2. **Linearity**: $\langle au + bv, w \rangle = a \langle u, w \rangle + b \langle v, w \rangle$
    3. **Conjugate Symmetry**: $\langle u, v \rangle = \overline{\langle v, u \rangle}$

!!! theorem "Theorem 08.1 (Cauchy-Schwarz Inequality)"
    For any $u, v \in V$:
    $$|\langle u, v \rangle| \le \|u\| \|v\|$$
    Equality holds if and only if $u$ and $v$ are linearly dependent.

---

## Exercises

1. **[Fundamentals] Verify that the standard dot product in $\mathbb{R}^n$ is an inner product.**
   ??? success "Solution"
       $\langle x, y \rangle = \sum x_i y_i$ satisfies linearity, symmetry (in $\mathbb{R}$), and positivity ($\sum x_i^2 \ge 0$). Thus it is an inner product.

2. **[Norm] Define the norm induced by an inner product.**
   ??? success "Solution"
       $\|v\| = \sqrt{\langle v, v \rangle}$. This generalizes the Euclidean notion of distance.

3. **[Triangle Inequality] Use Cauchy-Schwarz to prove the Triangle Inequality: $\|u+v\| \le \|u\| + \|v\|$.**
   ??? success "Solution"
       $\|u+v\|^2 = \langle u+v, u+v \rangle = \|u\|^2 + 2\operatorname{Re}\langle u, v \rangle + \|v\|^2 \le \|u\|^2 + 2\|u\|\|v\| + \|v\|^2 = (\|u\| + \|v\|)^2$.

4. **[Angle] Define the angle $\theta$ between two non-zero vectors in a real inner product space.**
   ??? success "Solution"
       $\cos \theta = \frac{\langle u, v \rangle}{\|u\| \|v\|}$. The Cauchy-Schwarz inequality ensures that this ratio is always in $[-1, 1]$.

5. **[Orthogonal Complement] If $W$ is a subspace, define $W^\perp$.**
   ??? success "Solution"
       $W^\perp = \{ v \in V : \langle v, w \rangle = 0, \; \forall w \in W \}$. It consists of all vectors perpendicular to the entire subspace.

6. **[Adjoint] What is the adjoint $T^*$ of a linear operator $T$?**
   ??? success "Solution"
       The unique operator $T^*$ such that $\langle Tu, v \rangle = \langle u, T^*v \rangle$ for all $u, v$. For matrices, $T^* = \overline{T}^T$.

7. **[Self-Adjoint] Why are self-adjoint operators ($T^* = T$) important?**
   ??? success "Solution"
       They have real eigenvalues and their eigenvectors form an orthonormal basis (Spectral Theorem). They correspond to symmetric or Hermitian matrices.

8. **[Unitary] Define a unitary operator.**
   ??? success "Solution"
       An operator $U$ such that $U^* U = I$, or equivalently $\langle Uu, Uv \rangle = \langle u, v \rangle$. These operators preserve distances and angles (isometries).

9. **[Projection] State the Projection Theorem for a subspace $W \subseteq V$.**
   ??? success "Solution"
       Every $v \in V$ has a unique decomposition $v = w + w^\perp$ where $w \in W$ and $w^\perp \in W^\perp$.

10. **[Complex Numbers] How does conjugate symmetry differ from standard symmetry?**
    ??? success "Solution"
        In complex spaces, $\langle u, v \rangle$ is the complex conjugate of $\langle v, u \rangle$. This ensures that $\langle v, v \rangle$ is always a real number, allowing for a well-defined norm.

## Chapter Summary

This chapter provides the geometric foundation for vector spaces:

1. **Metric Definition**: Defined inner products as the source of length, distance, and angle in abstract spaces.
2. **Fundamental Bounds**: Established the Cauchy-Schwarz inequality as the universal constraint on vector interactions.
3. **Operator Duality**: Developed the theory of adjoints, providing a link between an operator and its "mirror" image in the inner product.
4. **Geometric Decomposition**: Leveraged orthogonal complements and projections to solve best-approximation problems in high dimensions.
