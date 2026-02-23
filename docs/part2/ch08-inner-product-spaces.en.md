# Chapter 08: Inner Product Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Orthogonality (Ch7)

**Chapter Outline**: Definition of Abstract Inner Product → Norms and Metrics → Cauchy-Schwarz Inequality → Angle and Orthogonality → Completeness and Hilbert Spaces Introduction → Adjoint Operators $T^*$ → Normal and Unitary Operators → Spectral Theorem

**Extension**: Inner product spaces enhance vector spaces by introducing continuity and geometry; they are finite-dimensional precursors to functional analysis.

</div>

An inner product space is a vector space reinforced with an inner product. By defining this operation, we can talk not just about linear combinations, but also about lengths, angles, and approximations. It is the mathematical foundation for modern signal processing, quantum mechanics, and optimization theory.

---

## 08.1 Definitions and Core Properties

!!! definition "Definition 08.1 (Inner Product)"
    A vector space $V$ over field $F$ is equipped with an inner product $\langle \cdot, \cdot \rangle$ if it satisfies:
    1. **Positive Definiteness**: $\langle v, v \rangle \ge 0$, and equals 0 iff $v=0$.
    2. **Conjugate Symmetry**: $\langle u, v \rangle = \overline{\langle v, u \rangle}$.
    3. **Additivity/Homogeneity**: $\langle au+bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle$.

---

## Exercises

1. **[Axioms] Verify if $\langle u, v \rangle = u_1 v_1 + 2u_2 v_2$ on $\mathbb{R}^2$ constitutes an inner product.**
   ??? success "Solution"
       Yes. It is a weighted version of the standard dot product. Since the coefficients 1 and 2 are positive, it satisfies positive definiteness. Linearity and symmetry follow trivially. This is called a weighted inner product.

2. **[Cauchy-Schwarz] Prove: $|\langle u, v \rangle| \le \|u\| \|v\|$.**
   ??? success "Solution"
       Consider $f(t) = \|u - tv\|^2 = \langle u-tv, u-tv \rangle \ge 0$. Expanding gives $t^2 \|v\|^2 - 2t \operatorname{Re}\langle u, v \rangle + \|u\|^2 \ge 0$. For this quadratic to be non-negative, its discriminant must be $\le 0$, which yields the Cauchy-Schwarz inequality.

3. **[Adjoint] Let $T$ have matrix $A$ in the standard basis. Prove its adjoint $T^*$ has matrix $A^*$.**
   ??? success "Solution"
       By definition: $\langle Tv, w \rangle = \langle v, T^* w \rangle$.
       Left side $= (Av)^* w = v^* A^* w$.
       Right side $= v^* (T^* w)$.
       Since this holds for all $v, w$, the matrix for $T^*$ must be $A^*$ (conjugate transpose).

4. **[Unitary] Prove that unitary operators preserve the inner product: $\langle Uu, Uv \rangle = \langle u, v \rangle$.**
   ??? success "Solution"
       $\langle Uu, Uv \rangle = \langle u, U^* U v \rangle$.
       Since $U$ is unitary, $U^* U = I$.
       Thus $\langle u, Iv \rangle = \langle u, v \rangle$. This explains why rotations and reflections do not change angles or distances.

5. **[Normal Operators] Determine if $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ is a normal matrix.**
   ??? success "Solution"
       $A^* = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}$.
       $AA^* = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$.
       $A^*A = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$.
       $AA^* = A^*A$, so $A$ is normal.

6. **[Calculation] Compute the inner product of $f(x)=1$ and $g(x)=x$ under $\langle f, g \rangle = \int_0^1 f(x)g(x) dx$.**
   ??? success "Solution"
       $\langle 1, x \rangle = \int_0^1 x dx = [\frac{1}{2}x^2]_0^1 = 1/2$.

7. **[Polarization Identity] Write the formula for recovering the inner product from the norm in a real space.**
   ??? success "Solution"
       $\langle u, v \rangle = \frac{1}{4}(\|u+v\|^2 - \|u-v\|^2)$. This shows that in real spaces, length information is sufficient to reconstruct angular information.

8. **[Projection Theorem] Let $W$ be a subspace of $V$. Prove $V = W \oplus W^\perp$.**
   ??? success "Solution"
       For any $v \in V$, let $p = \operatorname{proj}_W v$. Let $e = v - p$. Clearly $p \in W$ and $\langle e, w \rangle = 0$ for all $w \in W$, so $e \in W^\perp$. This decomposition is unique.

9. **[Schur's Theorem App] Is every complex square matrix unitarily similar to an upper triangular matrix?**
   ??? success "Solution"
       Yes. This is Schur Decomposition. It guarantees that through unitary transformations (numerically stable), we can reveal the eigenvalues of any matrix on the diagonal.

10. **[Spectral Theorem] Prove: If $T$ is self-adjoint, its eigenvalues must be real.**
    ??? success "Solution"
        Let $Tv = \lambda v$, then $\langle Tv, v \rangle = \langle \lambda v, v \rangle = \lambda \|v\|^2$.
        Also $\langle Tv, v \rangle = \langle v, T^* v \rangle = \langle v, Tv \rangle = \bar{\lambda} \|v\|^2$.
        Since $v \neq 0$, $\lambda = \bar{\lambda}$.

## Chapter Summary

Inner product spaces establish the weights and measures of linear algebra:

1. **Geometrization**: Introduced norms, angles, and projections, giving abstract spaces physical intuition.
2. **Operator Symmetry**: The classification of self-adjoint, normal, and unitary operators forms the core of spectral theory.
3. **Essence of Approximation**: Orthogonal projection is the mathematical destination for all linear approximation and error minimization problems.
