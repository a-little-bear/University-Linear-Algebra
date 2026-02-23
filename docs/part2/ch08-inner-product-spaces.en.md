# Chapter 08: Inner Product Spaces

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Orthogonality Basics (Ch07)

**Chapter Outline**: Definition of Abstract Inner Products → Axioms → Cauchy-Schwarz Inequality & Geometric Proof → Generalization of Norms and Distance → Orthogonal Complement $W^\perp$ → Algebraic Essence of the Adjoint Operator $T^*$ → Normal, Self-adjoint (Hermitian), and Unitary Operators → Deep Dive into the Spectral Theorem → Positive Operators and Matrix Square Roots → Introduction to Completeness and Finite-dimensional Hilbert Spaces

**Extension**: Inner product spaces introduce continuity and geometry to linear spaces, serving as a finite-dimensional microcosm of functional analysis; the spectral theorem is the mathematical core of "observables" in quantum mechanics, revealing the purest structures of operators.

</div>

Inner product spaces are a reinforcement of vector spaces. By defining an abstract inner product, we can talk not only about linear combinations but also about lengths, angles, and approximations. It is the mathematical foundation for modern signal processing, quantum mechanics, and optimization theory. This chapter moves beyond the dot product in $\mathbb{R}^n$ to establish a universal framework for metric geometry.

---

## 8.1 Definition of Abstract Inner Products

!!! definition "Definition 08.1 (Inner Product)"
    A vector space $V$ over a field $F$ ($\mathbb{R}$ or $\mathbb{C}$) is equipped with an **inner product** $\langle \cdot, \cdot \rangle$ if it satisfies the following axioms:
    1.  **Positivity**: $\langle v, v \rangle \ge 0$, and $\langle v, v \rangle = 0$ iff $v=0$.
    2.  **Conjugate Symmetry**: $\langle u, v \rangle = \overline{\langle v, u \rangle}$.
    3.  **Additivity and Homogeneity (Linearity in first slot)**: $\langle au+bv, w \rangle = a\langle u, w \rangle + b\langle v, w \rangle$.

!!! note "Induced Norm"
    Every inner product naturally induces a norm: $\|v\| = \sqrt{\langle v, v \rangle}$.

---

## 8.2 Fundamental Inequalities and Complements

!!! theorem "Theorem 08.1 (Cauchy-Schwarz Inequality)"
    For any vectors $u, v$ in an inner product space:
    $$|\langle u, v \rangle| \le \|u\| \|v\|$$
    Equality holds iff $u$ and $v$ are linearly dependent. This inequality is the basis for defining the "angle" between vectors.

!!! definition "Definition 08.2 (Orthogonal Complement $W^\perp$)"
    Let $W$ be a subspace of $V$. The **orthogonal complement** of $W$ is the set of vectors orthogonal to every vector in $W$:
    $$W^\perp = \{ v \in V : \langle v, w \rangle = 0, \forall w \in W \}$$
    **Property**: $V = W \oplus W^\perp$ (The Projection Theorem).

---

## 8.3 Adjoint Operators and Classification

!!! definition "Definition 08.3 (Adjoint Operator $T^*$)"
    Let $T$ be a linear operator on $V$. If there exists an operator $T^*$ such that for all $u, v \in V$:
    $$\langle Tu, v \rangle = \langle u, T^* v \rangle$$
    then $T^*$ is called the **adjoint** of $T$. In an orthonormal basis, its matrix representation is the conjugate transpose $A^*$.

!!! theorem "Theorem 08.2 (Operator Classification)"
    1.  **Self-adjoint (Hermitian)**: $T^* = T$. All its eigenvalues are real.
    2.  **Unitary**: $T^* T = TT^* = I$. Preserves lengths and angles.
    3.  **Normal**: $T^* T = TT^*$.

---

## 8.4 The Spectral Theorem

!!! theorem "Theorem 08.3 (Complex Spectral Theorem)"
    A linear operator $T \in \mathcal{L}(V)$ is normal iff there exists an orthonormal basis of $V$ consisting of eigenvectors of $T$.
    **Physical Meaning**: This means that a normal operator can be completely decoupled into independent directional scalings by rotating the coordinate axes.

---

## Exercises

**1. [Axioms] Verify if $\langle u, v \rangle = u_1 v_1 + 2u_2 v_2$ on $\mathbb{R}^2$ is an inner product.**

??? success "Solution"
    **Verification Steps:**
    1. **Positivity**: $\langle u, u \rangle = u_1^2 + 2u_2^2$. Since squares are non-negative and coefficients are positive, the result is $\ge 0$. It is 0 iff $u_1=u_2=0$.
    2. **Symmetry**: $\langle u, v \rangle = u_1 v_1 + 2u_2 v_2 = v_1 u_1 + 2v_2 u_2 = \langle v, u \rangle$.
    3. **Linearity**: It is clearly linear with respect to the components of $u$.
    **Conclusion**: This is a weighted inner product and satisfies all axioms. It defines a geometry where the "unit sphere" is an ellipse with axes scaled by $\sqrt{2}$.

**2. [Cauchy-Schwarz] Prove: $|\langle u, v \rangle| \le \|u\| \|v\|$.**

??? success "Solution"
    **Proof:**
    1. If $v=0$, both sides are 0, and it holds.
    2. If $v \neq 0$, consider the norm squared of $u - \lambda v$:
       $\|u - \lambda v\|^2 = \langle u - \lambda v, u - \lambda v \rangle \ge 0$.
    3. Expand: $\langle u, u \rangle - \bar{\lambda}\langle u, v \rangle - \lambda\langle v, u \rangle + |\lambda|^2 \langle v, v \rangle \ge 0$.
    4. Let $\lambda = \frac{\langle u, v \rangle}{\langle v, v \rangle}$ (the projection coefficient).
    5. Substitute: $\|u\|^2 - \frac{|\langle u, v \rangle|^2}{\|v\|^2} \ge 0$.
    6. Multiply by $\|v\|^2$ and take the square root to get $|\langle u, v \rangle| \le \|u\| \|v\|$.

**3. [Adjoint] Let $T$ have matrix $A$ in an orthonormal basis. Prove the matrix of $T^*$ is $A^*$.**

??? success "Solution"
    **Proof:**
    1. Matrix expression of inner product: $\langle Tu, v \rangle = (Au)^* v = u^* A^* v$.
    2. By adjoint definition: $\langle u, T^* v \rangle = u^* [T^*]_B v$.
    3. Since this holds for all $u, v$, the matrices must be equal.
    **Conclusion**: The matrix representation of the adjoint is the conjugate transpose of the original matrix.

**4. [Unitary] Prove unitary operators preserve the inner product: $\langle Uu, Uv \rangle = \langle u, v \rangle$.**

??? success "Solution"
    **Derivation:**
    1. $\langle Uu, Uv \rangle = \langle u, U^* U v \rangle$ (using the definition of the adjoint).
    2. By unitary definition: $U^* U = I$.
    3. $= \langle u, Iv \rangle = \langle u, v \rangle$.
    **Geometric Insight**: Unitary transforms (or orthogonal transforms in $\mathbb{R}$) are rigid motions that do not change lengths or angles in space.

**5. [Normal] Determine if $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$ is a normal matrix.**

??? success "Solution"
    **Computation:**
    1. $A^* = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}$.
    2. $AA^* = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = 2I$.
    3. $A^*A = \begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = 2I$.
    **Conclusion**: Since $AA^* = A^*A$, the matrix is normal. By the spectral theorem, it is unitarily diagonalizable over $\mathbb{C}$.

**6. [Calculation] In $P_1$ with $\langle f, g \rangle = \int_0^1 f(x)g(x) dx$, find the angle between $1$ and $x$.**

??? success "Solution"
    **Steps:**
    1. Inner product: $\langle 1, x \rangle = \int_0^1 x dx = 0.5$.
    2. Norms: $\|1\|^2 = 1$; $\|x\|^2 = \int_0^1 x^2 dx = 1/3$.
    3. Cosine: $\cos\theta = \frac{0.5}{1 \cdot \sqrt{1/3}} = \frac{\sqrt{3}}{2}$.
    **Conclusion**: The angle $\theta = 30^\circ$ or $\pi/6$.

**7. [Polarization] Write the formula for the inner product using only norms in a real space.**

??? success "Solution"
    **Formula:**
    $$\langle u, v \rangle = \frac{1}{4}(\|u+v\|^2 - \|u-v\|^2)$$
    **Proof:**
    Expanding the RHS: $\frac{1}{4}(\langle u+v, u+v \rangle - \langle u-v, u-v \rangle) = \frac{1}{4}( (\|u\|^2 + 2\langle u, v \rangle + \|v\|^2) - (\|u\|^2 - 2\langle u, v \rangle + \|v\|^2) ) = \langle u, v \rangle$.
    This shows the inner product information is entirely contained within the norm (length) information.

**8. [Projection] Prove $V = W \oplus W^\perp$.**

??? success "Solution"
    **Proof Sketch:**
    1. **Existence**: For any $v$, let $p = \sum \langle v, e_i \rangle e_i$ where $\{e_i\}$ is an orthonormal basis for $W$. Clearly $p \in W$.
    2. Let $z = v - p$. One can check $\langle z, e_i \rangle = 0$ for all $i$, so $z \in W^\perp$. Thus $v = p + z$.
    3. **Uniqueness**: If $x \in W \cap W^\perp$, then $\langle x, x \rangle = 0 \implies x = 0$.
    **Conclusion**: The space is the direct sum of any subspace and its orthogonal complement.

**9. [Schur] Every complex square matrix is unitarily similar to which type of matrix?**

??? success "Solution"
    **Conclusion:**
    Every complex square matrix is **unitarily similar to an upper triangular matrix**.
    The diagonal entries are exactly the **eigenvalues** of the matrix. This is the foundation of spectral theory, proving that by rotating the basis, we can make any operator appear hierarchical.

**10. [Self-adjoint] Prove eigenvalues of a self-adjoint operator are real.**

??? success "Solution"
    **Proof:**
    1. Let $Tv = \lambda v$ with $v \neq 0$.
    2. Consider $\langle Tv, v \rangle = \lambda \|v\|^2$.
    3. Also, $\langle Tv, v \rangle = \langle v, T^* v \rangle = \langle v, Tv \rangle$ (since $T^*=T$).
    4. Thus $\langle v, \lambda v \rangle = \bar{\lambda} \|v\|^2$.
    5. So $\lambda \|v\|^2 = \bar{\lambda} \|v\|^2 \implies \lambda = \bar{\lambda}$ (since $\|v\|^2 > 0$).
    **Conclusion**: The eigenvalues must be real.

## Chapter Summary

Inner product spaces establish the metrics of linear algebra:

1.  **Geometrization**: The introduction of norms, angles, and projections gives abstract spaces physical intuition and approximation power.
2.  **Operator Symmetry**: The classification of self-adjoint, normal, and unitary operators forms the core of spectral theory, revealing the most stable structures of linear operators.
3.  **Approximation Essence**: Orthogonal projection is the mathematical endpoint for solving all linear approximation and error-minimization problems.
