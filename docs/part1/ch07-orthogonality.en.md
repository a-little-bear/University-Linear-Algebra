# Chapter 07: Orthogonality and the QR Decomposition

<div class="context-flow" markdown>

**Prerequisites**: Inner Product (Ch8) · Vector Spaces (Ch4) · Matrix Factorization (Ch10)

**Chapter Outline**: Inner Product and Norms → Orthogonality and Orthonormality → Orthogonal Projections → Gram-Schmidt Process → QR Decomposition ($A = QR$) → Orthogonal Matrices and their Properties → Best Approximation and Least Squares → Householder Reflections → Givens Rotations

**Extension**: Orthogonality is the geometry of "perpendicularity"; it provides the most numerically stable basis for computation and allows for projecting high-dimensional data onto lower-dimensional subspaces.

</div>

Orthogonality is the concept of "being at right angles," generalized to any dimension. In linear algebra, orthogonal bases are preferred because they minimize numerical error and simplify calculations (e.g., coordinates are found via simple projections). This chapter details the **Gram-Schmidt process**, which converts any basis into an orthonormal one, leading directly to the **QR decomposition**—a fundamental tool for solving least squares problems and finding eigenvalues.

---

## 07.1 Orthogonality and Projections

!!! definition "Definition 07.1 (Orthogonality)"
    Two vectors $u, v$ are **orthogonal** if their inner product is zero: $\langle u, v \rangle = 0$. A set of vectors is **orthonormal** if they are pairwise orthogonal and each has unit length ($\|v\| = 1$).

!!! theorem "Theorem 07.1 (Orthogonal Projection)"
    The orthogonal projection of $v$ onto a subspace $W$ with orthonormal basis $\{q_1, \dots, q_k\}$ is given by:
    $$\operatorname{proj}_W(v) = \sum_{i=1}^k \langle v, q_i \rangle q_i$$
    This is the vector in $W$ that is closest to $v$.

---

## Exercises

1. **[Fundamentals] Are $(1, 1)^T$ and $(1, -1)^T$ orthogonal in $\mathbb{R}^2$?**
   ??? success "Solution"
       $\langle u, v \rangle = 1(1) + 1(-1) = 0$. Yes, they are orthogonal.

2. **[Norm] Normalize the vector $v = (3, 4)^T$.**
   ??? success "Solution"
       $\|v\| = \sqrt{3^2 + 4^2} = 5$. The normalized vector is $\hat{v} = (0.6, 0.8)^T$.

3. **[Gram-Schmidt] Apply Gram-Schmidt to $\{(1, 1, 0), (1, 0, 1)\}$ to find an orthogonal basis.**
   ??? success "Solution"
       $q_1 = (1, 1, 0)$. Then $v_2' = (1, 0, 1) - \frac{(1,0,1) \cdot (1,1,0)}{2} (1,1,0) = (1, 0, 1) - (0.5, 0.5, 0) = (0.5, -0.5, 1)$.

4. **[Orthogonal Matrix] State the defining property of an $n \times n$ orthogonal matrix $Q$.**
   ??? success "Solution"
       $Q^T Q = I$, which implies $Q^{-1} = Q^T$. The columns of $Q$ form an orthonormal basis for $\mathbb{R}^n$.

5. **[Projection] Find the projection of $(1, 2, 3)^T$ onto the $xy$-plane.**
   ??? success "Solution"
       Using the orthonormal basis $\{e_1, e_2\}$, the projection is $1e_1 + 2e_2 = (1, 2, 0)^T$.

6. **[QR Decomposition] If $A = QR$ is the QR decomposition, what is $R$?**
   ??? success "Solution"
       $R$ is an upper triangular matrix with positive diagonal entries. Its entries are the coefficients of the Gram-Schmidt process.

7. **[Least Squares] How is QR used to solve $Ax = b$ in the least squares sense?**
   ??? success "Solution"
       $Ax = b \implies QRx = b \implies Rx = Q^T b$. Since $R$ is triangular, this is solved efficiently by back-substitution.

8. **[Invariance] Does an orthogonal transformation preserve the length of vectors?**
   ??? success "Solution"
       Yes. $\|Qv\|^2 = (Qv)^T (Qv) = v^T Q^T Q v = v^T I v = \|v\|^2$. Orthogonal matrices represent isometries (rotations and reflections).

9. **[Householder] What is a Householder reflection matrix?**
   ??? success "Solution"
       $H = I - 2uu^T$ where $\|u\|=1$. It reflects space across the hyperplane orthogonal to $u$. It is symmetric and orthogonal.

10. **[Pythagoras] State the Pythagorean Theorem for orthogonal vectors.**
    ??? success "Solution"
        If $u \perp v$, then $\|u+v\|^2 = \|u\|^2 + \|v\|^2$. This is the foundation for the best approximation property of projections.

## Chapter Summary

This chapter explores the geometric and numerical power of perpendicularity:

1. **Inner Product Space**: Defined orthogonality as the vanishing of the inner product, generalizing Euclidean angles.
2. **Structural Orthogonalization**: Developed the Gram-Schmidt process as the algorithm for constructing orthonormal bases.
3. **Factorization Stable**: Formulated the QR decomposition as the stable alternative to Gaussian elimination for solving linear systems.
4. **Distance Minimization**: Linked orthogonal projections to the best approximation problem, providing the mathematical basis for least squares.
