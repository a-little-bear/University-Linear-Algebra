# Chapter 03: Determinants

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Polynomials (Ch00) · Permutations

**Chapter Outline**: Definition of the Determinant → Leibniz Formula (Sum over Permutations) → Properties of Determinants (Linearity, Alternating, Multiplicativity) → Calculation via Cofactor Expansion → Calculation via Row Reduction → Cramer's Rule → Determinants and Invertibility → Geometric Interpretation (Oriented Volume)

**Extension**: The determinant is the scalar "signature" of a linear map; it is the unique alternating multilinear form that measures how a transformation scales space.

</div>

The determinant is a scalar value that encapsulates the essence of a square matrix. Geometrically, it measures the factor by which a transformation scales $n$-dimensional volume. Algebraically, it provides a definitive test for invertibility: a matrix is invertible if and only if its determinant is non-zero. From the combinatorial Leibniz formula to the practical cofactor expansion, this chapter explores the properties that make the determinant an essential tool in spectral theory and geometry.

---

## 03.1 Definitions and Volume Scaling

!!! definition "Definition 03.1 (Leibniz Formula)"
    For an $n \times n$ matrix $A$, the determinant is:
    $$\det A = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{i, \sigma(i)}$$
    This sum over all permutations captures the alternating nature of the volume form.

!!! theorem "Theorem 03.1 (Multiplicativity)"
    For any two square matrices $A$ and $B$, $\det(AB) = \det A \cdot \det B$. This mirrors the property that the volume scaling of a composed transformation is the product of the individual scalings.

---

## Exercises

1. **[Fundamentals] Compute the determinant of $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$.**
   ??? success "Solution"
       $\det A = ad - bc$. This represents the oriented area of the parallelogram spanned by the columns.

2. **[Triangular] Show that the determinant of a triangular matrix is the product of its diagonal entries.**
   ??? success "Solution"
       In the Leibniz sum, any permutation that picks an entry from below (or above) the diagonal for a triangular matrix will eventually be forced to pick a zero entry. The only non-zero product corresponds to the identity permutation $\sigma = id$.

3. **[Row Ops] How does a row swap affect the determinant? What about scaling a row by $k$?**
   ??? success "Solution"
       A row swap multiplies the determinant by $-1$ (alternating property). Scaling a single row by $k$ multiplies the determinant by $k$ (linearity).

4. **[Cramer] Solve $x+y=3, x+2y=5$ using Cramer's Rule.**
   ??? success "Solution"
       $D = \det \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix} = 1$. $D_x = \det \begin{pmatrix} 3 & 1 \\ 5 & 2 \end{pmatrix} = 1$. $D_y = \det \begin{pmatrix} 1 & 3 \\ 1 & 5 \end{pmatrix} = 2$. Solution: $x = D_x/D = 1, y = D_y/D = 2$.

5. **[Invertibility] Prove that $\det A \neq 0$ iff $A$ is invertible.**
   ??? success "Solution"
       If $A$ is invertible, $A A^{-1} = I \implies \det A \det A^{-1} = 1$, so $\det A \neq 0$. If $\det A \neq 0$, the adjoint matrix formula $A^{-1} = \frac{1}{\det A} \operatorname{adj}(A)$ provides an explicit inverse.

6. **[Transpose] Prove that $\det(A^T) = \det A$.**
   ??? success "Solution"
       The Leibniz formula can be summed over columns instead of rows. Since the sign of a permutation is the same as the sign of its inverse, the sum remains invariant.

7. **[Vandermonde] Compute the determinant of $\begin{pmatrix} 1 & a & a^2 \\ 1 & b & b^2 \\ 1 & c & c^2 \end{pmatrix}$.**
   ??? success "Solution"
       The Vandermonde determinant is $(b-a)(c-a)(c-b)$. It is non-zero iff the nodes $a, b, c$ are distinct.

8. **[Nilpotency] If $A^k = 0$, what is $\det A$?**
   ??? success "Solution"
       $(\det A)^k = \det(A^k) = \det(0) = 0$. Thus $\det A = 0$. Nilpotent matrices are always singular.

9. **[Orthogonality] If $Q$ is an orthogonal matrix, what are the possible values for $\det Q$?**
   ??? success "Solution"
       $Q^T Q = I \implies \det(Q^T) \det Q = 1 \implies (\det Q)^2 = 1$. Thus $\det Q = \pm 1$. $+1$ corresponds to rotations, $-1$ to reflections.

10. **[Block Matrix] Compute the determinant of $\begin{pmatrix} A & B \\ 0 & D \end{pmatrix}$.**
    ??? success "Solution"
        The determinant is $\det A \cdot \det D$. This follows from the block LU decomposition or by observing that any permutation contributing a non-zero term must stay within the diagonal blocks.

## Chapter Summary

This chapter explores the scalar invariant that governs solvability and geometry:

1. **Combinatorial Definition**: Formulated the determinant via permutations, establishing its alternating multilinear structure.
2. **Computational Paths**: Contrasted cofactor expansion with the efficiency of row reduction.
3. **Analytic Criterion**: Defined the non-zero determinant as the unique requirement for operator invertibility.
4. **Geometric Scaling**: Linked the determinant to the magnification of volumes under linear maps.
