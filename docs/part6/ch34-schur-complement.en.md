# Chapter 34: The Schur Complement

<div class="context-flow" markdown>

**Prerequisites**: Matrix Inversion (Ch2) · Determinants (Ch3) · Block Matrices · Positive Definiteness (Ch16)

**Chapter Outline**: Definition of Schur Complement → Block LU Decomposition → Determinant of Block Matrices → Matrix Inversion Lemma (Woodbury) → Positive Definiteness via Schur Complement → Inertia and Haynsworth's Law → Application in Gaussian Processes

**Extension**: The Schur complement is the algebraic engine for partial Gaussian elimination and the reduction of high-dimensional covariance matrices.

</div>

The **Schur complement** arises naturally when performing Gaussian elimination on block matrices. Given a partitioned matrix, the Schur complement of one block captures the "residual" information of the other block after accounting for the cross-correlations. It is a fundamental tool for computing inverses and determinants of large systems, and it provides a definitive criterion for the positive definiteness of block matrices.

---

## 34.1 Definitions and Block Factorization

!!! definition "Definition 34.1 (Schur Complement)"
    Let $M = \begin{pmatrix} A & B \ C & D \end{pmatrix}$ be a block matrix. If $A$ is invertible, the Schur complement of $A$ in $M$ is:
    $$M/A = D - C A^{-1} B$$

!!! theorem "Theorem 34.1 (Determinant Formula)"
    The determinant of $M$ is the product of the determinant of the block and the determinant of its Schur complement:
    $$\det M = \det A \cdot \det(D - C A^{-1} B)$$

---

## Exercises

1. **[Fundamentals] Compute the Schur complement of $A = (2)$ in $M = \begin{pmatrix} 2 & 1 \ 1 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       $M/A = 2 - (1)(2)^{-1}(1) = 2 - 0.5 = 1.5$. Note $\det M = 2(1.5) = 3$.

2. **[Block Inversion] State the formula for $M^{-1}$ using the Schur complement.**
   ??? success "Solution"
       $M^{-1} = \begin{pmatrix} A^{-1} + A^{-1}B(M/A)^{-1}CA^{-1} & -A^{-1}B(M/A)^{-1} \ -(M/A)^{-1}CA^{-1} & (M/A)^{-1} \end{pmatrix}$. This allows inverting a large matrix by inverting smaller sub-blocks.

3. **[Positive Definiteness] Show that $M = \begin{pmatrix} A & B \ B^* & D \end{pmatrix} \succ 0$ iff $A \succ 0$ and $M/A \succ 0$.**
   ??? success "Solution"
       Using the congruence $M = \begin{pmatrix} I & 0 \ B^* A^{-1} & I \end{pmatrix} \begin{pmatrix} A & 0 \ 0 & M/A \end{pmatrix} \begin{pmatrix} I & A^{-1}B \ 0 & I \end{pmatrix}$. Sylvester's Law of Inertia implies $M$ has the same number of positive eigenvalues as the block-diagonal matrix. Thus $M \succ 0 \iff A \succ 0$ and $M/A \succ 0$.

4. **[Woodbury Identity] State the Woodbury matrix identity (Matrix Inversion Lemma).**
   ??? success "Solution"
       $(A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1} + VA^{-1}U)^{-1}VA^{-1}$. This relates the inverse of a perturbed matrix to the Schur complement of the perturbation.

5. **[Variational] Express the Schur complement as a result of a minimization problem.**
   ??? success "Solution"
       For $M \succ 0$, the Schur complement $M/A$ appears in the minimization of the quadratic form: $\min_x \begin{pmatrix} x \ y \end{pmatrix}^T \begin{pmatrix} A & B \ B^T & D \end{pmatrix} \begin{pmatrix} x \ y \end{pmatrix} = y^T (M/A) y$. The optimal $x = -A^{-1}By$ eliminates the cross-term.

6. **[Haynsworth's Law] How does the inertia of $M$ relate to the inertia of $A$ and $M/A$?**
   ??? success "Solution"
       $\operatorname{In}(M) = \operatorname{In}(A) + \operatorname{In}(M/A)$. This additive law for the number of positive, negative, and zero eigenvalues is vital for stability analysis.

7. **[Orthogonal Projection] Interpret the Schur complement in terms of projections.**
   ??? success "Solution"
       The Schur complement $D - C A^{-1} B$ represents the portion of $D$ that is orthogonal to the subspace spanned by the columns of $B$ (in a weighted inner product sense).

8. **[Probability] In a Gaussian distribution $\begin{pmatrix} X_1 \ X_2 \end{pmatrix} \sim N(\mu, \Sigma)$, what is the conditional covariance of $X_2$ given $X_1$?**
   ??? success "Solution"
       The conditional covariance is exactly the Schur complement $\Sigma_{22} - \Sigma_{21} \Sigma_{11}^{-1} \Sigma_{12}$. This capture the "remaining uncertainty" in $X_2$ after $X_1$ is observed.

9. **[Monotonicity] Prove that if $M \succeq 0$, then $M/A \preceq D$.**
   ??? success "Solution"
       $M/A = D - C A^{-1} B$. Since $A \succeq 0 \implies A^{-1} \succeq 0$, the term $C A^{-1} B$ (which is $B^* A^{-1} B$ in the symmetric case) is positive semi-definite. Subtracting a PSD matrix from $D$ yields $M/A \preceq D$.

10. **[Rank] Relate the rank of $M$ to the ranks of $A$ and $M/A$.**
    ??? success "Solution"
        $\operatorname{rank}(M) = \operatorname{rank}(A) + \operatorname{rank}(M/A)$. This holds whenever the column space of $B$ is contained in the column space of $A$ and the row space of $C$ is contained in the row space of $A$.

## Chapter Summary

This chapter explores the reduction of block matrices via the Schur complement:

1. **Factorization Engine**: Positioned the Schur complement as the core of block LU decomposition and Gaussian elimination.
2. **Spectral Criteria**: Established definitive conditions for the positive definiteness and inertia of partitioned matrices.
3. **Inversion Lemmas**: Formulated the Woodbury identity as a powerful tool for updating inverses under low-rank perturbations.
4. **Statistical Depth**: Demonstrated its fundamental role in Gaussian conditioning and variance reduction in probability.
