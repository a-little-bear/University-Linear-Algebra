# Chapter 33: Generalized Inverses

<div class="context-flow" markdown>

**Prerequisites**: Matrix Rank (Ch2) · Orthogonal Projection (Ch5) · SVD (Ch11) · Least Squares (Ch72B)

**Chapter Outline**: Definition of Generalized Inverses → The Four Penrose Conditions → Moore-Penrose Inverse $A^\dagger$ → Existence and Uniqueness → Calculation via SVD and QR → Properties ($A^\dagger A$ as Projection) → Solving Overdetermined Systems → Minimum Norm Solutions

**Extension**: Generalized inverses are the engine of numerical linear algebra, providing stable solutions to rank-deficient and overdetermined systems.

</div>

Generalized inverses extend the concept of a matrix inverse to rectangular and singular matrices. While a standard inverse $A^{-1}$ only exists for square, non-singular matrices, the **Moore-Penrose inverse** $A^\dagger$ is unique and exists for *any* matrix. It provides the "best possible" solution to linear systems in the sense of least squares and minimum norm, mapping the column space of $A^T$ to the column space of $A$ bijectively.

---

## 33.1 The Moore-Penrose Conditions

!!! definition "Definition 33.1 (Moore-Penrose Inverse)"
    For any $m \times n$ matrix $A$, the Moore-Penrose inverse $A^\dagger$ is the unique $n \times m$ matrix satisfying the following four conditions:
    1. $A A^\dagger A = A$
    2. $A^\dagger A A^\dagger = A^\dagger$
    3. $(A A^\dagger)^* = A A^\dagger$ (Hermitian)
    4. $(A^\dagger A)^* = A^\dagger A$ (Hermitian)

!!! theorem "Theorem 33.1 (Construction via SVD)"
    If $A = U \Sigma V^*$ is the SVD of $A$, then $A^\dagger = V \Sigma^\dagger U^*$, where $\Sigma^\dagger$ is obtained by inverting the non-zero singular values and transposing.

---

## Exercises

1. **[Fundamentals] Does every matrix (including non-square and singular ones) have a Moore-Penrose inverse? Is it unique?**
   ??? success "Solution"
       Yes, every matrix $A$ over $\mathbb{C}$ has a unique Moore-Penrose inverse $A^\dagger$. While other generalized inverses (like $\{1\}$-inverses satisfying $AGA=A$) may exist and not be unique, $A^\dagger$ is uniquely determined by the four Penrose conditions.

2. **[Projections] Show that $P = A A^\dagger$ is the orthogonal projection onto the column space $\operatorname{Im}(A)$.**
   ??? success "Solution"
       Conditions 1 and 3 imply $P^2 = (A A^\dagger A) A^\dagger = A A^\dagger = P$ and $P^* = P$. Thus $P$ is an orthogonal projection. Since $P A = A$, its image contains $\operatorname{Im}(A)$, and since $P = A(A^\dagger)$, its image is contained in $\operatorname{Im}(A)$.

3. **[Full Rank] Derive $A^\dagger$ for a matrix $A$ with full column rank.**
   ??? success "Solution"
       If $A$ has full column rank, $A^* A$ is invertible. $A^\dagger = (A^* A)^{-1} A^*$. This is the standard "left inverse" used in least squares.

4. **[Linear Systems] Prove that $x = A^\dagger b$ is the minimum norm solution to the least squares problem $\min \|Ax - b\|^2$.**
   ??? success "Solution"
       Least squares requires $Ax = P_{\operatorname{Im}(A)} b = A A^\dagger b$. $x = A^\dagger b$ satisfies this. Any other solution $x'$ differs by an element in $\ker(A)$. Since $A^\dagger b \in \operatorname{Im}(A^*) = (\ker A)^\perp$, the Pythagorean theorem ensures $\|A^\dagger b\|$ is the minimum norm.

5. **[Properties] Is $(AB)^\dagger = B^\dagger A^\dagger$ always true?**
   ??? success "Solution"
       No. It holds if $A$ has full column rank and $B$ has full row rank, or if $A^* A B B^*$ is Hermitian. Generally, the identity fails due to the interaction of the column and row spaces.

6. **[Rank-1] Find the pseudoinverse of the rank-1 matrix $A = uv^*$.**
   ??? success "Solution"
       $A^\dagger = \frac{1}{\|u\|^2 \|v\|^2} v u^*$.

7. **[Inversion] Show that $(A^\dagger)^\dagger = A$.**
   ??? success "Solution"
       The conditions for $A^\dagger$ being the pseudoinverse of $A$ are symmetric to those for $A$ being the pseudoinverse of $A^\dagger$.

8. **[Calculation] Compute $A^\dagger$ for $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.**
   ??? success "Solution"
       $A$ has full column rank. $A^* A = (2)$. $(A^* A)^{-1} = (1/2)$. $A^\dagger = (1/2) \begin{pmatrix} 1 & 1 \end{pmatrix} = \begin{pmatrix} 0.5 & 0.5 \end{pmatrix}$.

9. **[Kernel] Relate the kernel of $A^\dagger$ to the adjoint of $A$.**
   ??? success "Solution"
       $\ker(A^\dagger) = \ker(A^*)$. The pseudoinverse maps the image of $A$ back to the row space and sends the noise (orthogonal to the image) to zero.

10. **[Continuity] Is the mapping $A \mapsto A^\dagger$ continuous?**
    ??? success "Solution"
        No. It is discontinuous at any matrix where the rank changes. A small perturbation can create a tiny non-zero singular value, making its inverse $1/\sigma$ explode. This is why truncated SVD is used in practice.

## Chapter Summary

This chapter establishes the Moore-Penrose inverse as the universal generalized inverse:

1. **Analytical Uniqueness**: Defined $A^\dagger$ via the four Penrose conditions, ensuring a unique "best" inverse for any linear operator.
2. **Geometric Mapping**: Demonstrated that $A^\dagger$ effectively solves the fundamental problem of projecting onto the range and inverting on the row space.
3. **Optimal Solving**: Positioned the pseudoinverse as the mathematical tool for minimum-norm least-squares solutions.
4. **Numerical Stability**: Linked the pseudoinverse to the SVD, highlighting the challenges of rank-deficiency and sensitivity.
