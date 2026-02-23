# Chapter 19: Kronecker Product and Vec Operator

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Matrix Equations (Ch20) · Tensor Products (Ch21)

**Chapter Outline**: Definition of the Kronecker Product → Basic Algebraic Properties (Product, Transpose, Inverse) → Spectral Properties (Eigenvalues and Trace) → The Vec Operator and its Properties → The Fundamental Identity $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ → The Kronecker Sum ($A \oplus B$) → Applications: Vectorization of Sylvester and Lyapunov Matrix Equations → Tensor Structures in Signal Processing

**Extension**: The Kronecker product is the concrete realization of the Tensor Product in the realm of matrices; it is the fundamental algebraic language for handling multi-variable systems, multi-dimensional signals, and quantum entangled states (Ch28).

</div>

When dealing with complex equations involving the interaction of multiple matrices (such as $AX + XB = C$), traditional matrix arithmetic often fails to provide a direct closed-form solution. The Kronecker product and the Vec operator provide a standard toolkit for transforming "matrix equations" into "vector equations," allowing us to tackle high-dimensional operator interactions using classic linear systems theory.

---

## 19.1 The Kronecker Product

!!! definition "Definition 19.1 (Kronecker Product)"
    Let $A$ be an $m \times n$ matrix and $B$ be a $p \times q$ matrix. Their **Kronecker product** $A \otimes B$ is an $mp \times nq$ block matrix:
    $$A \otimes B = \begin{pmatrix} a_{11}B & \cdots & a_{1n}B \\ \vdots & \ddots & \vdots \\ a_{m1}B & \cdots & a_{mn}B \end{pmatrix}$$

!!! theorem "Theorem 19.1 (Core Properties)"
    1.  **Mixed-product Property**: $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$.
    2.  **Spectral Property**: If $A$ has eigenvalues $\{\lambda_i\}$ and $B$ has eigenvalues $\{\mu_j\}$, then $A \otimes B$ has eigenvalues $\{\lambda_i \mu_j\}$.
    3.  **Trace and Determinant**:
        - $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A)\operatorname{tr}(B)$.
        - $\det(A \otimes B) = (\det A)^p (\det B)^m$.

---

## 19.2 The Vec Operator and Vectorization

!!! definition "Definition 19.2 (Vec Operator)"
    For an $m \times n$ matrix $X$, $\operatorname{vec}(X)$ is the $mn \times 1$ vector obtained by stacking the columns of $X$ one on top of the other.

!!! theorem "Theorem 19.2 (Matrix Equation Vectorization Identity)"
    For matrices $A, X, B$ of appropriate dimensions:
    $$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$$
    **Significance**: This is the most important formula in this chapter. It extracts the matrix variable $X$, transforming the equation into the standard form $My = f$.

---

## 19.3 The Kronecker Sum

!!! definition "Definition 19.3 (Kronecker Sum)"
    For square matrices $A$ (order $n$) and $B$ (order $m$):
    $$A \oplus B = A \otimes I_m + I_n \otimes B$$
    **Property**: The eigenvalues of $A \oplus B$ are $\{\lambda_i + \mu_j\}$.
    **Application**: Solving the Sylvester equation $AX + XB = C$ is equivalent to solving the linear system $(I \otimes A + B^T \otimes I) \operatorname{vec}(X) = \operatorname{vec}(C)$.

---

## Exercises

1. **[Calculation] Compute $\begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix} \otimes \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $\begin{pmatrix} 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 2 \\ 0 & 0 & 2 & 0 \end{pmatrix}$.

2. **[Eigenvalues] If $A$ has eigenvalues 1, 2 and $B$ has eigenvalues 3, 4, find the eigenvalues of $A \otimes B$.**
   ??? success "Solution"
       $\{1\cdot 3, 1\cdot 4, 2\cdot 3, 2\cdot 4\} = \{3, 4, 6, 8\}$.

3. **[Kronecker Sum] Find the eigenvalues of $A \oplus B$ for the matrices above.**
   ??? success "Solution"
       $\{1+3, 1+4, 2+3, 2+4\} = \{4, 5, 5, 6\}$.

4. **[Vectorization] Transform $AX=B$ into vectorized form.**
   ??? success "Solution"
       $(I \otimes A) \operatorname{vec}(X) = \operatorname{vec}(B)$.

5. **[Trace] Prove $\operatorname{tr}(A \otimes B) = \operatorname{tr}(B \otimes A)$.**
   ??? success "Solution"
       Both equal $\operatorname{tr}(A)\operatorname{tr}(B)$.

6. **[Inverse] If $A$ and $B$ are invertible, find $(A \otimes B)^{-1}$.**
   ??? success "Solution"
       $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$.

7. **[Rank] Prove $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A)\operatorname{rank}(B)$.**
   ??? success "Solution"
       Can be derived via SVD or eigenvalue multiplicity. The number of non-zero singular values is the product of the individual counts.

8. **[Vec] For $X = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, write $\operatorname{vec}(X)$.**
   ??? success "Solution"
       $(a, c, b, d)^T$.

9. **[Lyapunov] Vectorize the Lyapunov equation $AX + XA^T = Q$.**
   ??? success "Solution"
       $(I \otimes A + A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(Q)$.

10. **[Quantum] Why are joint states of two particles represented by a tensor product (Kronecker product)?**
    ??? success "Solution"
        Because the dimension of a composite system is the product of the dimensions of its subsystems. The Kronecker product perfectly characterizes the combination of degrees of freedom and the possibility of quantum entanglement.

## Chapter Summary

The Kronecker product and Vec operator provide a scheme for "dimension elevation and reduction" in matrix algebra:

1.  **Dimension Multiplication**: The Kronecker product integrates the actions of two independent operators into a single massive composite operator via a tiling-and-nesting approach—the only language for multi-body interactions.
2.  **Operator Deconstruction**: The Vec operator eliminates the two-dimensional topology of a matrix, reducing it to its most basic vector form in linear space, thereby unleashing the full power of classical linear solvers.
3.  **Equation Unification**: The vectorization identity bridges matrix equation theory and numerical linear algebra, proving that all linear matrix equations are essentially the same linear system viewed under different bases.
