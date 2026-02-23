# Chapter 02: Matrices and Matrix Operations

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01)

**Chapter Outline**: Definition and Notation → Basic Operations (Addition, Scalar Multiplication, Multiplication) → Non-commutativity of Multiplication → Transpose and its Properties → Special Matrices (Identity, Diagonal, Triangular, Symmetric) → Elementary Matrices & Row Operations → Inverse Matrices: Definition and Properties → Gauss-Jordan Method for Inversion → Block Matrix Operations → Trace of a Matrix

**Extension**: Matrices are not just containers for data but representation of linear operators; the definition of matrix multiplication reflects the composition of linear transformations (Ch05).

</div>

If systems of linear equations are the language of linear algebra, then matrices are its notation system. A matrix condenses complex linear relationships into concise rectangular arrays and endows them with a set of sophisticated algebraic rules. This chapter establishes the standard axioms of matrix algebra and explores the essential tool of the inverse matrix.

---

## 02.1 Basic Definitions and Operations

!!! definition "Definition 02.1 (Matrix)"
    An $m \times n$ **matrix** is a rectangular array of $m \cdot n$ elements arranged in $m$ rows and $n$ columns. Usually denoted by uppercase letters $A, B$.

!!! definition "Definition 02.2 (Matrix Multiplication)"
    If $A$ is an $m \times n$ matrix and $B$ is an $n \times p$ matrix, their product $C = AB$ is an $m \times p$ matrix with entries:
    $$c_{ij} = \sum_{k=1}^n a_{ik} b_{kj}$$
    **WARNING**: Matrix multiplication is generally not commutative, i.e., $AB \neq BA$.

---

## 02.2 Special Matrix Classes

!!! definition "Definition 02.3 (Special Matrices)"
    1.  **Identity Matrix $I$**: Diagonal entries are 1, all others 0. Satisfies $AI = IA = A$.
    2.  **Symmetric Matrix**: Satisfies $A^T = A$.
    3.  **Skew-symmetric Matrix**: Satisfies $A^T = -A$.
    4.  **Triangular Matrix**: Upper triangular (all entries below the main diagonal are 0) or Lower triangular.

---

## 02.3 Elementary Matrices and Inverses

!!! definition "Definition 02.4 (Inverse Matrix)"
    For a square matrix $A$, if there exists a matrix $B$ such that $AB = BA = I$, then $A$ is **invertible** (or non-singular). $B$ is called the **inverse** of $A$, denoted $A^{-1}$.

!!! theorem "Theorem 02.1 (Properties of Inverses)"
    1.  $(A^{-1})^{-1} = A$
    2.  $(AB)^{-1} = B^{-1} A^{-1}$ (Reversal law)
    3.  $(A^T)^{-1} = (A^{-1})^T$

!!! algorithm "Algorithm 02.1 (Gauss-Jordan Inversion)"
    Construct the block matrix $[A | I]$. Apply elementary row operations to transform the left side into $I$. The resulting right side is $A^{-1}$:
    $$[A | I] \xrightarrow{\text{row operations}} [I | A^{-1}]$$

---

## 02.4 Block Matrices

!!! technique "Technique: Block Operations"
    For large-scale matrices, it is useful to partition them into smaller sub-blocks. If block sizes are compatible, addition and multiplication rules are identical to those for standard matrices. This is vital for distributed computing and sparse matrix handling.

---

## Exercises

1. **[Fundamentals] Let $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}, B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. Calculate $AB$ and $BA$.**
   ??? success "Solution"
       $AB = \begin{pmatrix} 2 & 1 \\ 4 & 3 \end{pmatrix}, BA = \begin{pmatrix} 3 & 4 \\ 1 & 2 \end{pmatrix}$. Note $AB \neq BA$.

2. **[Identity] Prove that for any $n \times n$ matrix $A$, $AI = IA = A$.**
   ??? success "Solution"
       Using the multiplication definition: $(AI)_{ij} = \sum a_{ik} \delta_{kj}$. Since $\delta_{kj}$ is 1 only when $k=j$, the result is $a_{ij}$.

3. **[Transpose] Given $(AB)^T = B^T A^T$, find $(A^T B)^T$.**
   ??? success "Solution"
       $(A^T B)^T = B^T (A^T)^T = B^T A$.

4. **[Symmetry] If $A$ is symmetric, prove $A^2$ is also symmetric.**
   ??? success "Solution"
       $(A^2)^T = (AA)^T = A^T A^T = AA = A^2$.

5. **[Inversion] Find the inverse of $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**
   ??? success "Solution"
       $\det(A) = 4-6 = -2$. $A^{-1} = \frac{1}{-2} \begin{pmatrix} 4 & -2 \\ -3 & 1 \end{pmatrix} = \begin{pmatrix} -2 & 1 \\ 1.5 & -0.5 \end{pmatrix}$.

6. **[Powers] Calculate $A^k$ for $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $A^2 = \begin{pmatrix} 1 & 2 \\ 0 & 1 \end{pmatrix}, A^3 = \begin{pmatrix} 1 & 3 \\ 0 & 1 \end{pmatrix}$. By induction, $A^k = \begin{pmatrix} 1 & k \\ 0 & 1 \end{pmatrix}$.

7. **[Trace] Prove $\operatorname{tr}(AB) = \operatorname{tr}(BA)$.**
   ??? success "Solution"
       $\operatorname{tr}(AB) = \sum_i \sum_k a_{ik}b_{ki} = \sum_k \sum_i b_{ki}a_{ik} = \operatorname{tr}(BA)$.

8. **[Elementary] What is the effect of left-multiplying $A$ by an elementary matrix $E$?**
   ??? success "Solution"
       It performs the corresponding elementary row operation on $A$.

9. **[Block] Calculate $\begin{pmatrix} I & A \\ 0 & I \end{pmatrix} \begin{pmatrix} I & -A \\ 0 & I \end{pmatrix}$.**
   ??? success "Solution"
       $\begin{pmatrix} I & -A+A \\ 0 & I \end{pmatrix} = \begin{pmatrix} I & 0 \\ 0 & I \end{pmatrix}$. This shows the inverse is found by negating the top-right block.

10. **[Rank] What is the rank of $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$?**
    ??? success "Solution"
        The rank is 1. It has only one non-zero row.

## Chapter Summary

Matrix operations construct the computational syntax of linear algebra:

1.  **Non-commutativity**: The most fundamental difference between matrix and scalar multiplication, dictating that the order of operators cannot be arbitrarily swapped.
2.  **Structure Preservation**: Transpose and inverse operations maintain internal logical consistency, providing the basis for solving operator equations.
3.  **Computational Simplification**: Special matrices (Identity, Diagonal, Block) greatly simplify the analysis of complex systems and are key entry points for numerical algorithm optimization.
