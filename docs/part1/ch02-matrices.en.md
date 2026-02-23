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

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

Matrix operations construct the computational syntax of linear algebra:

1.  **Non-commutativity**: The most fundamental difference between matrix and scalar multiplication, dictating that the order of operators cannot be arbitrarily swapped.
2.  **Structure Preservation**: Transpose and inverse operations maintain internal logical consistency, providing the basis for solving operator equations.
3.  **Computational Simplification**: Special matrices (Identity, Diagonal, Block) greatly simplify the analysis of complex systems and are key entry points for numerical algorithm optimization.
