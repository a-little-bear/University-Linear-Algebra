# Chapter 19: Kronecker Product and Vec Operator

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Linear Transformations (Ch5) · Multilinear Algebra Intro (Ch21)

**Chapter Outline**: Definition of Kronecker Product $A \otimes B$ → Algebraic Properties (Distributivity, Associativity, Transpose, Inverse) → Spectral Properties (Eigenvalues and Singular Values) → Definition of Vec Operator → Core Identity $\operatorname{vec}(AXB) = (B^T \otimes A)\operatorname{vec}(X)$ → Applications (Matrix Equations, System Identification, Tensor Computation)

**Extension**: Kronecker product is a "linearization" tool for matrix equations, transforming complex operator operations into standard vector mappings.

</div>

The Kronecker product (a manifestation of the tensor product on matrices) provides a systematic way to combine two small matrices into one large matrix. Combined with the Vec operator, it can transform linear equations in matrix form into standard vector forms, allowing for direct application of existing linear system theory.

---

## 19.1 Definitions and Core Identities

!!! definition "Definition 19.1 (Kronecker Product)"
    Let $A$ be an $m \times n$ matrix and $B$ be a $p \times q$ matrix. $A \otimes B$ is an $mp \times nq$ block matrix:
    $$A \otimes B = \begin{pmatrix} a_{11}B & \dots & a_{1n}B \\ \vdots & \ddots & \dots \\ a_{m1}B & \dots & a_{mn}B \end{pmatrix}$$

!!! theorem "Theorem 19.3 (Vec Operator Identity)"
    For the matrix equation $AXB = C$, its vectorized form is:
    $$(B^T \otimes A) \operatorname{vec}(X) = \operatorname{vec}(C)$$

---

## Exercises

1. **[Basic Calculation] Calculate $\begin{pmatrix} 1 & 2 \end{pmatrix} \otimes \begin{pmatrix} 3 \\ 4 \end{pmatrix}$.**
   ??? success "Solution"
       According to the definition: $\begin{pmatrix} 1\begin{pmatrix} 3 \\ 4 \end{pmatrix} & 2\begin{pmatrix} 3 \\ 4 \end{pmatrix} \end{pmatrix} = \begin{pmatrix} 3 & 6 \\ 4 & 8 \end{pmatrix}$.

2. **[Determinant] If $A$ is $n \times n$ and $B$ is $m \times m$, what is the formula for $\det(A \otimes B)$?**
   ??? success "Solution"
       $\det(A \otimes B) = (\det A)^m (\det B)^n$.

3. **[Eigenvalues] If $A$ has eigenvalues $\lambda_i$ and $B$ has eigenvalues $\mu_j$, prove $A \otimes B$ has eigenvalues $\lambda_i \mu_j$.**
   ??? success "Solution"
       Let $Ax = \lambda x$ and $By = \mu y$.
       Then $(A \otimes B)(x \otimes y) = Ax \otimes By = (\lambda x) \otimes (\mu y) = (\lambda \mu)(x \otimes y)$.
       Since there are $nm$ such combinations, they form the complete spectrum of $A \otimes B$.

4. **[Inverse] Prove: $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$ (assuming inverses exist).**
   ??? success "Solution"
       $(A \otimes B)(A^{-1} \otimes B^{-1}) = (AA^{-1}) \otimes (BB^{-1}) = I \otimes I = I$.

5. **[Vec Application] Transform the matrix equation $AX + XA^T = C$ into vector form.**
   ??? success "Solution"
       $\operatorname{vec}(AX) + \operatorname{vec}(XA^T) = \operatorname{vec}(C)$.
       Using the identity: $(I \otimes A) \operatorname{vec}(X) + (A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(C)$.
       Thus $((I \otimes A) + (A \otimes I)) \operatorname{vec}(X) = \operatorname{vec}(C)$. This is the linearized form of the famous Lyapunov equation.

6. **[Rank] Prove $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A) \operatorname{rank}(B)$.**
   ??? success "Solution"
       Using SVD: the singular values of $A \otimes B$ are the products of the singular values of $A$ and $B$. The number of non-zero singular values is clearly the product of the ranks.

7. **[Trace] Prove $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \operatorname{tr}(B)$.**
   ??? success "Solution"
       $\operatorname{tr}(A \otimes B) = \sum \lambda_{ij} = \sum_i \sum_j \lambda_i \mu_j = (\sum \lambda_i) (\sum \mu_j) = \operatorname{tr}(A) \operatorname{tr}(B)$.

8. **[Commutation Matrix] Does there exist a permutation matrix $K$ such that $B \otimes A = K (A \otimes B) K^T$?**
   ??? success "Solution"
       Yes. This matrix $K$ is called the **Commutation Matrix**. It implements the non-commutative swap of the Kronecker product by rearranging entries.

9. **[Mixed-Product Property] Prove $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$.**
   ??? success "Solution"
       This is the most essential algebraic property of the Kronecker product, which can be verified directly via block matrix multiplication.

10. **[Application] Why is the Kronecker product frequently used in quantum information?**
    ??? success "Solution"
        Because the state space (Hilbert space) of a composite quantum system is formed by the tensor product of the subspaces. If system A has $n$ states and system B has $m$ states, the composite system A+B is described by an $n \times m$ dimensional vector (or a Kronecker product of density matrices).

## Chapter Summary

Kronecker product is the bridge for high-dimensional linear algebra:

1. **Dimensional Explosion and Reduction**: It provides a systematic way to increase dimensionality while maintaining linearity via the Vec operator.
2. **Spectral Inheritance**: The product structure of eigenvalues and singular values reveals the physical essence of composite systems.
3. **Equation Solving**: It is the standard route for transforming complex operator equations (like Lyapunov, Sylvester) into standard linear problems.
