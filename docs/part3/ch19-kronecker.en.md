# Chapter 19: The Kronecker Product

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Vector Spaces (Ch4) · Multilinear Algebra (Ch21) · Matrix Equations (Ch20)

**Chapter Outline**: Definition of Kronecker Product $A \otimes B$ → Mixed-product Property → Eigenvalues and Eigenvectors of $A \otimes B$ → Vec Operator and the Identity $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ → Kronecker Sum → Application in Sylvester and Lyapunov Equations → Tensor Products of Vector Spaces

**Extension**: The Kronecker product is the matrix representation of the tensor product; it is the essential tool for linearizing matrix equations and describing multi-particle systems in quantum mechanics.

</div>

The **Kronecker product** (or tensor product of matrices) $A \otimes B$ is a block matrix that captures all possible products between the entries of $A$ and $B$. While standard multiplication represents the composition of maps, the Kronecker product represents the action of two independent maps on a joint space. Its most powerful application is the **vec-Kronecker identity**, which converts matrix equations like $AX - XB = C$ into standard vector equations. This chapter details the properties of this product and its role in solving high-dimensional linear systems.

---

## 19.1 Definition and the Vec-Kronecker Identity

!!! definition "Definition 19.1 (Kronecker Product)"
    If $A \in M_{m \times n}$ and $B \in M_{p \times q}$, the Kronecker product $A \otimes B$ is the $mp \times nq$ block matrix:
    $$A \otimes B = \begin{pmatrix} a_{11}B & \dots & a_{1n}B \\ \vdots & \ddots & \vdots \\ a_{m1}B & \dots & a_{mn}B \end{pmatrix}$$

!!! theorem "Theorem 19.1 (The Vec Identity)"
    Let $\operatorname{vec}(X)$ be the vectorization of $X$ (stacking columns). Then:
    $$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$$
    This identity allows for "unrolling" matrix equations into linear systems.

---

## Exercises

1. **[Fundamentals] Compute $I_2 \otimes \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.**
   ??? success "Solution"
       $I_2 \otimes A = \begin{pmatrix} A & 0 \\ 0 & A \end{pmatrix} = \begin{pmatrix} 1 & 2 & 0 & 0 \\ 3 & 4 & 0 & 0 \\ 0 & 0 & 1 & 2 \\ 0 & 0 & 3 & 4 \end{pmatrix}$.

2. **[Eigenvalues] If $\sigma(A) = \{\lambda_i\}$ and $\sigma(B) = \{\mu_j\}$, what are the eigenvalues of $A \otimes B$?**
   ??? success "Solution"
       The eigenvalues are all possible products $\{\lambda_i \mu_j\}$. The eigenvectors are the Kronecker products of the individual eigenvectors $v_i \otimes w_j$.

3. **[Inversion] Prove that $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$.**
   ??? success "Solution"
       Using the mixed-product property $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$: $(A \otimes B)(A^{-1} \otimes B^{-1}) = (AA^{-1}) \otimes (BB^{-1}) = I \otimes I = I$.

4. **[Determinant] Show that $\det(A \otimes B) = (\det A)^n (\det B)^m$ for $A_{m \times m}$ and $B_{n \times n}$.**
   ??? success "Solution"
       Since the eigenvalues are $\{\lambda_i \mu_j\}$, the determinant is the product: $\prod_{i,j} \lambda_i \mu_j = (\prod \lambda_i)^n (\prod \mu_j)^m = (\det A)^n (\det B)^m$.

5. **[Symmetry] If $A$ and $B$ are symmetric, is $A \otimes B$ symmetric?**
   ??? success "Solution"
       Yes. $(A \otimes B)^T = A^T \otimes B^T = A \otimes B$. The Kronecker product preserves symmetry and positive definiteness.

6. **[Vec Operator] Express the equation $AX + XB = C$ as a vector equation.**
   ??? success "Solution"
       $\operatorname{vec}(AXI) + \operatorname{vec}(IXB) = \operatorname{vec}(C) \implies (I \otimes A + B^T \otimes I) \operatorname{vec}(X) = \operatorname{vec}(C)$. The matrix $(I \otimes A + B^T \otimes I)$ is called the **Kronecker sum**.

7. **[Trace] Prove $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \operatorname{tr}(B)$.**
   ??? success "Solution"
       The diagonal entries are $a_{ii}b_{jj}$. The sum is $\sum_i \sum_j a_{ii}b_{jj} = (\sum a_{ii})(\sum b_{jj}) = \operatorname{tr}(A) \operatorname{tr}(B)$.

8. **[Commutativity] Is $A \otimes B = B \otimes A$?**
   ??? success "Solution"
       No, but they are **permutation similar**. There exists a commutation matrix $P$ such that $B \otimes A = P (A \otimes B) P^T$.

9. **[Rank] What is $\operatorname{rank}(A \otimes B)$?**
   ??? success "Solution"
       $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A) \cdot \operatorname{rank}(B)$.

10. **[Quantum] How is the Kronecker product used to describe two-qubit systems?**
    ??? success "Solution"
        The state space of two qubits is the tensor product of their individual Hilbert spaces. If the qubits are in states $|\psi_1\rangle$ and $|\psi_2\rangle$, the joint state is $|\psi_1\rangle \otimes |\psi_2\rangle$. Entangled states are those that *cannot* be written as a Kronecker product.

## Chapter Summary

This chapter explores the arithmetic of joint operators:

1. **Multilinear Scaling**: Defined the Kronecker product as the matrix representation of tensor-product actions.
2. **Linearization Power**: Leveraged the vec-Kronecker identity to reduce matrix-valued equations to standard linear systems.
3. **Spectral Composition**: Established the eigenvalue-multiplication law, providing the spectral theory for multi-particle and joint systems.
4. **Structural Preservation**: Validated the preservation of symmetry, positivity, and rank under Kronecker multiplication.
