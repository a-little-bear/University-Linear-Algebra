# Chapter 44: Weyr Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Jordan Form (Ch12) · Eigenvalues (Ch6)

**Chapter Outline**: Limitations of Jordan Form → Weyr Characteristics → Definition of Weyr Form → Existence and Uniqueness → Duality with Jordan Form → Commuting Matrices → Centralizer Algebras

**Extension**: The Weyr form is more natural than the Jordan form for studying commuting families of matrices and matrix equations like $AX=XB$.

</div>

While the Jordan form is the most famous matrix canonical form, it is often cumbersome when studying commuting matrices. The **Weyr Canonical Form** (named after Eduard Weyr, 1885) provides an elegant alternative. It carries the same structural information as the Jordan form but organizes it "layer by layer" rather than "chain by chain," leading to a simpler description of the centralizer algebra.

---

## 44.1 Weyr Characteristics

!!! definition "Definition 44.2 (Weyr Characteristic)"
    For an eigenvalue $\lambda$, the **Weyr characteristic** $(w_1, w_2, \dots, w_s)$ is the sequence of dimensions of the increments of the kernels:
    $$w_j = \dim \ker(A - \lambda I)^j - \dim \ker(A - \lambda I)^{j-1}.$$
    It is the conjugate partition (transpose of the Young diagram) of the Jordan partition.

!!! theorem "Theorem 44.7 (Commuting Matrices)"
    A matrix $X$ commutes with a Weyr matrix $W$ if and only if $X$ is block upper triangular with respect to the Weyr partition, where the diagonal blocks are arbitrary.

---

## Exercises

****

??? success "Solution"
    
       $w_1 = |\{3, 2, 1\} \ge 1| = 3$.
       $w_2 = |\{3, 2, 1\} \ge 2| = 2$.
       $w_3 = |\{3, 2, 1\} \ge 3| = 1$.
       The Weyr characteristic is $(3, 2, 1)$ (Self-conjugate).

****

??? success "Solution"
    
       $n_1 = |\{4, 2, 2, 1\} \ge 1| = 4$.
       $n_2 = |\{4, 2, 2, 1\} \ge 2| = 3$.
       $n_3 = |\{4, 2, 2, 1\} \ge 3| = 1$.
       $n_4 = |\{4, 2, 2, 1\} \ge 4| = 1$.
       Jordan partition: $(4, 3, 1, 1)$.

****

??? success "Solution"
    
       $W = \begin{pmatrix} 0 & 0 & 1 \ 0 & 0 & 0 \ 0 & 0 & 0 \end{pmatrix}$.
       Here $\lambda I_2 = \begin{pmatrix} 0 & 0 \ 0 & 0 \end{pmatrix}$ and $F_1 = \begin{pmatrix} 1 \ 0 \end{pmatrix}$.

****

??? success "Solution"
    
       $\dim \mathcal{C}(A) = \sum w_j^2 = 2^2 + 2^2 + 1^2 = 4 + 4 + 1 = 9$.

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    
       $w_2 = 5 - 3 = 2$.
       $w_3 = 6 - 5 = 1$.
       Weyr characteristic is $(3, 2, 1)$.

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

The Weyr Canonical Form provides a powerful dual perspective to the Jordan form:

1. **Layered Organization**: Shifted from the "chain-based" view of Jordan to a "kernel-increment" view.
2. **Algebraic Simplicity**: Simplified the study of commuting matrices by inducing a block upper triangular structure on the centralizer.
3. **Partition Duality**: Leveraged the combinatorial property of conjugate partitions to classify similarity classes.
4. **Computational Advantage**: Provided a direct formula ($\sum w_j^2$) for the dimension of the space of commuting matrices.
