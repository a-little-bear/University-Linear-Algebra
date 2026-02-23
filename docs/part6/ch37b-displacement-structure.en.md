# Chapter 37B: Vandermonde, Cauchy Matrices and Displacement Structure

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Determinants (Ch3) · Decompositions (Ch10) · Numerical Linear Algebra (Ch22) · Toeplitz/Hankel (Ch37A)

**Chapter Outline**: Vandermonde Matrices & Interpolation → Cauchy Matrices & Displacement Rank 1 → Resultant & Krylov Matrices → Sylvester/Stein Displacement Operators → Calculation of Displacement Rank → Generalized Schur Algorithm → Gohberg-Semencul Formula → Hierarchical Matrices

**Extension**: Displacement structure theory (Kailath-Kung-Morf, 1979) unifies the fast algorithm frameworks for all classical structured matrix classes.

</div>

While Ch 37A focused on Toeplitz and Hankel matrices, this chapter explores Vandermonde and Cauchy matrices, characterized by power and kernel structures. We introduce the **Displacement Structure** framework, which reveals that while these matrices aren't low-rank, they become low-rank under specific "displacement operators." This insight allows for $O(n^2)$ or even $O(n \log^2 n)$ solvers.

---

## 37B.1 Core Concepts

!!! definition "Definition 37B.1 (Vandermonde Matrix)"
    Given nodes $x_1, \ldots, x_n \in \mathbb{C}$, the Vandermonde matrix $V$ is defined by $V_{ij} = x_i^{j-1}$. It is non-singular if and only if all $x_i$ are distinct.

!!! definition "Definition 37B.8 (Sylvester Displacement)"
    The displacement of matrix $A$ with respect to operators $(F, F')$ is $
abla_{F,F'}(A) = FA - AF'$. The **displacement rank** is the rank of this resulting matrix.

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
    
eq 1$.

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

Displacement structure is the unifying theory of fast linear algebra:

1. **Low-rank Shadow**: Structured matrices are characterized by having a low-rank image under displacement operators.
2. **Inheritance**: This structure is preserved under inversion and Schur complementation.
3. **Algorithmic Efficiency**: By operating on generators, we break the $O(n^3)$ barrier for wide classes of engineering problems.
