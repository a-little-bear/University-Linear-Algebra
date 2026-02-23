# Chapter 49: Exterior Algebra

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Determinants (Ch3) · Dual Spaces (Ch13a) · Multilinear Algebra (Ch21)

**Chapter Outline**: Exterior Product $\wedge$ → Antisymmetry and Alternating Tensors → Graded Algebra → Exterior Power $\Lambda^k(V)$ → Relation to Determinants → Plücker Embedding → Grassmannians $\operatorname{Gr}(k, n)$ → Hodge Star Operator → Interior Product

**Extension**: Exterior algebra is the mathematical language of differential forms in calculus on manifolds and electromagnetism (Maxwell's equations).

</div>

Exterior algebra, also known as **Grassmann algebra**, generalizes the concept of cross products and determinants to higher dimensions. By introducing the **exterior product** (wedge product $\wedge$), it provides a coordinate-independent way to describe oriented areas, volumes, and linear subspaces. This algebraic structure is the foundation for modern differential geometry and the study of the Grassmannian manifold.

---

## 49.1 The Exterior Product and Alternating Tensors

!!! definition "Definition 49.1 (Exterior Product)"
    The exterior product $u \wedge v$ of two vectors satisfies the **antisymmetry** property:
    $$u \wedge v = -v \wedge u$$
    Consequently, $v \wedge v = 0$ for any vector $v$.

!!! theorem "Theorem 49.1 (Linear Independence and $\wedge$)"
    A set of vectors $\{v_1, \dots, v_k\}$ is linearly independent if and only if their exterior product is non-zero:
    $$v_1 \wedge v_2 \wedge \dots \wedge v_k 
eq 0$$

---

## Exercises

1. **[Antisymmetry] Compute $(e_1 + e_2) \wedge (e_1 - e_2)$ in $\mathbb{R}^2$.**
   ??? success "Solution"
       Using linearity and antisymmetry: $(e_1 + e_2) \wedge (e_1 - e_2) = e_1 \wedge e_1 - e_1 \wedge e_2 + e_2 \wedge e_1 - e_2 \wedge e_2 = 0 - e_1 \wedge e_2 - e_1 \wedge e_2 - 0 = -2(e_1 \wedge e_2)$.

2. **[Dimension] Calculate the dimension of the $k$-th exterior power $\Lambda^k(V)$ for $\dim V = n$.**
   ??? success "Solution"
       The dimension is given by the binomial coefficient $\binom{n}{k}$. This is the number of ways to choose $k$ basis vectors to form a wedge product.

3. **[Determinant] Show that for a linear operator $T$, the action on the top exterior power $\Lambda^n(V)$ is multiplication by $\det T$.**
   ??? success "Solution"
       Let $\{e_1, \dots, e_n\}$ be a basis. $T(e_1) \wedge \dots \wedge T(e_n) = (\det T) e_1 \wedge \dots \wedge e_n$. This coordinate-independent definition serves as a modern foundation for determinant theory.

4. **[Linear Independence] Prove that $v \wedge w = 0$ implies $v$ and $w$ are linearly dependent.**
   ??? success "Solution"
       If $v, w$ are independent, they can be extended to a basis. Then $v \wedge w$ is a basis element of $\Lambda^2(V)$, which is non-zero. By contrapositive, if $v \wedge w = 0$, they must be dependent.

5. **[Plücker Coordinates] Define the Plücker coordinates of a $k$-dimensional subspace spanned by $\{v_1, \dots, v_k\}$.**
   ??? success "Solution"
       The Plücker coordinates are the coefficients of the wedge product $v_1 \wedge \dots \wedge v_k$ expressed in terms of the basis of $\Lambda^k(V)$. They uniquely identify the subspace up to a scalar multiple.

6. **[Grassmannian] Explain how the Grassmannian $\operatorname{Gr}(k, n)$ is embedded into the projective space $\mathbb{P}(\Lambda^k(V))$.**
   ??? success "Solution"
       This is the Plücker embedding. It maps each $k$-dimensional subspace to the line in $\Lambda^k(V)$ spanned by the wedge product of its basis vectors. The image is a projective algebraic variety defined by the Plücker relations.

7. **[Interior Product] Define the interior product (contraction) $i_v \omega$ and its geometric meaning.**
   ??? success "Solution"
       $i_v$ is an antiderivation of degree -1 that "inserts" a vector into an alternating form. Geometrically, it represents the reduction of a volume form to a lower-dimensional area form by fixing one direction.

8. **[Hodge Star] Describe the Hodge star operator $*: \Lambda^k(V) 	o \Lambda^{n-k}(V)$ in an oriented inner product space.**
   ??? success "Solution"
       The Hodge star maps a $k$-vector to its orthogonal $(n-k)$-complement. It satisfies $\alpha \wedge *\beta = \langle \alpha, \beta angle \omega$, where $\omega$ is the volume form.

9. **[Composition] If $T: V 	o V$ has eigenvalues $\lambda_1, \dots, \lambda_n$, what are the eigenvalues of the induced operator $\Lambda^k(T)$?**
   ??? success "Solution"
       The eigenvalues are all possible products of $k$ distinct eigenvalues of $T$: $\{\lambda_{i_1} \dots \lambda_{i_k} : 1 \le i_1 < \dots < i_k \le n\}$.

10. **[Decomposability] Define a "decomposable" $k$-vector and provide a non-decomposable example in $\Lambda^2(\mathbb{R}^4)$.**
    ??? success "Solution"
        A $k$-vector is decomposable if it can be written as $v_1 \wedge \dots \wedge v_k$. In $\mathbb{R}^4$, $\omega = e_1 \wedge e_2 + e_3 \wedge e_4$ is non-decomposable because $\omega \wedge \omega = 2 e_1 \wedge e_2 \wedge e_3 \wedge e_4 
eq 0$ (Plücker relation violation).

## Chapter Summary

This chapter formalizes the geometry of oriented volumes and subspaces through exterior algebra:

1. **Wedge Calculus**: Defined the exterior product as the fundamental operation for creating alternating multilinear forms.
2. **Structural Duality**: Explored the relationship between subspaces and decomposable $k$-vectors, leading to the Plücker embedding.
3. **Manifold Foundations**: Established the algebraic prerequisites for differential forms, Hodge duals, and the general Stokes' theorem.
4. **Spectral Mapping**: Demonstrated how linear operators lift to exterior powers, providing a deeper perspective on determinants and characteristic polynomials.
