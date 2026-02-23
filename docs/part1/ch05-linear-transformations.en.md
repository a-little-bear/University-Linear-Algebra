# Chapter 05: Linear Transformations

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04)

**Chapter Outline**: Definition of Linear Transformations → Linearity Criteria → Matrix Representation (Choice of Bases) → Kernel and Range → Injectivity, Surjectivity, and Isomorphisms → Operations on Transformations (Sum, Scalar Multiple, Composition) → Inverse Transformations → Matrix Representation under Change of Basis (Similarity) → Geometric Operators (Rotation, Reflection, Projection, Shear)

**Extension**: Linear transformations connect static vector spaces and are the concrete manifestation of "morphisms" in the category of vector spaces; they are the core for understanding Diagonalization (Ch06) and SVD (Ch11).

</div>

Vector spaces provide the stage, and linear transformations are the primary actors. A linear transformation is a mapping that preserves the addition and scalar multiplication structures of a vector space, allowing us to study the interconnections between different spaces. This chapter reveals a central fact: once bases are chosen, every linear transformation can be uniquely represented by a matrix.

---

## 05.1 Definition of Linear Transformations

!!! definition "Definition 05.1 (Linear Transformation)"
    A mapping $T: V \to W$ is called a **linear transformation** if it satisfies:
    1.  **Additivity**: $T(\mathbf{u} + \mathbf{v}) = T(\mathbf{u}) + T(\mathbf{v})$.
    2.  **Homogeneity**: $T(c\mathbf{v}) = cT(\mathbf{v})$.

!!! note "Properties"
    A linear transformation always maps the zero vector to the zero vector: $T(\mathbf{0}_V) = \mathbf{0}_W$.

---

## 05.2 Matrix Representation

!!! theorem "Theorem 05.1 (Matrix Representation)"
    Let $B = \{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$ be a basis for $V$ and $C$ be a basis for $W$. $T$ is completely determined by its action on the basis vectors $T(\mathbf{v}_j)$.
    The **matrix representation** of $T$ is $[T]_{C \leftarrow B} = [ [T(\mathbf{v}_1)]_C \ \cdots \ [T(\mathbf{v}_n)]_C ]$.
    Then: $[T(\mathbf{x})]_C = [T]_{C \leftarrow B} [\mathbf{x}]_B$.

---

## 05.3 Kernel and Range

!!! definition "Definition 05.2 (Kernel and Range)"
    1.  **Kernel $\ker(T)$**: The set of all vectors mapped to the zero vector: $\{\mathbf{v} \in V : T(\mathbf{v}) = \mathbf{0}\}$.
    2.  **Range $\operatorname{im}(T)$**: The set of all images of vectors in $V$: $\{T(\mathbf{v}) : \mathbf{v} \in V\}$.

!!! theorem "Theorem 05.2 (Rank-Nullity Theorem for Transformations)"
    $\dim \ker(T) + \dim \operatorname{im}(T) = \dim V$. This is essentially identical to the rank-nullity theorem for matrices.

---

## 05.4 Isomorphisms and Inverses

!!! definition "Definition 05.3 (Isomorphism)"
    If $T: V \to W$ is both injective (one-to-one) and surjective (onto), it is an **isomorphism**. In this case, $V$ and $W$ are algebraically equivalent.
    **Corollary**: Two finite-dimensional vector spaces are isomorphic if and only if they have the same dimension.

---

## Exercises

1. **[Criteria] Determine if $T(x, y) = (x+1, y)$ is a linear transformation.**

   ??? success "Solution"
       No. $T(0, 0) = (1, 0) \neq (0, 0)$. A linear transformation must preserve the zero vector.

2. **[Matrix] Let the matrix of $T$ in the standard basis of $\mathbb{R}^2$ be $\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$. Describe its geometry.**

   ??? success "Solution"
       It is a counter-clockwise rotation of 90 degrees about the origin.

3. **[Kernel] Find the kernel of the differentiation operator $D: P_2 \to P_1$ where $D(p) = p'$.**

   ??? success "Solution"
       $\ker(D)$ is the set of polynomials with zero derivative, which is the space of constant polynomials $\operatorname{span}\{1\}$. Its dimension is 1.

4. **[Range] Find the range of $D$ from the previous exercise.**

   ??? success "Solution"
       The basis $\{1, x, x^2\}$ is mapped to $\{0, 1, 2x\}$. Thus the range is $\operatorname{span}\{1, x\} = P_1$.

5. **[Composition] Prove the composition of two linear transformations is linear.**

   ??? success "Solution"
       $S(T(\mathbf{u}+\mathbf{v})) = S(T\mathbf{u} + T\mathbf{v}) = S(T\mathbf{u}) + S(T\mathbf{v})$. Homogeneity follows similarly.

6. **[Inversion] If the matrix of $T$ is $A$, what is the condition for $T^{-1}$ to exist?**

   ??? success "Solution"
       $A$ must be an invertible matrix (i.e., $\det A \neq 0$).

7. **[Projection] Write the matrix for the orthogonal projection onto the $x$-axis in $\mathbb{R}^2$.**

   ??? success "Solution"
       $\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$.

8. **[Isomorphism] Are $\mathbb{R}^n$ and $P_{n-1}$ isomorphic?**

   ??? success "Solution"
       Yes, because they both have dimension $n$.

9. **[Similarity] If $A$ and $B$ represent the same transformation in different bases, how are they related?**

   ??? success "Solution"
       $B = P^{-1}AP$; they are similar matrices.

10. **[Trace] Prove the trace of a linear transformation is independent of the basis.**

   ??? success "Solution"
        Trace is a similarity invariant. Since matrices in different bases are similar, their traces must be equal.

## Chapter Summary

Linear transformations build "bridges" between spaces:

1.  **Structure Preservation**: The essence of a linear transformation is maintaining the harmony of vector addition and scalar multiplication.
2.  **Matrix Realization**: Choosing a basis translates abstract mappings into matrix arithmetic, enabling the use of numerical algorithms for abstract logic.
3.  **Dimension Conservation**: The relationship between kernel and range (Rank-Nullity Theorem) reveals the patterns of information retention and loss during transformation.
