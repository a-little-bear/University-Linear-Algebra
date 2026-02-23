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

Linear transformations build "bridges" between spaces:

1.  **Structure Preservation**: The essence of a linear transformation is maintaining the harmony of vector addition and scalar multiplication.
2.  **Matrix Realization**: Choosing a basis translates abstract mappings into matrix arithmetic, enabling the use of numerical algorithms for abstract logic.
3.  **Dimension Conservation**: The relationship between kernel and range (Rank-Nullity Theorem) reveals the patterns of information retention and loss during transformation.
