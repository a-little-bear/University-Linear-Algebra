# Chapter 27: Linear Algebra in Graph Theory and Networks

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Eigenvalues (Ch06) · Combinatorial Matrix Structure (Ch65B)

**Chapter Outline**: Mapping from Graphs to Matrices → Adjacency Matrix ($A$) and its Properties → Degree Matrix ($D$) → The Core Object: The Laplacian Matrix $L = D - A$ → The Matrix-Tree Theorem → Introduction to Spectral Graph Theory: Eigenvalues and Connectivity → Algebraic Connectivity (Fiedler Vector) → Incidence Matrix ($M$) → Applications: Google’s PageRank (Eigenvector Centrality), Community Detection (Spectral Clustering), and Resistance Networks

**Extension**: Graph theory is discrete linear algebra; it treats rows of a matrix as nodes and non-zero entries as communication channels. It proves that the topological properties of a network are completely encoded in its spectral structure—the algebraic soul of modern social network analysis and web search.

</div>

Graphs are discrete structures composed of vertices and edges. **Spectral Graph Theory** studies the eigenvalues and eigenvectors of matrices associated with graphs, revealing their connectivity, bottlenecks, and symmetries. By transforming "path-finding" problems into "eigenvalue" problems, linear algebra provides extremely efficient global analysis tools for large-scale complex networks. This chapter introduces how to use matrix algebra to "read" the geometric blueprint of a graph.

---

## 27.1 Core Matrix Definitions

!!! definition "Definition 27.1 (Adjacency Matrix $A$)"
    For a graph $G$ with $n$ vertices, $A$ is an $n \times n$ matrix where:
    $$a_{ij} = \begin{cases} 1 & \text{if vertices } i \text{ and } j \text{ are connected} \\ 0 & \text{otherwise} \end{cases}$$
    **Property**: The $(i,j)$ entry of $A^k$ represents the number of paths of length $k$ from $i$ to $j$.

!!! definition "Definition 27.2 (Laplacian Matrix $L$)"
    Defined as $L = D - A$, where $D = \operatorname{diag}(\text{degrees of vertices})$.
    **Property**: $L$ is always positive semi-definite, and the multiplicity of its zero eigenvalue equals the number of connected components in the graph.

---

## 27.2 The Matrix-Tree Theorem

!!! theorem "Theorem 27.1 (Matrix-Tree Theorem)"
    The number of spanning trees of a graph $G$ is equal to the determinant of any cofactor of the Laplacian matrix $L$.
    **Significance**: This result enables the transition from massive combinatorial enumeration to simple determinant calculation.

---

## 27.3 Spectral Clustering and Partitioning

!!! technique "Technique: The Fiedler Vector"
    The second smallest eigenvalue $\lambda_2$ of the Laplacian matrix is known as the **algebraic connectivity**. Its corresponding eigenvector (the Fiedler vector) can be used to partition a graph into two communities by the signs of its components. This is the standard algorithm for image segmentation and social group detection.

---

## Exercises

**1. [Basics] Write the adjacency matrix $A$ for a triangle graph (3 vertices, all pairwise connected).**

??? success "Solution"
    **Construction:**
    All off-diagonal entries are 1, diagonal is 0.
    $A = \begin{pmatrix} 0 & 1 & 1 \\ 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix}$.

**2. [Calculation] Compute the Laplacian matrix $L$ for the triangle graph.**

??? success "Solution"
    **Steps:**
    1. Degree matrix: Each vertex has degree 2. $D = \operatorname{diag}(2, 2, 2)$.
    2. $L = D - A = \begin{pmatrix} 2 & -1 & -1 \\ -1 & 2 & -1 \\ -1 & -1 & 2 \end{pmatrix}$.
    **Verification**: The sum of each row is 0, a universal trait of Laplacian matrices.

**3. [Path Counting] If the $(1, 3)$ entry of $A^2$ is 2, what does this imply?**

??? success "Solution"
    **Conclusion:**
    It means there are **exactly 2 paths of length 2** from vertex 1 to vertex 3. For example, $1 \to 2 \to 3$ and $1 \to 4 \to 3$.

**4. [Matrix-Tree] Use the Matrix-Tree Theorem to find the number of spanning trees for the triangle graph.**

??? success "Solution"
    **Steps:**
    1. Take a principal minor of $L$ (delete row 1 and column 1): $\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$.
    2. Calculate the determinant: $2(2) - (-1)(-1) = 4 - 1 = 3$.
    **Conclusion**: There are 3 spanning trees. In a complete 3-node graph, removing any of the 3 edges results in a spanning tree.

**5. [Eigenvalues] Prove that the smallest eigenvalue of the Laplacian $L$ is always 0.**

??? success "Solution"
    **Proof:**
    1. Observe that each row of $L$ sums to 0.
    2. Let $\mathbf{1} = (1, 1, \ldots, 1)^T$.
    3. Then $L\mathbf{1} = \mathbf{0} = 0 \cdot \mathbf{1}$.
    **Conclusion**: $\mathbf{1}$ is the eigenvector corresponding to $\lambda = 0$.

**6. [Connectivity] If $L$ has two zero eigenvalues, what is the topology of the graph?**

??? success "Solution"
    **Conclusion: The graph consists of two disconnected subgraphs.**
    **Reasoning**: The geometric multiplicity of eigenvalue 0 equals the number of connected components. This implies the system has two independent "zero-energy" modes that do not interfere.

**7. [Bipartite] What symmetry exists in the spectrum of the adjacency matrix of a bipartite graph?**

??? success "Solution"
    **Conclusion: The spectrum is symmetric about the origin.**
    **Reasoning**: If $\lambda$ is an eigenvalue, then $-\lambda$ is also an eigenvalue. This occurs because the adjacency matrix of a bipartite graph can be permuted into a block anti-diagonal form $\begin{pmatrix} 0 & B \\ B^T & 0 \end{pmatrix}$.

**8. [Application] Is the core of Google's PageRank an eigenvalue problem?**

??? success "Solution"
    **Yes.**
    PageRank seeks the **principal eigenvector** (corresponding to $\lambda=1$) of a stochastic transition matrix $P$ constructed from the web link graph. The magnitudes of the eigenvector components $x_i$ represent the importance rankings of the pages.

**9. [Incidence] What is the incidence matrix $M$, and how does it relate to $L$?**

??? success "Solution"
    **Definition and Relation:**
    1. $M$ has rows for vertices and columns for edges. For a directed edge $k$ from $i$ to $j$, $M_{ik}=1$ and $M_{jk}=-1$.
    2. **Key Identity**: $L = MM^T$.
    This proves that the Laplacian is essentially a Gram matrix, explaining its positive semi-definiteness.

**10. [Application] Briefly describe the idea behind Spectral Clustering in image processing.**

??? success "Solution"
    1. Treat pixels as nodes and pixel similarity as edge weights.
    2. Construct the weighted Laplacian matrix.
    3. Compute the eigenvectors corresponding to the $k$ smallest eigenvalues.
    4. Use these eigenvectors as new coordinates to cluster the pixels.
    **Advantage**: It identifies complex, non-convex shapes that traditional K-means cannot handle.

## Chapter Summary

Linear algebra is the "universal microscope" for parsing network topology:

1.  **Spectral Representation of Topology**: It proves that macro-connectivity, clustering, and tree structures can be concentrated into micro-eigenvalue sequences, uniting discrete geometry and continuous algebra.
2.  **Algebraic Essence of Flow**: From PageRank to resistor networks, linear algebra reveals that diffusion and equilibrium on networks are essentially energy minimization under operator action.
3.  **Computational Leap**: Results like the Matrix-Tree theorem demonstrate how algebraic forms can transform impossible combinatorial tasks into standard matrix operations—the mathematical bedrock of modern search and social analysis.
