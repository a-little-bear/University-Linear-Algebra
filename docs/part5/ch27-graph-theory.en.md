# Chapter 27: Linear Algebra in Graph Theory and Networks

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Spectral Theorem (Ch08) · Non-negative Matrices (Ch17)

**Chapter Outline**: Matrix Representations of Graphs (Adjacency, Incidence, Laplacian) → Spectral Graph Theory → Laplacian and Connectivity (Algebraic Connectivity) → Kirchhoff's Matrix-Tree Theorem → Random Walks & PageRank → Cheeger Inequality & Graph Partitioning → Expander Graphs → Graph Coloring & Eigenvalues → Network Flows & Linear Programming

**Extension**: Graph theory and linear algebra intersect in Spectral Graph Theory, which uses matrix properties to reveal global topological features of a graph; PageRank is the most famous industrial application of this synergy.

</div>

Graphs are discrete structures consisting of vertices and edges, yet they can be completely encoded into matrices. By translating combinatorial problems into eigenvalue problems, we can "listen" to the shape of a network. This chapter explores how the spectrum of a graph dictates its connectivity, clusterability, and information flow.

---

## 27.1 Matrix Representations of Graphs

<div class="context-flow" markdown>

**Three Pillars**: The Adjacency matrix $A$ (neighbor relations), the Incidence matrix $B$ (vertex-edge relations), and the Laplacian $L = D - A$ (the most vital operator).

</div>

!!! definition "Definition 27.1 (Adjacency Matrix)"
    For a graph $G$ with $n$ vertices, the **Adjacency Matrix** $A$ is an $n \times n$ symmetric matrix where $A_{ij} = 1$ if vertices $i$ and $j$ are connected, and 0 otherwise.

!!! definition "Definition 27.2 (Laplacian Matrix)"
    The **Laplacian Matrix** $L$ is defined as $L = D - A$, where $D = \operatorname{diag}(d_1, \ldots, d_n)$ is the degree matrix. 
    **Key Property**: $L$ is positive semi-definite, and $L\mathbf{1} = \mathbf{0}$.

---

## 27.2 Spectral Graph Theory

<div class="context-flow" markdown>

**The Spectrum**: The set of eigenvalues of $A$ or $L$ carries information about the graph's regularity, diameter, and bipartite nature.

</div>

!!! theorem "Theorem 27.1 (Spectral Properties)"
    1.  The number of edges $|E| = \frac{1}{2} \operatorname{tr}(A^2)$.
    2.  The number of triangles in $G$ is $\frac{1}{6} \operatorname{tr}(A^3)$.
    3.  $G$ is bipartite iff its adjacency spectrum is symmetric about 0.

---

## 27.3 Laplacian and Connectivity

!!! theorem "Theorem 27.2 (Connectivity)"
    The multiplicity of the eigenvalue 0 in the Laplacian $L$ equals the number of connected components in the graph.
    - $G$ is connected iff $\lambda_2(L) > 0$.
    - $\lambda_2(L)$ is called the **Algebraic Connectivity** or the **Fiedler value**.

!!! theorem "Theorem 27.3 (Kirchhoff's Matrix-Tree Theorem)"
    The number of spanning trees in a graph $G$ is equal to any cofactor of the Laplacian matrix $L$. For a connected graph, this is $\frac{1}{n} \lambda_2 \lambda_3 \cdots \lambda_n$.

---

## 27.4 Random Walks and PageRank

<div class="context-flow" markdown>

**Markovian Flow**: A random walk on a graph is a Markov chain where the transition matrix is $P = D^{-1}A$.

</div>

!!! algorithm "Algorithm 27.1 (PageRank)"
    The PageRank of a web page is determined by the dominant eigenvector of the **Google Matrix**:
    $$G = \alpha P + (1-\alpha) \frac{1}{n} \mathbf{1}\mathbf{1}^T$$
    where $\alpha$ is the damping factor (usually 0.85). This ensures the matrix is strictly positive, guaranteeing a unique steady-state distribution by the Perron-Frobenius theorem.

---

## 27.5 Graph Partitioning and Cheeger Inequality

!!! theorem "Theorem 27.4 (Cheeger Inequality)"
    The spectral gap $\lambda_2$ of the normalized Laplacian provides bounds on the **conductance** $h(G)$ (the cost of the best cut):
    $$\frac{h(G)^2}{2} \le \lambda_2 \le 2h(G)$$
    This justifies **Spectral Clustering**, where we use the Fiedler vector to partition a network into communities.

---

## Exercises


****
??? success "Solution"
     It means each vertex has exactly one neighbor, so the graph is a collection of disjoint edges (a perfect matching).


****
??? success "Solution"
     $\operatorname{tr}(A^3) = 2^3 + 0^3 + (-2)^3 = 8 - 8 = 0$. The graph has zero triangles (it is a bipartite graph).


****
??? success "Solution"
     $D = \operatorname{diag}(2, 2, 2)$, $A = \begin{pmatrix} 0 & 1 & 1 \\ 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix} \implies L = \begin{pmatrix} 2 & -1 & -1 \\ -1 & 2 & -1 \\ -1 & -1 & 2 \end{pmatrix}$.


****
??? success "Solution"
     Because each row of $L$ sums to $d_i - \sum A_{ij} = d_i - d_i = 0$. Thus $L\mathbf{1} = \mathbf{0}$.


****
??? success "Solution"
     The multiplicity of 0 is 2, so the graph has 2 connected components.


****
??? success "Solution"
     To ensure the transition matrix is primitive and irreducible (it makes the graph strongly connected), allowing the power method to converge to a unique solution.


****
??? success "Solution"
     $L$ eigenvalues are $0, 3, 3$. Number of trees $= \frac{1}{3}(3 \cdot 3) = 3$.


****
??? success "Solution"
     Yes, the spectrum is perfectly symmetric about 0.


****
??? success "Solution"
     If $f$ is a vector of flows on edges, then $Bf$ is a vector of net flows at each vertex. Flow conservation (Kirchhoff's Current Law) is expressed as $Bf = 0$ for all internal nodes.

****
??? success "Solution"
    ## Chapter Summary

Linear algebra provides the "X-ray" for complex networks:


****: Transformed combinatorial objects (nodes/edges) into analytic objects (matrices), enabling the use of spectral tools.

****: Established that the eigenvalues of a graph are not just numbers but descriptors of connectivity, bipartite structure, and density.

****: Linked random walks and diffuse processes to the steady-state properties of stochastic matrices.

****: Used the Cheeger inequality to prove that the "physics" of the graph (eigenvalues) can solve the "logic" of the graph (min-cut problems).
