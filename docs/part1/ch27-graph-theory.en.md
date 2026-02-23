# Chapter 27: Applications of Linear Algebra in Graph Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Eigenvalues (Ch6) · Vector Spaces (Ch4)

**Chapter Outline**: Matrix Representations of Graphs (Adjacency, Incidence) → Graph Laplacian $L = D - A$ → Matrix Tree Theorem → Spectral Graph Theory (Spectral Radius, Algebraic Connectivity) → Random Walks and Eigenvalues → Network Flows → Applications (Centrality, Image Segmentation, Community Detection)

**Extension**: Spectral graph theory bridges discrete topology and continuous algebra; it provides insights into network connectivity and vulnerability through eigenvalue analysis.

</div>

Graph theory studies vertices and their relations. Linear algebra provides powerful analytical tools for complex topologies by "matrixizing" graphs. The spectrum (set of eigenvalues) of these matrices acts like a fingerprint, encoding the key structural information of the graph.

---

## 27.1 Core Matrices and Theorems

!!! definition "Definition 27.1 (Graph Laplacian)"
    For a simple undirected graph $G$, the Laplacian matrix is $L = D - A$, where $D$ is the degree matrix (diagonal) and $A$ is the adjacency matrix.

!!! theorem "Theorem 27.3 (Matrix Tree Theorem)"
    The number of spanning trees of graph $G$ equals the value of any cofactor of its Laplacian matrix $L$.

---

## Exercises

1. **[Fundamentals] Write the adjacency matrix $A$ for the triangle graph $K_3$.**
   ??? success "Solution"
       $A = \begin{pmatrix} 0 & 1 & 1 \ 1 & 0 & 1 \ 1 & 1 & 0 \end{pmatrix}$. Note the main diagonal is all zeros.

2. **[Laplacian Property] Prove: Every row sum of the Laplacian matrix $L$ is 0.**
   ??? success "Solution"
       $L_{ii} = d_i$ (node degree) and $L_{ij} = -1$ if $i, j$ are connected.
       Row sum $= d_i + \sum_{j 
eq i} (-1) = d_i - d_i = 0$.
       This means $\mathbf{1} = (1, \dots, 1)^T$ is always an eigenvector of $L$ for eigenvalue 0.

3. **[Connectivity] What does the multiplicity of eigenvalue 0 in the Laplacian represent?**
   ??? success "Solution"
       It equals the number of **connected components** of the graph. A connected graph has exactly one 0 eigenvalue.

4. **[Matrix Tree Theorem Calculation] Calculate the number of spanning trees for a path graph $P_3$ with 3 vertices.**
   ??? success "Solution"
       $L = \begin{pmatrix} 1 & -1 & 0 \ -1 & 2 & -1 \ 0 & -1 & 1 \end{pmatrix}$.
       Removing the last row and column, the minor is $\begin{vmatrix} 1 & -1 \ -1 & 2 \end{vmatrix} = 2 - 1 = 1$.
       The number of spanning trees is 1 (the path itself).

5. **[Algebraic Connectivity] What is the Fiedler vector? How is it used in community detection?**
   ??? success "Solution"
       The Fiedler vector is the eigenvector corresponding to the second-smallest eigenvalue $\lambda_2$ of the Laplacian. Since $\lambda_2$ measures the expander property, the signs of the Fiedler vector components can partition the graph into two loosely coupled parts (spectral clustering).

6. **[Incidence Matrix] Prove the undirected Laplacian satisfies $L = M M^T$, where $M$ is the incidence matrix.**
   ??? success "Solution"
       By expanding the terms of $M M^T$, one finds the diagonal entries are vertex degrees and off-diagonals correspond to negatives of the adjacency matrix. This decomposition proves $L$ is positive semi-definite.

7. **[Regular Graphs] If $G$ is a $k$-regular graph, what is the largest eigenvalue of its adjacency matrix $A$?**
   ??? success "Solution"
       The largest eigenvalue is $k$, with $\mathbf{1}$ as the corresponding eigenvector.

8. **[Random Walk] How are the eigenvalues of the transition matrix $P = D^{-1} A$ related to the Laplacian eigenvalues?**
   ??? success "Solution"
       $P = I - D^{-1} L$. If $G$ is regular, $P = I - \frac{1}{k} L$. Their eigenvalues $\mu_i$ and $\lambda_i$ satisfy $\mu_i = 1 - \lambda_i / k$.

9. **[Bipartite] Prove: A graph is bipartite if and only if its adjacency spectrum is symmetric about the origin.**
   ??? success "Solution"
       If bipartite, $A = \begin{pmatrix} 0 & B \ B^T & 0 \end{pmatrix}$. If $(\lambda, [u; v])$ is an eigenpair, then $(-\lambda, [u; -v])$ is also one. Hence the spectrum is symmetric.

10. **[Application] In PageRank, why is the transpose of the adjacency matrix more important than the matrix itself?**
    ??? success "Solution"
        Page importance is determined by in-links. Adjacency columns correspond to out-links, and rows to in-links. To sum the weights pointing to each page, we need row sums of the transpose (or columns of the original).

## Chapter Summary

The intersection of graph theory and linear algebra defines modern network science:

1. **Structural Encoding**: Matrices transform non-numerical connectivity into computable linear structures.
2. **Spectral Fingerprints**: Eigenvalue distributions reveal robustness, expansion, and hierarchy in networks.
3. **Algebra of Flows**: Incidence matrices provide conservation laws for studying electricity, water, and info flows in social media.
