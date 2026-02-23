# Chapter 27: Linear Algebra in Graph Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Eigenvalues (Ch6) · Combinatorial Structures (Ch65B)

**Chapter Outline**: Adjacency Matrix $A$ → Paths and Powers of $A$ → Degree Matrix $D$ → Laplacian Matrix $L = D - A$ → Matrix-Tree Theorem → Spectral Clustering → Graph Partitioning → Fiedler Value and Fiedler Vector → Incidence Matrix → Relation to Random Walks

**Extension**: Spectral graph theory uses the eigenvalues of the Laplacian to reveal the connectivity and "bottleneck" structure of a network.

</div>

Graph theory and linear algebra are deeply intertwined. A graph can be perfectly represented by its **Adjacency Matrix**, where matrix powers count the number of paths between nodes. Even more profound is the **Laplacian Matrix** $L = D - A$, whose eigenvalues provide information about the graph's connectivity, diameter, and partitioning. The second smallest eigenvalue, the **Fiedler value**, is a definitive measure of how easily a graph can be cut into two pieces. This chapter explores the "spectral" properties of networks, bridging the gap between discrete connections and continuous algebra.

---

## 27.1 Adjacency and Laplacian Matrices

!!! definition "Definition 27.1 (Adjacency Matrix)"
    For a graph $G$ with $n$ vertices, the adjacency matrix $A$ has $a_{ij} = 1$ if there is an edge between $i$ and $j$, and 0 otherwise.

!!! definition "Definition 27.2 (Laplacian Matrix)"
    The Laplacian is $L = D - A$, where $D = \operatorname{diag}(\text{degrees})$. $L$ is always symmetric and positive semi-definite.

---

## Exercises

1. **[Fundamentals] If $A$ is the adjacency matrix, what does $(A^k)_{ij}$ represent?**
   ??? success "Solution"
       It represents the number of paths of length exactly $k$ between vertex $i$ and vertex $j$.

2. **[Laplacian] Show that $\mathbf{1} = (1, \dots, 1)^T$ is always an eigenvector of $L$ with eigenvalue 0.**
   ??? success "Solution"
       The row sums of $L$ are $\sum_j (d_{ii}\delta_{ij} - a_{ij}) = d_i - d_i = 0$. Thus $L\mathbf{1} = 0$.

3. **[Connectivity] How does the multiplicity of $\lambda=0$ in the Laplacian relate to graph components?**
   ??? success "Solution"
       The number of connected components in the graph is exactly equal to the multiplicity of the eigenvalue 0 in $L$.

4. **[Fiedler Value] What is the Fiedler value $\lambda_2$?**
   ??? success "Solution"
       It is the second smallest eigenvalue of $L$. It is strictly greater than 0 if and only if the graph is connected. It measures the "algebraic connectivity."

5. **[Matrix-Tree] State the Matrix-Tree Theorem.**
   ??? success "Solution"
       The number of spanning trees in a graph is equal to any cofactor of the Laplacian matrix $L$.

6. **[Clustering] Describe Spectral Clustering using the Fiedler vector.**
   ??? success "Solution"
       To partition a graph, calculate the eigenvector $v_2$ corresponding to $\lambda_2$. Group vertices into two sets based on the sign of the entries in $v_2$. This minimizes the "cut" while keeping the groups balanced.

7. **[Regular Graphs] If $G$ is a $d$-regular graph, relate the spectra of $A$ and $L$.**
   ??? success "Solution"
       $L = dI - A$. Thus $\lambda_i(L) = d - \lambda_{n-i+1}(A)$. The Laplacian spectrum is a shifted and flipped version of the adjacency spectrum.

8. **[Incidence Matrix] Define the incidence matrix $B$ and show $L = B B^T$.**
   ??? success "Solution"
       For a directed graph, $B_{ve} = 1$ if vertex $v$ is the head of edge $e$, $-1$ if it is the tail, and 0 otherwise. $B B^T$ results in $D - A$, the Laplacian of the underlying undirected graph.

9. **[Quadratic Form] Express $x^T L x$ in terms of edges.**
   ??? success "Solution"
       $x^T L x = \sum_{\{i,j\} \in E} (x_i - x_j)^2$. This confirms that $L \succeq 0$ and shows that the Laplacian measures the "smoothness" of a signal $x$ on the graph.

10. **[Gershgorin] Use Gershgorin's theorem to bound the eigenvalues of $L$.**
    ??? success "Solution"
        Each disc is centered at $d_i$ with radius $d_i$. Thus all eigenvalues lie in $[0, 2 \max d_i]$.

## Chapter Summary

This chapter explores the spectral signature of discrete networks:

1. **Algebraic Encoding**: Represented graph topology through adjacency and degree matrices, linking path-counting to matrix power.
2. **Connectivity Spectrum**: Established the Laplacian as the definitive operator for analyzing graph components and algebraic connectivity.
3. **Combinatorial Laws**: Leveraged the Matrix-Tree theorem to link determinants to the number of spanning structures.
4. **Spectral Partitioning**: Positioned the Fiedler vector as the optimal tool for balanced graph cuts and data clustering.
