# Chapter 70B: Linear Algebra in Molecular Biology and Genomics

<div class="context-flow" markdown>

**Prerequisites**: SVD (Ch11) · Matrix Factorization (Ch10) · Graph Theory (Ch27) · Statistics (Ch29)

**Chapter Outline**: Sequence Alignment via Dynamic Programming (Matrix Trace) → Phylogeny and Distance Matrices → Substitution Matrices ($PAM, BLOSUM$) → Gene Expression Analysis (PCA/SVD) → Gene Regulatory Networks (Adjacency Matrices) → Protein-Protein Interaction (PPI) Networks → Structural Biology and Distance Geometry

**Extension**: SVD is the core of bio-informatics dimensionality reduction; substitution matrices derive from Markov models of amino acid evolution.

</div>

Life's microscopic code is essentially high-dimensional data. Linear algebra provides the tools to extract patterns from gene expression and reconstruct the history of evolution from sequence distances.

---

## 70B.1 Sequence Evolution and Invariant Forms

!!! definition "Definition 70B.1 (Substitution Matrix)"
    A substitution matrix $S$ (like PAM or BLOSUM) describes the probability or log-odds of one amino acid being replaced by another over evolutionary time. These matrices are derived from Markov chains on the space of 20 amino acids.

!!! theorem "Theorem 70B.3 (Phylogenetic Distance Consistency)"
    A distance matrix $D$ derived from DNA sequences represents a tree structure if and only if it satisfies the **four-point condition**: for any four taxa $i, j, k, l$, the two largest of $\{D_{ij}+D_{kl}, D_{ik}+D_{jl}, D_{il}+D_{jk}\}$ are equal.

---

## Exercises

1. **[Substitution Entropy] In a PAM-1 matrix, how does the trace $\operatorname{tr}(P)$ relate to the rate of sequence conservation?**
   ??? success "Solution"
       The trace $\operatorname{tr}(P) = \sum p_{ii}$ is the sum of the probabilities that each amino acid remains unchanged. A larger trace indicates higher sequence conservation over the given evolutionary distance (1 PAM unit).

2. **[PCA in Genomics] Explain why PCA is used to identify population stratification in genomic studies.**
   ??? success "Solution"
       Genotypes are represented as a large matrix $X$. The first few principal components often capture geographical or ethnic variations. By projecting individuals onto the PC subspace, researchers can detect clusters that correspond to ancestral lineages.

3. **[SVD Interpretation] If $A$ is a gene expression matrix (genes $\times$ samples), what is the biological meaning of the left and right singular vectors?**
   ??? success "Solution"
       The left singular vectors $u_i$ (eigen-genes) represent fundamental expression patterns. The right singular vectors $v_i$ (eigen-samples) represent the prevalence of these patterns across different biological conditions or tissues.

4. **[Calculation] Given a Jukes-Cantor distance matrix $D$, why must its entries satisfy $D_{ii} = 0$ and $D_{ij} = D_{ji}$?**
   ??? success "Solution"
       $D$ represents a metric on the sequence space. $D_{ii}=0$ reflects the identity axiom (distance to self is zero), and $D_{ij}=D_{ji}$ reflects the symmetry of the evolutionary distance between sequences $i$ and $j$.

5. **[Regulatory Networks] Represent a 3-gene feedback loop as an adjacency matrix and find its eigenvalues.**
   ??? success "Solution"
       $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$. The eigenvalues are the cube roots of unity $\{1, \omega, \omega^2\}$. The imaginary components reflect the oscillatory nature of the feedback system.

6. **[Distance Geometry] How is the SVD used to reconstruct 3D protein structures from nuclear magnetic resonance (NMR) distance constraints?**
   ??? success "Solution"
       The distance matrix $D$ is converted into a Gram matrix $G$ using $G_{ij} = \frac{1}{2}(D_{i0}^2 + D_{j0}^2 - D_{ij}^2)$. The 3D coordinates are obtained from the first 3 eigenvectors of $G$ scaled by the square root of their eigenvalues.

7. **[Network Robustness] In a Protein-Protein Interaction (PPI) network, how does the spectral gap of the Laplacian relate to network connectivity?**
   ??? success "Solution"
       The second smallest eigenvalue $\lambda_2$ (algebraic connectivity) measures how easily the network can be partitioned. A small $\lambda_2$ suggests the existence of weakly connected functional modules.

8. **[Phylogenetic Trees] Verify if the distance matrix $D = \begin{pmatrix} 0 & 3 & 7 \\ 3 & 0 & 6 \\ 7 & 6 & 0 \end{pmatrix}$ can be exactly represented by a rooted tree.**
   ??? success "Solution"
       For a tree, the distances must satisfy specific ultrametric or additive properties. Here, $D_{13} = 7$ and $D_{12}+D_{23} = 3+6=9$. Since $7 \neq 9$ and other additive checks must be performed, one evaluates the three-point or four-point conditions to determine consistency.

9. **[Mutational Equilibrium] For an amino acid substitution Markov chain with transition matrix $Q$, what does the stationary distribution $\pi$ represent?**
   ??? success "Solution"
       $\pi$ represents the equilibrium frequencies of the 20 amino acids across all proteins in the limit of infinite evolutionary time, assuming the selective pressures remain constant.

10. **[Complexity] Why is matrix factorization often preferred over direct graph algorithms for large-scale gene regulatory network analysis?**
    ??? success "Solution"
        Biological networks are often noisy and partially observed. Matrix factorization (like NMF or SVD) can extract robust latent factors and patterns while naturally handling missing data, whereas direct graph traversal is sensitive to every individual edge error.

## Chapter Summary

This chapter explores the linear algebraic framework of molecular life:

1. **Evolutionary Calculus**: Established substitution matrices as Markov transition operators for amino acid sequences.
2. **Spectral Genomics**: Utilized PCA and SVD to reduce the dimensionality of high-throughput gene expression data.
3. **Biological Networks**: Formulated the connectivity and robustness of regulatory and interaction systems using adjacency and Laplacian matrices.
4. **Distance Geometry**: Linked sequence metrics to phylogenetic trees and 3D structural reconstruction.
