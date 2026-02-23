# Chapter 70B: Matrix Applications in Molecular Biology

<div class="context-flow" markdown>

**Prerequisites**: Inner Product (Ch8) · Distance Metrics · Spectral Decomposition (Ch6) · Probability Matrices

**Chapter Outline**: Sequence Alignment Matrices (PAM/BLOSUM) → Evolutionary Distance Metrics → Jukes-Cantor Model → Kimura Model → Algebraic Reconstruction of Phylogenetic Trees → Microarray Data Analysis (PCA Application) → Contact Matrices in Protein Folding

**Extension**: Bioinformatics is a discipline built on the processing and pattern matching of large-scale sparse matrices; evolution matrices define the metric evolution within the genomic sequence space.

</div>

At the molecular scale, life is modeled as a stochastic replacement process of nucleotide or amino acid sequences in a state space. By introducing probability transition matrices and distance metrics, linear algebra provides rigorous computational means for reconstructing phylogenetic trees and identifying protein folding patterns.

---

## 70B.1 Substitution Matrices and Evolutionary Stochastic Processes

!!! definition "Definition 70B.1 (PAM Substitution Model)"
    The PAM (Point Accepted Mutation) matrix family is based on Markov chain properties. The PAM-1 matrix describes amino acid substitution probabilities over a unit of evolutionary time. Long-term evolution matrices are obtained via matrix power operations $M_k = (M_1)^k$.

!!! theorem "Theorem 70B.1 (Jukes-Cantor Distance Formula)"
    Under a four-state Markov model, the physical distance (average number of mutations per site) $d$ between two DNA sequences and the observed difference frequency $p$ satisfy a logarithmic mapping:
    $$d = -\frac{3}{4} \ln(1 - \frac{4}{3}p)$$
    This formula eliminates biases caused by "overlapping mutations" in random substitution processes via matrix logarithm operations.

---

## Exercises

1. **[Markov Powers] Explain why the PAM-250 matrix can be obtained through the 250th power operation of the PAM-1 matrix.**
   ??? success "Solution"
       This model assumes that sequence evolution is a time-homogeneous first-order Markov process. According to the Chapman-Kolmogorov equations, the transition probability matrix across multiple time steps equals the cumulative product of single-step transition matrices.

2. **[Score Matrix Symmetry] Analyze the algebraic prerequisites for amino acid substitution score matrices (e.g., BLOSUM) to typically possess symmetry.**
   ??? success "Solution"
       This is based on the Detailed Balance hypothesis, assuming the evolution process is reversible at equilibrium. The mutation flow $i 	o j$ and $j 	o i$ have equal probabilities, resulting in a symmetric log-likelihood ratio matrix.

3. **[Calculation] Given an observed difference frequency $p=0.3$ between two sequences, calculate their true evolutionary distance $d$ according to the Jukes-Cantor model.**
   ??? success "Solution"
       $d = -0.75 \ln(1 - 1.333 	imes 0.3) = -0.75 \ln(0.6) \approx 0.383$. Note that the true distance is greater than the observed difference, reflecting that multiple mutations have masked original differences.

4. **[Spectral Analysis] Analyze the spectral structure of the Jukes-Cantor transition matrix and state the physical state corresponding to the eigenvalue 1.**
   ??? success "Solution"
       The transition matrix has the form $M = (1-p)I + (p/3)(J-I)$. The eigenvalues are 1 (corresponding to the steady state $\pi = [0.25, \dots, 0.25]$) and $1-4p/3$ (triple, corresponding to the decay of differences). The eigenvalue 1 guarantees the balance of nucleotide frequencies after long-term evolution.

5. **[PCA Clustering] Describe how to utilize Principal Component Analysis (PCA) for dimensionality reduction and phylogenetic classification of large-scale genomic frequency matrices.**
   ??? success "Solution"
       Construct an $m 	imes n$ species-gene frequency matrix $X$. Perform spectral decomposition on the covariance matrix $X^T X$ to project species into a principal component subspace. Euclidean distances on this plane reflect similarities in genetic features among species.

6. **[Contact Matrix] Explain the definition of a protein Contact Matrix and its role in predicting protein tertiary structure.**
   ??? success "Solution"
       A contact matrix is a 0-1 adjacency matrix defined between residues. It captures the spatial topological constraints of a folded protein. Using matrix rank properties and Distance Geometry, one can partially reconstruct the 3D coordinates of residues from the contact matrix.

7. **[Sparse Reconstruction] Analyze how the low-rank sparse characteristics of gene expression matrices allow for the use of SVD to extract "eigen-genes".**
   ??? success "Solution"
       Gene expression is coordinated, leading to matrices with significant low-rank structures. Through truncated SVD (Ch11), the first few singular vectors define the eigen-genes that dominate changes in cell states, achieving massive information compression.

8. **[Kimura Model] Contrast the differences in transition matrix structure between the Kimura two-parameter model and the Jukes-Cantor model.**
   ??? success "Solution"
       The Kimura model distinguishes between transitions and transversions, giving the transition matrix a more complex block-cyclic structure. Its eigenvalues reflect the independent evolution rates of two different biochemical processes.

9. **[Phylogenetic Algebra] Briefly describe how to reconstruct a phylogenetic tree by minimizing the Frobenius norm error of matrices in the Neighbor-Joining (NJ) method.**
   ??? success "Solution"
       Given a distance matrix $D$ between species, find a path distance matrix $T$ induced by a tree topology such that $\|D - T\|_F$ is minimized. This is equivalent to finding the optimal approximation within a subspace satisfying Ultrametric constraints.

10. **[Information Entropy] Explain how the spectral entropy of an evolution matrix quantifies the randomization of genomic information sequences over time.**
    ??? success "Solution"
        As the matrix power increases, the transition matrix tends towards a rank-1 stationary matrix. The exponential decay of eigenvalues quantifies the growth of sequence entropy, reflecting the irreversible loss of original genetic information under mutation noise interference.

## Chapter Summary

This chapter discusses the quantitative role of linear algebra in parsing molecular information of life:

1. **Probability Transition Framework**: Established Markov transition matrices as the universal mathematical language for molecular evolution modeling.
2. **Distance Metrics**: Derived evolutionary distance correction models based on matrix logarithmic mapping, restoring true mutational history.
3. **Structural Parsing**: Established an algebraic bridge from topological connectivity to physical conformation via contact matrices and distance geometry.
4. **Population Inference**: Demonstrated feature clustering analysis of high-dimensional biological data on low-dimensional manifolds using PCA and SVD.
