# Technical variability in Spatial Transcriptomics technologies with sub-cellular resolution

<!-- badges: start -->

<!-- badges: end -->

This project aims to assess the technical variability, or conversely the replicability, of the feature counts measured in CosMx, MERSCOPE, and Xenium, based on the publicly available datasets initially released.

Adjacent slices were lined up using the [shinySTRegister package](https://github.com/JonSulc/shinySTRegister). Transcript counts were then binned to create pseudo-spots which could be compared across slices.

Correlation was assessed in two ways:

-   Spot-wise correlation is the correlation of feature counts in a spot in one slice with those in the same spot on the adjacent slice,

```math
cor_{spot}(S_{i,j}) = cor(\begin{bmatrix}
f_{1,i,j} \\
f_{2,i,j} \\
... \\
f_{n,i,j}
\end{bmatrix},
\begin{bmatrix}
f'_{1,i,j} \\
f'_{2,i,j} \\
... \\
f'_{n,i,j}
\end{bmatrix}
)
```
  where $f_{1,i,j}$ is the number of counts for the first feature in spot $i,j$ in the first slice and $f'_{1,i,j}$ is the corresponding measure in the adjacent slice.

-   Feature-wise correlation is the correlation of a given feature across all spots with the corresponding counts in the adjacent slice,

```math
cor_{feature}(f_i) = cor(\begin{bmatrix}
f_{i,1,1} \\
f_{i,1,2} \\
... \\
f_{i,2,1} \\
... \\
f_{i,l,m}
\end{bmatrix},
\begin{bmatrix}
f'_{i,1,1} \\
f'_{i,1,2} \\
... \\
f'_{i,2,1} \\
... \\
f'_{i,l,m}
\end{bmatrix}
)
```
## Setup

1.  Clone the code from github:

    `git clone git@github.com:JonSulc/subcell-correlation.git`

2.  Change directory to the project root:

    `cd subcell-correlation`

3.  Download the appropriate data folder from the cluster (using your username):

    `rsync -azP <user>@curnagl.dcsr.unil.ch:/work/PRTNR/CHUV/DIR/rgottar1/spatial/subcell-correlation/data .`

## Execution

The interface relies on a mix of Quarto and shiny. To run it, you can use the *Run document* button from Rstudio or run `quarto render && quarto preview` from within the project root directory.
