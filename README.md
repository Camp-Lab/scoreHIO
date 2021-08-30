# scoreHIO

scoreHIO leverages a published human developing multi-endodermal-organ cell atlas ([Yu et al, Cell, 2021](https://doi.org/10.1016/j.cell.2021.04.028)) to quantify the fidelity and epithelial stem cell maturity of human pluripotent stem cell derived intestine organoids

## Introduction
Due to the pluripotency of embryonic stem cells (ESCs) and induced pluripotent stem cells (iPSCs), even being provided with signaling molecules to steer intestinal cell fate, off-target cell types could be generated in the human pluripotent stem cell derived intestine organoids (HIOs). By comparing the organoid cells and the developing human multi-endodermal-organ cells, this package provides quantitative estimates about intestinal cell fate specification of the HIOs. For visualization, it provides the option to project the organoid cells to the developing atlas. In addition to assessment of fidelity, this package also supports intestinal or uncommitted epithelial stem cell maturity estimation by comparing to the develping or adult duodenum intestinal stem cells.

<img src="man/figures/art.jpeg" align="center" />

## Installation
scoreHIO reads Seurat object as input. The function of fidelity quantification by scoreHIO relies on R packages uwot and RANN. To write results to query Seurat object, Seurat package (>=3.0) is required. scoreHIO supports multiple ways to quantify stem cell maturity. If you would like to try the quadratic-programming-based method, you would need to install package quadprog.    

'uwot', 'RANN', 'Seurat', 'quadprog' could be installed from CRAN
```r
install.packages('Seurat')
install.packages('uwot')
install.packages('RANN')
install.packages('quadprog')
```

After installation of dependencies, you could install scoreHIO from github
```r
devtools::install_github('Camp-Lab/scoreHIO')
```

In addition to the R package, you'll also need to download the developing human multi-endodermal organ reference from [here](https://doi.org/10.17632/x53tts3zfr), and provide the path to the folder with these downloaded data to the parameter `organ_ref_dir` of the function `score_fidelity`. By default, `score_fidelity` assumes this folder is named as `Ref_data_for_projection_to_fetal_atlas` which is in parallel to the working folder.

# Quick start
```r
# Load Packages
library(scoreHIO)

# A `test_seurat_object` consisting of 1000 intestinal or non-fully specified stem cells of HIOs of two groups is automatically loaded for your test. It contains normalized expression data.

# You could estimate the organ fidelity of these query HIOs by running this:
seurat_object <- score_fidelity(
	que_obj = test_seurat_object,
	organ_ref_dir = "../Ref_data_for_projection_to_fetal_atlas/",
	group_by = "orig.ident"
	)

# You could estimate the maturity of the seleted cells with a simialrity based method by running this:
maturity_res <- score_maturity(
	que_obj = test_seurat_object,
	group_by = "orig.ident",
	score_method = "similarity_based"
	)

# or with a deconvolution-based method by running this:
maturity_res <- score_maturity(
	que_obj = test_seurat_object,
	group_by = "orig.ident",
	score_method = "deconvolution_based"
	)
```

## Citation
If you find this tool helpful, please consider citing our original paper [Yu et al, Cell, 2021](https://doi.org/10.1016/j.cell.2021.04.028)
