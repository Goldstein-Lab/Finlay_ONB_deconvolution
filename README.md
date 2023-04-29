# Finlay_ONB_deconvolution


# Deconstructing olfactory epithelium developmental pathways in olfactory neuroblastoma

## Authors
John B. Finlay<sup>1,2,3</sup>, Ralph Abi Hachem<sup>4</sup>, David W. Jang<sup>4</sup>, Nosayaba Osazuwa-Peters<sup>4</sup>, Bradley J. Goldstein<sup>2,3,4,*</sup>

<sup>1</sup>Medical Scientist Training Program, Duke University School of Medicine, Durham, NC
27710\
<sup>2</sup>Department of Head and Neck Surgery & Communication Sciences, Duke University
School of Medicine, Durham, NC 27710\
<sup>3</sup>Department of Cell and Molecular Biology, Duke University School of Medicine, Durham, NC
27710\
<sup>4</sup>Department of Neurobiology, Duke University School of Medicine, Durham, NC 27710\
<sup>*</sup>Corresponding author. Address: 201 Research Dr, Room 485, Durham, NC 27710
Email: bradley.goldstein@duke.edu

## Abstract
Olfactory neuroblastoma is a rare tumor arising from the olfactory cleft region of the nasal cavity. Due to the low incidence of this tumor, as well as an absence of established cell lines and murine models, understanding the mechanisms driving olfactory neuroblastoma pathobiology has been challenging. Here, we sought to apply advances from research on the human olfactory epithelial neurogenic niche, along with new biocomputational approaches, to better understand the cellular and molecular factors in low- and high-grade olfactory neuroblastoma and how specific transcriptomic markers may predict prognosis. We analyzed a total of 19 olfactory neuroblastoma samples with available bulk RNA-Sequencing and survival data, along with 10 samples from normal olfactory epithelium. A bulk RNA-Sequencing deconvolution model identified a significant increase in globose basal cell (GBC) and CD8 T cell identities in high grade tumors (GBC from approximately 0% to 8%, CD8 T cell from 0.7% to 2.2%), and significant decreases in mature neuronal, Bowman’s gland, and olfactory ensheathing programs, in high grade tumors (mature neuronal from 3.7% to approximately 0%, Bowman’s gland from 18.6% to 10.5%, olfactory ensheathing from 3.4% to 1.1%). Trajectory analysis identified potential regulatory pathways in proliferative olfactory neuroblastoma cells, including PRC2, which was validated by immunofluorescence staining. Survival analysis guided by gene expression in bulk RNA-Sequencing data identified favorable prognostic markers such as SOX9, S100B, and PLP1 expression. Our analyses provide a basis for additional research on olfactory neuroblastoma management, as well as identification of potential new prognostic markers.


## Manuscript
For more details, please see our full manuscript at (to be posted upon final acceptance)

# Data
All datasets used in this manuscript are available on NCBI GEO at accession numbers [GSE118995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118995), [GSE80249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80249), [GSE139522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139522), and [GSE184117](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184117). 

  
# Example Code
Code to replicate analyses performed in [Finlay et al. 2023 Cancer Research Communications]

1. Download data from NCBI GEO
2. See Jupyter Notebooks above for manuscript specific dataset deconvolution, and R scripts for deSeq2 and survival analysis performed in this study.
3. Please refer to [Scanpy](https://scanpy.readthedocs.io/en/stable/), [scvi-tools](https://docs.scvi-tools.org/en/stable/tutorials/index.html), [RNA-Sieve](https://github.com/songlab-cal/rna-sieve), and [deSeq2](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) tutorials for package-specific code that was adapted for datasets in this manuscript. 

# Contact
Please consult methods described in our manuscript for more details or [contact](bradley.goldstein@duke.edu) corresponding author for specific requests.
