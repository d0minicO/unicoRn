# UnicoRn
## _R multiple sequence alignment tool for protein conservation analysis_

UnicoRn is a custom R function that takes gene names and outputs a pdf for each gene showing a pretty alignment of the amino acid sequence for that gene across multiple predefined species ([![example](https://ibb.co/RYNdk5r)](https://github.com/d0minicO/unicoRn/blob/main/example_output_files/example_output.PNG))

---
# Install
### Part 1 -- install package dependencies if you have not already got them
```sh
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("msa")
BiocManager::install("Biostrings")
install.packages("seqinr")
install.packages("tidyverse")
install.packages("data.table")
install.packages("tools")
install.packages("devtools")
install.packages("tinytex")
```
### Part 2 -- source unicoRn from github
```sh
devtools::source_url("https://github.com/d0minicO/unicoRn/blob/main/unicoRn.R?raw=TRUE")
```
### Part 3 (optional) -- download the deleted uniprot IDs database
- download this Rds file (https://1drv.ms/u/s!Ah6q8jTg5ESfgaM4owHIdjhCGqcTwg?e=MTklQ4) of deleted uniprot IDs locally to your computer somewhere accessible to your Rstudio session
-- you will pass the location of this downloaded file as an argument to unicoRn (see below)
---


# Arguments
- base = route / working directory _(eg. "C:/user/baloons/")_
- subs name = name of this list of genes to save in their own folder _(eg. "UnicoRn_analysis")_
- genes = character string or vector eg "Gene", or c("Gene1", "Gene2")
- len = the length of the first how many AAs to plot (or if you just want the whole length make len="whole" or any characters will do)
- speciesToUse = character string of the species to use
- del_data = location where the databse of IDs deleted from uniprot are _(eg. "C:/user/baloons/del_data.Rds")_
 -- these must be formatted in the way that ensembl recognises them and separated with a boolean **or** _(eg. "hsapiens|mmusculus|clfamiliaris")_ will get sequences for human, mouse, and dog
- check_delID = TRUE or FALSE. Some Uniprot IDs get deleted by them but are still found in the biomart database. So you might need to filter these out before trying to get sequences from uniprot otherwise you will receive an error and no output. This takes a long time (several Gb of RAM to load the database...) so default should be to try running without checking for deleted IDs _(eg FALSE)_

---

# Output
[![Example image](https://ibb.co/RYNdk5r)](https://github.com/d0minicO/unicoRn/blob/main/example_output_files/example_output.PNG)

See https://github.com/d0minicO/unicoRn/blob/main/example_output_files/ for examples of the .fasta alignment, .tex file, and alignment.pdf

---

# How does it work?

UnicoRn uses the powerful biomaRt package (https://bioconductor.org/packages/release/bioc/html/biomaRt.html) to locate gene names across species and link them to a uniprot ID. Then, the curated uniprot database (https://www.uniprot.org/) is queried to extract the amino acid sequences for those proteins. The package msa (https://bioconductor.org/packages/release/bioc/html/msa.html) is used to generate a multiple sequence alignment. The alignment is aved (.fasta format) and a custom .tex file is created utilising the texshade package (https://www.ctan.org/pkg/texshade) inspired by this example (https://www.overleaf.com/latex/templates/standalone-msa-figure/rbgrxrmctccc). Finally, the tinytex implimentation of latex (https://www.rdocumentation.org/packages/tinytex/versions/0.32) is used to compile and save a pdf of the alignment.

---

# Troubleshooting

## Tinytex
Latex must be installed and communicating properly with R in order to compile and save the pdf. Tinytex (https://www.rdocumentation.org/packages/tinytex/versions/0.32) is a fairly straightforward way to get Rstudio and latex to communicate. If the pdf output step fails then it is likely that Rstudio is not able to communicate properly with tinytex/latex. You can try testing out the following code to see where the problem lies

```sh
## See if tinytex can compile this basic test pdf
writeLines(c(
   '\\documentclass{article}',
   '\\begin{document}',
   "Help, get it to work!",
   '\\end{document}'
 ), "test.tex")
tinytex::pdflatex("test.tex")

# now test if tinytex can properly use texshade
## first make an example alignment file
writeLines(c(
  ">DDX21_Human",
  "MPGKLRSDAG--",
  ">DDX21_Mouse",
  "MPGKLRSGAK--",
  ">DDX21_Dog",
  "MPGKLLSDAG--",
  ">DDX21_Zebrafish",
  "--EKWQDSRRWT"),
  "example.fasta")

## then write a minimal .tex file
writeLines(c(
  '\\documentclass{article}',
  "\\usepackage{texshade}",
  '\\begin{document}',
  "\\begin{texshade}",
  "example.fasta}",
  "\\end{texshade}",
  '\\end{document}'
), "example.tex")

# test out compiling
tinytex::pdflatex("example.tex")

## if this doesn't work try this method to compile instead
## This is what unicoRn actually uses as tinytex::pdflatex did not work for me
tools::texi2pdf("example.tex", clean=TRUE)
```

If these didn't both work, then you will have to try make sure tinytex and texshade are properly installed

```sh
# reinstall the tinytex package
install.packages("tinytex")
library(tinytex)

# install the full version of the package
tinytex:::install_prebuilt('TinyTeX')

# list installed packages available to tinytex
# check that "texshade" is listed
tl_pkgs()

# install a lot of the additional common packages
tinytex:::install_yihui_pkgs()

# force installation of texshade
tlmgr_install("texshade")

# If having problems use this to try to get more informative error messages
options(tinytex.verbose = TRUE)

## Alternative way to install tinytex with texshade
install_tinytex(
  force = FALSE,
  dir = "auto",
  version = "daily",
  repository = "ctan",
  extra_packages = "texshade",
  add_path = TRUE
)
```

For more help with tinytex/latex issues check out the developer's debugging instructions (https://yihui.org/tinytex/r/#debugging)


## Deleted Uniprot IDs

It is annoying but Uniprot is a curated database meaning that they deleted old / spurious IDs. Many of these are still associated with genes in biomaRt. If you receive an error when trying to get sequences from uniprot and you have been running the default mode of check_delID=FALSE, then it is likely that you will need to download the database (see Part 3 of installation instructions above) and run Unicorn with check_delID=TRUE, and provide the location of the downl deleted ids database Rds file to unicoRn as the del_data argument. This step takes a long time as the database is millions of entries long and takes a while to load and search.

When running with check_delID=TRUE if you receive an error like "cannot allocate vector of size...x" then it is likely that your computer ran out of RAM while loading the database of deleted IDs. Here are some things that might help:
```sh
gc() # garbage cleanup
```
-  restarting a fresh R session without any variables saved 
-  closing open programs on your computer
-  running on a machine with more RAM (16Gb Windows 10 works for me)