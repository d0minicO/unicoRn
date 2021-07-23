# reinstall the tinytex package
#install.packages("tinytex")
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
  "{example.fasta}",
  "\\end{texshade}",
  '\\end{document}'
), "example.tex")

# test out compiling
tinytex::pdflatex("example.tex")

## if this doesn't work try this method to compile instead
## This is what unicoRn actually uses as tinytex::pdflatex did not work for me
tools::texi2pdf("example.tex", clean=TRUE)

## Alternative way to install tinytex with texshade
install_tinytex(
  force = FALSE,
  dir = "auto",
  version = "daily",
  repository = "ctan",
  extra_packages = "texshade",
  add_path = TRUE
)