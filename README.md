# BIGDAWGv2

## Install

### 1.) Install BIGDAWGv2 from GitHub using devtools

```{r install}

library(devtools)
if(!system.file(package='BIGDAWGv2')){
  devtools::install_github('mmariani123/BIDAWGv2')
}

```
## Run

### 2.) Run the Labrador (Lab) test data making sure to select
### 'dla' as the species parameter

```{r run}

library(BIGDAWGv2)

#Select your species:
Species <- 'dla'

#Select your input data:
Data <- BIGDAWGv2::labData
#Data <- BIGDAWGv2::akitaData  

#Choose your loci of interest:
LociSet <- list('DRB1','DQA1')

#Set your output directory:
resultsDir <- '<path/to/output/dir>'

#Set your other parameters and run BIGDAWGv2:
bigdawg.results <- BIGDAWGv2(
  Data         = Data,
  Species      = Species,
  HLA          = FALSE,
  Run.Tests    = c('HWE', 'A'),
  Loci.Set     = LociSet, #Must be a list or will throw an error
  Exon         = c(1),
  All.Pairwise = FALSE,
  Trim         = FALSE,
  Res          = 2,
  EVS.rm       = FALSE,
  Missing      = 2,
  Strict.Bin   = FALSE,
  Cores.Lim    = 1L,
  Results.Dir  = ResultsDir,
  Return       = TRUE,
  Output       = TRUE,
  Merge.Output = FALSE,
  Verbose      = TRUE)

```

## Notes

### 3.) Ongoing notes

############################## NOTES ###################################
########################################################################
########################################################################
########################################################################
########################################################################

BIGDAWGv2 Software in progress with Professor Richard Single PhD of the 
University of Vermont and Professor Steven Mack of the University of 
California San Francisco. R package for statistical analyses of highly 
polymorphic genetic systems (e.g. HLA genes), now incorporating various 
non-human species.  

Original BIGDAWG software:

Data sets and functions for chi-squared Hardy-Weinberg and case-control
association tests of highly polymorphic genetic data [e.g., human 
leukocyte antigen (HLA) data]. Performs association tests at multiple 
levels of polymorphism (haplotype, locus and HLA amino-acids) as 
described in Pappas DJ, Marin W, Hollenbach JA, Mack SJ (2016) 
<doi:10.1016/j.humimm.2015.12.006>. Combines rare variants to a 
common class to account for sparse cells in tables as described by 
Hollenbach JA, Mack SJ, Thomson G, Gourraud PA (2012) 
<doi:10.1007/978-1-61779-842-9_14>.
    
For more details, please visit http://tools.immunogenomics.org.

BIGDAWG version2 in progress with additions by Michael P. Mariani PhD
and Richard Single PhD.

Updates for BIGDAWGv2:

Code in the main BIGDAWG(v2) function has been refactored somewhat to 
allow for more species. New, smaller, sub, functions have been created
from the code in this file. The main new paramtere is "Species" which
can be set to "hla" - human, 'dla' - dog, 'bla' - cow, and 'gla' - 
chicken.

The Hardy Weinburg test will be updated to run without specifying 
controls in the input data set as this should run regardless.

The amino anid analaysis functionality needed to be updated to allow for
species selection, and specifically, the uption to create a new 
ExonProtAlign object for a new species has to be incorporated and 
updateable Just like the 'hla' Exon object (.rda) is updateable

.rda objects will be stored in the top-level 'data' folder. Files used
internally by the program are stored in the 'extdata' folder. 

Thus the new code will generate an 'UpdatePtnAlign.RData' object and 
store it for use when 'dla' (or another species) is selected for 
anlaysis.

Additionally functionality was created in the create_files_functions.R 
file in order to create prot (_prot.txt) and pGroup
(_nom_p.txt) files for the dla and genes and dla overall respectively. 
This functionality is not included as part of BIGDAWGv2 but had to be 
created for both testing and to create the Exon object that was needed, 
e.g. for use with the AA test. 

The general AA test update workflow stacked as follows:
A-Exon_Ptn_Align_functions_dog.R --> A_wrapper.R --> A.R
These scripts call on the updated A_support functions file. 

To create/update the new dog dla exon object the following files had to 
be modified both to allow for the creation of a new .rda object (not 
just updating) and to allow for updating  of said object for new 
species. This new/updated Exon .rda object is stored in the top-level 
data folder as mentioned above

ExonPtnAlign.Create is found in A_ExonPtnAlign_functions_dog.R
just as it was found in the original A_ExonPthAlign_functions.R
file. The new AlignObj.Create and former AlignObj.Update functions are 
also found in this file.

Labrador and Akita (dog/dla) test data sets will be included in 
BIGDAWGv2 internally. 

Specifically, when going through the AminoAcid Analysis component for 
dog,

debug(A.wrapper)

the Check.Params function was modified to allow for multiple
species conditions

the run hla_checks was reorganized given the need for multiple
species now and a new run_dla_checks function was created

The more generic check_functions.R file was modified in some 
cases again to add additional cases for additional species, e.g. dog.

mult_set_dup_check function was put in its own file
and this function was commented out when dealing with 
dla data in the main BIGDAWGv2.R file/function. 

CheckAlleles commnented out in run_dla_checks.R
can consider something similar for dog

The run_amino_acid_analysis function has been adapated for
multiple species now

A_wrapper.R runs a variety of functions, e.g., ExonFilter(),
that are found in A_support_functions.R

Exon.Filter() is found in A_support_functions.R,
AlignmentFilter, Condense.EPL is found here too.

The other functions in A_support_functions.R are run 
in the A.R file. 

Thus run_amino_acid_analysis() will run A_wrapper(), 
which will run ExonFilter(), Condense.EPL(), and A().

#####################################################

end_analysis() function in BIGDAWGv2 main:

Do we want to set a specific path here?

For the ExonPtnAlign_functions_dog.R file, added the Species
parameter and output path to extdata + species for the 
functions contained therein. 

Note system.file vs. package.file

Probably could clean up ExonPtnAlign.Create() function more
after adjustments

Remember had to use truncated allele for dog right now, see 
the ALignMatrix in A_ExonPtnAlign_functions_dog.R

Is the release.txt file being stored in the right location?

A_ExonPtnAlign_functions_dog.R could probably just be kept as
A_ExonPtnAlign_functions.R with the additional species 
switching functionality included

Consider Bengstrom's future batch tools as well. 

Still need to remove 'browser' commands where needed

Need to finish dog unit test and vignette
