# postMUT

postMUT (posterior probability of a given mutation being deleterious) 
  is a tool written in Perl and R which combines the functional predictions from the *in silico* 
  algorithms SIFT, PolyPhen-2 and Xvar to provide a POSTerior probability of a given 
  MUTation being deleterious. This tool requires the installation of the pre-computed
  predictions from the *in silico* algorithms.  This tool will extract all missense 
  mutations from a given VCF file, extract the functional predictions from SIFT, 
  PolyPhen-2 and Xvar, and compute the postMUT posterior probability for each mutation 
  that has a prediction from all three *in silico* algorithms. 

For help with postMUT, there is an EXAMPLE provided in the directory. 
  
### Installation

 See the INSTALL file for installation instructions. 
 

