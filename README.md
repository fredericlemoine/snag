# SNAG (SequeNce simulAtion in Golang) simulates sequences given a set of input trees.

## Models
Evolutionary model can be chosen with -model <name> option.

### Available nucleotide models are:
- jc
- k2p : parameter kappa can be set with -parameters <kappa>
- f81 : parameters can be set with -parameters <piA,piC,piG,piT>
- gtr : parameters can be set with -parameters <d,f,b,e,a,c,piA,piC,piG,piT>, with gtr matrix being:
       ⌈ *  d  f  b ⌉
       | d  *  e  a |
       | f  e  *  c |
       ⌊ b  a  c  * ⌋

### Available protein models are:
- jtt (https://doi.org/10.1093/bioinformatics/8.3.275)
- wag (https://doi.org/10.1093/oxfordjournals.molbev.a003851)
- lg (https://doi.org/10.1093/molbev/msn067)
- hivb (https://doi.org/10.1371/journal.pone.0000503)
- ab (https://doi.org/10.1093/molbev/msu340)
- mtrev (https://doi.org/10.1007/bf02498640)

## Site evolutionary rates
By default, site rates follow a discrete gamma distribution with a shape parameter (alpha) of 1.0 and 4 categories, but it is possible to :

- disable the gamma distribution with: -gamma=false
- use a continuous gamma distribution with: -discrete=false
- change the number of categories with: -gamma-cat=6
- change the alpha parameter with: -alpha=0.8
- provide an input file with site rates, with `-rates <file>`. In this case, other rate related options are ignored (`-gamma`, `-gamma-cat`, `-alpha`, `-discrete`), and site rates are taken from the file. There must be exactly one rate per site. Otherwise, the other options are taken into account.

## SNAG Usage:
```
Usage of ./snag:
  -alpha float
    	gamma shape parameter (default 1)
  -ancestral
    	If true, then write ancestral sequences as internal nodes comments in the output tree file
  -discrete
    	discrete gamma distribution (default true)
  -gamma
    	enable gamma distribution of site rates (default true)
  -gamma-cat int
    	number of gamma categories (default 4)
  -help
    	help
  -intree string
    	Input tree (default "stdin")
  -length int
    	Simulated alignment length (default 100)
  -model string
    	Evolutionary model (for dna: jc, k2p, f81, gtr; for aa: jtt, wag, lg, hivb, mtrev, ab) (default "k2p")
  -num-aligns int
    	number of alignments to simulate per input tree (default 1)
  -out-align string
    	Output alignment file (default "stdout")
  -out-rates string
    	Output site rates file (default "stdout")
  -out-trees string
    	Output tree file with real nb mutations as branch lengths and potentially ancestral sequences at internal nodes (default "stdout")
  -parameters string
    	Model parameters: k2p: 'kappa'; f81: 'piA,piC,piG,piT'; gtr: 'd,f,b,e,a,c,piA,piC,piG,piT, amino-acid models: piA, piR, piN, piD, piC, piQ, piE, piG, piH, piI, piL, piK, piM, piF, piP, piS, piT, piW, piY, piV (otherwise: take model frequencies)'
  -rates string
    	Input site rate file (one rate per line), if given, -gamma, -alpha, -discrete and -gamma-cat are ignored (default "none")
  -root-seq string
    	Fasta file with sequence to take as root (invalidates -length)
  -seed int
    	Random Seed parameter (default 1699004613724063000)
  -version
    	Display version of snag
```
