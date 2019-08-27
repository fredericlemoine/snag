# SNAG
SNAG (SequeNce simulAtion in Golang) simulates sequences given a set of input trees.

So far only LG (with `-aa`) and K2P (default) models are available on the command line.

By default, site rates follows a discrete gamma distribution with a shape parameter (alpha) of 1.0 and 4 categories, but it is possible to :

- disable the gamma distribution with: `-gamma=false`
- use a continuous gamma distribution with: `-discrete=false`
- change the number of categories with: `-gamma-cat=6`
- change the alpha parameter with: `-alpha=0.8`

## Usage

```
Usage of snag:
  -aa
    	Simulate protein sequence
  -alpha float
    	gamma shape parameter (default 1)
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
  -kappa float
    	Kappa parameter of K2P model (default 4)
  -length int
    	Simulated alignment length (default 100)
  -num-aligns int
    	number of alignments to simulate per input tree (default 1)
  -out-align string
    	Output alignment file (default "stdout")
  -out-rates string
    	Output site rates file (default "stdout")
  -out-trees string
    	Output tree file (default "stdout")
  -seed int
    	Random Seed parameter (default )
```
