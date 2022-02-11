package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/phylip"
	"github.com/evolbioinfo/goalign/models"
	"github.com/evolbioinfo/goalign/models/dna"
	"github.com/evolbioinfo/goalign/models/protein"
	"github.com/evolbioinfo/gotree/io/utils"
	"github.com/evolbioinfo/gotree/tree"
)

var Version string = "Unknown"

// Simulator Interface
type Snag interface {
	genrootseq() (root []int)
	Simulate(t *tree.Tree) <-chan SnagSimu
	SetRoot(rootseq []uint8) error
	SetAncestral(b bool)
	SetRates(r []float64)
}

// Simulator implementation
type snagImpl struct {
	pi       []float64    // Character frequencies, used for root and simulation
	ns       int          // number of character states (4 or 20)
	discrete bool         // discrete gamma or not
	alpha    float64      // alpha parameter of gamma distribution
	gamma    bool         // use a gamma distribution
	ncat     int          // number of categories of discrete gamma distribution
	naligns  int          // number of alignments
	rates    []float64    // Use given rates
	seed     int64        // random seed
	aa       bool         // simate aa or nt
	m        models.Model // simulation model
	l        int          // length of the alignment to generate
	rootseq  []int        // rootsequence to take into account (if nil, then it is randomly generated for each new simulated alignment)
	ancest   bool         // write ancestral sequences at internal node comments of output trees
}

// Result of a simulation
type SnagSimu struct {
	t *tree.Tree      // Tree with branch lengths corresponding to number of mutations
	a align.Alignment // Simulated alignment
	r []float64       // rates per sites
}

// Wether or not ancestral sequences are
// written at internal node of the output
// tree or not.
func (s *snagImpl) SetAncestral(b bool) {
	s.ancest = b
}

// SetRates will give the site rates to the simulator
// Will then desactivate gamma distribution, etc.
func (s *snagImpl) SetRates(r []float64) {
	s.rates = r
}

func (s *snagImpl) index2Seq(ind []int) (seq []uint8, err error) {
	seq = make([]uint8, len(ind))
	for i, c := range ind {
		if s.aa {
			seq[i], err = align.Index2AA(c)
		} else {
			seq[i], err = align.Index2Nt(c)
		}
	}
	return
}

func (s *snagImpl) seq2Index(seq []uint8) (ind []int, err error) {
	ind = make([]int, len(seq))
	for i, c := range seq {
		if s.aa {
			ind[i], err = align.AA2Index(c)
		} else {
			ind[i], err = align.Nt2Index(c)
		}
	}
	return
}

// Sets the root sequence to take to start simulations
// If not set, then root sequence is randomly generated
// for each new simulated alignment, by taking into account
// state frequencies s.pi.
//
// if rootseq is nil or with a size==0 or with unknown characters,
// then returns an error.
func (s *snagImpl) SetRoot(rootseq []uint8) (err error) {
	// Check that the given rootseq is initialized
	if rootseq == nil || len(rootseq) == 0 {
		err = fmt.Errorf("Nil or 0-length root sequence, cannot use it")
		return
	}
	s.rootseq = make([]int, len(rootseq))
	s.l = len(rootseq)
	// Convert it to indices
	if s.rootseq, err = s.seq2Index(rootseq); err != nil {
		err = fmt.Errorf("Error while setting root sequence, non existent character: %v", err)
	}
	return
}

// If s.rootseq is not nil, then returns s.rootseq
// and does not take into account l.
//
// Otherwise, will randomly generate the rootsequence
// given the frequencies of states given in s.pi
func (s *snagImpl) genrootseq() (root []int) {
	if s.rootseq != nil {
		root = s.rootseq
	} else {
		root = make([]int, s.l)
		for c := 0; c < s.l; c++ {
			tmp := 0.0
			rc := rand.Float64()
			nt := -1
			for rc > tmp {
				nt++
				tmp += s.pi[nt]
			}
			root[c] = nt
		}
	}
	return
}

// Checks that pi sums to 1.0.
// If not, rescales the values
func checkPI(pi []float64) {
	sum := .0
	for _, v := range pi {
		sum += v
	}
	if sum != 1.0 {
		fmt.Fprintf(os.Stderr, "[Warning] State frequencies do not add up to 1.0, rescaling\n")
		for i, v := range pi {
			pi[i] = v / sum
			fmt.Fprintf(os.Stderr, "pi[%d]=%f\n", i, pi[i])
		}
	}
}

func NewSnag(ns, l int, gamma, discrete bool, alpha float64, ncat int,
	params []float64, naligns int, seed int64, model string) (s *snagImpl, err error) {

	if ns != 4 && ns != 20 {
		err = fmt.Errorf("Number of character state can only be 4 or 20")
		return
	}

	s = &snagImpl{
		pi:       make([]float64, ns),
		ns:       ns,
		discrete: discrete,
		alpha:    alpha,
		gamma:    gamma,
		ncat:     ncat,
		naligns:  naligns,
		seed:     seed,
		aa:       false,
		m:        nil,
		l:        l,
		rootseq:  nil,
		ancest:   false,
	}

	switch model {
	case "jc":
		for i := 0; i < ns; i++ {
			s.pi[i] = 1. / float64(ns)
		}
		dm := dna.NewJCModel()
		dm.InitModel()
		s.m = dm
	case "k2p":
		if len(params) != 1 {
			err = fmt.Errorf("Wrong parameters for k2p model: %v", params)
			return
		}
		for i := 0; i < ns; i++ {
			s.pi[i] = 1. / float64(ns)
		}
		dm := dna.NewK2PModel()
		dm.InitModel(params[0])
		s.m = dm
	case "f81":
		if len(params) != 4 {
			err = fmt.Errorf("Wrong parameters for f81 model: %v", params)
			return
		}
		for i := 0; i < ns; i++ {
			s.pi[i] = params[i]
		}
		checkPI(s.pi)
		dm := dna.NewF81Model()
		dm.InitModel(s.pi[0], s.pi[1], s.pi[2], s.pi[3])
		s.m = dm
	case "gtr":
		if len(params) != 10 {
			err = fmt.Errorf("Wrong parameters for gtr model: %v", params)
			return
		}
		for i := 0; i < ns; i++ {
			s.pi[i] = params[i+6]
		}
		checkPI(s.pi)
		dm := dna.NewGTRModel()
		dm.InitModel(params[0], params[1], params[2], params[3], params[4], params[5], s.pi[0], s.pi[1], s.pi[2], s.pi[3])
		s.m = dm
	default:
		var modelint int
		var pm *protein.ProtModel
		s.aa = true

		switch model {
		case "jtt":
			modelint = protein.MODEL_JTT
		case "wag":
			modelint = protein.MODEL_WAG
		case "lg":
			modelint = protein.MODEL_LG
		case "hivb":
			modelint = protein.MODEL_HIVB
		default:
			err = fmt.Errorf("Wrong model: %s", model)
			return
		}

		if pm, err = protein.NewProtModel(modelint, gamma, alpha); err != nil {
			return
		}

		// Initialize aa frequencies as defined in the model
		for i := 0; i < ns; i++ {
			s.pi[i] = pm.Pi(i)
		}
		pm.InitModel(nil)
		s.m = pm
	}
	return
}

func (s *snagImpl) Simulate(t *tree.Tree) (simuChan <-chan SnagSimu) {
	// all simulated sequences
	var seqs [][]int
	// num nodes
	var nnodes int = 0
	// num edges
	var nedges int = 0
	// all br lengths
	var lengths []float64 = make([]float64, 0, 0)
	// alignment
	var a align.Alignment
	// proba matrix per branch and per gamma category
	var pijs [][]*models.Pij
	// root sequence
	var rootseq []int
	// gamma rates
	var rates []float64
	// gamma cats (1 -> ncat)
	var cats []int
	// temp probability
	var tmp, rc float64
	// single model pij
	var pij *models.Pij

	// Index variables
	var ia, cati, c int
	var sitecat, nt int

	// Output Channel
	var tmpChan chan SnagSimu = make(chan SnagSimu, 1)
	simuChan = tmpChan

	// Set Id of all nodes and edges
	t.PostOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) (keep bool) {
		cur.SetId(nnodes)
		if prev != nil {
			e.SetId(nedges)
			lengths = append(lengths, e.Length())
			nedges++
		}
		nnodes++
		return true
	})

	// One sequence per tree node
	seqs = make([][]int, nnodes)

	// Pij per rate category (if any) and per branch
	pijs = make([][]*models.Pij, s.ncat)
	for cati = 0; cati < s.ncat; cati++ {
		pijs[cati] = make([]*models.Pij, nnodes)
	}

	go func() {
		for ia = 0; ia < s.naligns; ia++ {
			rootseq = s.genrootseq()
			if s.rates != nil && len(s.rates) == len(rootseq) {
				rates = s.rates
				cats = make([]int, len(rootseq))
			} else {
				if s.rates != nil {
					fmt.Fprintf(os.Stderr, "Warning: A rate file is given but it does not have the same length as the root sequence, we then generate rates according to other given options")
				}
				rates, cats = models.GenerateRates(s.l, s.gamma, s.alpha, s.ncat, s.discrete)
			}
			if s.aa {
				a = align.NewAlign(align.AMINOACIDS)
			} else {
				a = align.NewAlign(align.NUCLEOTIDS)
			}
			t.PreOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) (keep bool) {
				// We keep the information about the number of
				// realized mutations in the branch length
				var nbmuts = 0
				if prev == nil {
					seqs[cur.Id()] = rootseq
					// We write the sequence at the root
					if s.ancest {
						var aaseq []uint8
						aaseq, _ = s.index2Seq(rootseq)
						cur.AddComment(string(aaseq))
					}
				} else {
					prevseq := seqs[prev.Id()]
					curseq := make([]int, s.l)
					seqs[cur.Id()] = curseq
					for c = 0; c < s.l; c++ {
						sitecat = 0
						if s.rates == nil && s.gamma && s.discrete {
							sitecat = cats[c]
						}
						// We initialize the pij for that branch and that category
						if pijs[sitecat][cur.Id()] == nil {
							pijs[sitecat][cur.Id()], _ = models.NewPij(s.m, lengths[e.Id()]*rates[c])
						}
						pij = pijs[sitecat][cur.Id()]
						// If not discrete, we must recompute the pijs
						// with given most probably unique (rate,brlen) pair
						if s.rates != nil || (s.gamma && !s.discrete) {
							pij.SetLength(lengths[e.Id()] * rates[c])
						}
						// we simulate the next character state along the branch
						tmp = 0.0
						rc = rand.Float64()
						nt = -1
						for rc > tmp {
							nt++
							tmp += pij.Pij(prevseq[c], nt)
						}
						curseq[c] = nt
						if prevseq[c] != curseq[c] {
							nbmuts++
						}
					}
					e.SetLength(float64(nbmuts))
					if cur.Tip() || s.ancest {
						var aaseq []uint8
						aaseq, _ = s.index2Seq(curseq)
						if cur.Tip() {
							a.AddSequenceChar(cur.Name(), aaseq, "")
						}
						if s.ancest {
							cur.AddComment(string(aaseq))
						}
					}
				}
				return true
			})
			tmpChan <- SnagSimu{t.Clone(), a, rates}
		}
		close(tmpChan)
	}()
	return
}

func snagMain() (exitcode int) {

	var treeChan <-chan tree.Trees
	var rootsequence []uint8
	var treeReader *bufio.Reader
	var treeFile io.Closer
	var err error
	var paramslice []float64
	var rates []float64
	var tmpf float64
	var outalignfile, outtreefile, outratefile *os.File
	var snag Snag

	exitcode = 0

	// Arguments

	helpmessage := `SNAG (SequeNce simulAtion in Golang) simulates sequences given a set of input trees.

# Models
Evolutionary model can be chosen with -model <name> option. 

## Available nucleotide models are:
- jc
- k2p : parameter kappa can be set with -parameters <kappa>
- f81 : parameters can be set with -parameters <piA,piC,piG,piT>
- gtr : parameters can be set with -parameters <d,f,b,e,a,c,piA,piC,piG,piT>, with gtr matrix being:
       ⌈ *  d  f  b ⌉
       | d  *  e  a |
       | f  e  *  c |
       ⌊ b  a  c  * ⌋

## Available protein models are:
- jtt
- wag
- lg
- hivb

# Site evolutionary rates
By default, site rates follow a discrete gamma distribution with a shape parameter (alpha) of 1.0 and 4 categories, but it is possible to :

- disable the gamma distribution with: -gamma=false
- use a continuous gamma distribution with: -discrete=false
- change the number of categories with: -gamma-cat=6
- change the alpha parameter with: -alpha=0.8
- use site rates given by the user -rates <file> (-gamma and -discrete are useless in this case)

`
	ns := 4
	discrete := flag.Bool("discrete", true, "discrete gamma distribution")
	alpha := flag.Float64("alpha", 1.0, "gamma shape parameter")
	gamma := flag.Bool("gamma", true, "enable gamma distribution of site rates")
	ncat := flag.Int("gamma-cat", 4, "number of gamma categories")
	intree := flag.String("intree", "stdin", "Input tree")
	inrates := flag.String("rates", "none", "Input site rate file (one rate per line), if given, -gamma, -alpha, -discrete and -gamma-cat are ignored")
	naligns := flag.Int("num-aligns", 1, "number of alignments to simulate per input tree")
	seed := flag.Int64("seed", time.Now().UTC().UnixNano(), "Random Seed parameter")
	l := flag.Int("length", 100, "Simulated alignment length")
	model := flag.String("model", "k2p", "Evolutionary model (for dna: jc, k2p, f81, gtr; for aa: jtt, wag, lg, hivb)")
	parameters := flag.String("parameters", "", "Model parameters: k2p: 'kappa'; f81: 'piA,piC,piG,piT'; gtr: 'd,f,b,e,a,c,piA,piC,piG,piT'")
	rootseq := flag.String("root-seq", "", "Fasta file with sequence to take as root (invalidates -length)")
	outalign := flag.String("out-align", "stdout", "Output alignment file")
	ancestral := flag.Bool("ancestral", false, "If true, then write ancestral sequences as internal nodes comments in the output tree file")
	outtrees := flag.String("out-trees", "stdout", "Output tree file with real nb mutations as branch lengths and potentially ancestral sequences at internal nodes")
	outrates := flag.String("out-rates", "stdout", "Output site rates file")
	version := flag.Bool("version", false, "Display version of snag")
	help := flag.Bool("help", false, "help")
	aa := false
	nt := true
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, helpmessage)
		fmt.Fprintf(os.Stderr, "Usage of %s:\n", os.Args[0])
		flag.PrintDefaults()
	}

	flag.Parse()

	if *help {
		flag.Usage()
		exitcode = 1
		return
	}

	if *version {
		fmt.Fprintf(os.Stderr, "%s version %s\n", os.Args[0], Version)
		exitcode = 0
		return
	}

	aa = (*model == "jtt" || *model == "wag" || *model == "lg" || *model == "hivb")
	nt = (*model == "jc" || *model == "k2p" || *model == "f81" || *model == "gtr")

	if !aa && !nt {
		fmt.Fprintf(os.Stderr, "Wrong model: %s\n", *model)
		flag.Usage()
		exitcode = 1
		return
	}

	// Create output files
	if *outalign == "stdout" {
		outalignfile = os.Stdout
	} else {
		if outalignfile, err = os.Create(*outalign); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
		defer outalignfile.Close()
	}

	if *inrates != "none" {
		var fi io.Closer
		var r *bufio.Reader
		var b float64
		rates = make([]float64, 0)

		if fi, r, err = utils.GetReader(*inrates); err != nil {
			fmt.Fprintf(os.Stderr, "Error while parsing site rate file: %s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
		scanner := bufio.NewScanner(r)
		// optionally, resize scanner's capacity for lines over 64K, see next example
		for scanner.Scan() {
			if b, err = strconv.ParseFloat(scanner.Text(), 64); err != nil {
				fmt.Fprintf(os.Stderr, "Error while parsing site rate file: %s\n", err)
				flag.Usage()
				exitcode = 1
				return
			}
			rates = append(rates, b)
		}
		fi.Close()
	}

	if *outtrees == "stdout" {
		outtreefile = os.Stdout
	} else {
		if outtreefile, err = os.Create(*outtrees); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
		defer outtreefile.Close()
	}
	if *outrates == "stdout" {
		outratefile = os.Stdout
	} else {
		if outratefile, err = os.Create(*outrates); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
		defer outratefile.Close()
	}

	if aa {
		ns = 20
	}

	rand.Seed(*seed)

	if !*gamma || !*discrete {
		*ncat = 1
	}

	// Parse parameters
	paramslice = make([]float64, 0, 0)
	if *parameters != "" {
		tmpslice := strings.Split(*parameters, ",")
		for _, v := range tmpslice {
			if tmpf, err = strconv.ParseFloat(v, 64); err != nil {
				fmt.Fprintf(os.Stderr, "Wrong parameters arguments: %s\n", *parameters)
				flag.Usage()
				exitcode = 1
				return
			}
			paramslice = append(paramslice, tmpf)
		}
	}

	// Parse root sequence if one is given
	if *rootseq != "" {
		var fi io.Closer
		var r *bufio.Reader

		var tmpali align.Alignment

		if fi, r, err = utils.GetReader(*rootseq); err != nil {
			fmt.Fprintf(os.Stderr, "Error while parsing root sequence: %s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
		if tmpali, err = fasta.NewParser(r).Parse(); err != nil {
			fmt.Fprintf(os.Stderr, "Error while parsing root sequence: %s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
		fi.Close()
		if tmpali.NbSequences() != 1 {
			fmt.Fprintf(os.Stderr, "Root Sequence file must contain exactly one sequence\n")
			flag.Usage()
			exitcode = 1
			return
		}

		if (aa && tmpali.Alphabet() != align.AMINOACIDS) ||
			(!aa && tmpali.Alphabet() != align.NUCLEOTIDS) {
			fmt.Fprintf(os.Stderr, "Root Sequence alphabet is not compatible with simulation model\n")
			flag.Usage()
			exitcode = 1
			return
		}
		rootsequence, _ = tmpali.GetSequenceCharById(0)
	}

	// Parse Trees
	if treeFile, treeReader, err = utils.GetReader(*intree); err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err)
		flag.Usage()
		exitcode = 1
		return
	}
	defer treeFile.Close()
	treeChan = utils.ReadMultiTrees(treeReader, utils.FORMAT_NEWICK)

	// Simulate alignments
	if snag, err = NewSnag(ns, *l, *gamma, *discrete, *alpha, *ncat, paramslice, *naligns, *seed, *model); err != nil {
		fmt.Fprintf(os.Stderr, "%s\n", err)
		flag.Usage()
		exitcode = 1
		return
	}
	snag.SetAncestral(*ancestral)
	snag.SetRates(rates)

	// If a root sequence is given we give it to the simulator
	if rootsequence != nil {
		if err = snag.SetRoot(rootsequence); err != nil {
			fmt.Fprintf(os.Stderr, "%s\n", err)
			flag.Usage()
			exitcode = 1
			return
		}
	}
	for trees := range treeChan {
		for si := range snag.Simulate(trees.Tree) {
			fmt.Fprintln(outratefile, si.r)
			fmt.Fprintln(outtreefile, si.t.Newick())
			fmt.Fprintln(outalignfile, phylip.WriteAlignment(si.a, false, false, false))
		}
	}

	return
}

func main() {
	os.Exit(snagMain())
}
