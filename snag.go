package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"math/rand"
	"os"
	"time"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/phylip"
	"github.com/evolbioinfo/goalign/models"
	"github.com/evolbioinfo/goalign/models/dna"
	"github.com/evolbioinfo/goalign/models/protein"
	"github.com/evolbioinfo/gotree/io/utils"
	"github.com/evolbioinfo/gotree/tree"
)

// Simulator Interface
type Snag interface {
	genrootseq(l int) (root []int)
	Simulate(t *tree.Tree, l int) <-chan SnagSimu
}

// Simulator implementation
type snagImpl struct {
	pi       []float64    // Character frequencies, used for root and simulation
	ns       int          // number of character states (4 or 20)
	discrete bool         // discrete gamma or not
	alpha    float64      // alpha parameter of gamma distribution
	gamma    bool         // use a gamma distribution
	ncat     int          // number of categories of discrete gamma distribution
	kappa    float64      // kappa parameter if k2p model is chosen
	naligns  int          // number of alignments
	seed     int64        // random seed
	aa       bool         // simate aa or nt
	m        models.Model // simulation model
}

// Result of a simulation
type SnagSimu struct {
	t *tree.Tree      // Tree with branch lengths corresponding to number of mutations
	a align.Alignment // Simulated alignment
	r []float64       // rates per sites
}

func (s *snagImpl) genrootseq(l int) (root []int) {
	root = make([]int, l)
	for c := 0; c < l; c++ {
		tmp := 0.0
		rc := rand.Float64()
		nt := -1
		for rc > tmp {
			nt++
			tmp += s.pi[nt]
		}
		root[c] = nt
	}
	return
}

func NewSnag(ns int, gamma, discrete bool, alpha float64, ncat int,
	kappa float64, naligns int, seed int64, aa bool) (s *snagImpl, err error) {

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
		kappa:    kappa,
		naligns:  naligns,
		seed:     seed,
		aa:       aa,
		m:        nil,
	}

	if !aa {
		// Initialize base frequencies at 1/4
		for i := 0; i < ns; i++ {
			s.pi[i] = 1. / float64(ns)
		}
		dm := dna.NewK2PModel()
		dm.InitModel(kappa)
		s.m = dm
	} else {
		if pm, err2 := protein.NewProtModel(protein.MODEL_LG, gamma, alpha); err2 != nil {
			return
		} else {
			// Initialize aa frequencies as defined in the model
			for i := 0; i < ns; i++ {
				s.pi[i] = pm.Pi(i)
			}
			//pm.PrintMat()
			pm.InitModel(nil)
			//pm.PrintMat()
			//pm.PrintFreqs()
			s.m = pm
		}
	}
	return
}

func (s *snagImpl) Simulate(t *tree.Tree, l int) (simuChan <-chan SnagSimu) {
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
			rootseq = s.genrootseq(l)
			rates, cats = models.GenerateRates(l, s.gamma, s.alpha, s.ncat, s.discrete)
			if s.aa {
				a = align.NewAlign(align.AMINOACIDS)
			} else {
				a = align.NewAlign(align.NUCLEOTIDS)
			}
			t.PreOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) (keep bool) {
				// We keep the information about the number of
				// realized mutations in the branch length
				nbmuts := 0
				if prev == nil {
					seqs[cur.Id()] = rootseq
				} else {
					prevseq := seqs[prev.Id()]
					curseq := make([]int, l)
					seqs[cur.Id()] = curseq
					for c = 0; c < l; c++ {
						sitecat = 0
						if s.gamma && s.discrete {
							sitecat = cats[c]
						}
						// We initialize the pij for that branch and that category
						if pijs[sitecat][cur.Id()] == nil {
							pijs[sitecat][cur.Id()], _ = models.NewPij(s.m, lengths[e.Id()]*rates[c])
						}
						pij = pijs[sitecat][cur.Id()]
						// If not discrete, we must recompute the pijs
						// with given most probably unique (rate,brlen) pair
						if s.gamma && !s.discrete {
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

					if cur.Tip() {
						aaseq := make([]rune, l)
						for c = 0; c < l; c++ {
							if s.aa {
								aaseq[c], _ = align.Index2AA(curseq[c])
							} else {
								aaseq[c], _ = align.Index2Nt(curseq[c])
							}
						}
						a.AddSequenceChar(cur.Name(), aaseq, "")
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

func main() {
	var treeChan <-chan tree.Trees
	var treeReader *bufio.Reader
	var treeFile io.Closer
	var err error
	var outalignfile, outtreefile, outratefile *os.File
	var snag Snag

	// Arguments
	ns := 4
	discrete := flag.Bool("discrete", true, "discrete gamma distribution")
	alpha := flag.Float64("alpha", 1.0, "gamma shape parameter")
	gamma := flag.Bool("gamma", true, "enable gamma distribution of site rates")
	ncat := flag.Int("gamma-cat", 4, "number of gamma categories")
	kappa := flag.Float64("kappa", 4.0, "Kappa parameter of K2P model")
	intree := flag.String("intree", "stdin", "Input tree")
	naligns := flag.Int("num-aligns", 1, "number of alignments to simulate per input tree")
	seed := flag.Int64("seed", time.Now().UTC().UnixNano(), "Random Seed parameter")
	l := flag.Int("length", 100, "Simulated alignment length")
	aa := flag.Bool("aa", false, "Simulate protein sequence")
	outalign := flag.String("out-align", "stdout", "Output alignment file")
	outtrees := flag.String("out-trees", "stdout", "Output tree file")
	outrates := flag.String("out-rates", "stdout", "Output site rates file")
	help := flag.Bool("help", false, "help")

	flag.Parse()

	// fmt.Println(*discrete)
	// fmt.Println(*alpha)
	// fmt.Println(*gamma)
	// fmt.Println(*ncat)
	// fmt.Println(*kappa)
	// fmt.Println(*intree)
	// fmt.Println(*naligns)
	// fmt.Println(*l)
	// fmt.Println(*aa)
	// fmt.Println(*outalign)
	// fmt.Println(*outtrees)
	// fmt.Println(*help)

	if *help {
		flag.Usage()
		return
	}

	// Create output files
	if *outalign == "stdout" {
		outalignfile = os.Stdout
	} else {
		if outalignfile, err = os.Create(*outalign); err != nil {
			return
		}
		defer outalignfile.Close()
	}
	if *outtrees == "stdout" {
		outtreefile = os.Stdout
	} else {
		if outtreefile, err = os.Create(*outtrees); err != nil {
			return
		}
		defer outtreefile.Close()
	}
	if *outrates == "stdout" {
		outratefile = os.Stdout
	} else {
		if outratefile, err = os.Create(*outrates); err != nil {
			return
		}
		defer outratefile.Close()
	}

	if *aa {
		ns = 20
	}

	rand.Seed(*seed)

	if !*gamma || !*discrete {
		*ncat = 1
	}

	// Parse Trees
	if treeFile, treeReader, err = utils.GetReader(*intree); err != nil {
		return
	}
	defer treeFile.Close()
	treeChan = utils.ReadMultiTrees(treeReader, utils.FORMAT_NEWICK)

	// Simulate alignments
	if snag, err = NewSnag(ns, *gamma, *discrete, *alpha, *ncat, *kappa, *naligns, *seed, *aa); err != nil {
		panic(err)
	}
	for trees := range treeChan {
		for si := range snag.Simulate(trees.Tree, *l) {
			fmt.Fprintln(outratefile, si.r)
			fmt.Fprintln(outtreefile, si.t.Newick())
			fmt.Fprintln(outalignfile, phylip.WriteAlignment(si.a, false, false, false))
		}
	}
}
