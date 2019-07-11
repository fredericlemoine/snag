package main

import (
	"bufio"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/phylip"
	"github.com/evolbioinfo/goalign/models"
	"github.com/evolbioinfo/goalign/models/dna"
	"github.com/evolbioinfo/goalign/models/protein"
	"github.com/evolbioinfo/gotree/io/utils"
	"github.com/evolbioinfo/gotree/tree"
)

func genrootseq(pi []float64, l int) (root []int) {
	root = make([]int, l)
	for c := 0; c < l; c++ {
		tmp := 0.0
		rc := rand.Float64()
		nt := -1
		for rc > tmp {
			nt++
			tmp += pi[nt]
		}
		root[c] = nt
	}
	return
}

func main() {
	var treeChan <-chan tree.Trees
	var treeReader *bufio.Reader
	var treeFile io.Closer
	var rootseq []int
	var pi []float64
	var err error
	var pij *models.Pij
	var pijs [][]*models.Pij
	var proba float64
	var seqs [][]int
	var nnodes int
	var cati int
	var sitecat int
	var tmp, rc float64
	var nt int
	var c int
	var ia int
	var rates []float64
	var cats []int
	var a align.Alignment
	var outalignfile, outtreefile *os.File

	ns := 4
	discrete := true
	alpha := 1.0
	gamma := true
	ncat := 4
	kappa := 4.0
	intree := os.Args[1]
	naligns, _ := strconv.Atoi(os.Args[2])
	seed, _ := strconv.ParseInt(os.Args[3], 10, 64)
	l, _ := strconv.Atoi(os.Args[4])
	aa, _ := strconv.ParseBool(os.Args[5])
	outalign := os.Args[6]
	outtrees := os.Args[7]

	if outalignfile, err = os.Create(outalign); err != nil {
		return
	}
	defer outalignfile.Close()
	if outtreefile, err = os.Create(outtrees); err != nil {
		return
	}
	defer outtreefile.Close()

	if aa {
		ns = 20
	}

	rand.Seed(seed)

	if !gamma || !discrete {
		ncat = 1
	}
	//fmt.Println(rates)

	if treeFile, treeReader, err = utils.GetReader(intree); err != nil {
		return
	}
	defer treeFile.Close()
	treeChan = utils.ReadMultiTrees(treeReader, utils.FORMAT_NEWICK)

	var m models.Model
	// Init frequencies (aa or nt)
	pi = make([]float64, ns)
	if !aa {
		// Initialize base frequencies at 1/4
		for i := 0; i < ns; i++ {
			pi[i] = 1. / float64(ns)
		}
		dm := dna.NewK2PModel()
		dm.InitModel(kappa)
		m = dm
	} else {
		if pm, err2 := protein.NewProtModel(protein.MODEL_LG, gamma, alpha); err2 != nil {
			return
		} else {
			// Initialize aa frequencies as defined in the model
			for i := 0; i < ns; i++ {
				pi[i] = pm.Pi(i)
			}
			pm.PrintMat()
			pm.InitModel(nil)
			pm.PrintMat()
			pm.PrintFreqs()
			m = pm
		}
	}

	for trees := range treeChan {
		nnodes = 0
		nedges := 0
		lengths := make([]float64, len(trees.Tree.Edges()))
		trees.Tree.PostOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) (keep bool) {
			cur.SetId(nnodes)
			if prev != nil {
				e.SetId(nedges)
				lengths[nedges] = e.Length()
				nedges++
			}
			nnodes++
			return true
		})
		seqs = make([][]int, nnodes)
		// Pij per rate category (if any) and per branch
		pijs = make([][]*models.Pij, ncat)

		for cati = 0; cati < ncat; cati++ {
			pijs[cati] = make([]*models.Pij, nnodes)
		}

		for ia = 0; ia < naligns; ia++ {
			rootseq = genrootseq(pi, l)
			rates, cats = models.GenerateRates(l, gamma, alpha, ncat, discrete)
			if aa {
				a = align.NewAlign(align.AMINOACIDS)
			} else {
				a = align.NewAlign(align.NUCLEOTIDS)
			}
			trees.Tree.PreOrder(func(cur *tree.Node, prev *tree.Node, e *tree.Edge) (keep bool) {
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
						if gamma && discrete {
							sitecat = cats[c]
						}
						// We initialize the pij for that branch and that category
						if pijs[sitecat][cur.Id()] == nil {
							pijs[sitecat][cur.Id()], _ = models.NewPij(m, lengths[e.Id()]*rates[c])
						}
						pij = pijs[sitecat][cur.Id()]
						// If not discrete, we must recompute the pijs
						// with given unique rate/brlen combo
						if gamma && !discrete {
							pij.SetLength(lengths[e.Id()] * rates[c])
						}
						tmp = 0.0
						rc = rand.Float64()
						nt = -1
						for rc > tmp {
							nt++
							proba = pij.Pij(prevseq[c], nt)
							tmp += proba
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
							if aa {
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
			fmt.Fprintln(outtreefile, trees.Tree.Newick())
			fmt.Fprintln(outalignfile, phylip.WriteAlignment(a, false, false, false))
		}
	}
}
