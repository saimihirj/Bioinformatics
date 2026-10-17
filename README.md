# Bioinformatics

My solutions to *Biology Meets Programming: Bioinformatics for Beginners*, the University of California San Diego course on Coursera taught by Pavel Pevzner and Phillip Compeau (completed August 2021; [certificate](Certificate.pdf)).

The course asks two biological questions and answers them with short algorithms written from scratch in Python.

| Lab | Question | What it implements |
|---|---|---|
| [Week 1](Labs/week1.py) | Where in the genome does replication begin? | Pattern counting, most frequent k-mers, reverse complements |
| [Week 2](Labs/week2.py) | Where in the genome does replication begin? (continued) | Symbol arrays and the GC skew, whose minimum points to the replication origin |
| [Week 3](Labs/week3.py) | Which DNA patterns act as molecular clocks? | Motif counts, profile matrices, consensus strings and scores, greedy motif search |
| [Week 4](Labs/week4.py) | Which DNA patterns act as molecular clocks? (continued) | Pseudocounts, randomised motif search and Gibbs sampling |

## Data

| File | Contents |
|---|---|
| `Datasets/E_coli.txt` | The *E. coli* genome (4.6 million bases) |
| `Datasets/DosR.txt` | Upstream regions of genes regulated by DosR in *Mycobacterium tuberculosis* |

## Trying it

The labs use only the Python standard library. For example, the minimum GC skew of the *E. coli* genome falls at position 3,923,620, close to its known origin of replication:

```python
import sys; sys.path.insert(0, "Labs")
import week2
genome = open("Datasets/E_coli.txt").read().strip()
print(week2.MinimumSkew(genome)[:3])   # [3923620, 3923621, 3923622]
```

## License

MIT. See [LICENSE](LICENSE).
