Errors can be introduced to sequenced molecules in multiple stages during sequencing experiments. Typically, the primary stages where errors typically occur include PCR amplification, sequencing, and bead synthesis. Inaccuracies are

* **Substitution errors**: Incorrectly identifying one nucleotide as another. Mistakes made by DNA polymerase during amplification (PCR). Base calling failure (sequencing)
* **Insertion errors**: Adding extra nucleotides that are not present in the original sequence. 
* **Deletion errors**: Deleting nucleotides that are present in the original sequence.

To simulate reads as more genuinely as possible, Tresor implements a powerful function to allow any types of errors to bo added during PCR amplification, sequencing, and bead synthesis. Users can pass certain error rates onto the following attributes. Then, Tresor will dispose of the error-adding process.

1. Bead synthesis
``` py
# substitution error
bead_mutation=True, # True False
bead_mut_rate=1e-4, # 0.016 0.00004
# deletion
bead_deletion=True, # True False
bead_del_rate=0.1,  # 0.016 0.00004, 2.4e-7
# insertion
bead_insertion=True, # True False
bead_ins_rate=7.1e-7,  # 0.011 0.00001, 7.1e-7
```

2. PCR
``` py
# substitution error
pcr_err=0.00001
# deletion
pcr_deletion=True, # False True
pcr_del_rate=2.4 * 10e-6,
# insertion
pcr_insertion=True, # False True
pcr_ins_rate=7.1 * 10e-7,
```

3. Sequencing
``` py
# substitution error
seq_err=0.001,
# deletion
seq_deletion=True, # False True
seq_del_rate=2.4 * 10e-6,
# insertion
seq_insertion=True, # False True
seq_ins_rate=7.1 * 10e-7,
```