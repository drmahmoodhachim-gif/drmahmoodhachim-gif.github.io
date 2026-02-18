# BLAST Activity – Sequence Similarity Search

**What you'll learn:** Use BLAST to find similar sequences in databases. No software to install – we use the free NCBI BLAST website.

---

## What Is BLAST?

BLAST (Basic Local Alignment Search Tool) compares your sequence (DNA or protein) against huge databases to find similar sequences. Use it to:

- Identify unknown sequences
- Find homologs in other species
- Check if a sequence matches a known gene

---

## Activity: Run BLAST on Example Sequences

### Step 1: Open NCBI BLAST

Go to: **https://blast.ncbi.nlm.nih.gov/Blast.cgi**

---

### Step 2: Choose Nucleotide BLAST

Click **"Nucleotide BLAST"** (blastn) – we'll search DNA/RNA sequences.

---

### Step 3: Enter Your Sequence

Open the file **`sequences_for_blast.fa`** in this folder. Copy one of the sequences (the line of A, C, G, T – not the `>` line) and paste it into the **"Enter Query Sequence"** box.

**Or use this example** (copy from the second line):

```
>Example_gene_fragment
ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCC
```

Copy only: `ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCCC`

---

### Step 4: Set Search Options

- **Database:** Nucleotide collection (nr/nt)
- **Organism:** Optional – type "Homo sapiens" to search only human sequences
- Leave other options as default

---

### Step 5: Click BLAST

Click the blue **"BLAST"** button. Wait 10–30 seconds.

---

### Step 6: Interpret Results

You'll see:

- **Description** – Name of matching sequences (e.g. gene names, species)
- **Alignments** – How your sequence aligns to matches; **Identities** = % similarity
- **E-value** – Lower = better match (e.g. 1e-50 is very strong)

**Things to check:**
- Top hit – which gene does it match?
- Percent identity – how similar?
- E-value – is it significant? (&lt; 0.01 is usually good)

---

## Try Sequences Related to Your Project

We've included example sequences in **`sequences_for_blast.fa`**:

| Sequence        | Related to                          |
|-----------------|-------------------------------------|
| SCN1A fragment  | Epilepsy (Mutho)                   |
| INS fragment    | Diabetes (Maryam)                  |
| NKX2-5 fragment | Cardiac genetics (Tuleen)           |
| IL4 fragment    | Asthma (Yara)                       |

Blast each and see what they match in the database.

---

## Brief Report to Add to Your Submission

After running BLAST, add 1–2 sentences to your submission:

- *"I BLASTed the SCN1A sequence and it matched sodium channel genes in humans and mice, with 98% identity. This fits my epilepsy project because SCN1A is a known epilepsy gene."*

---

## Link

**NCBI BLAST:** https://blast.ncbi.nlm.nih.gov/Blast.cgi
