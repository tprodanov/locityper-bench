# Locityper benchmarking scripts

Let us assume that [base repository](https://github.com/tprodanov/locityper) is installed at `~/Code/locityper`
and this repository is cloned to `~/Code/locityper-bench`.

## Accuracy evaluation

First, predictions need to be converted to CSV format:
```bash
~/Code/locityper/extra/into_csv.py -i analysis/./* -o summary.csv.gz
```
Haplotypes alignments are constructed using `locityper align`
```bash
locityper align -i loci/$locus/haplotypes.fa.gz -o loci/$locus/haplotypes.paf.gz -A -D 1
```

Next, `eval_accuracy` script is called:
```bash
~/Code/locityper/extra/eval_accuracy.py \
    -i summary.csv.gz -o eval.csv.gz \
    -a db/loci/{}/haplotypes.paf.gz -d dbs/327 [--loo]
```

## Trio concordance calculation

To calculate trio concordance, we use `trio_conc.py` and 1KGP pedigree file.

```bash
~/Code/locityper/extra/trio_conc.py -i summary.csv.gz -p 1KGP.ped \
    -a db/loci/{}/haplotypes.paf.gz -o trio_conc.csv.gz
```

## Comparing VCF files

In order to convert Locityper predictions into VCF, we use the following command:
```bash
~/Code/locityper/extra/into_vcf.py -i summary.csv.gz -d db \
    -v hprc-v1.1-grch38.vcf.gz -g GRCh38 -o locityper-vcfs
```

To extract 1KGP call set haplotypes, we used the following commands:
```bash
bcftools view -s $sample -ac1 -R locus.bed -e 'ALT~"<.*>" && ALT != "<DEL>"' 1KGP.vcf.gz
```
and reconstructed haplotypes using `bcftools consensus`.
Next, accuracy evaluation, similar to Locityper, can be performed for the reconstructed haplotypes.

In order to compare VCF files, we first decomposed HPRC, Pangenie and Locityper VCF files:
```bash
for i in 0 1; do
    # Start with `vt decompose_blocksub -a`, then try without `-a`.
    if [[ i -eq 0 ]]; then
        args=(-a)
    else
        args=()
    fi
    if (bcftools norm -m-any -f "$genome/genome.fa" "$invcf" | \
        vt decompose_blocksub ${args[@]+"${args[@]}"} - | \
        vt normalize -r "$genome/genome.fa" - | \
        bcftools sort - | \
        vt uniq -o "$outvcf" -)
    then
        tabix -p vcf "$outvcf"
        return 0
    fi
done
```
1KGP call sets were not decomposed, as non-decomposed variants showed higher accuracy.
Then, RTG-tools were used to compare two call sets:
```bash
rtg vcfeval -b base.vcf.gz -c calls.vcf.gz -e locus.bed -t genome.fa.sdf -o out
```
Finally, we combined `vcfeval` results using:
```bash
~/Code/locityper-bench/combine_vcfevals.py \
    -i loci/{locus}/evals/{tool}.{sample} -o combined_evals.csv.gz \
    -b loci/{locus}/decompose/baseline.{sample}.vcf.gz -c loci/{locus}/decompose/{tool}.{sample}.vcf.gz \
    -t 0,1,10
```

## Evaluating MHC/KIR-calling

We used existing Immuannot annotations for the HPRC samples, combined them with
```bash
~/Code/locityper-bench/hla/convert_annot.py -a Immuannot/{H,N,CHM13,GRCh38}*.gtf.gz -o Immuannot/annot.csv
```

Then, we converted Locityper and T1K predictions to the same format:
```bash
~/Code/locityper-bench/hla/convert_locityper.py -i summary.csv.gz -l loci.txt -a Immuannot/annot.csv \
    -o locityper.csv
~/Code/locityper-bench/hla/convert_locityper.py -i summary_loo.csv.gz -l loci.txt -a Immuannot/annot.csv \
    -o locityper_loo.csv
~/Code/locityper-bench/hla/convert_t1k.py -t t1k/analysis/*/./* -o t1k.csv
```
where `loci.txt` represented two-column file with HLA/KIR genes and the corresponding Locityper loci.

Finally, we evaluated accuracy with
```bash
~/Code/locityper-bench/hla/accuracy.py -b Immuannot/annot.csv \
    -p t1k.csv -o eval/t1k.csv
~/Code/locityper-bench/hla/accuracy.py -b Immuannot/annot.csv \
    -p locityper.csv -o eval/locityper.csv -a Immuannot/annot.csv
~/Code/locityper-bench/hla/accuracy.py -b Immuannot/annot.csv \
    -p locityper_loo.csv -o eval/locityper_loo.csv -a Immuannot/annot.csv --loo
```
