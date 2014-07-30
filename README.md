rbwa
====

R bindings for bwa ( https://github.com/lh3/bwa )

Code is inspired from Heng Li's https://github.com/lh3/bwa/blob/master/example.c

## Author

Pierre Lindenbaum PhD

@yokofakun


## Install and Compilation

```bash
$ git clone https://github.com/lindenb/rbwa.git
$ cd rbwa
```

download the submodule for bwa:

```
$ git submodule foreach git pull
```
then compile the bindings for R

```bash
$ make
```
or you can pass the path to the R include files by specifying **RCFLAGS**:

```bash
$ make -B librbwa.so RCFLAGS=" -I/commun/data/packages/R/R-3.1.1/R-3.1.1/include/ "
```


## Example

```bash
 make test
```


```R
> source("rbwa.R")
> 
> bwt <- bwa.open("test_files/chrM.fa")
> 
> for(s in c(
+ 		"GCATGTGTAATCTTACTAAGAGCTAATAGAAAGGCTAGGACCAAACCTAT",
+ 		"GCATGTGTAATCTTACTAAGCTAATAGAAAGGCTAGGACCAAACCTAT",
+ 		"CTATCTACTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTAT",
+ 		"TACTAAACCC",
+ 		"GCGAACCCAACTTCGATTCCCTCGCCGATCTCCGACGGAGCCGTGTGCAT"	
+ 		))
+ 	{
+ 	print(paste("QUERY:",s));
+ 	hits<-bwa.map(bwt,s)
+ 	print(hits)
+ 	}
[1] "QUERY: GCATGTGTAATCTTACTAAGAGCTAATAGAAAGGCTAGGACCAAACCTAT"
  chrom pos strand mapq NM secondary
1  chrM 650      1   60  0         0
[1] "QUERY: GCATGTGTAATCTTACTAAGCTAATAGAAAGGCTAGGACCAAACCTAT"
  chrom pos strand mapq NM secondary
1  chrM 650      1   60  2         0
[1] "QUERY: CTATCTACTTCAAATTCCTCCCTGTACGAAAGGACAAGAGAAATAAGGCCTCACAAAGCGCCTTCCCCCGTAAATGATATCATCTCAACTTAGTAT"
  chrom  pos strand mapq NM secondary
1  chrM 3100      0   60  4         0
[1] "QUERY: TACTAAACCC"
[1] chrom     pos       strand    mapq      NM        secondary
<0 rows> (or 0-length row.names)
[1] "QUERY: GCGAACCCAACTTCGATTCCCTCGCCGATCTCCGACGGAGCCGTGTGCAT"
[1] chrom     pos       strand    mapq      NM        secondary
<0 rows> (or 0-length row.names)
> 
> 
> bwa.close(bwt);
[1] TRUE
> 
> 
> 
```

