R exercise
========================================================
css: Rpress.css
author: Leonie Roos and Nevena Cvetesic
date: 11. 01. 2017.
autosize: true
width: 1640
height: 1200
css: Rpress.css
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css: style.css

Overview
========================================================
type: section

Overview
========================================================
- IRanges
- Genomic Ranges

IRanges 
========================================================
type: section

IRanges
========================================================
- before learning to manipulate CAGE data, we will repeat some very important data structures and packages necessary for manipulation of genomic data

- first we will do some exercises using IRanges:  

```r
library(IRanges)
rng <- IRanges(start=1, end=15)
rng
```

```
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        15        15
```
- IRanges and GenomicRanges, are 1-based and use closed intervals, as most R packages 

IRanges
========================================================
- IRanges can also be created using start and width position:  

```r
rng <- IRanges(start=5, width=10)
```

- you can also create an IRanges object using a vector of arguments to specify start, end or width elements:  

```r
rng <- IRanges(start=c(2,3,7,15), end=c(15, 25, 46, 78))
rng
```

```
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         2        15        14
  [2]         3        25        23
  [3]         7        46        40
  [4]        15        78        64
```

IRanges
========================================================
- in addition, like other objects, ranges can be named:  

```r
names(rng) <- letters[1:4]
rng
```

```
IRanges object with 4 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         2        15        14
  b         3        25        23
  c         7        46        40
  d        15        78        64
```

IRanges
========================================================
- keep in mind that this type of object is not a dataframe and that the object's class determines its behaviour in R. (test its class)

```r
class(rng)
```

```
[1] "IRanges"
attr(,"package")
[1] "IRanges"
```

- to access the start, end and width positions, we can use accessor functions

```r
start(rng)
```

```
[1]  2  3  7 15
```

```r
end(rng)
```

```
[1] 15 25 46 78
```

```r
width(rng)
```

```
[1] 14 23 40 64
```

IRanges
========================================================
- subsetting works equally as for any other R object
- lets first create an IRange object wiith several ranges 

```r
rng <- IRanges(start=c(1,5,10), width=c(3,5,9))
rng
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1         3         3
  [2]         5         9         5
  [3]        10        18         9
```

```r
rng[2:3]
```

```
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         5         9         5
  [2]        10        18         9
```

IRanges
========================================================
- we can create an index for subsetting

```r
index <- start(rng) <=5
rng[index]
```

```
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1         3         3
  [2]         5         9         5
```

IRanges
========================================================
- the advantage of IRanges is the ability to transform our data. How addition or substraction work?  

```r
rng + 1
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         0         4         5
  [2]         4        10         7
  [3]         9        19        11
```

```r
rng - 1
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         2         2         1
  [2]         6         8         3
  [3]        11        17         7
```
- the added or substracted values have been applied to both the start and the end.  

IRanges
========================================================
- there are also special functions such as **flank**, lets see what that does:  

```r
rng
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1         3         3
  [2]         5         9         5
  [3]        10        18         9
```

```r
flank(rng, width = 2)
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        -1         0         2
  [2]         3         4         2
  [3]         8         9         2
```
- **flank** returns the regions that are on the side of each range in the IRanges object with the specified width, and by default it operates to create a range with 'width' position upstream of each start position

IRanges
========================================================
- to *flank* the regions downstream, set *start = false*:  

```r
flank(rng, width = 2, start = FALSE)
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4         5         2
  [2]        10        11         2
  [3]        19        20         2
```

IRanges
========================================================
- next important function is **reduce**
- this allows us to transform overlapping ranges into a set of non-overlapping ranges that cover the same positions 
- lets first create some overlapping ranges:  

```r
rng <- IRanges(start = c(45, 18, 27, 14), width = rep(5,4))
rng
```

```
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        45        49         5
  [2]        18        22         5
  [3]        27        31         5
  [4]        14        18         5
```

```r
reduce(rng)
```

```
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        14        22         9
  [2]        27        31         5
  [3]        45        49         5
```
- the function *reduce* is useful when we have a set of genomic regions and we are interested only in which regions are "covered", and not in the exact structure of the regions

IRanges
========================================================
- s similar function is *gap*, which returns the regions that are not covered 

```r
gaps(rng)
```

```
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        23        26         4
  [2]        32        44        13
```
- notice that the *gaps* function returns only the gaps between the ranges, and not from beginning of the sequence to the first range
- this is because we have not specified the beginning or the end of the sequence, if we wish to include those gaps also, we have to specify it:   

```r
gaps(rng, start = 1, end = 100)
```

```
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        13        13
  [2]        23        26         4
  [3]        32        44        13
  [4]        50       100        51
```

IRanges
========================================================
- further useful functions are **set operations**
- IRanges objects can be percieved as set of integers (4-7 is 4,5,6,7) and therefore difference (setdiff()), intersection (intersect()) and union (union()) are useful operations.  
- note that in some cases, the order of IRanges objects matters.  


```r
rng1 <- IRanges(start = 3, end = 15)
rng2 <- IRanges(start = 13, end = 19)
  
union(rng1, rng2)
```

```
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         3        19        17
```

```r
intersect(rng1, rng2)
```

```
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        13        15         3
```

IRanges
========================================================
- note that in some cases, the order of IRanges objects matters.  

```r
setdiff(rng1, rng2)
```

```
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         3        12        10
```

```r
setdiff(rng2, rng1)
```

```
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        16        19         4
```

IRanges
========================================================
- lastly, one of the most important functions is **findOverlaps** 
- the basic task of this function is to find overlaps between two sets of IRanges objects. Lets see how it works:    

```r
query <- IRanges(start = c(1, 27, 18, 12, 22, 8), end = c(15, 29, 18, 15, 23, 8), names = letters[1:6])
query
```

```
IRanges object with 6 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        15        15
  b        27        29         3
  c        18        18         1
  d        12        15         4
  e        22        23         2
  f         8         8         1
```

```r
subject <- IRanges(start = c(1,18,9), end = c(4, 28, 15), names = letters[15:17])
subject
```

```
IRanges object with 3 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  o         1         4         4
  p        18        28        11
  q         9        15         7
```

```r
hits <- findOverlaps(query, subject)
```

IRanges
========================================================

```r
hits
```

```
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         1           1
  [2]         1           3
  [3]         2           2
  [4]         3           2
  [5]         4           3
  [6]         5           2
  -------
  queryLength: 6 / subjectLength: 3
```

- hits object represents a mapping between query and subject - which query ranges overlap wich subject range/ranges
- first column represents the index of the query object, while the second the index of the subject

IRanges
========================================================
- hits object can be used to subset ovelapping ranges:  

```r
query[queryHits(hits)]
```

```
IRanges object with 6 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        15        15
  a         1        15        15
  b        27        29         3
  c        18        18         1
  d        12        15         4
  e        22        23         2
```

```r
subject[subjectHits(hits)]
```

```
IRanges object with 6 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  o         1         4         4
  q         9        15         7
  p        18        28        11
  p        18        28        11
  q         9        15         7
  p        18        28        11
```

IRanges
========================================================
- similarly, we can also use *subsetByOverlaps*:  

```r
subsetByOverlaps(query, subject)
```

```
IRanges object with 5 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        15        15
  b        27        29         3
  c        18        18         1
  d        12        15         4
  e        22        23         2
```

IRanges
========================================================
- what is the difference between using *subsetByOverlaps* and the following:  

```r
query[queryHits(hits)]
```

```
IRanges object with 6 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        15        15
  a         1        15        15
  b        27        29         3
  c        18        18         1
  d        12        15         4
  e        22        23         2
```

```r
query[unique(queryHits(hits))]
```

```
IRanges object with 5 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        15        15
  b        27        29         3
  c        18        18         1
  d        12        15         4
  e        22        23         2
```

IRanges
========================================================
- from the output we can see that *findOverlaps* by default works to find any overlaping ranges, even if just a part overlaps
- we can also specify that a query is overlapping only if it completely falls within the subject range
  

```r
hits_within <- findOverlaps(query, subject, type = "within")
hits_within
```

```
Hits object with 3 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         3           2
  [2]         4           3
  [3]         5           2
  -------
  queryLength: 6 / subjectLength: 3
```

```r
query[queryHits(hits_within)]
```

```
IRanges object with 3 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  c        18        18         1
  d        12        15         4
  e        22        23         2
```

IRanges
========================================================
- last set of functions we will examine is to find nearest, preceding or following ranges and calculate distances.  

```r
query <- IRanges(start = 2, end = 7, name = "query")
query
```

```
IRanges object with 1 range and 0 metadata columns:
            start       end     width
        <integer> <integer> <integer>
  query         2         7         6
```

```r
subject
```

```
IRanges object with 3 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  o         1         4         4
  p        18        28        11
  q         9        15         7
```

```r
nearest(query, subject)
```

```
[1] 1
```

IRanges
========================================================

```r
query[nearest(query, subject)]
```

```
IRanges object with 1 range and 0 metadata columns:
            start       end     width
        <integer> <integer> <integer>
  query         2         7         6
```

```r
precede(query, subject)
```

```
[1] 3
```

```r
follow(query, subject)
```

```
[1] NA
```
- nearest returns the nearest subject range even if its overlapping the query range
- precede and follow find ranges that query precedes or follows (remember these functions are relative to the query)

IRanges
========================================================
- to calculate some distances between ranges, we can use the following:  

```r
query
```

```
IRanges object with 1 range and 0 metadata columns:
            start       end     width
        <integer> <integer> <integer>
  query         2         7         6
```

```r
subject
```

```
IRanges object with 3 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  o         1         4         4
  p        18        28        11
  q         9        15         7
```

```r
distanceToNearest(query, subject)
```

```
Hits object with 1 hit and 1 metadata column:
      queryHits subjectHits |  distance
      <integer>   <integer> | <integer>
  [1]         1           1 |         0
  -------
  queryLength: 1 / subjectLength: 3
```

IRanges
========================================================

```r
distance(query, subject)
```

```
[1]  0 10  1
```
- we can see that **distanceToNearest** works similarly as **findOverlaps** - it returns a mapping hits object and an additional column containing calculated distances 
- distances returns pairwise calculated distances between each query and subject ranges

Genomic Ranges - GRanges
========================================================
type: section

Genomic Ranges - GRanges
========================================================
- GRanges is build on IRanges, with additional functionalities such as storing sequence name (chromosome names), strand name (+ or -) and any other information in the form of metadata.  
- lets see how to create a GRanges object.  

```r
library(GenomicRanges)
gr <- GRanges(seqname = c("chr1", "chr1", "chr2"), 
              ranges = IRanges(start = 10:12, width = 10),
              strand = c("+", "-", "+"),
              gc_perc = round(runif(3), 2))
gr
```

```
GRanges object with 3 ranges and 1 metadata column:
      seqnames    ranges strand |   gc_perc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1  [10, 19]      + |      0.99
  [2]     chr1  [11, 20]      - |       0.3
  [3]     chr2  [12, 21]      + |      0.34
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRanges
========================================================
- very useful way to create a GRanges object is from a dataframe
- first we will read data from a file *ranges.txt* into a dataframe and then create a GRanges object

```r
rng <- read.table(file = "../Tutorials/ranges.txt", sep = "\t", header = TRUE)
rng
```

```
   chr start end strand gc_perc
1 chr1     0  16      -    0.78
2 chr2    28  35      -    0.54
3 chr2    12  22      +    0.34
4 chr3     7  14      +    0.23
5 chr3     5  15      -    0.66
6 chr7     7  23      +    0.65
```

Genomic Ranges - GRanges
========================================================
- we can specify everything manually

```r
rng.gr <- GRanges(seqnames = rng$chr,
                  ranges = IRanges(start = rng$start, end = rng$end),
                  strand = rng$strand,
                  gc_perc = rng$gc_perc)
rng.gr
```

```
GRanges object with 6 ranges and 1 metadata column:
      seqnames    ranges strand |   gc_perc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1  [ 0, 16]      - |      0.78
  [2]     chr2  [28, 35]      - |      0.54
  [3]     chr2  [12, 22]      + |      0.34
  [4]     chr3  [ 7, 14]      + |      0.23
  [5]     chr3  [ 5, 15]      - |      0.66
  [6]     chr7  [ 7, 23]      + |      0.65
  -------
  seqinfo: 4 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRanges
========================================================
- there is also an option to convert the dataframe into a GRanges object directly using makeGRangesFromDataFrame

```r
rng.gr <- makeGRangesFromDataFrame(rng, keep.extra.columns = TRUE)
rng.gr
```

```
GRanges object with 6 ranges and 1 metadata column:
      seqnames    ranges strand |   gc_perc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1  [ 0, 16]      - |      0.78
  [2]     chr2  [28, 35]      - |      0.54
  [3]     chr2  [12, 22]      + |      0.34
  [4]     chr3  [ 7, 14]      + |      0.23
  [5]     chr3  [ 5, 15]      - |      0.66
  [6]     chr7  [ 7, 23]      + |      0.65
  -------
  seqinfo: 4 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRanges
========================================================
- similar to IRanges, we can use accessor functions (everything that worked on IRanges will work on GRanges objects):

```r
start(rng.gr)
```

```
[1]  0 28 12  7  5  7
```

```r
end(rng.gr)
```

```
[1] 16 35 22 14 15 23
```

```r
width(rng.gr)
```

```
[1] 17  8 11  8 11 17
```

```r
seqnames(rng.gr)
```

```
factor-Rle of length 6 with 4 runs
  Lengths:    1    2    2    1
  Values : chr1 chr2 chr3 chr7
Levels(4): chr1 chr2 chr3 chr7
```

Genomic Ranges - GRanges
========================================================

```r
ranges(rng.gr)
```

```
IRanges object with 6 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         0        16        17
  [2]        28        35         8
  [3]        12        22        11
  [4]         7        14         8
  [5]         5        15        11
  [6]         7        23        17
```

```r
strand(rng.gr)
```

```
factor-Rle of length 6 with 4 runs
  Lengths: 2 2 1 1
  Values : - + - +
Levels(3): + - *
```

Genomic Ranges - GRanges
========================================================
- we can also see the number of ranges in our GRanges object and name each range:  

```r
length(rng.gr)
```

```
[1] 6
```

```r
names(rng.gr) <- letters[1:6]
rng.gr
```

```
GRanges object with 6 ranges and 1 metadata column:
    seqnames    ranges strand |   gc_perc
       <Rle> <IRanges>  <Rle> | <numeric>
  a     chr1  [ 0, 16]      - |      0.78
  b     chr2  [28, 35]      - |      0.54
  c     chr2  [12, 22]      + |      0.34
  d     chr3  [ 7, 14]      + |      0.23
  e     chr3  [ 5, 15]      - |      0.66
  f     chr7  [ 7, 23]      + |      0.65
  -------
  seqinfo: 4 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRanges
========================================================
- to access metadata use:  

```r
mcols(rng.gr)
```

```
DataFrame with 6 rows and 1 column
    gc_perc
  <numeric>
1      0.78
2      0.54
3      0.34
4      0.23
5      0.66
6      0.65
```

Genomic Ranges - GRanges
========================================================
- important function that will be often used with genomics data is coverage  
- coverage calculates how many individual positions are overlapped with the ranges within the .gr object. Lets see how it works:    

```r
# first check you .gr object
rng.gr
```

```
GRanges object with 6 ranges and 1 metadata column:
    seqnames    ranges strand |   gc_perc
       <Rle> <IRanges>  <Rle> | <numeric>
  a     chr1  [ 0, 16]      - |      0.78
  b     chr2  [28, 35]      - |      0.54
  c     chr2  [12, 22]      + |      0.34
  d     chr3  [ 7, 14]      + |      0.23
  e     chr3  [ 5, 15]      - |      0.66
  f     chr7  [ 7, 23]      + |      0.65
  -------
  seqinfo: 4 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRanges
========================================================

```r
coverage(rng.gr)
```

```
RleList of length 4
$chr1
integer-Rle of length 16 with 1 run
  Lengths: 16
  Values :  1

$chr2
integer-Rle of length 35 with 4 runs
  Lengths: 11 11  5  8
  Values :  0  1  0  1

$chr3
integer-Rle of length 15 with 4 runs
  Lengths: 4 2 8 1
  Values : 0 1 2 1

$chr7
integer-Rle of length 23 with 2 runs
  Lengths:  6 17
  Values :  0  1
```

Genomic Ranges - GRangesList
========================================================
type: sub-section


Genomic Ranges - GRangesList
========================================================
- a very usefull GRanges structure is a GRangesList that can be used to group the data together - this will be frequently used with genomic data.  


```r
gr1 <- GRanges(seqname = c("chr1", "chr2", "chr2"), 
               ranges = IRanges(start = 2:4, width = 10),
               strand = c("+", "+", "-"))
gr1
```

```
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [2, 11]      +
  [2]     chr2   [3, 12]      +
  [3]     chr2   [4, 13]      -
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRangesList
========================================================

```r
gr2 <- GRanges(seqname = c("chr2", "chr2", "chr3", "chr1"), 
               ranges = IRanges(start = 6:9, width = 12),
               strand = c("+", "+", "-","+"))
gr2
```

```
GRanges object with 4 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr2   [6, 17]      +
  [2]     chr2   [7, 18]      +
  [3]     chr3   [8, 19]      -
  [4]     chr1   [9, 20]      +
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRangesList
========================================================

```r
grl <- GRangesList(gr1, gr2)
grl
```

```
GRangesList object of length 2:
[[1]] 
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [2, 11]      +
  [2]     chr2   [3, 12]      +
  [3]     chr2   [4, 13]      -

[[2]] 
GRanges object with 4 ranges and 0 metadata columns:
      seqnames  ranges strand
  [1]     chr2 [6, 17]      +
  [2]     chr2 [7, 18]      +
  [3]     chr3 [8, 19]      -
  [4]     chr1 [9, 20]      +

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRangesList
========================================================
- GRangesList behaves as a regular R list

```r
length(grl)
```

```
[1] 2
```

```r
unlist(grl)
```

```
GRanges object with 7 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [2, 11]      +
  [2]     chr2   [3, 12]      +
  [3]     chr2   [4, 13]      -
  [4]     chr2   [6, 17]      +
  [5]     chr2   [7, 18]      +
  [6]     chr3   [8, 19]      -
  [7]     chr1   [9, 20]      +
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRangesList
========================================================

```r
grl[1]
```

```
GRangesList object of length 1:
[[1]] 
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [2, 11]      +
  [2]     chr2   [3, 12]      +
  [3]     chr2   [4, 13]      -

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

```r
grl[[1]]
```

```
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [2, 11]      +
  [2]     chr2   [3, 12]      +
  [3]     chr2   [4, 13]      -
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRangesList
========================================================
- to work with lists as data structures, we will often need to use *lapply* or *sapply*
- these functions enable us to iterate through all the elements of a list and apply a function to each element
- they work with R's regular lists, and also with GRangesLists. Lets test how it works:  

```r
# first check your GRangesList object
grl
```

```
GRangesList object of length 2:
[[1]] 
GRanges object with 3 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1   [2, 11]      +
  [2]     chr2   [3, 12]      +
  [3]     chr2   [4, 13]      -

[[2]] 
GRanges object with 4 ranges and 0 metadata columns:
      seqnames  ranges strand
  [1]     chr2 [6, 17]      +
  [2]     chr2 [7, 18]      +
  [3]     chr3 [8, 19]      -
  [4]     chr1 [9, 20]      +

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRangesList
========================================================

```r
# extract start positions for each GRangesList element using lapply
lapply(grl, function(x) start(x))
```

```
[[1]]
[1] 2 3 4

[[2]]
[1] 6 7 8 9
```

Genomic Ranges - GRangesList
========================================================

```r
# extract width of each element
lapply(grl, function(x) width(x))
```

```
[[1]]
[1] 10 10 10

[[2]]
[1] 12 12 12 12
```
What is the difference between lapply and sapply? Test the following:  

```r
# extract the number of ranges in each list elements
lapply(grl, length)
```

```
[[1]]
[1] 3

[[2]]
[1] 4
```

```r
sapply(grl, length)
```

```
[1] 3 4
```

Genomic Ranges - Exercise
========================================================
type: prompt

Genomic Ranges - Exercise
========================================================
  1. Import *Danio rerio* gene coordinates from *danRer7genes.txt* file and create a GRanges object.
  2. Create a GRangesList object with genes split according to chromosome name (each list element contains all genes located on the corresponding chromosome).
  3. Create a new GRanges object containing promoters of genes located on chromosome 11 (see *promoters* function).

Genomic Ranges - Solutions
========================================================
type: sub-section

Genomic Ranges - Solutions
========================================================

1) Import *Danio rerio* gene coordinates from *danRer7genes.txt* file and create a GRanges object.  
  

```r
# read the data from a file into a dataframe
genes.df <- read.table("../Tutorials/danRer7genes.txt", header = TRUE) 
head(genes.df)
```

```
  seqnames    start      end strand            gene_id
1     chr9 35083371 35093143      - ENSDARG00000000001
2     chr9 35060460 35084513      + ENSDARG00000000002
3     chr4 14146777 14169064      - ENSDARG00000000018
4     chr4 14076733 14125268      + ENSDARG00000000019
5    chr12 34979366 35032034      + ENSDARG00000000068
6    chr24 22678238 22711533      - ENSDARG00000000069
```

Genomic Ranges - Solutions
========================================================
1) Import *Danio rerio* gene coordinates from *danRer7genes.txt* file and create a GRanges object.  
  

```r
# create a GRanges object from a dataframe
genes.gr <- makeGRangesFromDataFrame(genes.df, keep.extra.columns = TRUE) 
head(genes.gr)
```

```
GRanges object with 6 ranges and 1 metadata column:
    seqnames               ranges strand |            gene_id
       <Rle>            <IRanges>  <Rle> |           <factor>
  1     chr9 [35083371, 35093143]      - | ENSDARG00000000001
  2     chr9 [35060460, 35084513]      + | ENSDARG00000000002
  3     chr4 [14146777, 14169064]      - | ENSDARG00000000018
  4     chr4 [14076733, 14125268]      + | ENSDARG00000000019
  5    chr12 [34979366, 35032034]      + | ENSDARG00000000068
  6    chr24 [22678238, 22711533]      - | ENSDARG00000000069
  -------
  seqinfo: 499 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - Solutions
========================================================
2) Create a GRangesList object with genes split according to chromosome name (each list element contains all genes located on the corresponding chromosome).  
 

```r
genes.grl <- split(genes.gr, seqnames(genes.gr)) 
genes.grl
```

```
GRangesList object of length 499:
$chr1 
GRanges object with 1494 ranges and 1 metadata column:
        seqnames               ranges strand |            gene_id
           <Rle>            <IRanges>  <Rle> |           <factor>
      7     chr1 [47492921, 47586925]      - | ENSDARG00000000086
     43     chr1 [20837616, 21057496]      + | ENSDARG00000000606
     69     chr1 [ 1318859,  1331666]      - | ENSDARG00000001009
     96     chr1 [47613228, 47619236]      - | ENSDARG00000001470
    139     chr1 [ 1446561,  1462945]      + | ENSDARG00000001870
    ...      ...                  ...    ... .                ...
  32356     chr1 [25633941, 25636188]      + | ENSDARG00000096667
  32357     chr1 [29774503, 29779515]      - | ENSDARG00000096668
  32359     chr1 [57591358, 57593975]      + | ENSDARG00000096670
  32360     chr1 [57874379, 57886230]      - | ENSDARG00000096671
  32887     chr1 [29886254, 29887453]      - | ENSDARG00000097202

...
<498 more elements>
-------
seqinfo: 499 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - Solutions
========================================================
3) Create a new GRanges object containing promoters of genes located on chromosome 11 (see *promoters* function).  

```r
promoters.gr <- promoters(genes.grl[["chr11"]], upstream = 2000, downstream = 500, strand = TRUE)
head(promoters.gr)
```

```
GRanges object with 6 ranges and 1 metadata column:
      seqnames               ranges strand |            gene_id
         <Rle>            <IRanges>  <Rle> |           <factor>
   29    chr11 [24878281, 24880780]      + | ENSDARG00000000472
  111    chr11 [ 8401538,  8404037]      - | ENSDARG00000001712
  124    chr11 [ 8273738,  8276237]      + | ENSDARG00000001782
  149    chr11 [32100656, 32103155]      - | ENSDARG00000001897
  166    chr11 [12859556, 12862055]      - | ENSDARG00000001969
  168    chr11 [21594167, 21596666]      - | ENSDARG00000001972
  -------
  seqinfo: 499 sequences from an unspecified genome; no seqlengths
```
