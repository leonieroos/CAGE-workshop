R exercise
========================================================
css: Rpress.css
author: Leonie Roos and Nevena Cvetesic
date: 11. 01. 2017.
autosize: true
width: 1440
height: 1100
css: stylesheet.css
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css

Overview
========================================================
- IRanges
- Genomic Ranges


  
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

```r
# we can create an index for subsetting
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
- this can be used to subset ovelapping ranges:  

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

```r
subject[subjectHits(hits_within)]
```

```
IRanges object with 3 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  p        18        28        11
  q         9        15         7
  p        18        28        11
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
  [1]     chr1  [10, 19]      + |      0.18
  [2]     chr1  [11, 20]      - |      0.42
  [3]     chr2  [12, 21]      + |       0.9
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

Genomic Ranges - GRanges
========================================================
- very useful way to create a GRanges object is from a dataframe
- first we will read data from a file *ranges.txt* into a dataframe and then create a GRanges object
























```
Error in file(file, "rt") : cannot open the connection
```
