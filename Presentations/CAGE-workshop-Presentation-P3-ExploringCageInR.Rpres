Exploring CAGE data in R
========================================================
css: style/Rpress.css
author: Leonie Roos and Nevena Cvetesic
date: 11. 01. 2017.
autosize: true
width: 1440
height: 1100
css: stylesheet.css
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css: style/style.css

Overview
========================================================
type: sub-section
<br>
- Data formats for downstream manipulation
- Annotation for genomic features
- Promoter classification
- Oligonucleotide heatmaps

What can we get out of a CAGEset?
========================================================
<br>
__dataframes__ & __matrices__ returned by _CAGEr_ functions:

- Per individual sample ( _raw_ or _normalized_ tag counts)
  
  - CTSS 
  - Tag clusters & interquantile widths
  - Consensus clusters & interquantile widths
<br>    
- Across samples (one output for multiple samples)

  - Consensus clusters & interquantile widths
<br>

These can be easily manipulated by the user in R.
<br>
> We'll go through this in more detail in the tutorial


CAGE data per sample
========================================================
<br>
Other than consensus clusters most data is extracted per sample as they are __sample specific__ 

- genomic coordinates of CTSS
- genomic coordinates of Tag clusters

<br>
Likely you want to perform the same analyses, checks, plots, etc for all the samples.
<br> <br>
_Solution?_ <br><br>
Automate what you're doing! Functions are great to avoid mistakes

- create a list for example that contain all the info of samples per slot

> We'll create a few of these in the tutorial


CAGE data - Genomic Features
========================================================
<br>
Where does most of the signal fall? <br>
![](images/genomic_feat.png)

<br>
Nepal, C., et al. (2013). Dynamic regulation of coding and non-coding transcription initiation landscape at single nucleotide resolution during vertebrate embryogenesis. Genome Research, 23(11):1938-1950. 


CAGE data - Genomic Features
========================================================
<br>
we count the overlap our TCs per sample for each feature and count the occurances per sample:
<br>
![](images/genome-feat-our.png)

* promoter
* 5kb upstream of promoter
* exons
* introns
* gene


Promoter Width
========================================================
<br>
__Promoter interquantile widths__

- determined earlier in _CAGEr_
- set between q0.1 - 0.9

<br>
__ggplot2__ to create other types of graphs than _CAGEr_ offers

![](images/dens-tc-width.png)

Sequence Features of Core Promoters
========================================================
<br>
__Distinguishing features of sharp and broad promoters__

* Genomic position of dominant TSS
* Add a window (200 bp for example) around it
* Plot the average oligonucleotide profiles

![](images/seqpattern.png)
Based on interquantile widths: __sharp < 10__ & __broad >= 10__

Taken from R package: [seqPattern](http://bioconductor.org/packages/release/bioc/html/seqPattern.html).


Heatmaps
========================================================
<br>
Another great visualization tool is density heatmaps: <br>

__These plot the  density  of  oligonucleotides__ 
such as TA, CG, WW, SS

![](images/heatmaps-ordered-iq.png)

__Striped line indicates dominant CTSS)__

- 400 bp up and downstream
- Sequences are ordered by interquantile width of Tag cluster

Produced by [heatmaps](https://github.com/mgperry/heatmaps) R package


Heatmaps - Order is important
========================================================

The order of sequences is important when looking for architectures:

![](images/clus-shift-heatmap.png)
<br>
__Similar sample but ordered differently__


Tutorial
========================================================
type: section
<br>
The tutorial that is linked with this presentation:
<br>

__Tutorials dir__ 

CAGE-workshop-Tutorial-P3-ExploringCageInR


