# lcQTH

Rapid quantitative trait mapping through tracing parental haplotype with ultra-low-coverage sequencing

View homepage at [https://esctrionsit.github.io/lcQTH](https://esctrionsit.github.io/lcQTH)

## Recommanded environments with version tested:

+ Linux system with 64 bit CPU
+ Bcftools v1.11
+ Bedtools v2.29.2
+ python v3.7.6
+ python packages:
  + hmmlearn v0.2.4
  + scikit_learn v1.3.2
  + openpyxl v3.0.7
+ R v3.6
+ R packages:
  + ComplexHeatmap v2.0.0
  + qtl v1.46-2
+ Golang 1.20.3 linux/amd64
+ gcc version 11.2.1 20220127 (Red Hat 11.2.1-9) (GCC)
+ Node v14.17.5
+ Node packages:
  + d3@4.13.0
  + jsdom@15.2.1
  + jspdf@2.5.1
  + svg2pdf.js@2.2.1

## Installation

Clone this git repository, and run `python3 lcQTH init --config=[path]` with a configuration file.

## Usage

```text
Usage: python3 lcQTH.py [options] COMMAND

Toolkit for tracing parental haplotypes and perform simple QTL mapping with aligned ultra-low-coverage sequencing data.

Options:
  --help            Show this message and exit.

Commands:
  init              Initialize this lcQTH project, cleaning previous results.
  run               Run the analysis program continually (with no options) or from a specific step.
  exporthmp         Export parental haplotypes of RILs in hapmap format.
  exportxlsx        Export parental haplotypes of RILs in hapmap format and saved in an Excel spreadsheet

Run 'python3 lcQTH.py COMMAND --help' for more information on a command.
```

## Input file format

Example files are in `ExampleFiles/`.

**Note:** Some files in `ExampleFiles/` are empty to reduce the size of the whole repo. A runnable example dataset can be downloaded at [https://doi.org/10.57760/sciencedb.18595](https://doi.org/10.57760/sciencedb.18595). (The value of `ResolutionInKb` is recommanded to be changed as `500` in practical)

## Citation

Wenxi Wang, Zhe Chen, Zhengzhao Yang, Zihao Wang, Jilu Liu, Jie Liu, Huiru Peng, Zhenqi Su, Zhongfu Ni, Qixin Sun, Weilong Guo. lcQTH: rapid quantitative trait mapping through tracing parental haplotype with ultra-low-coverage sequencing, Plant Communications, 2024. https://doi.org/10.1016/j.xplc.2024.101008.
