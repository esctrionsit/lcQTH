# lcQTH

Rapid quantitative trait mapping through tracing parental haplotype with ultra-low-coverage sequencing

## Recommanded environments with version tested:

+ Linux system with a 64 bit CPU
+ Bcftools v1.11
+ Bedtools v2.29.2
+ python v3.7.6
+ python packages:
  + hmmlearn v0.2.4
  + openpyxl
+ R v3.6
+ R packages:
  + ComplexHeatmap
  + qtl
+ Golang 1.20.3 linux/amd64
+ gcc version 11.2.1 20220127 (Red Hat 11.2.1-9) (GCC)
+ Node v14.17.5
+ Node packages:
  + d3@4.13.0
  + jsdom@15.2.1
  + jspdf
  + svg2pdf.js

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

Example files are in `ExampleFiles`.

## Citation

Wenxi Wang, Zhe Chen, Zhengzhao Yang, Zihao Wang, Jilu Liu, Jie Liu, Huiru Peng, Zhenqi Su, Zhongfu Ni, Qixin Sun, Weilong Guo. lcQTH: rapid quantitative trait mapping through tracing parental haplotype with ultra-low-coverage sequencing, Plant Communications, 2024. https://doi.org/10.1016/j.xplc.2024.101008.
