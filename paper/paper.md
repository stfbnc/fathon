---
​---
title: 'fathon: A Python package for a fast computation of detrendend fluctuation analysis and related algorithms'
tags:
  - Python
  - Cython
  - Detrended fluctuation analysis
  - Multifractal spectrum
authors:
  - name: Stefano Bianchi
  - orcid: 0000-0001-6729-8655
  - affiliation: 1
affiliations:
  - Department of Mathematics and Physics, Roma Tre University
  - index: 1
date: 13 October 2019
bibliography: paper.bib
​---
---

# Summary

Detrended fluctuation analysis (DFA) was first developed by Peng to study memory effects in sequences of DNA [@peng_1993]. Since then, it found various applications in other fields, like geophysics or economy [@varotsos_2007; @koscielny_1998; @talkner_2000; @fraedrich_2003; @matsoukas_2000; @kavasseri_2004; @liu_1999; @janosi_1999; @grau_2000; @ivanov_2004]. DFA is useful in recognising long memory processes by means of the scaling properties of the fluctuation function $F (n) \sim n^H$ . $H$ is the Hurst exponents, and its value can tell if a process (most likely a time series) is persistent or anti-persistent:

- $H \in [0.0, 0.5) \rightarrow$ antipersistency;
- $H = 0.5 \rightarrow$ uncorrelated process;
- $H \in (0.5, 1.0] \rightarrow$ persistency;
- $H > 1.0$ nonstationary process, stronger long-range correlations are present.

 The algorithm has later evolved to multifractal detrended fluctuation analysis (MFDFA) [@kantelhardt_2002] to take into account the variability of time series and the possible multiple scaling properties of the fluctuation function. Other algorithms have been further developed, such as the time dependent Hurst exponent [@ihlen_2012] that examines changes of persistency within time intervals, and the detrended cross-correlation analysis (DCCA) [@podobnik_2008] that studies cross-correlations in terms of persistency between non-linear time series. Recently, a cross-correlation index $\rho$ has been introduced in order to make the interpretation of DCCA clearer [@zebende_2011].



``fathon`` is a Python package for DFA (Detrended Fluctuation Analysis) and related algorithms. It aims at providing a simple and interactive interface to perform multiple analyses on multiple time series. Its purpose is to gather all the algorithms together, in order to have a single package including different methodologies. They are optimised, easy to use in their basic or advanced versions, and the user can easily switch between them. The Python interface is class-based and easy to use, while the underlying code is written in Cython and C allowing a fast computation of the algorithms. ``fathon`` implements common operations for these kind of algorithms, such as the fluctuation function and its fit in different ranges, the multifractal spectrum, and the detrended cross-correlation index and its confidence intervals.

# References

