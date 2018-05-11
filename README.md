# SWOTAprimeCalcs

This repository contains code to compute fits on height-width data, using an Errors in Variables model. Description in "ProposalForAreaCalculations.pdf".

* TestAprimeCalcsPepsiLive.mlx: Live Script that does Aprime calculations on Pepsi Rivers
* TestAprimeCalcsPepsi.m: Script that does Aprime calculations on Pepsi Rivers (built from TestAprimeCalcsPepsiLive)
* CalcWHFitsEIV.m: Function that interfaces and bookkeeps for the main optimization call to fit data (WH_EIV.m)
* WH_EIV.m: Function that does the fits on the width-height data as a constrained non-linear optimization problem. 
* CalculatedAEIV.m: Function that computes Aprime given the data and computed fits on H-W.
