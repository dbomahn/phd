#!/bin/bash
#PBS -N momip
#PBS -l nodes=1:ppn=1 
#PBS -l walltime=20:00:00
module() { eval `/usr/bin/modulecmd bash $*`; }
module add julia
module add cplex
julia /home/k2g00/k2g3475/multiobjective/solvers/ep+FP/fpep.jl $data $presol $pf
