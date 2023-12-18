#!/usr/bin/env python3

##################################################
# Program Name:         TCFreeze
# Original Author:      Mark A. Hix
# Date Last Updated:    2022 October 11
# Version:              1.2
# Software Required:    Python3,TeraChem
# Contributors/Testers: Ashlyn Murphy,
#                       Solomon Yamoah Effah
##################################################
################# CHANGE LOG #####################
# 2022-08-16: Added "Finalize()" function to complete optimization without restraints. (v1.2)
##################################################
# library imports
import glob
import subprocess
import os
import numpy as np
import argparse
from datetime import datetime

def DateTime():
    now = datetime.now()
    return now.strftime("%Y/%m/%d %H:%M:%S")

def PadHeading(string: str):
    if (len(string) == len(string.strip())) and string != "":
        string = " " + string + " "
    start = (80-len(string))//2
    end = 80 - start - len(string)
    return ''.join(["#" for x in range(start)]) + string + ''.join(["#" for x in range(end)]) + "\n"

def Print_Frozen_To_Input_Log(frozen_indices: list,step_number: int,JobSettings: dict,logfile: str,**kwargs):
    with open(logfile,"a") as f:
        if step_number == 0:
            f.write(PadHeading("Getting Initial Gradient"))
            f.write(PadHeading(f"Started at {DateTime()}"))
            for key,value in kwargs.items():
                f.write(PadHeading(f"{key:<20}: {value:>10}"))
        if step_number <= 1:
            f.write(PadHeading("Beginning TCFreeze Optimizations"))
            f.write(PadHeading(f"Started at {DateTime()}"))
            for key,value in JobSettings.items():
                f.write(f"{key:<30}{value:30}\n")
        
        f.write("\n")
        if step_number>0:
            frozen_indices.sort()
            num_chars = 0
            if len(frozen_indices) > 0:
                num_chars=int(len(str(frozen_indices[-1])))
            num_indices_before_line = 80//(num_chars+1)
            f.write(PadHeading(f"Iteration {step_number}"))
            f.write(PadHeading(f"Started at {DateTime()}"))
            f.write(PadHeading(f"Frozen Atom Indices: {len(frozen_indices)}"))
            for i in range(len(frozen_indices)):
                f.write(f"{frozen_indices[i]:<{num_chars}} ")
                if i%num_indices_before_line == num_indices_before_line-1:
                    f.write("\n")
                elif i+1 == len(frozen_indices):
                    f.write("\n")
        f.write(PadHeading(""))
    return

def ProcessOutputFile(outputfile):
    subcycle= 0
    E_total = 0
    E_rms   = 0
    S_max   = 0
    S_rms   = 0
    G_max   = 0
    G_rms   = 0
    outputstring = ""
    for line in open(outputfile).readlines():
        if "Testing convergence  in cycle" in line:
            subcycle = line.split()[-1]
        elif "Energy calculation finished, energy" in line:
            E_total  = float(line.split()[-1])
        elif all(["Energy" in line,"Target" in line, "converged" in line]):
            E_rms    = float(line.split()[1])
        elif all(["Max step" in line,"Target" in line, "converged" in line]):
            S_max    = float(line.split()[2])
        elif all(["RMS step" in line,"Target" in line, "converged" in line]):
            S_rms    = float(line.split()[2])
        elif all(["Max grad" in line,"Target" in line, "converged" in line]):
            G_max    = float(line.split()[2])
        elif all(["RMS grad" in line,"Target" in line, "converged" in line]):
            G_rms    = float(line.split()[2])
            outputstring += f"Subcycle Step # {subcycle:>20}\n"
            outputstring += f"Total Energy:   {E_total:>20.5f}\n"
            outputstring += f"RMS Energy:     {E_rms:>20.5e}\n"
            outputstring += f"Max Step Size:  {S_max:>20.5e}\n"
            outputstring += f"RMS Step Size:  {S_rms:>20.5e}\n"
            outputstring += f"Max Gradient:   {G_max:>20.5e}\n"
            outputstring += f"RMS Gradient:   {G_rms:>20.5e}\n"
            outputstring += f"####################################\n\n"
            E_total = 0
            E_rms   = 0
            S_max   = 0
            S_rms   = 0
            G_max   = 0
            G_rms   = 0
    return outputstring

# Primary class definition
class System():
    def __init__(self,input_file,threshold,steps_per_cycle):
        self.steps_per_cycle = steps_per_cycle
        if "%" in threshold:
           self.threshold = threshold
        elif str(threshold).upper() != "UNRESTRAINED":
            self.threshold = float(threshold)
        else:
            self.threshold = threshold
        self._parse_job_settings(input_file)
        os.chdir(os.path.dirname(os.path.abspath(input_file)))
        self.coords  = self.job_settings["coordinates"]
        self.regions = self.job_settings["qmindices"]
        self.prmtop  = self.job_settings["prmtop"]
        self.work_dir= f"{self.coords.split('.')[0]}_{str(self.threshold).replace('%','pct')}thresh_{self.steps_per_cycle}spc"
        self.maininpfile=os.path.join(os.path.dirname(os.path.abspath(input_file)),"maininpfile.txt")
        self.mainlogfile=os.path.join(os.path.dirname(os.path.abspath(input_file)),"mainlogfile.txt")
        self.mainerrfile=os.path.join(os.path.dirname(os.path.abspath(input_file)),"mainerrfile.txt")
        if steps_per_cycle == 1:
            self.steps_per_cycle = 2
        self.converged = False
        self.frozen_atoms=[]
        self.iterations = len(glob.glob(self.work_dir+"/optim_*.dcd"))
        self._get_natoms_prmtop()
        self._get_atom_indices()
        self.gradient=np.zeros(int(self.natoms),dtype=float)
        return

    def _parse_job_settings(self,input_file):
        self.job_settings={}
        for line in open(input_file,"r").readlines():
            strline = line.strip()
            if not strline:
                continue
            if strline[0] == "$":
                break
            if strline[0] == "#":
                continue
            key=strline.split()[0]
            value=strline.split()[1]
            self.job_settings[key]=value
        if self.job_settings["run"] != "conical":
            self.job_settings["run"] = "minimize"
        self.job_settings["new_minimizer"] = "no"
        self.job_settings["min_coordinates"] = "cartesian"
        self.job_settings["nstep"] = self.steps_per_cycle
        self.job_settings["scrdir"] = "./"
        self.job_settings["mincheck"] = "false"
        return
    
    def _get_natoms_prmtop(self):
        with open(self.prmtop,"r") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if "%FORMAT(10I8)" in lines[i]:
                    self.natoms = int(lines[i+1].split()[0])
                    break
        return
    
    def _get_atom_indices(self):
        with open(self.regions,"r") as f:
            qmatoms = [int(x) for x in f.read().split()]
        mmatoms = list(range(self.natoms))
        for atom in qmatoms:
            mmatoms.pop(atom)
        self.indices = list(qmatoms) + list(mmatoms)
        del mmatoms
        del qmatoms
        return

    def _get_gradient(self,output_file):
        self.gradient=np.zeros(int(self.natoms),dtype=float)
        with open(output_file,"r") as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if "Gradient units are Hartree/Bohr" in lines[i]:
                last_grad_start = i+3
        newlines = lines[last_grad_start:]
        i = 0
        for line in newlines:
            if "---------------------------------------------------" in line:
                break
            if not line:
                break
            if (line == "\n") or (line == "------- MM / Point charge part -------\n") or ("0.0000000000     0.0000000000     0.0000000000" in line):
                continue
            [x,y,z] = line.split()
            j = self.indices[i]
            self.gradient[j] = np.sqrt(float(x)**2 + float(y)**2 + float(z)**2)
            i+=1
        return
    
    def _freeze_atoms(self):
        self.frozen_atoms=[]
        if type(self.threshold) == str:
            if str(self.threshold).upper() == "UNRESTRAINED":
                curr_thresh = 0
            elif "%" in self.threshold:
                curr_thresh = float(self.threshold.replace("%",""))*0.01*max(self.gradient)
        elif type(self.threshold) == float:
            curr_thresh=self.threshold
        for i in range(len(self.gradient)):
            if self.gradient[i] < curr_thresh:
                self.frozen_atoms.append(i)
        return
    
    def _write_input_file(self,com_file):
        with open(com_file,"w+") as of:
            for key,value in self.job_settings.items():
                if key == "nstep" and self.steps_per_cycle == 0:
                    self.steps_per_cycle = 10
#                    continue
                of.write(f"{key:<30}{value:<50}\n")
            if self.frozen_atoms != []:
                of.write("$constraints\n")
                of.write("".join([f"atom {x}\n" for x in self.frozen_atoms]))
                of.write("$end\n")
        of.close()
        return

    def Check_If_Complete(self,output_file):
        self.converged = False
        with open(output_file,"r",encoding='utf-8') as f:
            for line in f.readlines():
                if not line:
                    break
                if "Converged!" in line:
                    self.converged=True
                if "NOT CONVERGED" in line:
                    self.converged=False
                if "DIE called" in line:
                    print("TCFREEZE JOB FAILURE.  REVIEW OUTPUTS.")
                    self.converged=True
        return
  
    def Macroiteration(self):
        ### get frozen atoms
        ### write input file
        comfile = f"minim_{self.iterations}.in"
        outfile = f"minim_{self.iterations}_output"
        errfile = f"minim_{self.iterations}_error"
        if self.iterations == 0:
            self.job_settings["run"] = "gradient"
        if self.iterations > 0:
            self._freeze_atoms()
            self.job_settings["run"] = "minimize"
        if self.iterations > 1:
            self.job_settings["coordinates"] = "optim.rst7"
        self._write_input_file(comfile)
        Print_Frozen_To_Input_Log(self.frozen_atoms,self.iterations,self.job_settings,self.maininpfile,Threshold=self.threshold,Subcycles=self.steps_per_cycle)
        ### run input file with terachem
        self.silent_shell(f"terachem {comfile} 1> {outfile} 2> {errfile}")
        #### check for convergence
        if self.iterations > 0:
            self.Check_If_Complete(outfile)
            self.silent_shell(f"mv optim.dcd optim_{self.iterations}.dcd")
        #### Process input, output, and error files to respective main logs.
        with open(self.mainlogfile,"a") as of:
            of.write(PadHeading(f"Macroiteration {self.iterations}"))
            of.write(ProcessOutputFile(outfile))
        with open(errfile,"r") as f:
            errorfile = f.read().strip()
            if errorfile:
                self.silent_shell(f"cat {errfile} >> {self.mainerrfile}")
                self.silent_shell(f"echo '\n\n#####################################\n\n' >> {self.mainerrfile}")
        #### update gradient
        self._get_gradient(outfile)
        self.iterations +=1

    def Optimize(self,max_steps=0):
        self.silent_shell(f"mkdir -p {self.work_dir}")
        self.silent_shell(f"cp * {self.work_dir}")
        os.chdir(self.work_dir)
        while not self.converged:
            self.Macroiteration()
            if max_steps!=0:
                if self.iterations >= max_steps:
                    print("Maximum cycles reached.  Optimization Terminating.")
                    self.converged=True
        os.chdir("../")
        return
    def Finalize(self):
        self.threshold = "unrestrained"
        self.converged = False
        os.chdir(self.work_dir)
        while not self.converged:
            self.Macroiteration()
        os.chdir("../")
        return
    def silent_shell(self,com_line):
        devnull = open(os.devnull,"w")
        subprocess.call(com_line,shell=True,stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        devnull.close()
        return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input-file",dest="infile",help="/path/to/Terachem/input",required=True)
    parser.add_argument("-t","--threshold",dest="thresh",help="Threshold below which atomic coordinates are frozen for optimization.  Terachem has a default MaxGrad convergence value of 4.5E-4.  May also be given as a percentage (20%%, 10%%, etc.).  Using a percentage will freeze atomic coordinates if the gradient at that atom is lower than the given percentage of the CURRENT MaxGrad.",required=True)
    parser.add_argument("-n","--nsteps-per-opt",dest="nsteps",help="Number of optimization steps before recalculating which atoms to freeze.",required=True)
    parser.add_argument("-m","--max-opt-steps",dest="maxsteps",help="Maximum total number of optimization/freezing cycles (The maximum number of times the system will recalculate frozen atoms).",default=0)
    args = parser.parse_args()
    job = System(args.infile,args.thresh,int(args.nsteps))
    job.Optimize(int(args.maxsteps))
    job.Finalize()
