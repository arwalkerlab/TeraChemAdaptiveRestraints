#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <experimental/filesystem>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <locale>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <iomanip>
#include "tcfreeze.h"

bool DEBUG_VERBOSE = false;

namespace fs = std::experimental::filesystem;

void Print_Frozen_To_Input_Log(std::vector<int> frozen_indices, int step_number, std::map<std::string,std::string> JobSettings, std::string currfile, std::string logfile)
{   
    std::ofstream outFile;
    
    outFile.open(logfile,std::ios_base::app);
    if (step_number == 0)
    {   
        std::string dummy;
        dummy = "Getting Initial Gradient";
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
        dummy = "Started at ";
        dummy += utils::DateTime();
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
    }
    if (step_number <= 1)
    {
        std::string dummy;
        dummy = "Beginning TCFreeze Optimizations";
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
        dummy = "Started at ";
        dummy += utils::DateTime();
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
        utils::catFiles(currfile,logfile);

    }
    outFile << std::endl;
    if (step_number > 1)
    {
        std::string dummy;
        int num_chars = 0;
        int n_frozen = frozen_indices.size();
        if (n_frozen > 0)
            num_chars = std::to_string(frozen_indices[n_frozen-1]).length();
        int num_indices_before_line = 80 / ( num_chars + 1 );
        utils::Debug("Number of frozen indices per line: ");
        utils::Debug(std::to_string(num_indices_before_line));
        dummy = " Iteration " + std::to_string(step_number);
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
        dummy = "Started at ";
        dummy += utils::DateTime();
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
        dummy = "Frozen Atom Indices: " + std::to_string(n_frozen);
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
        for (int i=0; i < n_frozen; i++)
        {
            std::string tmp_idx = std::to_string(frozen_indices[i]);
            for (int j=0; j < num_chars - tmp_idx.length()+1; j++)
            {
                outFile << " ";
            }
            outFile << tmp_idx;
            if (i % num_indices_before_line == num_indices_before_line - 1)
                outFile << std::endl;
        }
        outFile << std::endl;
        dummy = "";
        utils::PadHeading(dummy);
        outFile << dummy << std::endl;
    }
    outFile.close();
}

std::string ProcessOutputFile(std::string outputfile)
{
    utils::Debug("Processing Output File:");
    utils::Debug(outputfile);
    int subcycle   = 0;
    double E_total = 0;
    double E_rms   = 0;
    double S_max   = 0;
    double S_rms   = 0;
    double G_max   = 0;
    double G_rms   = 0;
    std::string outputstring = "";
    std::string line;
    std::ifstream inFile;
    std::string dummy;
    inFile.open(outputfile);
    utils::Debug("Opened Output File for reading.");
    while (getline(inFile,line))
    {
        if (utils::isSubstring(line,"Testing convergence  in cycle"))
        {
            utils::Debug("Found cycle number.");
            subcycle = stoi(line.substr(line.find_last_of(" ")+1));
        }
        else if (utils::isSubstring(line,"Energy calculation finished, energy"))
        {
            utils::Debug("Found E_total line.");
            E_total = stof(line.substr(line.find_last_of(" ")+1));
            utils::Debug("Found E_total");
        }
        else if (utils::isSubstring(line,"Energy") && utils::isSubstring(line,"Target") && utils::isSubstring(line,"converged"))
        {
            utils::Debug("Found E_rms line.");
            dummy = utils::trim(line.substr(10,22));
            E_rms = stof(dummy);
            utils::Debug("Found E_rms");
        } 
        else if (utils::isSubstring(line,"Max step") && utils::isSubstring(line,"Target") && utils::isSubstring(line,"converged"))
        {
            utils::Debug("Found S_max line.");
            dummy = utils::trim(line.substr(10,22));
            S_max = stof(dummy);
            utils::Debug("Found S_max");
        }
        else if (utils::isSubstring(line,"RMS step") && utils::isSubstring(line,"Target") && utils::isSubstring(line,"converged"))
        {
            utils::Debug("Found S_rms line.");
            dummy = utils::trim(line.substr(10,22));
            S_rms = stof(dummy);
            utils::Debug("Found S_rms");
        }
        else if (utils::isSubstring(line,"Max grad") && utils::isSubstring(line,"Target") && utils::isSubstring(line,"converged"))
        {
            utils::Debug("Found G_max line");
            dummy = utils::trim(line.substr(10,22));
            G_max = stof(dummy);
            utils::Debug("Found G_max");
        }
        else if (utils::isSubstring(line,"RMS grad") && utils::isSubstring(line,"Target") && utils::isSubstring(line,"converged"))
        {
            utils::Debug("Found G_rms line.");
            dummy = utils::trim(line.substr(10,22));
            G_rms = stof(dummy);
            utils::Debug("Found G_rms");
            std::stringstream buffer;

            buffer << "Subcycle Step # ";
            buffer.precision(0);
            buffer << std::setfill(' ') << std::setw(20) << subcycle;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            buffer << "Total Energy:   ";
            buffer.precision(5);
            buffer << std::setfill(' ') << std::setw(20) << E_total;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            buffer << "RMS Energy:     ";
            buffer.precision(5);
            buffer << std::setfill(' ') << std::setw(20) << E_rms;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            buffer << "Max Step Size:  ";
            buffer.precision(5);
            buffer << std::setfill(' ') << std::setw(20) << S_max;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            buffer << "RMS Step Size:  ";
            buffer.precision(5);
            buffer << std::setfill(' ') << std::setw(20) << S_rms;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            buffer << "Max Gradient:   ";
            buffer.precision(5);
            buffer << std::setfill(' ') << std::setw(20) << G_max;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            buffer << "RMS Gradient:   ";
            buffer.precision(5);
            buffer << std::setfill(' ') << std::setw(20) << G_rms;
            buffer << std::endl;
            outputstring += buffer.str();
            buffer.str("");

            outputstring += "####################################\n\n";
            utils::Debug("Resetting step result values.");
            E_total = 0;
            E_rms   = 0;
            S_max   = 0;
            S_rms   = 0;
            G_max   = 0;
            G_rms   = 0;
        }
    }
    inFile.close();
    utils::Debug("Closed Output File.");
    return outputstring;
}

namespace fs = std::experimental::filesystem;

namespace utils
{
    bool isSubstring(std::string s1, std::string s2)
    {
        return (s1.find(s2) != std::string::npos);
    }
    
    std::string trim(std::string s)
    {
        int start = s.find_first_not_of(' ');
        int end   = s.find_last_not_of(' ');
        if (start == -1 )
        {
            return {};
        }
        return s.substr(start,end-start+1);
    }

    bool isComment(std::string s)
    {
        return (trim(s)[0] == '#');
    }

    bool isEmpty(std::string s)
    {
        return (trim(s).empty());
    }

    void silent_shell(const char* cmd)
    {
        std::array<char, 128> buffer;
        std::string result;
        std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
        if (!pipe) throw std::runtime_error("popen() failed!");
        while (!feof(pipe.get())) {
            if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
                result += buffer.data();
        }
    }

    bool isInArray(std::vector<int> array, int value)
    {
        for (int i=0; i < array.size(); i++)
        {
            if (array[i] == value)
                return true;
        }
        return false;
    }

    std::string DateTime()
    {
        auto t = std::time(nullptr);
        auto tm = *std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(&tm, "%Y/%m/%d %H:%M:%S");
        auto str = oss.str();
        return str;
    }

    void PadHeading(std::string& line)
    {
        line = trim(line);
        if (isEmpty(line))
            line = "################################################################################";
        else
        {
            line = " " + line + " ";
            int str_len = line.size();
            int start = (80 - str_len)/2;
            int end = 80 - start - str_len;
            std::string newline = "";
            for (int i=0; i<start; i++)
            {
                newline += "#";
            }
            newline += line;
            for (int i=0; i<end; i++)
            {
                newline += "#";
            }
            line = newline;
        }
    }

    void catFiles(std::string fromfile,std::string tofile)
    {
        std::ifstream inFile;
        std::ofstream outFile;
        std::string line;
        inFile.open(fromfile,std::ios::in);
        outFile.open(tofile,std::ios::app);
        while (getline(inFile,line))
        {
            outFile << line << std::endl;
        }
        outFile.close();
        inFile.close();
    }
 
    bool IsFlag(char* bigstring) //checks if a string is a command line flag or a value.
    {
        if (bigstring[0] == '-')
            return true;
        return false;
    }

    void ReadArgs(int argc,
                char** argv,
                std::vector<std::vector<std::string>>& flags)
    {
        int j = 0;
        for (int i = 1; i < argc; i++)
        {
            if (IsFlag(argv[i]))
            {
                j++;
                flags.push_back({argv[i]});
            }
            else
                flags[j].push_back(argv[i]);

        }
    }

    int FindFlag(std::vector<std::vector<std::string>>& flags, char* target)
    {
        for (int i = 1; i<flags.size(); i++)
        {
            if (flags[i][0] == target)
                return i;
        }
        return 0;
    }

    void ReportFailure(std::string failure)
    {
        std::cout << failure << "  Exiting." << std::endl;
        // close all open files?

        // exit program with error.
        exit(EXIT_FAILURE);
    }

    void Debug(std::string debug)
    {
        if (DEBUG_VERBOSE)
            std::cout << debug << std::endl;
    }

    void CheckForExternalPrograms(std::vector<char*> programs)
    {
        bool allprogsavail;
        allprogsavail=true;
        for (int i=0; i<programs.size(); i++)
        {
            bool progavail;
            progavail = CheckProgAvailable(programs[i]);
            if (!progavail)
                allprogsavail=false;
        }
        if (!allprogsavail)
        {
            ReportFailure("Required programs were missing.");
        }
    }
 
    bool CheckProgAvailable(char* program)
    {
            std::string result;
            std::string cmd;
            cmd = "which ";
            cmd += program;
            result=GetSysResponse(cmd.c_str());
            if (result.empty())
            {
                std::cout << "Missing program: " << program << std::endl;
                return false;
            }
            return true;
    }

    std::string GetSysResponse(const char* cmd) 
    {
        std::array<char, 128> buffer;
        std::string result;
        std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
        if (!pipe) throw std::runtime_error("popen() failed!");
        while (!feof(pipe.get())) {
            if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
                result += buffer.data();
        }
        return result;
    }

}

void System::_parse_job_settings()
{
    std::ifstream inFile;
    std::string line;
    inFile.open(SubmittedInputFile,std::ios::in);
    while (getline(inFile,line))
    {

        if (utils::isComment(line) || utils::isEmpty(line))
           continue;
        std::stringstream dummy;
        std::string key;
        std::string value;
        dummy << line;
        dummy >> key;
        dummy >> value;
        if (utils::isEmpty(key)||utils::isEmpty(value))
            continue;
        JobSettings[key] = value.c_str();
    }
    JobSettings["run"] = "minimize";
    JobSettings["new_minimizer"] = "no";
    JobSettings["min_coordinates"] = "cartesian";
    JobSettings["nstep"] = std::to_string(steps_per_cycle);
    JobSettings["scrdir"] = "./";
    JobSettings["mincheck"] = "false";
    coords  = JobSettings["coordinates"];
    regions = JobSettings["qmindices"];
    prmtop  = JobSettings["prmtop"];
    inFile.close();
}

void System::_get_natoms_prmtop()
{
    std::ifstream inFile;
    inFile.open(prmtop,std::ios::in);
    std::string line;
    while(getline(inFile,line))
    {
        if (utils::isSubstring(line,"FORMAT(10I8)"))
        {
            getline(inFile,line);
            std::string dummy;
            dummy = line.substr(line.find_first_not_of(' '));
            dummy = dummy.substr(0,dummy.find_first_of(' '));
            natoms = stoi(dummy);
            break;
        }
    }
}

void System::_get_atom_indices()
{
    std::vector<int> qm_indices;
    std::vector<int> mm_indices;
    std::ifstream inFile;
    std::string line;
    inFile.open(regions,std::ios::in);
    while (getline(inFile,line)) 
    {
        std::stringstream dummy;
        dummy << line;
        int x;
        while (dummy >> x)
        {
            qm_indices.push_back(x);
        }
        if( inFile.eof() ) 
            break;    
    }
    inFile.close();
    for (int i=0; i < natoms; i++)
    {
        if (!utils::isInArray(qm_indices, i))
            mm_indices.push_back(i);
    }
    for (int i=0;i<qm_indices.size(); i++)
    {
        indices.push_back(qm_indices[i]);
    }
    for (int i=0;i<mm_indices.size(); i++)
    {
        indices.push_back(mm_indices[i]);
    }
    utils::Debug("QM Atoms:");
    utils::Debug(std::to_string(qm_indices.size()));
    utils::Debug("MM Atoms:");
    utils::Debug(std::to_string(mm_indices.size()));
}

void System::_get_gradient(std::string outputfile)
{
    gradient = std::vector<double> (natoms,0.0);
    utils::Debug("Initialized gradient.");
    std::ifstream inFile;
    std::string line;
    int i = 0;
    utils::Debug("Initialized several variables inside _get_gradient()");
    inFile.open(outputfile,std::ios::in);
    utils::Debug("Opened the logfile from the previous step.  Attempting to read in gradient information.");
    int grads_encountered = 0;
    while (getline(inFile,line))
    {
        if (utils::isSubstring(line,"Gradient units are Hartree/Bohr"))
        {
            getline(inFile,line);
            getline(inFile,line);
            grads_encountered++;
        }
        if (grads_encountered == stoi(JobSettings["nstep"]))
        {
            break;
        }
    }
    while (getline(inFile,line))
    {
        if (utils::isSubstring(line,"---------------------------------------------------"))
        {
            break;
        }
        if (utils::isEmpty(line) || utils::isSubstring(line,"------- MM / Point charge part -------\n") || utils::isSubstring(line,"0.0000000000     0.0000000000     0.0000000000"))
        {
            continue;
        }
        std::stringstream dummy;
        dummy << line;
        double x;
        double y;
        double z;
        dummy >> x >> y >> z;
        int j = indices[i];
        double gradnorm = pow(x,2) + pow(y,2) + pow(z,2);
        gradient[j] = gradnorm;
        i++;
    }
    inFile.close();
    max_grad = *std::max_element(gradient.begin(), gradient.end());
    utils::Debug("Obtained the maximum gradient value: ");
    utils::Debug(std::to_string(max_grad));
}

void System::_freeze_atoms()
{
    frozen_atoms={};
    if (no_threshold)
        current_threshold = 0.0;
    else if (dynamic_threshold)
        current_threshold = max_grad * threshold;
    else if (static_threshold)
        current_threshold = threshold;
    for (int i=0;i<gradient.size(); i++)
    {
        if (gradient[i] < current_threshold*current_threshold)
        {
            frozen_atoms.push_back(i);
        }
    }
}
    
void System::_write_input_file(std::string inputfilename)
{
    std::ofstream outFile;
    outFile.open(inputfilename,std::ios::out);
    for(std::map<std::string,std::string>::iterator iter = JobSettings.begin(); iter != JobSettings.end(); ++iter)
    {
        //iterate through JobSettings keys and values.
        std::string key =  iter->first;
        std::string value = iter->second;
        if (key == "nstep" && steps_per_cycle == 0)
        {
            continue;
        }
        std::stringstream buffer;
        buffer.str();
        buffer << key;
        for (int i=key.size(); i < 30; i++)
        {
            buffer << " ";
        }
        buffer << value;
        for (int i=value.size(); i < 50; i++)
        {
            buffer << " ";
        }
        buffer << std::endl;
        outFile << buffer.str();
    }   
    if (frozen_atoms.size() > 0)
    {
        outFile << "$constraints" << std::endl;
        for (int i = 0; i < frozen_atoms.size(); i++)
        {
            outFile << frozen_atoms[i] << std::endl;
        }
        outFile << "$end" << std::endl;
    }
    outFile.close();

}

void System::Check_If_Complete(std::string outputfile)
{
    std::ifstream inFile;
    inFile.open(outputfile,std::ios::in);
    std::string line;
    while (getline(inFile,line))
    {
        if (utils::isSubstring(line,"Converged!"))
        {
            converged=true;
        }
        if (utils::isSubstring(line,"NOT CONVERGED"))
        {
            converged=false;
        }
        if (utils::isSubstring(line,"DIE called"))
        {
            std::cout << "TCFREEZE JOB FAILURE.  REVIEW OUTPUTS." << std::endl;
            converged=true;
        }
    }
    inFile.close();
}

void System::Macroiteration()
{
    std::stringstream buffer;
    buffer.str("");
    buffer << "minim_";
    buffer << std::setfill('0') << std::setw(4) << iterations;
    buffer << ".in";
    curr_inp_file = buffer.str();

    buffer.str("");
    buffer << "minim_";
    buffer << std::setfill('0') << std::setw(4) << iterations;
    buffer << ".out";
    curr_out_file = buffer.str();

    buffer.str("");
    buffer << "minim_";
    buffer << std::setfill('0') << std::setw(4) << iterations;
    buffer << ".err";
    curr_err_file = buffer.str();

    utils::Debug("Current Input, Output, and Error Files:");
    utils::Debug(curr_inp_file);
    utils::Debug(curr_out_file);
    utils::Debug(curr_err_file);

    _freeze_atoms();
    utils::Debug("Atoms Frozen.");

    if (iterations == 0)
    {
        JobSettings["run"] = "gradient";   
    }
    if (iterations > 0)
    {
        _freeze_atoms();
        JobSettings["run"] = "minimize";   
    }
    if (iterations > 1)
    {
        JobSettings["coordinates"] = "optim.rst7";           
    }
    utils::Debug("Accounted for iteration number in JobSettings.");

    _write_input_file(curr_inp_file);
    utils::Debug("Input File written.");

    Print_Frozen_To_Input_Log(frozen_atoms,iterations,JobSettings,curr_inp_file,maininpfile);
    utils::Debug("Printed frozen atoms to main input log.");
    
    buffer.str("");
    buffer << "terachem " << curr_inp_file;
    buffer << " 1> " << curr_out_file;
    buffer << " 2> " << curr_err_file;

    utils::Debug("Calling terachem with:");
    utils::Debug(buffer.str().c_str());
    utils::silent_shell(buffer.str().c_str());
    
    utils::Debug("Terachem has completed its run.");
    if (iterations > 0)
    {
        utils::Debug("Checking for convergence.");
        Check_If_Complete(curr_out_file);
        utils::Debug("Moving optimization trajectory to unique-id file.");
        buffer.str("");
        buffer << "mv optim.dcd optim_";
        buffer << std::setfill('0') << std::setw(4) << iterations;
        buffer << ".dcd";
        utils::silent_shell(buffer.str().c_str());
    }
    utils::Debug("Processing output file into main output log.");
    std::ofstream outFile;
    outFile.open(mainlogfile,std::ios::app);
    buffer.str("");
    buffer << "Macroiteration ";
    buffer << std::setfill('0') << std::setw(4) << iterations;
    std::string dummy = buffer.str();
    utils::PadHeading(dummy);
    outFile << dummy << std::endl;
    utils::Debug("Processing current output file.");
    dummy = ProcessOutputFile(curr_out_file);
    outFile << dummy << std::endl;
    outFile.close();
    utils::Debug("Processing current error file.");
    utils::catFiles(curr_err_file,mainerrfile);
    utils::Debug("Output files processed.  Obtaining current gradient.");
    _get_gradient(curr_out_file);
    utils::Debug("gradient obtained.  Updating iteration counter.");
    iterations++;
}

System::System(std::string input_file,std::string thresh, int spc, int msteps)
{
    SubmittedDirectory = fs::absolute(fs::current_path());
    WorkingDirectory = SubmittedDirectory.string() + "/TCFreeze_scr/";
    utils::Debug("Identified Submitted and Working directories:");
    utils::Debug(SubmittedDirectory.string());
    utils::Debug(WorkingDirectory.string());
    fs::create_directory(WorkingDirectory);
    SubmittedInputFile = fs::absolute(input_file);
    utils::Debug("Found SubmittedInputFile:");
    utils::Debug(SubmittedInputFile.string());
    steps_per_cycle = spc;
    max_steps = msteps;
    // Boolean Flags
    converged = false;
    static_threshold = false;
    dynamic_threshold = false;
    no_threshold = false;
    utils::Debug("Initialized all internal starting boolean flags.");
    if (utils::isSubstring(thresh,"%"))
    {
        std::string dummy;
        dynamic_threshold = true;
        dummy = thresh.substr(0,thresh.find_first_of('%'));
        threshold = stof(dummy);
    }
    else if (utils::isSubstring(thresh,"unrestrained"))
    {
        no_threshold=true;
        threshold = 0.0;
    }
    else
    {
        static_threshold=true;
        threshold = stof(thresh);
    }
    utils::Debug("Processed provided threshold to set type (static, dynamic, none).");
    iterations = 0;
    _parse_job_settings();
    utils::Debug("Parsed Job Settings.");
    fs::current_path(fs::absolute(SubmittedInputFile).parent_path());
    utils::Debug("Moved into directory with SubmittedInputFile.");
    maininpfile = fs::path (SubmittedDirectory.string() + "/maininpfile.txt");
    mainerrfile = fs::path (SubmittedDirectory.string() + "/mainerrfile.txt");
    mainlogfile = fs::path (SubmittedDirectory.string() + "/mainlogfile.txt");
    utils::Debug("Main Logfiles:");
    utils::Debug(maininpfile);
    utils::Debug(mainlogfile);
    utils::Debug(mainerrfile);
    std::ofstream outFile;
    outFile.open(mainlogfile);
    outFile << "Initializing TCFreeze Protocols." << std::endl;
    outFile.close();
    if (steps_per_cycle == 1)
    {
        steps_per_cycle = 2;
    }
    frozen_atoms = {};
    utils::Debug("Initialized steps-per-cycle and frozen_atoms vectors.");
    _get_natoms_prmtop();
    utils::Debug("Got natoms from prmtop. ");
    _get_atom_indices();
    utils::Debug("Got QM and MM atom indices.");
    gradient = std::vector<double> (0.0,natoms);
    utils::Debug("Initialized the gradient to zeros.");
}

System::~System()
{
}

void System::Optimize()
{
    char buffer[200];
    std::stringstream cmdline;
    cmdline.str("");
    cmdline << "mkdir -p " << fs::absolute(WorkingDirectory).string();
    utils::silent_shell(cmdline.str().c_str());

    cmdline.str("");
    cmdline << "cp " << prmtop << " " << fs::absolute(WorkingDirectory).string();
    utils::silent_shell(cmdline.str().c_str());

    cmdline.str("");
    cmdline << "cp " << coords << " " << fs::absolute(WorkingDirectory).string();
    utils::silent_shell(cmdline.str().c_str());

    cmdline.str("");
    cmdline << "cp " << regions << " " << fs::absolute(WorkingDirectory).string();
    utils::silent_shell(cmdline.str().c_str());

    fs::path workpath {fs::absolute(WorkingDirectory)};
    fs::path startdir {fs::absolute(fs::current_path())};
    fs::current_path(workpath);
    while (!converged)
    {
        std::cout << "Moving into Macroiteration()" << std::endl;
        Macroiteration();
        if (max_steps != 0)
        {
            if (iterations >= max_steps)
            {
                std::cout << "Maximum cycles reached.  Optimization Terminating." << std::endl;
                converged=true;
            }
        }
    }
    fs::current_path(startdir);
}

void System::Finalize()
{
    threshold = 0.0;
    no_threshold = true;
    static_threshold = false;
    dynamic_threshold = false;
    converged = false;
    fs::path workpath {fs::absolute(WorkingDirectory)};
    fs::path startdir {fs::absolute(fs::current_path())};
    fs::current_path(workpath);
    JobSettings["nstep"] = "10000";
    while (!converged)
    {
        Macroiteration();
    }
    fs::current_path(startdir);
}

int main(int argc, char **argv)
{
    int max_steps = 0;
    std::vector<char*> programs = {"terachem"};
    utils::CheckForExternalPrograms(programs);
    
    std::vector<std::vector<std::string>> flags = {{}};
    utils::ReadArgs(argc,argv,flags);
    
    int foundindex;
    foundindex = utils::FindFlag(flags,"--debug");
    if (foundindex)
    {
        DEBUG_VERBOSE=true;
    }

    foundindex = utils::FindFlag(flags,"-i");
    if (foundindex == -1)
        utils::ReportFailure("Missing Input File.");
    std::string inputfilename = flags[foundindex][1];
    utils::Debug("Found inputfilename.");

    foundindex = utils::FindFlag(flags,"-t");
    if (foundindex == -1)
        utils::ReportFailure("Missing Threshold Value.");
    std::string threshold = flags[foundindex][1];
    utils::Debug("Found threshold value.");

    foundindex = utils::FindFlag(flags,"-n");
    if (foundindex == -1)
        utils::ReportFailure("Missing Steps/Cycle Value.");
    int steps_per_cycle = stoi(flags[foundindex][1]);
    utils::Debug("Found steps/cycle value.");

    foundindex = utils::FindFlag(flags,"-m");
    utils::Debug("foundindex for -m");
    std::cout << foundindex << std::endl;
    if (foundindex != 0)
        max_steps = stoi(flags[foundindex][1]);
    utils::Debug("Found maxsteps.");
    System job (inputfilename,threshold,steps_per_cycle,max_steps);
    utils::Debug("Completed Constructor, let's see how she flies!");
    
    job.Optimize();
    utils::Debug("Finished restrained optimization.");
    job.Finalize();
    utils::Debug("Finished final optimization without restraints.");
    return 0;
}
