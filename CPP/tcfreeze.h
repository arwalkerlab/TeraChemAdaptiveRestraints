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

namespace fs = std::experimental::filesystem;
void Print_Frozen_To_Input_Log(std::vector<int> frozen_indices, int step_number, std::map<std::string,std::string> JobSettings, std::string currfile, std::string logfile);
std::string ProcessOutputFile(std::string outputfile);
namespace fs = std::experimental::filesystem;
namespace utils
{
    bool isSubstring(std::string s1, std::string s2);
    
    std::string trim(std::string s);

    bool isComment(std::string s);

    bool isEmpty(std::string s);

    void silent_shell(const char* cmd);

    bool isInArray(std::vector<int> array, int value);

    std::string DateTime();

    void PadHeading(std::string& line);

    void catFiles(std::string fromfile,std::string tofile);

    void ReadArgs(int argc,char** argv,std::vector<std::vector<std::string>>& flags);

    int FindFlag(std::vector<std::vector<std::string>>& flags, char* target);

    void ReportFailure(std::string failure);
    
    void Debug(std::string debug);

    void CheckForExternalPrograms(std::vector<char*> programs);

    bool CheckProgAvailable(char* program);

    std::string GetSysResponse(const char* cmd);
}
class System
{
    private:
        //Internal Variables
        fs::path SubmittedDirectory;
        fs::path WorkingDirectory;
        fs::path SubmittedInputFile; //Obtained from constructor arguments
        double threshold;            //Obtained from constructor arguments
        int steps_per_cycle;         //Obtained from constructor arguments
        std::map<std::string,std::string> JobSettings; //mapped values obtained in _parse_job_settings()
        double current_threshold;
        fs::path maininpfile;
        fs::path mainlogfile;
        fs::path mainerrfile;
        std::string curr_inp_file;
        std::string curr_out_file;
        std::string curr_err_file;
        int iterations;
        int max_steps;

        // Chemical System Information
        int natoms;          //Obtained in _parse_job_settings()
        std::string coords;  //Obtained in _parse_job_settings()
        std::string regions; //Obtained in _parse_job_settings()
        std::string prmtop;  //Obtained in _parse_job_settings()
        std::vector<int> frozen_atoms;
        std::vector<double> gradient;
        std::vector<int> indices;
        double max_grad;

        // Boolean Flags
        bool converged;
        bool static_threshold;
        bool dynamic_threshold;
        bool no_threshold;

        //Internal Functions
        void _parse_job_settings();
        void _get_natoms_prmtop();
        void _get_atom_indices();
        void _get_gradient(std::string outputfile);
        void _freeze_atoms();
        void _write_input_file(std::string);
        void Check_If_Complete(std::string);
        void Macroiteration();

    public:
        System(std::string,std::string,int,int);
        ~System();
        void Optimize();
        void Finalize();
};
