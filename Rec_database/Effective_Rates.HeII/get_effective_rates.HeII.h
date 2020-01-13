//====================================================================================================================
// get_effective_rates.HeI.h
//====================================================================================================================

#ifndef GET_EFFECTIVE_RATES_HeII
#define GET_EFFECTIVE_RATES_HeII

#include <string>
#include <vector>

#include "Atom.h"

namespace Effective_rates_HeII
{
    void read_parameter_file(string fname, Gas_of_Atoms &HIA, vector<int> &res_Level_indices);
    
    void load_rates(string path, int nShell_max, Gas_of_Atoms &HIA);
    
    void get_rates(double Tg, vector<double> &Ai, vector<double> &Bi, vector<vector<double> > &RijVec);
    
    void get_rates_all(double Tg, vector<double> &Ai, vector<double> &Bi,
                       vector<double> &Bitot, vector<vector<double> > &RijVec);
    
    void get_Bitot(double Tg, vector<double> &Bitot);
    void get_Bitot(double Tg, double &Bitot, int ires);
    
    int get_res_index(int index_of_HI_state);
}

int get_HeII_index(int index_of_res_state);
int get_number_of_resolved_levels_HeII();

#endif

//====================================================================================================================
//====================================================================================================================
