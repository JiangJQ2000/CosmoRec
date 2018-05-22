//===========================================================================================================
void print_message(string mess) // for messages to user
{
    if(parameters.verbosity>=0)
    {
        cout << "\n " << setfill('=') << setw(100) << "=" << endl;
        cout << " || " << mess << endl;
        cout << " " << setfill('=') << setw(100) << "=" << endl << endl;
    }
    
    return;
}

void throw_error(string funcname, string message)
{
    cerr << " " << funcname << " : " << message << endl;
    exit(1);
}

//===========================================================================================================
void set_detailed_name(const Parameters_form_initfile &inp, string &dname);  // defined in output_methods.cpp

//===========================================================================================================
//===========================================================================================================
