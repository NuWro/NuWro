

#include "dirs.h"
#include <cstring>
#include <cstdlib>
#include <fstream>


static string data_dir;
static string bin_dir;

string get_data_dir ()
{
	return data_dir;
}


string get_bin_dir ()
{
	return bin_dir;
}

void set_dirs (char* exename)
{    
	 
	 
	 
     int n=strlen(exename);
	 do {
		 --n;
	 }while(n>=0&& exename[n]!='/');
	 if(exename[n]=='/') // path was given
	 { 	
	 	exename[n]=0;
	 	bin_dir=string(exename)+"/";
	 	data_dir=bin_dir+"../data/";
	 	exename[n]='/';
	 	return;
	 }
	 else // no path given - exemine PATH try to find which exename was used
	 { char *p=getenv("PATH");
	   if(p)
	   { 
               string duplicated_path = p;
               size_t idx_begin = 0;
               size_t idx_end = duplicated_path.find(':');
               while (true) {
                   std::string bin_path = duplicated_path.substr(idx_begin, idx_end-idx_begin);
                   if (ifstream((bin_path+"/"+exename).c_str())) {
                       bin_dir=bin_path+"/";
                       data_dir=bin_dir+"../data/";
                       return;
                   }
                   if (idx_end == string::npos) {
                       break;
                   }
                   idx_begin = idx_end+1;
                   idx_end = duplicated_path.find(':', idx_begin);
               }

   	   } 	
	 } 
	 // try to use the environment variable "NUWRO"
	 char *nuwro_dir=getenv("NUWRO");
	 if(nuwro_dir)
	   {bin_dir=string(nuwro_dir)+"/";
	   	data_dir=bin_dir+"data/";
	    return;
	   }
}

bool open_data_file(ifstream& incfile,string name)
{
		incfile.open(name.c_str());
		if(incfile.is_open()) 
			return true;
		incfile.open((data_dir+name).c_str());
		if(incfile.is_open()) 
			return true;
		else
		    return false;
}

