#ifndef PARTAGE_H
#define PARTAGE_H

#include <unistd.h>
#define GetCurrentDir getcwd
#include <string>

// Le mot-clé extern dit "la variable existe, mais elle est définie ailleurs"
extern int nxyz,ph_ncomp,ph_nspc,ph_gcomp,ph_nsolu;

std::string get_current_dir();
// On déclare la variable avec extern (SANS le '= ...')
extern std::string cur_dir;
#endif
