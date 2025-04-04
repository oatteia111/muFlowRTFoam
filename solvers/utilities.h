#include <dirent.h>
#include <unistd.h>
#include <stdio.h>

std::vector<int> indexC(labelList &cells, std::vector<float> &data)
{
    std::vector<int> c1(cells.size(),0);//Info<<"in index "<<cells.size()<<endl;
	for (int i=0; i<cells.size();i++)  // reads the first ncells lines
		{
		auto iter = std::find(cells.begin(), cells.end(), static_cast<int>(data[i*4+1]));
		int i1 = {std::distance(cells.begin(), iter)};  //Info << i << " icd "<< icd << " cll " << cells_[i] << " indx " << a << endl; // cell number in cellsData //index of cellsData in cells_
		c1[i] = cells[i1];
		}
    return c1;
}

inline bool fexists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

struct outData {float t; std::vector<float> d;};
 
// a function to get data from binary file
outData getCbuffer(string fname, int itime, int ncell) 
{
    std::vector<float> data(ncell*4);
	std::ifstream inputData{cur_dir+"/constant/options/"+fname, std::ios::binary}; //
	inputData.seekg(ncell*itime*4*sizeof(float)); //each line is composed of 4 numbers and there are two values at the beginning

    inputData.read(reinterpret_cast<char*>(&data[0]), ncell*4*sizeof(float));
	float time;	
	inputData.read(reinterpret_cast<char*>(&time), sizeof(float));
	outData output;
	output.t = time;std::cout<<"readbin "<<fname <<" "<<itime<<" "<<ncell;
	output.d = data;
    return output;
}

struct outTable {std::vector<std::string> headers; std::vector<float> data;int ncol;int nrow;};

//a function to read a table from file (1st line nrow ncol, then one header/row and data
outTable readTable(string fname)
{
	std::vector<std::string> headers;
	std::string h;
	std::ifstream if0(fname);
	int nrow, ncol;float d;
	if0 >> nrow; if0 >> ncol;std::cout<<"in readtable "<<nrow<<" "<<ncol<<"\n";
	headers.resize(nrow);
	std::vector<float> data;
	data.resize(nrow*ncol);
	for (int i=0;i<nrow;i++) {
		if0 >> h;headers[i]=h;//std::cout<<" "<<headers[i];
		for (int j=0;j<ncol;j++) {if0 >> d;data[i*ncol+j]=d;}//std::cout<<" "<<d;}
		}
	outTable output;
	output.headers = headers;
	output.ncol=ncol;output.nrow=nrow;
	output.data = data;
	return output;
}

//delete all files in a folder (from chatgpt)
void deleteFilesInDirectory(const std::string& path) {
    DIR *directory = opendir(path.c_str());
    struct dirent *entry;

    if (directory != NULL) {
        while ((entry = readdir(directory)) != NULL) {
            if (strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0) {
                std::string filepath = std::string(path) + "/" + std::string(entry->d_name);
                if (remove(filepath.c_str()) != 0) {
                    std::cerr << "Erreur lors de la suppression de " << filepath << std::endl;
                } else {
                    std::cout << "Fichier supprimé : " << filepath << std::endl;
                }
            }
        }
        closedir(directory);
    } else {
        std::cerr << "Impossible d'ouvrir le répertoire " << path << std::endl;
    }
}