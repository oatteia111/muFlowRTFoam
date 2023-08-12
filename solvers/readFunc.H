// file reader of myfunc
std::string fname;
fname = cur_dir+"/myFunc.txt";
//verify if the file exists
std::ifstream f(fname.c_str());
std::vector<double> parms(1.,0.);

if (f.good())
{
// open the file
std::ifstream if0(fname);
// read by keywords
std::string key;
std::string ftype;
int nb;
double aa;
if0 >> key;//Info<<" key "<<key<<endl;
if (key=="fDe_T")
{
	// read the type
	if0 >> ftype;
	//Info<<" type "<<ftype<<endl;
	fDe_T.setType(ftype);
	// read the number of parameters
	if0 >> nb;
	parms.resize(nb);
	for (i=0;i<nb;i++) {if0 >> aa;parms[i] = aa;}
	fDe_T.setParms(parms);
}
}
else
{
fDe_T.setType("linear");
fDe_T.setParms(parms);
}