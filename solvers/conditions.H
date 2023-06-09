
// -------------- conditions for head -------------
labelList cellsHfx; std::vector<float> cellsHfxData;std::vector<double> condHfx;
labelList cellsHw; std::vector<float> cellsHwData;std::vector<double> condHw;
labelList cellsHdrn; std::vector<float> cellsHdrnData;
labelList cellsHriv; std::vector<float> cellsHrivData;
labelList cellsHghb; std::vector<float> cellsHghbData;
labelList cellsHrch; std::vector<float> cellsHrchData;

// reading the data for hfix
std::ifstream inputHfxdata{cur_dir+"/constant/options/hfix"}; //
cellsHfxData = {std::istream_iterator<float>{inputHfxdata}, {}};
if (cellsHfxData.size()>0) { cellSet selectedHfx(mesh, "hfix"); cellsHfx = selectedHfx.toc();}
condHfx.resize(cellsHfx.size());

// reading the hwell data for wells
std::ifstream inputHwdata{cur_dir+"/constant/options/hwel"}; //
cellsHwData = {std::istream_iterator<float>{inputHwdata}, {}};
if (cellsHwData.size()>0) {	cellSet selectedHwel(mesh, "hwel");cellsHw = selectedHwel.toc();}
condHw.resize(cellsHw.size());

// reading the hdrn data for drains
std::ifstream inputHdrndata{cur_dir+"/constant/options/hdrn"}; //
cellsHdrnData = {std::istream_iterator<float>{inputHdrndata}, {}};
if (cellsHdrnData.size()>0){cellSet selectedHdrn(mesh, "hdrn");cellsHdrn = selectedHdrn.toc();}

// reading the hriv data for riv cells
std::ifstream inputHrivdata{cur_dir+"/constant/options/hriv"}; //
cellsHrivData = {std::istream_iterator<float>{inputHrivdata}, {}};
if (cellsHrivData.size()>0){cellSet selectedHriv(mesh, "hriv");cellsHriv = selectedHriv.toc();}
	
// reading the hghb data for GHb cells
std::ifstream inputHghbdata{cur_dir+"/constant/options/hghb"}; //
cellsHghbData = {std::istream_iterator<float>{inputHghbdata}, {}};
if (cellsHghbData.size()>0){cellSet selectedHghb(mesh, "hghb");cellsHghb = selectedHghb.toc();}

// reading the rch data for Recharge cells
std::ifstream inputHrchdata{cur_dir+"/constant/options/hrch"}; //
cellsHrchData = {std::istream_iterator<float>{inputHrchdata}, {}};
if (cellsHrchData.size()>0){cellSet selectedHrch(mesh, "hrch");cellsHrch = selectedHrch.toc();}
