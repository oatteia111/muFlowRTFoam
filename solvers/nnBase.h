#include "torch/torch.h"
#include <ATen/ATen.h>

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <iterator>
#include <random>

// this base is from https://kailaix.github.io/ADCME.jl/v0.3/pytorchnn/

struct Net : torch::nn::Module {
  Net(int n_in,int n_hidden,int n_out) {
    fc1 = register_module("fc1", torch::nn::Linear(n_in, n_hidden));
    fc2 = register_module("fc2", torch::nn::Linear(n_hidden, n_hidden));
    fc3 = register_module("fc3", torch::nn::Linear(n_hidden, n_out));
    //fc4 = register_module("fc4", torch::nn::Linear(n_hidden, n_out));
  }

  torch::Tensor forward(torch::Tensor x) {
    x = torch::tanh(fc1->forward(x));
    x = torch::tanh(fc2->forward(x));
    x = torch::tanh(fc3->forward(x)); // NB : relu produces only positive values
    //x = torch::tanh(fc4->forward(x));
    //x = torch::tanh(fc5->forward(x));
    return x;
  }

  torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr};//, fc4{nullptr};//, fc5{nullptr};
};

// defintion close to https://github.com/pytorch/pytorch/issues/82933
class CustomDataset : public torch::data::datasets::Dataset<CustomDataset> {
    using Example = torch::data::Example<>;

    public:
	torch::Tensor __data;
	torch::Tensor __target;
		
	CustomDataset(torch::Tensor& data, torch::Tensor& target):__data(data),__target(target) {}
	
    torch::optional<size_t> size() const {
        return __target.sizes()[0];
    };
	
	torch::data::Example<> get(size_t index) override {
		torch::Tensor sample_data = __data[index];
		torch::Tensor sample_target = __target[index];
		return {sample_data.clone(), sample_target.clone()};
	};
	
};

class my_NN {
	public:
	
	std::vector<float> data;
	std::vector<float> target;
	std::vector<int> parms {10,12,1}; // network
	std::vector<float> rparms {6.,16.,1e-3,0.75}; //run params (epoc,batch,learning rate,%train)

	int trained=0;
	torch::Tensor output;
	int nd,nvar;
	void setNNParms(std::vector<int> parms) {this->parms = parms;}
	void setRunParms(std::vector<float> rparms) {this->rparms = rparms;}
	void setData(std::vector<float> data){this->data = data;}
	void setTarget(std::vector<float> target){this->target = target;}
	void setNd(int nd){this->nd = nd;}
	//https://www.cppstories.com/2015/02/non-static-data-members-initialization/#the-case-with-auto	
	// valid only with c++17
	static inline auto nn = std::make_shared<Net>(10,12,9);
	
	void init() {
		this->nn = std::make_shared<Net>(this->parms[0],this->parms[1],this->parms[2]); //static inline 
		//torch::nn::init::xavier_uniform_(this->nn->fc1->weight, 1.0);
	}
	
	
	//------- train the model (also test on the last 25% of the dataset)
	float train() {	
		auto options = torch::TensorOptions().dtype(torch::kFloat);
		torch::manual_seed(0.5);
		//torch::optim::Adam optimizer(this->nn->parameters(), torch::optim::AdamOptions(this->rparms[2]));
		torch::optim::SGD optimizer(this->nn->parameters(), torch::optim::SGDOptions(this->rparms[2]));
		int n_in = this->parms[0];
		int n_out = this->parms[2];
		int n_epochs = static_cast<int>(this->rparms[0]); // Number of epochs
		float ls=0.;
		int batch_size = static_cast<int>(this->rparms[1]);
		int nd1 = this->nd*this->rparms[3]; //parm 3 is the prportion of maple to use for training
		std::vector<float> data1 = {this->data.begin(), this->data.begin()+nd1*n_in};
		std::vector<float> target1 = {this->target.begin(), this->target.begin()+nd1*n_out};

		torch::Tensor tdata = torch::from_blob(data1.data(), {nd1,n_in}, torch::TensorOptions(torch::kFloat));//.to(torch::kFloat64);
		torch::Tensor ttarget = torch::from_blob(target1.data(), {nd1,n_out},torch::TensorOptions(torch::kFloat));//.to(torch::kDouble);
		auto dataset = CustomDataset(tdata, ttarget).map(torch::data::transforms::Stack<>()); // map is needed even if nothing in
		auto dataloader = torch::data::make_data_loader<torch::data::samplers::SequentialSampler>(std::move(dataset),torch::data::DataLoaderOptions().batch_size(batch_size));
		this->nn->train(); 
		//std::cout<<"parms "<<this->nn->parameters()<<std::endl;
		namespace F = torch::nn::functional;
		for(int epoch=1; epoch<=n_epochs; epoch++) {
			//std::cout << "Train Epoch: "<< epoch <<std::endl;//<<"b size "<<batch.size()
			optimizer.zero_grad();int nb=0; //zero_grad must be here not in the loop below
			for (torch::data::Example<>& batch : *dataloader) {
				auto data = batch.data; 
				auto target = batch.target; 
				// Clear the optimizer parameters, create output and loss fct
				auto output = this->nn->forward(data);
				//auto loss = torch::mse_loss(output, target);
				auto loss = F::mse_loss(output, target, F::MSELossFuncOptions(torch::kSum));
				// Backpropagate the loss
				loss.backward();
				// Update the parameters
				optimizer.step();nb+=1;
				ls = loss.item<float>();
			  }
			if (epoch==n_epochs-1) {std::cout << " epoc "<<epoch<<" loss " << ls << " nb "<<nb<<std::endl;}
		}
		//std::cout<<"parm "<<nn->parameters()<<std::endl; // ok parameters vary
		// testing
		char sep ='/'; //don't understand why here sep is needed as char while in muFlow just work directly
		std::ofstream outNNoutput(cur_dir+sep+"NNoutput.txt");
		outNNoutput<<this->nn->forward(tdata);
		std::ofstream outNNdata1{cur_dir+sep+"NNdata1.txt"};
		std::ofstream outNNtarget1(cur_dir+sep+"NNtarget1.txt");
		std::ofstream outNNoutput1(cur_dir+sep+"NNoutput1.txt");
		int nd2 = nd -nd1;
		data1 = {this->data.begin()+nd1*n_in, this->data.end()};
		target1 = {this->target.begin()+nd1*n_out, this->target.end()};
		tdata = torch::from_blob(data1.data(), {nd2,n_in}, torch::TensorOptions(torch::kFloat));//.to(torch::kFloat64);
		outNNdata1<<tdata;
		ttarget = torch::from_blob(target1.data(), {nd2,n_out},torch::TensorOptions(torch::kFloat));//.to(torch::kDouble);
		outNNtarget1<<ttarget;
		this->nn->eval(); 
		output = this->nn->forward(tdata);
		outNNoutput1<<output; //torch::nn::MSELossOptions(
		auto loss = F::mse_loss(output, ttarget, F::MSELossFuncOptions(torch::kSum));std::cout<<"loss "<<loss<<"\n";//here loss is mean((ypred-y)**2)
		//float sqr = 0;
		//for (int i=0;i<nd2;i++) {sqr+=std::pow(target1[i],2.);}
		return loss.item<float>()/n_out;
	}

	void eval() {
		this->nn->eval(); 
		int n_in= this->parms[0];int n_out=this->parms[2];
		torch::Tensor tdata = torch::from_blob(this->data.data(), {this->nd,n_in}, torch::TensorOptions(torch::kFloat));//.to(torch::kFloat64);
		std::cout<<"nd "<<this->nd<<" data_eval size "<<tdata.sizes()<<std::endl;
		this->output = this->nn->forward(tdata);
		std::cout<<"output size "<<this->output.sizes()<<std::endl;
	}
};