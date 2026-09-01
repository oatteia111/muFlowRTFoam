struct Net : torch::nn::Module {
  Net(int n_in,int n_hidden,int n_out) {
    fc1 = register_module("fc1", torch::nn::Linear(n_in, n_hidden));
    fc2 = register_module("fc2", torch::nn::Linear(n_hidden, n_hidden));
    fc3 = register_module("fc3", torch::nn::Linear(n_hidden, n_out));
    //fc3 = register_module("fc3", torch::nn::Linear(n_hidden*2, n_hidden));
    //fc4 = register_module("fc4", torch::nn::Linear(n_hidden, n_out));
  }

  torch::Tensor forward(torch::Tensor x) {
    x = torch::tanh(fc1->forward(x));
    x = torch::tanh(fc2->forward(x));
    x = torch::relu(fc3->forward(x));
    //x = torch::relu(fc4->forward(x));
    return x;
  }

  torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr};
};
