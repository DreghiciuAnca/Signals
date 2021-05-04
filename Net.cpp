#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <math.h> 
using namespace std;
constexpr float FLOAT_MIN = -1;
constexpr float FLOAT_MAX = 1;

struct Connection{
double weight;
};
class Neuron;

typedef vector<Neuron> Layer;
class Neuron{
public:
    Neuron(unsigned numOutputs, unsigned myIndex);
    void setOutputVal(double val);
    double getOutputVal() const {return outputVal;};
    void feedForward(const Layer &prevLayer);
    void updateInputWeights(Layer &prevLayer);

private:
    static double randomWeight(void) {return FLOAT_MIN + (float)(rand()) / ((float)(RAND_MAX/(FLOAT_MAX - FLOAT_MIN)));};
    double outputVal;
    vector<Connection> m_outputWeights;
    unsigned m_myIndex;
    
};



void Neuron::setOutputVal(double val) {
    outputVal = val;
}

void Neuron::updateInputWeights(Layer &prevLayer){

    for(unsigned n=0; n<prevLayer.size(); ++n)
    {
        Neuron &neuron = prevLayer[n];
        double oldWeight = neuron.m_outputWeights[m_myIndex].weight;
        double newWeight =  neuron.getOutputVal() * oldWeight;        
        neuron.m_outputWeights[m_myIndex].weight += newWeight;

    }
}


Neuron::Neuron(unsigned numOutputs, unsigned myIndex){
    for (unsigned c =0 ;c< numOutputs ; ++c){
        m_outputWeights.push_back(Connection());
        m_outputWeights.back().weight = randomWeight();
    }

    m_myIndex =myIndex;
}
void Neuron::feedForward(const Layer &prevLayer){
    double sum = 0.0;

    for(unsigned n= 0; n< prevLayer.size(); ++n)
    {
        sum+= prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myIndex].weight;
    }

    outputVal = sum;
}



class Network{
public:
    Network(const vector<unsigned> &topology);
    void feedForward(const vector<double> &inputVals);
    void backProp();
    void getResults(vector<double>resultVals) const;

private:
    vector<Layer> layers;
    
};

void Network::getResults(vector<double>resultVals) const{
    resultVals.clear();
    for(unsigned n =0; n< layers.back().size() - 1; ++n)
    {
        
        resultVals.push_back(layers.back()[n].getOutputVal());
        cout<<resultVals[n];
    }
}
Network::Network(const vector<unsigned> &topology)
{
    unsigned numLayers = topology.size();
    for(unsigned layerNum = 0 ; layerNum <numLayers; ++layerNum){
        layers.push_back(Layer());
        unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

        for(unsigned neuronNum=0; neuronNum <= topology[layerNum]; ++neuronNum){
            layers.back().push_back(Neuron(numOutputs, neuronNum));
            cout<< "Made a neuron"<< endl;
        }
        layers.back().back().setOutputVal(1.0);
    }
}

void Network::feedForward(const vector<double> &inputVals)
{
    assert(inputVals.size() == layers[0].size() -1);
    for(unsigned i = 0; i < inputVals.size(); ++i)
    {
        layers[0][i].setOutputVal(inputVals[i]);
    }

    for(unsigned layerNum = 1; layerNum < layers.size() ;++layerNum)
    {
        Layer &prevLayer= layers[layerNum -1];
        for(unsigned n = 0; n < layers[layerNum].size() -1 ;++n)
        {
            layers[layerNum][n].feedForward(prevLayer);
        }
    }
}

void Network::backProp()
{
    
    for(unsigned layerNum = layers.size() -1; layerNum >0 ; --layerNum)
    {
        Layer &layer= layers[layerNum];
        Layer &prevLayer = layers[layerNum -1];

        for(unsigned n = 0; n< layer.size() -1; ++n)
        {
            layer[n].updateInputWeights(prevLayer);
        }
    }
}

int main()
{
    vector<unsigned> topology;
    topology.push_back(2); //number of inputs
    topology.push_back(2); //number of nodes in 2 layer
    topology.push_back(1); //number of outputs
    Network myNet(topology);
    
    vector<double> inputVals;
    //values of inputs
    inputVals.push_back(0.5); 
    inputVals.push_back(1.0);
    myNet.feedForward(inputVals);

    myNet.backProp();

    vector<double> resultsVals;
    myNet.getResults(resultsVals);
   
}