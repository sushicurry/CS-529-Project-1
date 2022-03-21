#ifndef NODE_H
#define NODE_H
/*Header file that contains a class defining nodes if our tree 
* for now this will apply only to the training data.
* Will look at implementing something maybe different for the testing data
*/

#include "Data.h"

struct best_split { //structure to hold information about the best split we will calculate
            double best_info_gain;
            string best_feature;
            int best_pos;
};
//Our Node class
class Node { 
    protected: //private variables
        //pointers to child nodes in the tree
        Node* p_leftChild;
        Node* p_rightChild;
        Data p_data; //our training data
        string p_label;
        

        //Will possibly add more variables as we go
    public:
    best_split p_bs;
    Node(); //default constructor
    Node(const Data& data); //constructor 
    ~Node(); //deconstructor 

    //getters
    Node* getLChild();
    Node* getRChild(); 
    Data  getData();
    string getLabel();
    //setters
    void setLChild(Node* node);
    void setRChild(Node* node);
    void setData(const Data& data); 

    best_split getBestSplit(); //calculates the info gained by splitting the data on each feature and returns the feature with max info gain
    vector<Data> splitData(string feature, int pos); //returns a vector of Data elements containing each split of data split using said feature at given position, one with one without
    void createLabel();

};

Node::Node() {p_leftChild = NULL; p_rightChild = NULL;};



Node::~Node() {};

Node* Node::getLChild() { return p_leftChild; };
Node* Node::getRChild() { return p_rightChild;};
Data Node::getData() {return p_data; };

void Node::setLChild(Node* node) { p_leftChild = node;};
void Node::setRChild(Node* node) {p_rightChild = node;};
void Node::setData(const Data& data) { p_data = data;};

vector<Data> Node::splitData(string feature, int pos) {
    //split the data set into two - one where the feature appears at the given position and one where it does not
    Data has, not_have;
    for(auto it = p_data.cbegin(); it!= p_data.cend(); it++) {
        string sequence = it-> second[0];
        string middle = data::nucleotideAtPos(sequence, pos);
        if(middle == feature) { 
            has[it->first] = it-> second;
        } else { 
            not_have[it->first] = it->second;
        }
    }

    vector<Data> split_data;
    split_data.push_back(has);
    split_data.push_back(not_have);

    return split_data;
};

best_split Node::getBestSplit() { 
    string b_split = ""; //string that holds feature that gives best split
    auto features = data::features; //get features in data
    double best_ig = 0; //best value for info gain
    int best_p = 0;
    best_split b;
    b.best_info_gain = 0;
    b.best_feature = "";
    b.best_pos = 0;
    for(auto feat: features) {
        for(int i = 0; i < 60; i++) {
            vector<Data> split_data = splitData(feat, i);
            auto ig = data::infoGainEnt(p_data, split_data);
            //cout << "Info gain using " << feat <<  " at position " << i<< " is " << ig << endl;
            if(ig > best_ig) {
                best_ig = ig;
                b_split = feat;
                best_p = i;
            }
        }
        
    }

   b.best_info_gain = best_ig;
   b.best_feature = b_split;
   b.best_pos = best_p;

   return b;
};

void Node::createLabel() {
    auto count = data::countClasses(p_data);
    //p_label = to_string(count[0]) + ", " + to_string(count[1]) + ", " + to_string(count[2]);
    if(count[0] > count[1] && count[0] > count[2]) {
        p_label = "N";
    } else if(count[1] >= count[0] && count[1] >= count[2] ) {
        p_label = "IE";
    } else if(count[2] >= count[0] && count[2] >= count[1]) {
        p_label = "EI";
    } else {
        p_label = "Undetermined";
    }
};

string Node::getLabel() {
    return p_label;
};

Node::Node(const Data& data) {
    p_data = data;
    p_leftChild = NULL;
    p_rightChild = NULL;
    p_bs = getBestSplit();
    createLabel();
};


#endif
