//
//  Tree.h
//  
//
//  Created by Sushyam Ravi on 2/26/22.
//

#ifndef Tree_h
#define Tree_h

#include "Node.h"

class Tree {
public:
    Node* root;
    //double alph;
    
    //constructor
    Tree(string path, double alpha);
    //deconstructor
    ~Tree();
    
    //Tree functions
    void traverse() const;
    static vector<Data> split(Data original, string feature, int pos);
    static void construct(Node* node, double alpha);
    void deleteChild(Node* node);
    //predict functions
    static string getPrediction(Node* n, string sequence);
    vector<string> predict(string file);
    vector<string> predict(Data d);
    
};

vector<Data> Tree::split(Data original, string feature, int pos) {
    Data has;
    Data does_not_have;
    
    for(auto it = original.cbegin(); it!= original.cend(); it++) {
        string sequence = it-> second[0];
        string atPos = data::nucleotideAtPos(sequence, pos);
        if(atPos == feature) {
            has[it->first] = it-> second;
        } else {
            does_not_have[it->first] = it->second;
        }
    }

    vector<Data> split_data;
    split_data.push_back(has);
    split_data.push_back(does_not_have);

    return split_data;
}

void Tree::construct(Node* node, double alpha) {
    //This is where we use our chi square test. If the best split calculated gives a significant change we create our child nodes. We do this recursively
    
    auto childData = split(node->getData(), node->p_bs.best_feature, node->p_bs.best_pos);
    
    if(data::isSignificant(childData, data::getProbs(node->getData()), 1 - alpha)) {
        //if true then create child nodes
        //cout << "Significant," << endl;
        if(childData[0].size() > 0) {
           // cout << "Left Child created." << endl;
            Node* lChild = new Node(childData[0]); //left child hss feature
            node->setLChild(lChild);
            construct(node->getLChild(), alpha);
        }
        
        if(childData[1].size() > 0) {
            //cout << "Right child created." << endl;
            Node* rChild = new Node(childData[1]);
            node->setRChild(rChild);
            construct(node->getRChild(),alpha);
        }
    } else {
        //cout << "Not significant" << endl;
    }
};

void Tree::deleteChild(Node* node) {
    if(node->getLChild() != NULL) {
        deleteChild(node->getLChild());
    }
    if(node->getRChild() != NULL) {
        deleteChild(node->getRChild());
    }
    
    delete node;
};

string Tree::getPrediction(Node* node, string sequence) {
    string prediction = "";
    string char_at_bs = sequence.substr(node->p_bs.best_pos, 1);
    
    if(node->getLChild() == NULL && node->getRChild() == NULL) {
        //if no child nodes then
       // cout << "No Child" << endl;
        prediction = node->getLabel();
    } else if(char_at_bs == node->p_bs.best_feature) {
        //if the sequence matches the nucleotide at the best split position then it is in the left child so recurse on left child
        //cout << "Searching left child." << endl;
        prediction = getPrediction(node->getLChild(), sequence);
    } else {
        //has to be on right child
       // cout << "Searching right child." << endl;
        prediction = getPrediction(node->getRChild(), sequence);
    }
    
    return prediction;
};

vector<string> Tree::predict(string file) {
    //load data from a file and then get predictions
    auto test_data = data::getDataFromFile(file);
    
    vector<string> predictions;
    
    for(auto it = test_data.cbegin(); it!=test_data.cend(); it++) {
        predictions.push_back(getPrediction(root, it->second[0]));
    }
    
    return predictions;
};

vector<string> Tree::predict(Data d) {
    vector<string> predictions;
    
    for(auto it = d.cbegin(); it!=d.cend(); it++) {
        predictions.push_back(getPrediction(root, it->second[0]));
    }
    
    return predictions;
}
Tree::Tree(string path, double alpha) {
    Data d = data::getDataFromFile(path);
    Node* r = new Node(d);
    //this->alph = alpha;
    construct(r, alpha);
    this->root = r;

};

Tree::~Tree() {
    //deconstructor just goes through and deletes the tree
    deleteChild(this->root);
}

#endif /* Tree_h */
