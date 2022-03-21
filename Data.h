#ifndef DATA_H
#define DATA_H
/*Header file containing a data struct. This will contain our analysis functions
*/


#include<string>
#include<vector>
#include <map>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;
using boost::math::chi_squared;
using boost::math::quantile;
//alias for our data container which is a map
using Data = map<int, vector<string>>;

namespace data {
        array<string, 4> features = {"A", "C", "G", "T"};
        array<string, 16> double_features = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "TA", "TC", "TG", "TT", "GA", "GC", "GT", "GG"};
        //functions to analyze the data
        array<int,3> countClasses(Data p_data); //returns an arry containing count of each outcome in our data
        double missclassError(Data p_data); //returns the missclassification error of our data
        double entropy(Data p_data) ;//returns the entropy of our data
        double giniIndex(Data p_data); //returns the calculated gini index of our data

        double infoGainME(Data p_data, vector<Data> split_data); //calculates info gain using p_data as original data and missclassification error
        double infoGainEnt(Data p_data, vector<Data> split_data); //calculates info gain using p_data as original data and entropy
        double infoGainGI(Data p_data, vector<Data> split_data); //using GI

        string nucleotideAtPos(string dna_sequence, int pos);
        
        vector<double> getProbs(Data data);
        double chiSquareCV(vector<Data> split, vector<double> probs);
        bool isSignificant(vector<Data> split, vector<double> probs, double alpha);

        //functions to read data from csv
        string readFileToString(const string path);
        Data getDataFromFile(const string path);
        void writePredictionsToFile(vector<string> predictions);

};

//Defining our functions

//Our analysis functions
array<int, 3> data::countClasses(Data p_data) {
    array<int, 3> count = {0, 0, 0}; //initial array. In order the numbers will represent N, IE and EI

    for(auto it = p_data.cbegin(); it!= p_data.cend(); it++) { 
        //for loop iterating over our data map
        //if the second value in our vector<string> matches one of our classes we add one to that element
        //of our array. 
        if(it->second[1] == "N") { 
            count[0]++;
        } else if(it -> second[1] == "IE") { 
            count[1]++;
        } else if(it -> second[1] == "EI") { 
            count[2]++;
        }
    }

    return count;
};

//function to calculate the missclassification error of our data
double data::missclassError(Data p_data) {
    auto count = countClasses(p_data); //count of our classes
    double err; //what we will return
    int total = 0; //total number of values in our data
    for(auto i = 0; i < count.size(); i++) { 
        total+=count[i];
    }

    double probs[count.size()]; //array to calculate the probability of each class in our data
    for(int i = 0; i < count.size(); i++) { 
        probs[i] = static_cast<double>(count[i])/static_cast<double>(total);
    }

    err = 1.0 - *(max_element(probs, probs+count.size())); //max_element returns a pointer to the max element in our array
    return err;
};

//function to calculate the entropy of our data
double data::entropy(Data p_data) {
    auto count = countClasses(p_data);
    double ent = 0;
    int total = 0;

    for(auto i = 0; i < count.size(); i++) { 
        total+=count[i];
    }
    double probs[count.size()]; //array to calculate the probability of each class in our data
    for(int i = 0; i < count.size(); i++) { 
        probs[i] = static_cast<double>(count[i])/static_cast<double>(total);
    }
    for(int i = 0; i < count.size(); i++) {
        if(probs[i] != 0) {
            //cout << probs[i] << ", " << log2(probs[i]) << ", " << -probs[i]*log2(probs[i]) << endl;
            ent+= log2(pow(probs[i], -probs[i]));
            //cout << ent << endl;
        }
    }
    //cout << ent << endl;
    return ent;
};

//function to calcualte the gini index of our data
double data::giniIndex(Data p_data) {
    auto count = countClasses(p_data);
    double gi = 1;
    int total = 0;

    for(auto i = 0; i < count.size(); i++) { 
        total+=count[i];
    }

    double probs[count.size()]; //array to calculate the probability of each class in our data
    for(int i = 0; i < count.size(); i++) { 
        probs[i] = static_cast<double>(count[i])/static_cast<double>(total);
    }

    for(int i = 0; i < count.size(); i++) { 
        gi -= probs[i]*probs[i];
    }

    return gi;
};

double data::infoGainME(Data p_data, vector<Data> split_data) {
    double info_gain; //info gain
    //calculating totals in each split and original
    size_t size = p_data.size();
    size_t size_of_splits[2] = {split_data[0].size(), split_data[1].size()};
    //counting classes
    auto classes = countClasses(p_data);
    auto sp1_classes = countClasses(split_data[0]);
    auto sp2_classes = countClasses(split_data[1]);
    //calculating missclass error
    auto err = missclassError(p_data);
    double split_error[2] = {missclassError(split_data[0]), missclassError(split_data[1])};

    info_gain = err; //initialize info gain
    for(int i = 0; i < 2; i++) { 
        info_gain -= (static_cast<double>(size_of_splits[i])/static_cast<double>(size))*split_error[i];
    }

    return info_gain;
};

double data::infoGainEnt(Data p_data, vector<Data> split_data) {
    double info_gain; //info gain
    //calculating totals in each split and original
    size_t size = p_data.size();
    size_t size_of_splits[2] = {split_data[0].size(), split_data[1].size()};
    //counting classes
    auto classes = countClasses(p_data);
    //cout << "Original class count - " << classes << endl;
    auto sp1_classes = countClasses(split_data[0]);
    //cout << "Split 1 class count - " << sp1_classes << endl;
    auto sp2_classes = countClasses(split_data[1]);
    //cout << "Split 2 class cout - " << sp2_classes << endl;
    //calculating missclass error
    auto err = entropy(p_data);
    //cout << "Original data entropy - " << err << endl;
    double split_error[2] = {entropy(split_data[0]), entropy(split_data[1])};
    //cout << "Split 1 entropy - " << split_error[0] << " , " << "Split 2 entropy - " << split_error[1] << endl;
    info_gain = err; //initialize info gain
    for(int i = 0; i < 2; i++) { 
        info_gain -= (static_cast<double>(size_of_splits[i])/static_cast<double>(size))*split_error[i];
    }

    return info_gain;
};

double data::infoGainGI(Data p_data, vector<Data> split_data) {
    double info_gain; //info gain
    //calculating totals in each split and original
    size_t size = p_data.size();
    size_t size_of_splits[2] = {split_data[0].size(), split_data[1].size()};
    //counting classes
    auto classes = countClasses(p_data);
    auto sp1_classes = countClasses(split_data[0]);
    auto sp2_classes = countClasses(split_data[1]);
    //calculating missclass error
    auto err = giniIndex(p_data);
    double split_error[2] = {giniIndex(split_data[0]), giniIndex(split_data[1])};

    info_gain = err; //initialize info gain
    for(int i = 0; i < 2; i++) { 
        info_gain -= (static_cast<double>(size_of_splits[i])/static_cast<double>(size))*split_error[i];
    }

    return info_gain;
};

string data::nucleotideAtPos(string dna_sequence, int pos) {
    //A function that returns the nucleotide at given position
    string n_split = "";
    //auto midpoint = (double)dna_sequence.size()/2.0;
    n_split = dna_sequence.substr(pos, 1);

    return n_split;
};

string data::readFileToString(const string path) {
    auto ss = ostringstream{};
    ifstream input_file(path);
    if(!input_file.is_open()) { //check to make sure file is open
        cerr << "Could not open the file - '" << path << "'" << endl;
        exit(EXIT_FAILURE);
    }

    ss << input_file.rdbuf(); //read file
    return ss.str(); //return string with file contents
};

Data data::getDataFromFile(const string path) {
    char delim = ','; //what seperates entries in a csv file
    string contents;
    Data data;
    contents = readFileToString(path); //read the file into a string
    
    istringstream sstream(contents); //input stream for string contents
    vector<string> items; //list of items in the csv file
    string record; //string to hold a line
    int counter = 0; //counter for number of entries
    while(getline(sstream, record)) { //while we can get a new line
        /*We know our data is of the form entry#, DNA sequence, class.
         So to remove extreneous information we remove the entry# portion so we are just left with the information that we want to store - DNA sequence and class.
         */
        const auto p = record.find_first_not_of("0123456789"); //find first position in line that is not a number
        if(p != string::npos) {
            record.erase(0, p + 1); //remove said numbers
        }
        istringstream line(record); //store the new line in record
        while(getline(line, record, delim)) { //get parts delimited by commas
            record.erase(remove_if(record.begin(), record.end(), ::isspace), record.end()); //deal with spaces, endline, tab etc
            items.push_back(record); //add items in the same line
        }
        
        data[counter] = items;
        items.clear();
        counter+=1;
    }
    
    return data;
};

double data::chiSquareCV(vector<Data> split, vector<double> probs) {
    int total = 0;
    double cv = 0.0;
    auto c1 = countClasses(split[0]);
    auto c2 = countClasses(split[1]);
    
    for(auto n: c1) {
        total+= n;
    }
    for(auto n: c2) {
        total+= n;
    }
    //cout<< split[0].size() << ", " << split[1].size() << endl;
    for(int i  = 0; i< probs.size(); i++) {
        double exp1 = probs[i]*(split[0].size());
        //cout << c1[i] << ", " << exp1 << endl;
        cv += ((c1[i] - exp1)*(c1[i] - exp1))/exp1;
    }
    for(int i = 0; i < probs.size(); i++) {
        double exp2 = probs[i]*(split[1].size());
        //cout << c2[i] << ", " << exp2 << endl;
        //cout << ((c2[i] - exp2)*(c2[i] - exp2))/exp2 << endl;
        cv += ((c2[i] - exp2)*(c2[i] - exp2))/exp2;
    }
    
    return cv;
};

vector<double> data::getProbs(Data data) {
    auto classes = countClasses(data);
    int total = data.size();
    
    vector<double> probs;
    for(int i = 0; i < classes.size(); i++) {
        probs.push_back(static_cast<double>(classes[i])/static_cast<double>(total));
        //cout<< static_cast<double>(classes[i])/static_cast<double>(total) << endl;
    }
    
    return probs;
};

bool data::isSignificant(vector<Data> split, vector<double> probs, double alpha) {
    chi_squared dist(3);
    bool sig = false;
    auto cv = chiSquareCV(split, probs);
    if(cv >= quantile(dist, 1-alpha)) {
        sig = true;
    }
    
    return sig;
}

void data::writePredictionsToFile(vector<string> predictions) {
        ofstream output;
        output.open("predictions.csv");
        output << "ID, Class, \n";
        for(auto i = 0; i < predictions.size(); i++) {
            output << i + 2001 << "," << predictions[i] << ",\n";
        }

        output.close();
}

#endif
