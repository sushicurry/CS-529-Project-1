#include "Tree.h"
#include <iostream>
#include <chrono>
#include <random>



int main(void) { 
    
    cout << "Cunstructing our Tree" << endl;
    
    Tree test("train.csv", 0.05);
    
    cout << "Getting our predictions" << endl;
    
    auto p = test.predict("test.csv");
    
    cout << "Now writing our predictions to a file. " << endl;
    
    data::writePredictionsToFile(p);
    
    /*//Want to double check against random sampling from training file.
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed);
    
    Data d = data::getDataFromFile("train.csv");
    Data ran;
    for(int i = 0; i < 1000; i++) {
        auto num = (gen() % (999 - 0 + 1));
        ran[i] = d[num];
    }
    
    auto ran_p = test.predict(ran);
    
    ofstream ranout("random_predictions.csv");
    ranout << "ID, Predicted Class, Actual Class, \n";
    auto c = 0;
    for(auto it = ran.cbegin(); it!=ran.cend(); it++) {
        ranout << c + 1 << ", " << ran_p[c] << ", " << it->second[1] << ",\n";
        c++;
    }
    
    ranout.close();*/
    
    

    return 0;
}
