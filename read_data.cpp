//
// Created by thinkcomputational on 25/11/21.
//
#include <string>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream
#include "read_data.h"

using namespace std;

vector<vector<double>> parse2DCsvFile(string inputFileName) {

    vector<vector<double>> data;
    ifstream inputFile(inputFileName);
    int l = 0;

    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> record;

            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }

            data.push_back(record);
        }
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return data;
}

double InterpolateWhetherData(int loc, int prop, int iter) {
    //[i][0]= Tamb, [i][1]= wind, [i][2]=Q
    return MetData[loc-1][prop] + (MetData[loc][prop] - MetData[loc-1][prop])*(iter * dt)/300;
}
