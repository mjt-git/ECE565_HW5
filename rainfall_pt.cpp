//
//  main.cpp
//  ece565hw5
//
//  Created by mac on 11/15/20.
//  Copyright © 2020 mac. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <sstream>
#include <utility>
#include <float.h>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <thread>
#include <mutex>

using namespace std;

mutex** mtxes;
mutex countMtx;

void initElevations(vector<vector<int>> & elevations, const char * filePath) {
    ifstream infile(filePath);
    string line;
    while(getline(infile, line)) {
        vector<int> vec;
        stringstream ss(line);
        string s = "";
        while(ss >> s) {
            vec.push_back(stoi(s));
        }
        elevations.push_back(vec);
    }
}

// each cell of lowests stores coordinates of neighbors that have same lowest elevation
void initLowests(vector<vector<vector<pair<int, int>>>> & lowests, vector<vector<int>> & elevations, const int N) {
    int dr[4] = {-1, 1, 0, 0};
    int dc[4] = {0, 0, -1, 1};
    for(int r = 0; r < N; ++r) {
        for(int c = 0; c < N; ++c) {
            vector<pair<int, pair<int, int>>> vec;
            for(int k = 0; k < 4; ++k) {
                int nextR = r + dr[k];
                int nextC = c + dc[k];
                if(nextR >= 0 && nextR < N && nextC >= 0 && nextC < N) {
                    vec.push_back(make_pair(elevations[nextR][nextC], make_pair(nextR, nextC)));
                }
            }
            sort(vec.begin(), vec.end(), [](const pair<int, pair<int, int>> p1, const pair<int, pair<int, int>> p2) -> bool {
                return p1.first < p2.first;
            });
            if(vec[0].first >= elevations[r][c]) {continue;}    // current cell's elevation is the lowest among its neighbors
            int minVal = vec[0].first;
            for(size_t m = 0; m < vec.size(); ++m) {
                if(vec[m].first > minVal) {break;}
                lowests[r][c].push_back(make_pair(vec[m].second.first, vec[m].second.second));
            }
        }
    }
}

bool isZero(vector<vector<float>> & states, int N) {
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            if(abs(states[i][j]) > FLT_EPSILON) {return false;}
        }
    }
    return true;
}

void processSimulation(vector<vector<vector<pair<int, int>>>> & lowests, vector<vector<float>> & absorbed, vector<vector<float>> & states, vector<vector<float>> & trickled, const int id, const int N, int timeStep, const int timeOfRain, const float absorbRate, const int threadNum, int & count) {
    int threadN = N / threadNum;    // threadN is the number of rows that one thread needs to deal with
    // Receive a new raindrop (if it is still raining) for each point
    // If there are raindrops on a point, absorb water into the point
    // Calculate the number of raindrops that will next trickle to the lowest neighbor(s)
    for(int i = id * threadN; i < (id+1) * threadN; ++i) {
        for(int j = 0; j < N; ++j) {
            if(timeStep <= timeOfRain) {
                states[i][j] += 1.0;
            }
            if(states[i][j] > 0.0) {
                float absorbedVal = states[i][j] > absorbRate ? absorbRate : states[i][j];
                states[i][j] -= absorbedVal;
                absorbed[i][j] += absorbedVal;
            }
            if(lowests[i][j].size() > 0 && states[i][j] > 0.0) {
                float trickledValTotal = states[i][j] >= 1.0 ? 1.0 : states[i][j];
                states[i][j] -= trickledValTotal;
                
                float trickledVal = trickledValTotal / lowests[i][j].size();
                
                for(size_t k = 0; k < lowests[i][j].size(); ++k) {
                    int neiR = lowests[i][j][k].first;
                    int neiC = lowests[i][j][k].second;
                    mtxes[neiR][neiC].lock();
                    trickled[neiR][neiC] += trickledVal;     // need to be locked
                    mtxes[neiR][neiC].unlock();
                }

            }
        }
    }
    
    countMtx.lock();
    count++;
    countMtx.unlock();
    while(count != threadNum) {}
    
    // update the number of raindrops at each lowest neighbor
    for(int i = id * threadN; i < (id+1) * threadN; ++i) {
        for(int j = 0; j < N; ++j) {
            states[i][j] += trickled[i][j];
            trickled[i][j] = 0.0;
        }
    }
}

int main(int argc, const char * argv[]) {
    if(argc != 6) {
        cout << "WRONG INPUT! Should be ./rainfall <P> <M> <A> <N> <elevation_file>" << endl;
    }
    
    // initialization & fetch arguments
    const int threadNum = atoi(argv[1]);
    const int timeOfRain = atoi(argv[2]);
    char * pEnd;
    const float absorbRate = strtof(argv[3], &pEnd);
    const int N = atoi(argv[4]);
    
    vector<vector<int>> elevations;
    initElevations(elevations, argv[5]);
    
    vector<vector<vector<pair<int, int>>>> lowests(N, vector<vector<pair<int, int>>>(N, vector<pair<int, int>>()));
    initLowests(lowests, elevations, N);

    mtxes = new mutex*[N];
    for(int i = 0; i < N; ++i) {
        mtxes[i] = new mutex[N];
    }
    
    // process the simulation
    vector<vector<float>> absorbed(N, vector<float>(N, 0.0));    // calculate how many rains has been absorbed to each cell
    vector<vector<float>> states(N, vector<float>(N, 0.0));      // represent state of each cell at certain time stamp
    vector<vector<float>> trickled(N, vector<float>(N, 0.0));    // calculated how many rains are trickled to this cell at thie time stamp
    
    const clock_t beginTime = clock();
    
    int timeStep = 1;
    while(true) {
        int count = 0;
        vector<thread> threads(threadNum);
        for(int i = 0; i < threadNum; ++i) {
            threads[i] = thread(processSimulation, ref(lowests), ref(absorbed), ref(states), ref(trickled), i, N, timeStep, timeOfRain, absorbRate, threadNum, ref(count));
        }
        for(int i = 0; i < threadNum; ++i) {
            threads[i].join();
        }
        
        // check if states of all points are zero, if so, break;
        if(isZero(states, N)) {
            break;
        }
        ++timeStep;
    }
    
    const clock_t endTime = clock();
    float timeUsed = float(endTime - beginTime) / CLOCKS_PER_SEC;
    
    cout << "Rainfall simulation completed in " << timeStep << " time steps" << endl;
    cout << "Runtime = " << timeUsed << " seconds" << endl;
    cout << "The following grid shows the number of raindrops absorbed at each point: \n" << endl;
    
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            cout << setw(8) << setprecision(6) << absorbed[i][j];
        }
        cout << endl;
    }

    for(int i = 0; i < N; ++i) {
        delete[] mtxes[i];
    }
    delete[] mtxes;
    
    return 0;
}
