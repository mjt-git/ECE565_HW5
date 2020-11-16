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
#include <atomic>
#include <condition_variable>
#include <pthread.h>
#include <omp.h>

using namespace std;

pthread_mutex_t** mtxes;
pthread_mutex_t globalFinishedMtx = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t printMtx = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barrier1;
pthread_barrier_t barrier2;
pthread_barrier_t barrier3;
pthread_barrier_t barrier4;

bool globalFinished = true;
atomic<int> roundNum(0);
int cnt = 0;

int threadNum;
int timeOfRain;
char * pEnd;
float absorbRate;
int N;

vector<vector<int>> elevations;
vector<vector<vector<pair<int, int>>>> lowests;
vector<vector<float>> absorbed;    // calculate how many rains has been absorbed to each cell
vector<vector<float>> states;      // represent state of each cell at certain time stamp
vector<vector<float>> trickled;    // calculated how many rains are trickled to this cell at thie time stamp


void initElevations(const char * filePath) {
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
void initLowests() {
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

void * processSimulation(void* val) {
    int id = *((int*) val);
    int threadN = N / threadNum;    // threadN is the number of rows that one thread needs to deal with
    // Receive a new raindrop (if it is still raining) for each point
    // If there are raindrops on a point, absorb water into the point
    // Calculate the number of raindrops that will next trickle to the lowest neighbor(s)
    int timeStep = 1;
    while(true) {
        // cout << "thread " << id << " timeStep: " << timeStep << endl;

        bool isFinished = true;

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
                        pthread_mutex_lock(&mtxes[neiR][neiC]);
                        trickled[neiR][neiC] += trickledVal;     // need to be locked
                        pthread_mutex_unlock(&mtxes[neiR][neiC]);
                    }
                }
            }
        }

        // barrier needed
        // cout << "thread " << id << " before barrier 1" << endl;
        pthread_barrier_wait(&barrier1);
        // cout << "thread " << id << " after barrier 1" << endl;

        // for test
        // pthread_mutex_lock(&printMtx);
        // cout << "thread " << id << " states matrix: " << endl;
        // for(int i = id * threadN; i < (id+1) * threadN; ++i) {
        //     for(int j = 0; j < N; ++j) {
        //         cout << setw(8) << setprecision(6) << states[i][j];
        //     }
        //     cout << endl;
        // }
        // cout << endl;
        // pthread_mutex_unlock(&printMtx);
        
        // update the number of raindrops at each lowest neighbor
        for(int i = id * threadN; i < (id+1) * threadN; ++i) {
            for(int j = 0; j < N; ++j) {
                states[i][j] += trickled[i][j];
                trickled[i][j] = 0.0;
                if(abs(states[i][j]) > FLT_EPSILON) {
                    isFinished = false;
                    // cout << "inside isFinished false" << endl;
                }
            }
        }

        pthread_mutex_lock(&globalFinishedMtx);
        // cout << "thread " << id << " globalFinished && isFinished:" << endl;
        // cout << "globalFinished: " << globalFinished << ", " << "isFinished" << isFinished << endl;
        globalFinished = globalFinished && isFinished;
        pthread_mutex_unlock(&globalFinishedMtx);

        // cout << "thread " << id << " before barrier 2" << endl;
        pthread_barrier_wait(&barrier2);
        // cout << "thread " << id << " after barrier 2" << endl;

        if(globalFinished) {
            roundNum.store(timeStep, memory_order_relaxed);
            // cout << "thread " << id << " inside globalFinished TRUE" << endl;
            // cout << "thread " << id << " inside globalFinished TRUE timeStep: " << timeStep << endl;
            break;
        } 
        // cout << "thread " << id << " before barrier 3" << endl;
        pthread_barrier_wait(&barrier3);
        // cout << "thread " << id << " after barrier 3" << endl;

        pthread_mutex_lock(&globalFinishedMtx);
        globalFinished = true;
        pthread_mutex_unlock(&globalFinishedMtx);
        timeStep++;

        pthread_barrier_wait(&barrier4);
    }
}

int main(int argc, const char * argv[]) {
    if(argc != 6) {
        cout << "WRONG INPUT! Should be ./rainfall <P> <M> <A> <N> <elevation_file>" << endl;
    }
    
    // initialization & fetch arguments
    threadNum = atoi(argv[1]);
    timeOfRain = atoi(argv[2]);
    absorbRate = strtof(argv[3], &pEnd);
    N = atoi(argv[4]);
    pthread_barrier_init (&barrier1, NULL, threadNum);
    pthread_barrier_init (&barrier2, NULL, threadNum);
    pthread_barrier_init (&barrier3, NULL, threadNum);
    pthread_barrier_init (&barrier4, NULL, threadNum);

    lowests = vector<vector<vector<pair<int, int>>>>(N, vector<vector<pair<int, int>>>(N, vector<pair<int, int>>()));
    absorbed = vector<vector<float>>(N, vector<float>(N, 0.0));    // calculate how many rains has been absorbed to each cell
    states = vector<vector<float>>(N, vector<float>(N, 0.0));      // represent state of each cell at certain time stamp
    trickled = vector<vector<float>>(N, vector<float>(N, 0.0));    // calculated how many rains are trickled to this cell at thie time stamp


    initElevations(argv[5]);
    initLowests();

    // for test
    // cout << "elevations: " << endl;
    // for(int i = 0; i < N; ++i) {
    //     for(int j = 0; j < N; ++j) {
    //         cout << setw(8) << setprecision(6) << elevations[i][j];
    //     }
    //     cout << endl;
    // }
    // cout << endl;


    mtxes = new pthread_mutex_t*[N];
    for(int i = 0; i < N; ++i) {
        mtxes[i] = new pthread_mutex_t[N];
        for(int j = 0; j < N; ++j) {
            mtxes[i][j] = PTHREAD_MUTEX_INITIALIZER;
        }
    }
    
    // process the simulation
    const double beginTime = omp_get_wtime();

    pthread_t** threads = new pthread_t*[threadNum];
    int* args = new int[threadNum];
    for(int i = 0; i < threadNum; ++i) {
        args[i] = i;
        threads[i] = new pthread_t;
        pthread_create(threads[i], NULL, processSimulation, (void*)&(args[i])); 
    }

    // cout << "before join" << endl;
    for(int i = 0; i < threadNum; ++i) {
        pthread_join(*threads[i], NULL);
    }
    // cout << "after join" << endl;
    
    const double endTime = omp_get_wtime();
    double timeUsed = endTime - beginTime;
    
    cout << "Rainfall simulation completed in " << roundNum.load(memory_order_relaxed) << " time steps" << endl;
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

    delete[] args;
    for(int i = 0; i < threadNum; ++i) {
        delete threads[i];
    }
    delete[] threads;
    
    return 0;
}
