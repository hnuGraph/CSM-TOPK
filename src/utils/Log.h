//
// Created by ¸ßÉ­É­ on 2022/12/7.
//

#ifndef BASELINE_LOG_H
#define BASELINE_LOG_H
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
class Log {
public:
    static ofstream* f_track1;
    static ofstream* f_track2;
    static ofstream* f_track3;

    static void init_track1(string track_path);
    static void init_track2(string track_path);
    static void init_track3(string track_path);
    static void track1(stringstream & ss);
    static void track1(string _s,string _lat="\n");
    static void track2(stringstream & ss);
    static void track2(string _s,string _lat="\n");
    static void track3(stringstream & ss);
    static void track3(string _s,string _lat="\n");
    static void close();
    static void finalize();
};


#endif //BASELINE_LOG_H
