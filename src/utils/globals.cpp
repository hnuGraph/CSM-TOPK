#include "globals.h"

std::atomic<bool> reach_time_limit = {false};
long long windowSize=500000;
uint ts=0;
uint k=100 ;
float mw=100001;
bool print_result= false;