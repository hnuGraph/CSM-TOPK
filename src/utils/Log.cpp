//
// Created by ��ɭɭ on 2022/12/7.
//

#include "Log.h"
ofstream* Log::f_track1 = NULL;
ofstream* Log::f_track2 = NULL;
ofstream* Log::f_track3 = NULL;

void Log::init_track1(std::string track_path) {
    if(Log::f_track1!=NULL){
        exit(-1);
    }
    Log::f_track1=new ofstream(track_path.c_str(),ios::out);
    if(!(*f_track1)){
        cout<<"err failed open"<<track_path<<endl;
        exit(-1);
    }
}

void Log::init_track2(std::string track_path) {
    if(Log::f_track2!=NULL){
        exit(-1);
    }
    Log::f_track2=new ofstream(track_path.c_str(),ios::out);
    if(!(*f_track2)){
        cout<<"err failed open"<<track_path<<endl;
        exit(-1);
    }
}

void Log::init_track3(std::string track_path) {
    if(Log::f_track3!=NULL){
        exit(-1);
    }
    Log::f_track3=new ofstream(track_path.c_str(),ios::app);
    if(!(*f_track3)){
        cout<<"err failed open"<<track_path<<endl;
        exit(-1);
    }
}

void Log::track1(std::stringstream &ss) {
    Log::track1(ss.str(),"");
}
void Log::track1(std::string _s, std::string _lat){
    if(Log::f_track1==NULL){
        cout<<"err NULL f_track"<<endl;
        exit(-1);
    }
    *(Log::f_track1)<<_s<<_lat;
    Log::f_track1->flush();
}

void Log::track2(std::stringstream &ss) {
    Log::track2(ss.str(),"");
}
void Log::track2(std::string _s, std::string _lat){
    if(Log::f_track2==NULL){
        cout<<"err NULL f_track"<<endl;
        exit(-1);
    }
    *(Log::f_track2)<<_s<<_lat;
    Log::f_track2->flush();
}

void Log::track3(std::stringstream &ss) {
    Log::track3(ss.str(),"");
}
void Log::track3(std::string _s, std::string _lat){
    if(Log::f_track3==NULL){
        cout<<"err NULL f_track"<<endl;
        exit(-1);
    }
    *(Log::f_track3)<<_s<<_lat;
    Log::f_track3->flush();
}



void Log::close() {
    Log::f_track1->close();
    Log::f_track2->close();
    Log::f_track3->close();
}
void Log::finalize() {

    if(Log::f_track1!=NULL){
        Log::f_track1->close();
        delete Log::f_track1;
        Log::f_track1=NULL;
    }
    if(Log::f_track2!=NULL){
        Log::f_track2->close();
        delete Log::f_track2;
        Log::f_track2=NULL;
    }
    if(Log::f_track3!=NULL){
        Log::f_track3->close();
        delete Log::f_track3;
        Log::f_track3=NULL;
    }

}