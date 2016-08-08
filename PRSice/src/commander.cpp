#include "commander.h"

Commander::Commander(const int argc, const char *argv[]){
    if(argc<=1){
        usage();
        throw "Please provide the required parameters";
    }
    static const char *optString = "b:t:p:mh?";
    static const struct option longOpts[]={
	    //Qt specific parameter
		{"base",required_argument,NULL,'b'},
		{"target",required_argument,NULL,'t'},
		{"thread",required_argument,NULL,'p'},
		{"multiPheno",no_argument,NULL,'m'},
		{"help",no_argument,NULL,'h'},
		{NULL, 0, 0, 0}
	};
}

Commander::Commander()
{
    //ctor
}

Commander::~Commander()
{
    //dtor
}

void Commander::usage(){
}
