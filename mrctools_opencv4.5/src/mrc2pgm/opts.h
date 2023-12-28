#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C" {
#include <getopt.h>
}
#include "util/exception.h"

struct options {
    std::string input;
    std::string idxs;
    std::string output;
};

inline int GetOpts(int argc, char **argv, options* opts_){

    static struct option longopts[] = {
        { "help",            no_argument,            NULL,              'h' },
        { "output",          required_argument,      NULL,              'o' },
        { "input",    	     required_argument,      NULL,              'i' },
		{ "pages",    	     required_argument,      NULL,              'p' },
        { NULL,              0,                      NULL,               0  }
    };
	
    if((argc != 7) && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1){
		EX_TRACE("[-i INPUT FILENAME][-p INDEX,INDEX][-o OUTPUT FOLDER]\n");
		return -1;
    }
    
    int ch;
    while ((ch = getopt_long(argc, argv, "hi:o:p:", longopts, NULL))!= -1) {
        switch (ch) {

        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;

        case 'h':
            EX_TRACE("[-i INPUT FILENAME][-p INDEX,INDEX][-o OUTPUT FOLDER]\n");
			return 0;

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
				return -1;
			}
        }
        break;
		
		case 'p': 
        {
            std::istringstream iss(optarg);
            iss >> opts_->idxs;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
				return -1;
			}
        }
        break;
        
        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
		
        case 0: 
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif
