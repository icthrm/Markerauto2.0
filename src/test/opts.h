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
    char input[255];
    char inputangle[255];
    char outputxf[255];
    char outputangle[255];
};

inline int GetOpts(int argc, char **argv, opts* opts_){

    static struct option longopts[] = {
        { "help",            no_argument,            NULL,              'h' },
        { "output",          required_argument,      NULL,              'o' },
        { "input",    	     required_argument,      NULL,              'i' },
        { "initangle",       required_argument,      NULL,              'a' },
        { "newangle",        required_argument,      NULL,              'n' },
        { NULL,              0,                      NULL,               0  }
    };
    
    if(argc != 9  && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1){
	EX_TRACE("[-i INPUT FILENAME][-a INIT ANGLE FILENAME][-o OUTPUT FILENAME][-n NEW ANGLE FILENAME]\n");
	return -1;
    }
    
    int ch;
    while ((ch = getopt_long(argc, argv, "ho:i:a:n:", longopts, NULL))!= -1) {
        switch (ch) {

        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;

        case 'h':
            EX_TRACE("  [--input(-i) INPUT FILENAME]\n  [--initangle(-a) INIT ANGLE FILENAME]\n  [--output(-o) OUTPUT FILENAME]\n  [--newangle(-n) NEW ANGLE FILENAME]\n");
            return 0;

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'a': 
        {
            std::istringstream iss(optarg);
            iss >> opts_->inputangle;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->outputxf;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'n':
        {
            std::istringstream iss(optarg);
            iss >> opts_->outputangle;
            if (iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
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