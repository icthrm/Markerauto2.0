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
    std::string xform;
    std::string output;
    int mode;
};

inline int GetOpts(int argc, char **argv, options* opts_){

    static struct option longopts[] = {
        { "help",            no_argument,            NULL,              'h' },
        { "output",          required_argument,      NULL,              'o' },
        { "xform",           required_argument,      NULL,              'x' },
        { "input",    	     required_argument,      NULL,              'i' },
        { "mode",    	     required_argument,      NULL,              'm' },
        { NULL,              0,                      NULL,               0  }
    };
	
    if((argc != 9)&&(argc != 7) && argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1){
		EX_TRACE("[-i INPUT FILENAME][-x TRANSFORM][-o OUTPUT FILENAME]([-m MODE])\n");
		EX_TRACE("  MODE is the interpolation modes to transform images.\n");
		EX_TRACE("  MODE can be 1 for bicubic interpolation over 4x4 pixel neighborhood, 2 for lanczos interpolation over 8x8 pixel neighborhood.\n");
          
		return -1;
    }
    
    int ch;
    while ((ch = getopt_long(argc, argv, "hi:o:x:m:", longopts, NULL))!= -1) {
        switch (ch) {

        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;

        case 'h':
		EX_TRACE("[-i INPUT FILENAME][-x TRANSFORM][-o OUTPUT FILENAME]([-m MODE])\n");
		EX_TRACE("  MODE is the interpolation modes to transform images.\n");
		EX_TRACE("  MODE can be 1 for bicubic interpolation over 4x4 pixel neighborhood, 2 for lanczos interpolation over 8x8 pixel neighborhood.\n");
			return 0;
        case 'm': 
        {
            std::istringstream iss(optarg);
            iss >> opts_->mode;
            if (iss.fail()){
                EX_TRACE("Invalid argument '%s'.", optarg);
                return -1;
            }
        }
        break;
	
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
		
        case 'x':
        {
            std::istringstream iss(optarg);
            iss >> opts_->xform;
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
