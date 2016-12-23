//
// Created by aaron on 12/23/16.
// nocontact.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"

#define CONTACT_VERSION_STRING "1.0" __DATE__ " " __TIME__

/* Contact Constants (set by configuration file). */

extern "C" void contact_init(job_t *job, size_t id);

extern "C" void resolve_contact(job_t *job, size_t id);

/*----------------------------------------------------------------------------*/
void contact_init(job_t *job, size_t id) {

    std::cout << "Contact properties (N/A).\n";

    /*printf("%s:%s: (contact version %s) done initializing contact.\n",
           __FILE__, __func__, CONTACT_VERSION_STRING);*/
    std::cout << "Done initializing contact {" << job->contacts[id].bodyIDs[0]  << ", " << job->contacts[id].bodyIDs[1] << "}.\n";
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void resolve_contact(job_t *job, size_t id) {
    return;
}
/*----------------------------------------------------------------------------*/