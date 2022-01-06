/**
 * Repoduction of Emsenhuber et al. 2020 Perfect Merger simulations
 *
 * Given an initial simulation file with ~155 planetessimals
 * Program will simulate 400 Myrs of evolution
 * Collisions will be treated as perfecrt Mergers with mass and momentum conserved (not energy)
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "collresolve.h"


void perfect_merge(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

double sec2year = 31536000.0;
double rad2Deg = 57.29577951308232;
double years = 1e9;
//double years = 2e7;
double au2m = 149597870700;
double V0D0 = 100*1; // 100m/s for a 1km obj
double obliquities[8] = {0};

int main(int argc, char* argv[]){
    if(argc != 6){
        printf("Must include 'path_to_sim', 'path_to_arch', 'run_num', 'p_size', and 'size_file_flag  as command line inputs!\nProgram Terminated\n");
        exit(0);
    }

    double tmax = years * sec2year;

    char* path_to_sim = argv[1];
    char* path_to_arch = argv[2];
    unsigned int  run_num = atoi(argv[3]);
    double p_size;
    double p_sizes[8];
    char* size_file_flag = argv[5];
    if (strcmp(size_file_flag, "1") == 0){
	char* file_name = argv[4];
	double all_sizes[20000];
	float tek;
	FILE* size_file = fopen(file_name, "r");
	if (size_file == NULL) {
            perror("Error: Failed to open size file.");
            exit(0);
	}
	// read values from file until EOF is returned by fscanf
	for (int i = 0; i < 20000; ++i) {
	    if (fscanf(size_file, "%f\n", &tek) == 1) {
		all_sizes[i] = tek;
	    } else {
		break;
	    }
	}
	fclose(size_file);
	for (int i = 0; i < 8; i++){
	    p_sizes[i] = all_sizes[4*run_num + i/2];
	}
    } else {
	p_size = atof(argv[4]); // in km
	for (int i = 0; i < 8; i++){
	    p_sizes[i] = p_size;
	}
    }
    printf("Run # %i\n", (unsigned int) run_num);
    
    srand((unsigned int)run_num); // set random seed based off ever changing run_num

    struct reb_simulationarchive* sa = reb_open_simulationarchive(path_to_sim);
    if (sa==NULL){
        printf("Can not open sim file.\nProgram Temrinating\n");
        exit(0);
    }
    // Get a simulation from the file (if possible, otherwise NULL is returned)
    struct reb_simulation* r_init = reb_create_simulation_from_simulationarchive(sa,-1);
    // Whenever you've opened a SimulationArchive and don't need it anymore, close it.
    reb_close_simulationarchive(sa);

    const double G = r_init->G;
    const double init_t = r_init->t;
    printf("INITIAL TIME: %f\n", init_t/sec2year);
    if (init_t < 1000000*sec2year){
	reb_remove_by_hash(r_init, reb_hash("Mercury"), 1);
    
    	struct reb_particle* floraptr = reb_get_particle_by_hash(r_init, reb_hash("Flora"));
    	double f_x = floraptr->x;    double f_y = floraptr->y;    double f_z = floraptr->z;
    	double f_vx = floraptr->vx;    double f_vy = floraptr->vy;    double f_vz = floraptr->vz;

    	reb_remove_by_hash(r_init, reb_hash("Flora"), 1);

    	// Initial Setup is complete. Now define new Sim.
    	struct reb_simulation* r = reb_create_simulation();
    	// Copy over particles from old sim
    	struct reb_particle* particles = r_init->particles;
    	const int N = r_init->N;
    	for (int i=0;i<N;i++){
	    const struct reb_particle p = particles[i];
	    reb_add(r, p);
    	}
    	reb_free_simulation(r_init); // delete connection to initFile

    	// Add 8 Florks with random ejection vector and specified size
    	for (int i=0;i<8;i++){
	    struct reb_particle p = {0};
	    p_size = p_sizes[i];
	    double V = (double)V0D0/p_size; // particles escape velcity scales with size
	    double evx = (double)rand() / (double)((unsigned)RAND_MAX + 1) - 0.5;
	    double evy = (double)rand() / (double)((unsigned)RAND_MAX + 1) - 0.5;
	    double evz = (double)rand() / (double)((unsigned)RAND_MAX + 1) - 0.5;
	    double norm = sqrt(evx*evx + evy*evy + evz*evz);
	    p.x = f_x;	p.y = f_y;	p.z = f_z;
	    p.vx = f_vx + V*evx/norm;	p.vy = f_vy + V*evy/norm;	p.vz = f_vz + V*evz/norm;
	    p.m = 0;
	    p.r = p_size/2*1e3;
	    reb_add(r, p);
	    double obl = cos((double)rand() / (double)((unsigned)RAND_MAX + 1)*M_PI);
	    obliquities[i] = obl;
	    printf("P%.0i: D=%.3f, OBL=%.3f, vx=%.5f, vy=%.5f, vz=%.5f\n", i+1, p_size, obl, V*evx/norm, V*evy/norm, V*evz/norm);
        }
	reb_simulationarchive_automate_interval(r,path_to_arch,10000*sec2year);
	r->force_is_velocity_dependent  = 1;
	r->integrator		    = REB_INTEGRATOR_WHFAST;
	r->additional_forces            = yarko_da;    // setup callback function for velocity dependent forces
	r->heartbeat                    = heartbeat;
    	r->N_active                     = 8;
    	r->dt			    = 3e6;
    	r->G			    = G;
    	reb_integrate(r, tmax);
    	reb_free_pointers(r);
    	exit(0);
    }else{
	printf("Long Simulation already exits! Continuing Simulation to final time!\n");
    	struct reb_simulation* r = reb_create_simulation();
    	// Copy over particles from old sim
    	struct reb_particle* particles = r_init->particles;
    	const int N = r_init->N;
    	for (int i=0;i<N;i++){
	    const struct reb_particle p = particles[i];
	    reb_add(r, p);
    	}
	reb_free_simulation(r_init);
	printf("START TIME: %f\n", init_t/sec2year);
	r->t = init_t;
	if (r->t >= tmax){
	    printf("Simulation has finished!\n");
	    reb_free_pointers(r);
	    exit(0);
        }else{
	    for (int i = 0; i < 8; i++){
		int scrap = (double)rand();
		scrap = (double)rand();
		scrap = (double)rand();
		double obl = ((double)rand() / (double)((unsigned)RAND_MAX + 1) - 0.5)*2;
		obliquities[i] = obl;
	    }
	    reb_simulationarchive_automate_interval(r,path_to_arch,10000*sec2year);
	    r->force_is_velocity_dependent  = 1;
	    r->integrator		    = REB_INTEGRATOR_WHFAST;
	    r->additional_forces            = yarko_da;    // setup callback function for velocity dependent forces
	    r->heartbeat                    = heartbeat;
	    r->N_active                     = 8;
	    r->dt			    = 3e6;
	    r->G			    = G;
	    reb_integrate(r, tmax);
	    reb_free_pointers(r);
	    exit(0);
	}
    }
}

void perfect_merge(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    const int N = r->N;
    const struct reb_particle star = particles[0];                // cache
    int j = 0;
#pragma omp parallel for
    for (int i=0;i<N;i++){
        const struct reb_particle p = particles[i];             // cache
        if (p.m!=0.) continue;                         // only test particles feel yarko

	// int          obliquity = -1;
	// if (i%2==0) {obliquity = 1;} // Alternating direction of drift
	double obliquity = obliquities[j];

        const double prx  = p.x-star.x;
        const double pry  = p.y-star.y;
        const double prz  = p.z-star.z;
        const double pr   = sqrt(prx*prx + pry*pry + prz*prz);         // distance relative to star
        
        const double prvx = p.vx-star.vx;
        const double prvy = p.vy-star.vy;
        const double prvz = p.vz-star.vz;
        const double v2   = prvx*prvx + prvy*prvy + prvz*prvz;

        const double gm         = r->G*star.m;
        const double energy     = 0.5*v2 - gm/pr;
        const double a          = -0.5*gm/energy;
        const double auPmyr2mPs = au2m/sec2year/1e6;
	const double D 		= (p.r*2)*1e-3; // Particle Diameter in km
        const double dadt       = 4e-4/D*auPmyr2mPs*obliquity; // drift rate in m/s
        const double k          = 0.5*dadt*gm/(a*a);

        particles[i].ax += k * prvx/v2;
        particles[i].ay += k * prvy/v2;
        particles[i].az += k * prvz/v2;
	j++;
    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 5e5*sec2year)){        // print some information to screen
        reb_output_timing(r, years * sec2year);
	printf("\n");
    }
}
