/**
 * Repoduction of Emsenhuber et al. 2020 C2019 simulations
 *
 * Given an initial simulation file with ~155 planetessimals
 * Program will simulate 400 Myrs of evolution
 * Collisions will be treated by the Cambioni et al. 2019 giant collision resolver
 * 
 * Units in AU, Msun, Days
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "collresolve.h"


void heartbeat(struct reb_simulation* r);
int C2019(struct reb_simulation* const r, struct reb_collision c);

double year_to_sec = 31557600.0;
double day_to_sec = 518400.0;
double year_to_day = 365.25;
double rad2Deg = 57.29577951308232;
double years = 5e2;
double au_to_m = 149597870700;
double m_earth = 3.003e-6; // Earth mass in Solar Units
char* path_to_arch;

int main(int argc, char* argv[]){
    if(argc != 3){
        printf("Must include 'path_to_sim' and 'path_to_arch'  as command line inputs!\nProgram Terminated\n");
        exit(0);
    }

    double tmax = years * year_to_day;

    char* path_to_sim = argv[1];
    path_to_arch = argv[2];
	remove(path_to_arch);

    struct reb_simulationarchive* sa = reb_open_simulationarchive(path_to_sim);
    if (sa==NULL){
        printf("Can not open sim file.\nProgram Temrinating\n");
        exit(0);
    }
    // Get a simulation from the file (if possible, otherwise NULL is returned)
    struct reb_simulation* r_init = reb_create_simulation_from_simulationarchive(sa,-1);
    // Whenever you've opened a SimulationArchive and don't need it anymore, close it.
    reb_close_simulationarchive(sa);

    const double G = r_init->G; // Should be in AU, Msun, Days (G ~ 4Pi^2*365.25)
    const double init_t = r_init->t;
    printf("INITIAL TIME: %f\n", init_t/year_to_day);

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

	// reb_simulationarchive_automate_interval(r,path_to_arch,10000*year_to_day);
	reb_simulationarchive_automate_interval(r,path_to_arch,50*year_to_day);
	// reb_simulationarchive_automate_step(r,path_to_arch,6088);
	r->integrator		    		= REB_INTEGRATOR_MERCURIUS;
	r->ri_mercurius.safe_mode		= 0;
	r->ri_whfast.corrector			= 11;
	// r->integrator		    		= REB_INTEGRATOR_IAS15;
	r->collision					= REB_COLLISION_DIRECT;
	r->collision_resolve			= C2019;
	r->heartbeat                    = heartbeat;
	r->N_active                     = r->N;
	r->G			    			= G;
	r->dt			    			= 1;
	reb_integrate(r, tmax);
	reb_free_pointers(r);
	exit(0);
}

int C2019(struct reb_simulation* const r, struct reb_collision c){
    printf("COLLISION TIME: %f yrs\n", r->t/year_to_day);
	struct reb_particle p1 = r->particles[c.p1];
	struct reb_particle p2 = r->particles[c.p2];
	printf("pid %d: m=%e, r=%e, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e\n",p1.hash,p1.m,p1.r,p1.x,p1.y,p1.z,p1.vx,p1.vy,p1.vz);
	printf("pid %d: m=%e, r=%e, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e\n",p2.hash,p2.m,p2.r,p2.x,p2.y,p2.z,p2.vx,p2.vy,p2.vz);

	struct collresolve_body big;
	struct collresolve_body small;
	int swap = 0;
	if (p1.m >= p2.m){
		big.mass = p1.m; big.radius = p1.r; big.pos[0] = p1.x; big.pos[1] = p1.y; big.pos[2] = p1.z; big.vel[0] = p1.vx; big.vel[1] = p1.vy; big.vel[2] = p1.vz;
		small.mass = p2.m; small.radius = p2.r; small.pos[0] = p2.x; small.pos[1] = p2.y; small.pos[2] = p2.z; small.vel[0] = p2.vx; small.vel[1] = p2.vy; small.vel[2] = p2.vz;
	}else{
		big.mass = p2.m; big.radius = p2.r; big.pos[0] = p2.x; big.pos[1] = p2.y; big.pos[2] = p2.z; big.vel[0] = p2.vx; big.vel[1] = p2.vy; big.vel[2] = p2.vz;
		small.mass = p1.m; small.radius = p1.r; small.pos[0] = p1.x; small.pos[1] = p1.y; small.pos[2] = p1.z; small.vel[0] = p1.vx; small.vel[1] = p1.vy; small.vel[2] = p1.vz;
		swap = 1;
	}

	double E0 = reb_tools_energy(r);
	struct reb_vec3d ang_vec0 = reb_tools_angular_momentum(r);

	struct collresolve_conf* conf = collresolve_conf_new();
	collresolve_conf_unit_merc(conf);
	collresolve_conf_model(conf, COLLRESOLVE_MODEL_C2019);
	collresolve_conf_sep_after(conf, 1.2);

	struct collresolve_body ret[3];
	collresolve_resolve(conf, big, small, 2, ret);
	collresolve_conf_free(conf);

	big = ret[0];
	small = ret[1];
	// struct collresolve_body debris = ret[2];
	if (small.mass > 1e-3*m_earth){
		if (p1.m >= p2.m){
			p1.m = big.mass; p1.r = big.radius; p1.x = big.pos[0]; p1.y = big.pos[1]; p1.z = big.pos[2]; p1.vx = big.vel[0]; p1.vy = big.vel[1]; p1.vz = big.vel[2];
			p2.m = small.mass; p2.r = small.radius; p2.x = small.pos[0]; p2.y = small.pos[1]; p2.z = small.pos[2]; p2.vx = small.vel[0]; p2.vy = small.vel[1]; p2.vz = small.vel[2];
		}else{
			p2.m = big.mass; p2.r = big.radius; p2.x = big.pos[0]; p2.y = big.pos[1]; p2.z = big.pos[2]; p2.vx = big.vel[0]; p2.vy = big.vel[1]; p2.vz = big.vel[2];
			p1.m = small.mass; p1.r = small.radius; p1.x = small.pos[0]; p1.y = small.pos[1]; p1.z = small.pos[2]; p1.vx = small.vel[0]; p1.vy = small.vel[1]; p1.vz = small.vel[2];
		}
	}else{
		if (p1.m >= p2.m){
			p1.m = big.mass; p1.r = big.radius; p1.x = big.pos[0]; p1.y = big.pos[1]; p1.z = big.pos[2]; p1.vx = big.vel[0]; p1.vy = big.vel[1]; p1.vz = big.vel[2];
			p2.m = 0;
		}else{
			p2.m = big.mass; p2.r = big.radius; p2.x = big.pos[0]; p2.y = big.pos[1]; p2.z = big.pos[2]; p2.vx = big.vel[0]; p2.vy = big.vel[1]; p2.vz = big.vel[2];
			p1.m = 0;
		}
	}
	r->particles[c.p1] = p1; // Update particle information in sim
	r->particles[c.p2] = p2; // Update particle information in sim

	double E1 = reb_tools_energy(r);
	struct reb_vec3d ang_vec1 = reb_tools_angular_momentum(r);
	
	printf("Debris Energy Lost: %e\n", E0-E1);
	printf("Angular Momentum Change: h0 - [%e, %e, %e] h1 - [%e, %e, %e]\n", ang_vec0.x, ang_vec0.y, ang_vec0.z, ang_vec1.x, ang_vec1.y, ang_vec1.z);
	
	r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
	r->ri_mercurius.recalculate_dcrit_this_timestep = 1;

	if (small.mass > 1e-3*m_earth){
		printf("pid %d: m=%e, r=%e, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e\n",p1.hash,p1.m,p1.r,p1.x,p1.y,p1.z,p1.vx,p1.vy,p1.vz);
		printf("pid %d: m=%e, r=%e, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e\n",p2.hash,p2.m,p2.r,p2.x,p2.y,p2.z,p2.vx,p2.vy,p2.vz);
		printf("COLLISION RESOLVED - Particles Remaining: %d\n", r->N);
		return 0;
	}else{
		if (!swap){
			printf("%d merged %d: m=%e, r=%e, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e\n",p1.hash,p2.hash,p1.m,p1.r,p1.x,p1.y,p1.z,p1.vx,p1.vy,p1.vz);
			printf("COLLISION RESOLVED - Particles Remaining: %d\n", r->N-1);
			return 2;
		}else{
			printf("%d merged %d: m=%e, r=%e, x=%e, y=%e, z=%e, vx=%e, vy=%e, vz=%e\n",p2.hash,p1.hash,p2.m,p2.r,p2.x,p2.y,p2.z,p2.vx,p2.vy,p2.vz);
			printf("COLLISION RESOLVED - Particles Remaining: %d\n", r->N-1);
			return 1;
		}
	}
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 50*year_to_day)){        // print some information to screen
        reb_output_timing(r, years * year_to_day);
		printf("\n");
    }
}
