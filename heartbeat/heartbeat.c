#include "rebound.h"

double min_distance_from_sun_squared = 0;
double max_distance_from_sun_squared = 0;

struct hb_event {
    uint32_t hash;
    double time;
    unsigned int new;
};

struct hb_event hb_escapes[500];
struct hb_event hb_sun_collisions[500];

int hb_escape_index = 0;
int hb_sun_collision_index = 0;

int needs_synchronize = 0

void heartbeat(struct reb_simulation *sim) {
    if ((sim->steps_done % 100) == 0) {
        const struct reb_particle *const particles = sim->particles;
        int N = sim->N - sim->N_var;
        for (int i = 1; i < N; i++) { // skip sun
            struct reb_particle p = particles[i];
            double distance_squared = p.x * p.x + p.y * p.y + p.z * p.z;
            if (distance_squared > max_distance_from_sun_squared) {
                printf("remove %u at t=%f (max)\n", p.hash, sim->t);
                reb_remove_by_hash(sim, p.hash, 1);
                hb_escapes[hb_escape_index].hash = p.hash;
                hb_escapes[hb_escape_index].time = sim->t;
                hb_escapes[hb_escape_index].new = 1;

                hb_escape_index++;
                needs_synchronize = 1
            } else if (distance_squared < min_distance_from_sun_squared) {
                printf("remove %u at t=%f (min)\n", p.hash, sim->t);
                reb_remove_by_hash(sim, p.hash, 1);
                hb_sun_collisions[hb_sun_collision_index].hash = p.hash;
                hb_sun_collisions[hb_sun_collision_index].time = sim->t;
                hb_sun_collisions[hb_sun_collision_index].new = 1;

                hb_sun_collision_index++;
                needs_synchronize = 1
            }
            if (needs_synchronize) {
                N--;
                reb_move_to_com(sim);
                reb_integrator_synchronize(sim);
                sim->ri_mercurius.recalculate_coordinates_this_timestep = 1;
                sim->ri_mercurius.recalculate_dcrit_this_timestep = 1;
            }
        }
    }
}
