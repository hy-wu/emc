#include <vector>
#include <deque>
#include <cmath>
// #include <cassert>
#include <random>
#include <iostream>

const double Boltzmann = 1.380649e-23;

struct Particle {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

struct Cell {
    std::deque<Particle> particles;
};

struct V3 {
    double x;
    double y;
    double z;
};

struct V33 {
    double xx;
    double xy;
    double xz;
    double yx;
    double yy;
    double yz;
    double zx;
    double zy;
    double zz;
};

class Cubic {
public:
    Cubic(int n_particles, double cell_len, int cubic_size, int box_size, double r, double m, double T, bool equ=true, std::vector<double> vs={}) {
        /* unit: nm, kg, K, ... */
        k_B = Boltzmann * 1e9;
        this->n_particles = n_particles;
        this->cell_len = cell_len;
        this->box_size = box_size;
        this->cell_volume = pow(cell_len, 3);
        this->cubic_size = cubic_size;
        // assert(cubic_size % box_size == 0);
        this->cubic_box_size = cubic_size / box_size;
        this->cubic_len = cubic_size * cell_len;
        this->box_len = box_size * cell_len;
        this->r = r;
        this->m = m;
        this->T0 = T;
        this->t = 0;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0, this->cubic_len);
        std::normal_distribution<> disN(0, sqrt(this->k_B * T / m));
        cells.resize(cubic_size);
        for (int i = 0; i < cubic_size; i++) {
            cells[i].resize(cubic_size);
            for (int j = 0; j < cubic_size; j++) {
                cells[i][j].resize(cubic_size);
            }
        }
        for (int i = 0; i < n_particles; i++) {
            Particle particle;
            particle.x = dis(gen);
            particle.y = dis(gen);
            particle.z = dis(gen);
            if (equ) {
                particle.vx = disN(gen);
                particle.vy = disN(gen);
                particle.vz = disN(gen);
            } else {
                particle.vx = vs[0];
                particle.vy = vs[1];
                particle.vz = vs[2];
            }
            particles.push_back(particle);
            cells[particle.x / cell_len][particle.y / cell_len][particle.z / cell_len].particles.push_back(particle);
        }
        std::cout << "nσ^3 = " << n_particles / pow(cubic_size, 3) * pow((2 * r / cell_len), 3) << std::endl;
    }

    void run(int n_steps, double dt, bool quiet=false, bool z_loop=true) {
        // dt to big (dt > 1e-1) may cause error
        for (int step = 0; step < n_steps; step++) {
            t += dt;
            ts.push_back(t);
            // clear cells
            for (int i = 0; i < cubic_size; i++) {
                for (int j = 0; j < cubic_size; j++) {
                    for (int k = 0; k < cubic_size; k++) {
                        cells[i][j][k].particles.clear();
                    }
                }
            }
            // free movement
            for (int i = 0; i < n_particles; i++) {
                particles[i].x += particles[i].vx * dt;
                particles[i].y += particles[i].vy * dt;
                particles[i].z += particles[i].vz * dt;
                if (particles[i].x < 0) {
                    particles[i].x = -particles[i].x;
                    particles[i].vx = -particles[i].vx;
                } else if (particles[i].x >= cubic_len) {
                    particles[i].x = 2 * cubic_len - particles[i].x;
                    particles[i].vx = -particles[i].vx;
                }
                if (particles[i].y < 0) {
                    particles[i].y = -particles[i].y;
                    particles[i].vy = -particles[i].vy;
                } else if (particles[i].y >= cubic_len) {
                    particles[i].y = 2 * cubic_len - particles[i].y;
                    particles[i].vy = -particles[i].vy;
                }
                if (z_loop) {
                    if (particles[i].z < 0) {
                        particles[i].z += cubic_len;
                    } else if (particles[i].z >= cubic_len) {
                        particles[i].z -= cubic_len;
                    }
                } else {
                    if (particles[i].z < 0) {
                        particles[i].z = -particles[i].z;
                        particles[i].vz = -particles[i].vz;
                    } else if (particles[i].z >= cubic_len) {
                        particles[i].z = 2 * cubic_len - particles[i].z;
                        particles[i].vz = -particles[i].vz;
                    }
                }
            }
            // update cells
            for (int i = 0; i < n_particles; i++) {
                cells[particles[i].x / cell_len][particles[i].y / cell_len][particles[i].z / cell_len].particles.push_back(particles[i]);
            }
        }
    }

private:

public:
    double k_B;
    int n_particles;
    double cell_len;
    int box_size;
    double cell_volume;
    int cubic_size;
    int cubic_box_size;
    double cubic_len;
    double box_len;
    double r;
    double m;
    double T0;
    double t;
    std::vector<Particle> particles;
    std::vector<double> ts;
    std::vector<double> Ts;
    std::vector<double> ωs;
    std::vector<V3> us;
    std::vector<V33> pks;
    std::vector<V3> qks;
    std::vector<std::vector<std::vector<Cell>>> cells;
};

int mc() {
    Cubic cubic(1e6, 2, 150, 5, 1, 1e-26, 300);
    cubic.run(10, 1e-8);
    return 0;
}
