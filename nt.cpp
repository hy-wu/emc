#include <vector>
#include <deque>
#include <cmath>
// #include <cassert>
#include <indicators/progress_bar.hpp>
#include <random>
#include <iostream>
#include <fstream>
// #include <torch/torch.h>

constexpr double Boltzmann = 1.380649e-23;

struct Particle {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

struct Cell {
    std::deque<Particle> ptcs;
};

struct Box {
    std::deque<Particle> ptcs;
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
    double k_B;
    int n_ptcs;
    double cell_l;
    int box_size;
    double cell_volume;
    int cubic_size;
    int cubic_box_size;
    double cubic_len;
    double box_len;
    double r;
    double m;
    double T0;
    double T;
    double t;
    std::vector<Particle> ptcs;
    std::vector<double> ts;
    std::vector<double> Ts;
    std::vector<double> ωs;
    std::vector<bool> accepts;
    // std::vector<V3> us;
    std::vector<std::vector<std::vector<std::vector<V3>>>> us;
    // std::vector<std::vector<std::vector<std::vector<torch::Tensor>>>> us_tensor;
    // std::vector<V33> pks;
    std::vector<std::vector<std::vector<std::vector<V33>>>> pks;
    // std::vector<std::vector<std::vector<std::vector<torch::Tensor>>>> pks_tensor;
    // std::vector<V3> qks;
    std::vector<std::vector<std::vector<std::vector<V3>>>> qks;
    // std::vector<std::vector<std::vector<std::vector<torch::Tensor>>>> qks_tensor;
    std::vector<std::vector<std::vector<Cell>>> cells;
    std::vector<std::vector<std::vector<Box>>> boxes;

    Cubic(const int n_ptcs, const double cell_l, const int cubic_size, const int box_size, const double r, const double m, const double T, const bool equ=true, const std::vector<double>& vs={}) {
        /* unit: nm, kg, K, ... */
        k_B = Boltzmann * 1e9;
        this->n_ptcs = n_ptcs;
        this->cell_l = cell_l;
        this->box_size = box_size;
        this->cell_volume = pow(cell_l, 3);
        this->cubic_size = cubic_size;
        // assert(cubic_size % box_size == 0);
        this->cubic_box_size = cubic_size / box_size;
        this->cubic_len = cubic_size * cell_l;
        this->box_len = box_size * cell_l;
        this->r = r;
        this->m = m;
        this->T0 = T;
        this->T = T;
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
        for (int i = 0; i < n_ptcs; i++) {
            Particle particle{};
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
            ptcs.push_back(particle);
            cells[static_cast<int>(particle.x / cell_l)][static_cast<int>(particle.y / cell_l)][static_cast<int>(particle.z / cell_l)].ptcs.push_back(particle);
            // cells[particle.x / cell_l][particle.y / cell_l][particle.z / cell_l].ptcs.push_back(particle);
        }
        std::cout << "nσ^3 = " << n_ptcs / pow(cubic_size, 3) * pow((2 * r / cell_l), 3) << std::endl;
    }

    static double bound(const double x, const double cubic_len) {
        if (x < 0) {
            return -x;
        } if (x >= cubic_len) {
            return 2 * cubic_len - x;
        } {
            return x;
        }
    }

    static double cycle_bound(const double x, const double cubic_len) {
        if (x < 0) {
            return x + cubic_len;
        } if (x >= cubic_len) {
            return x - cubic_len;
        } {
            return x;
        }
    }

    static std::tuple<double, double, double> bound(double x, double y, double z, const double l) {
        if (x < 0) {
            x = -x;
        } else if (x >= l) {
            x = 2 * l - x;
        }
        if (y < 0) {
            y = -y;
        } else if (y >= l) {
            y = 2 * l - y;
        }
        if (z < 0) {
            z = -z;
        } else if (z >= l) {
            z = 2 * l - z;
        }
        return std::make_tuple(x, y, z);
    }

    [[nodiscard]] std::tuple<double, double, double> z_cycle(double x, double y, double z) const {
        /* periodic boundary condition in z direction, elastic boundary condition in x and y direction
        */
        if (x < 0) {
            x = -x;
        } else if (x >= cubic_len) {
            x = 2 * cubic_len - x;
        }
        if (y < 0) {
            y = -y;
        } else if (y >= cubic_len) {
            y = 2 * cubic_len - y;
        }
        if (z < 0) {
            z += cubic_len;
        } else if (z >= cubic_len) {
            z -= cubic_len;
        }
        return std::make_tuple(x, y, z);
    }

    static V33 outer_add(const V33& v33, const V3& v3) {
        return {v33.xx + v3.x * v3.x, v33.xy + v3.x * v3.y, v33.xz + v3.x * v3.z,
                v33.yx + v3.y * v3.x, v33.yy + v3.y * v3.y, v33.yz + v3.y * v3.z,
                v33.zx + v3.z * v3.x, v33.zy + v3.z * v3.y, v33.zz + v3.z * v3.z};
    }

    void run(int n_steps, double dt, bool quiet=false, bool z_loop=true, int sample_step=10) {
        // dt to big (dt > 1e-1) may cause error
        // pks_tensor.resize(n_steps);
        // qks_tensor.resize(n_steps);
        // us_tensor.resize(n_steps);
        us.resize(n_steps);
        pks.resize(n_steps);
        qks.resize(n_steps);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dis(0, 1);
        cells.resize(cubic_size);
        for (int i = 0; i < cubic_size; i++) {
            cells[i].resize(cubic_size);
            for (int j = 0; j < cubic_size; j++) {
                cells[i][j].resize(cubic_size);
            }
        }
        boxes.resize(cubic_box_size);
        for (int i = 0; i < cubic_box_size; i++) {
            boxes[i].resize(cubic_box_size);
            for (int j = 0; j < cubic_box_size; j++) {
                boxes[i][j].resize(cubic_box_size);
            }
        }
        indicators::ProgressBar bar{
            indicators::option::PrefixText{"Simulation"},
            indicators::option::ShowElapsedTime{true},
            indicators::option::ShowRemainingTime{true},
            indicators::option::ShowPercentage{true},
        };
        double xJ, yJ, zJ, nl3, dv2;
        for (int step = 0; step < n_steps; step++) {
            std::vector<std::vector<double>> σhats(n_ptcs, std::vector<double>(3));
            std::vector<std::vector<double>> σs(n_ptcs, std::vector<double>(3));
            std::vector<double> accepteds(n_ptcs);
            // generate random directions
            for (int i = 0; i < n_ptcs; i++) {
                double θ = acos(2 * dis(gen) - 1);
                double φ = 2 * M_PI * dis(gen);
                σhats[i][0] = sin(θ) * cos(φ);
                σhats[i][1] = sin(θ) * sin(φ);
                σhats[i][2] = cos(θ);
                σs[i][0] = σhats[i][0] * r * 2;
                σs[i][1] = σhats[i][1] * r * 2;
                σs[i][2] = σhats[i][2] * r * 2;
                accepteds[i] = dis(gen);
            }
            if (!quiet) {
                bar.set_option(indicators::option::PostfixText{std::to_string(step + 1) + "/" + std::to_string(n_steps)});
                bar.set_progress((step + 1) * 100 / n_steps);
            }

            t += dt;
            ts.push_back(t);
            // clear cells
            for (int i = 0; i < cubic_size; i++) {
                for (int j = 0; j < cubic_size; j++) {
                    for (int k = 0; k < cubic_size; k++) {
                        cells[i][j][k].ptcs.clear();
                    }
                }
            }

            /* free movement */
            for (int i = 0; i < n_ptcs; i++) {
                ptcs[i].x += ptcs[i].vx * dt;
                ptcs[i].y += ptcs[i].vy * dt;
                ptcs[i].z += ptcs[i].vz * dt;
                if (ptcs[i].x < 0) {
                    ptcs[i].x = -ptcs[i].x;
                    ptcs[i].vx = -ptcs[i].vx;
                } else if (ptcs[i].x >= cubic_len) {
                    ptcs[i].x = 2 * cubic_len - ptcs[i].x;
                    ptcs[i].vx = -ptcs[i].vx;
                }
                if (ptcs[i].y < 0) {
                    ptcs[i].y = -ptcs[i].y;
                    ptcs[i].vy = -ptcs[i].vy;
                } else if (ptcs[i].y >= cubic_len) {
                    ptcs[i].y = 2 * cubic_len - ptcs[i].y;
                    ptcs[i].vy = -ptcs[i].vy;
                }
                if (z_loop) {
                    if (ptcs[i].z < 0) {
                        ptcs[i].z += cubic_len;
                    } else if (ptcs[i].z >= cubic_len) {
                        ptcs[i].z -= cubic_len;
                    }
                } else {
                    if (ptcs[i].z < 0) {
                        ptcs[i].z = -ptcs[i].z;
                        ptcs[i].vz = -ptcs[i].vz;
                    } else if (ptcs[i].z >= cubic_len) {
                        ptcs[i].z = 2 * cubic_len - ptcs[i].z;
                        ptcs[i].vz = -ptcs[i].vz;
                    }
                }
                cells[static_cast<int>(ptcs[i].x / cell_l)][static_cast<int>(ptcs[i].y / cell_l)][static_cast<int>(ptcs[i].z / cell_l)].ptcs.push_back(ptcs[i]);
            }
            std::vector<std::vector<std::vector<Cell>>> new_cells(cubic_size, std::vector<std::vector<Cell>>(cubic_size, std::vector<Cell>(cubic_size)));

            /* collision */
            for (int i = 0; i < n_ptcs; i++) {
                std::tie(xJ, yJ, zJ) = z_cycle(ptcs[i].x + σs[i][0], ptcs[i].y + σs[i][1], ptcs[i].z + σs[i][2]);
                Cell cell_j = cells[static_cast<int>(xJ / cell_l)][static_cast<int>(yJ / cell_l)][static_cast<int>(zJ / cell_l)];
                int len_j = static_cast<int>(cell_j.ptcs.size());
                if (len_j > 0) {
                    // int j = rand() % len_j;
                    std::uniform_int_distribution<> distribution(0, len_j - 1);
                    int j = distribution(gen);
                    double g_ij[3] = {cell_j.ptcs[j].vx - ptcs[i].vx, cell_j.ptcs[j].vy - ptcs[i].vy, cell_j.ptcs[j].vz - ptcs[i].vz};
                    double σhat_i_dot_g_ij = σhats[i][0] * g_ij[0] + σhats[i][1] * g_ij[1] + σhats[i][2] * g_ij[2];
                    if (σhat_i_dot_g_ij > 0) {
                        double ω_ij = 16 * M_PI * pow(r, 2) * σhat_i_dot_g_ij * len_j / cell_volume * dt;
                        ωs.push_back(ω_ij);
                        if (accepteds[i] < ω_ij) {
                            ptcs[i].vx -= σs[i][0], ptcs[i].vy -= σs[i][1], ptcs[i].vz -= σs[i][2];
                        }
                    }
                }
                // std::tie(ptcs[i].x, ptcs[i].y, ptcs[i].z) = z_cycle(ptcs[i].x, ptcs[i].y, ptcs[i].z);
                new_cells[static_cast<int>(ptcs[i].x / cell_l)][static_cast<int>(ptcs[i].y / cell_l)][static_cast<int>(ptcs[i].z / cell_l)].ptcs.push_back(ptcs[i]);
                // new_cells[x / cell_l][y / cell_l][z / cell_l].ptcs.push_back(ptcs[i]);
                // new_cells.at(x / cell_l).at(y / cell_l).at(z / cell_l).ptcs.push_back(ptcs[i]);
            }
            cells = new_cells;

            /* temperature */
            T = 0;
            for (int i = 0; i < n_ptcs; i++) {
                T += m * (pow(ptcs[i].vx, 2) + pow(ptcs[i].vy, 2) + pow(ptcs[i].vz, 2)) / (3 * k_B * n_ptcs);
            }
            Ts.push_back(T);

            /* velocity, pressure and heat flux */
            if (step % sample_step != 0) {
                continue;
            }
            for (int i = 0; i < cubic_box_size; i++) {
                for (int j = 0; j < cubic_box_size; j++) {
                    for (int k = 0; k < cubic_box_size; k++) {
                        boxes[i][j][k].ptcs.clear();
                    }
                }
            }
            for (int i = 0; i < n_ptcs; i++) {
                boxes[static_cast<int>(ptcs[i].x / box_len)][static_cast<int>(ptcs[i].y / box_len)][static_cast<int>(ptcs[i].z / box_len)].ptcs.push_back(ptcs[i]);
            }
            V3 u_mean = {0, 0, 0};
            V33 pk_mean = {0, 0, 0, 0, 0, 0, 0, 0, 0};
            V3 qk_mean = {0, 0, 0};
            V3 dv = {0, 0, 0};
            us[step].resize(cubic_box_size);
            pks[step].resize(cubic_box_size);
            qks[step].resize(cubic_box_size);
            for (int i = 0; i < cubic_box_size; i++) {
                us[step][i].resize(cubic_box_size);
                pks[step][i].resize(cubic_box_size);
                qks[step][i].resize(cubic_box_size);
                for (int j = 0; j < cubic_box_size; j++) {
                    us[step][i][j].resize(cubic_box_size);
                    pks[step][i][j].resize(cubic_box_size);
                    qks[step][i][j].resize(cubic_box_size);
                    for (int k = 0; k < cubic_box_size; k++) {
                        int len = static_cast<int>(boxes[i][j][k].ptcs.size());
                        V3 u = {0, 0, 0};
                        dv = {0, 0, 0};
                        if (len > 0) {
                            for (int l = 0; l < len; l++) {
                                u.x += boxes[i][j][k].ptcs[l].vx;
                                u.y += boxes[i][j][k].ptcs[l].vy;
                                u.z += boxes[i][j][k].ptcs[l].vz;
                            }
                            u.x /= len;
                            u.y /= len;
                            u.z /= len;
                            pk_mean = {0, 0, 0, 0, 0, 0, 0, 0, 0};
                            qk_mean = {0, 0, 0};
                            for (int l = 0; l < len; l++) {
                                dv.x = boxes[i][j][k].ptcs[l].vx - u.x;
                                dv.y = boxes[i][j][k].ptcs[l].vy - u.y;
                                dv.z = boxes[i][j][k].ptcs[l].vz - u.z;
                                pk_mean = outer_add(pk_mean, dv);
                                dv2 = pow(dv.x, 2) + pow(dv.y, 2) + pow(dv.z, 2);
                                qk_mean = {qk_mean.x + dv2 * dv.x, qk_mean.y + dv2 * dv.y, qk_mean.z + dv2 * dv.z};
                            }
                            nl3 = len * pow(cubic_box_size, 3);
                            pk_mean = {pk_mean.xx / nl3, pk_mean.xy / nl3, pk_mean.xz / nl3,
                                       pk_mean.yx / nl3, pk_mean.yy / nl3, pk_mean.yz / nl3,
                                       pk_mean.zx / nl3, pk_mean.zy / nl3, pk_mean.zz / nl3};
                            qk_mean = {qk_mean.x / nl3, qk_mean.y / nl3, qk_mean.z / nl3};
                            us[step][i][j][k] = u;
                            pks[step][i][j][k] = pk_mean;
                            qks[step][i][j][k] = qk_mean;
                        }
                    }
                }
            }
        }
    }

private:

};

Cubic nt() {  // ROOT entrance
    int num_steps;
    int sampling_step;
    std::cout << "num_steps: ";
    std::cin >> num_steps;
    std::cout << "sampling_step: ";
    std::cin >> sampling_step;
    Cubic cubic(1e5, 2, 50, 5, 1, 1e-26, 300);
    cubic.run(num_steps, 1e-8, false, true, sampling_step);
    std::ofstream file("data.csv");

    for (int step = 0; step < num_steps ; step += sampling_step) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    file << cubic.ts[step] << "," << i << "," << j << "," << k;
                    file << "," << cubic.us[step][i][j][k].x << "," << cubic.us[step][i][j][k].y << "," << cubic.us[step][i][j][k].z;
                    file << "," << cubic.pks[step][i][j][k].xx << "," << cubic.pks[step][i][j][k].xy << "," << cubic.pks[step][i][j][k].xz;
                    file << "," << cubic.pks[step][i][j][k].yx << "," << cubic.pks[step][i][j][k].yy << "," << cubic.pks[step][i][j][k].yz;
                    file << "," << cubic.pks[step][i][j][k].zx << "," << cubic.pks[step][i][j][k].zy << "," << cubic.pks[step][i][j][k].zz;
                    file << "," << cubic.qks[step][i][j][k].x << "," << cubic.qks[step][i][j][k].y << "," << cubic.qks[step][i][j][k].z;
                    file << "\n";
                }
            }
        }
    }
    file.close();
    file.open("data_Ts.csv");
    for (int i = 0; i < num_steps; ++i) {
        file << cubic.ts[i] << "," << std::setprecision(15) << cubic.Ts[i] << "\n";
    }
    file.close();
    return cubic;
}

int main() {  // g++ entrance
    int num_steps;
    constexpr int sampling_step = 10;
    std::cout << "num_steps: ";
    std::cin >> num_steps;
    Cubic cubic(1e5, 2, 50, 5, 1, 1e-26, 300);
    cubic.run(num_steps, 1e-8, false, true, sampling_step);
    std::ofstream file("data.csv");

    for (int step = 0; step < num_steps ; step += sampling_step) {
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                for (int k = 0; k < 10; ++k) {
                    file << step << "," << i << "," << j << "," << k;
                    file << "," << cubic.us[step][i][j][k].x << "," << cubic.us[step][i][j][k].y << "," << cubic.us[step][i][j][k].z;
                    file << "," << cubic.pks[step][i][j][k].xx << "," << cubic.pks[step][i][j][k].xy << "," << cubic.pks[step][i][j][k].xz;
                    file << "," << cubic.pks[step][i][j][k].yx << "," << cubic.pks[step][i][j][k].yy << "," << cubic.pks[step][i][j][k].yz;
                    file << "," << cubic.pks[step][i][j][k].zx << "," << cubic.pks[step][i][j][k].zy << "," << cubic.pks[step][i][j][k].zz;
                    file << "," << cubic.qks[step][i][j][k].x << "," << cubic.qks[step][i][j][k].y << "," << cubic.qks[step][i][j][k].z;
                    file << "\n";
                }
            }
        }
    }
    file.close();
    file.open("data_Ts.csv");
    for (int i = 0; i < num_steps; ++i) {
        file << cubic.ts[i] << "," << cubic.Ts[i] << "\n";
    }
    file.close();
    return 0;
}
