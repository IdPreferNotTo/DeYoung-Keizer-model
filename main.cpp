#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <pwd.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    // -------------------------- Noise parameters & random number generator -------------------------------------------

    /*
    std::string arg1 = argv[1];
    std::string arg2 = argv[2];
    try {
        std::size_t pos;
        int n = std::stoi(arg1, &pos);
        if (pos < arg1.size()) {
            std::cerr << "Trailing characters after number: " << arg1 << '\n';
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg1 << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg1 << '\n';
    }

    try {
        std::size_t pos;
        float ip3 = std::stof(arg2, &pos);
        if (pos < arg2.size()) {
            std::cerr << "Trailing characters after number: " << arg2 << '\n';
        }
    } catch (std::invalid_argument const &ex) {
        std::cerr << "Invalid number: " << arg2 << '\n';
    } catch (std::out_of_range const &ex) {
        std::cerr << "Number out of range: " << arg2 << '\n';
    }
    */

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> uniform_dist(0, 1);


    const int N_total = std::stoi(argv[1]);
    const double v1 = 6.;
    const double v2 = 0.11;
    const double v3 = 0.9;

    const double c0 = 2.;
    const double c1 = 0.185;

    const double k3 = 0.1;

    const double a1 = 400;    // binding of any IP3+
    const double a2 = 0.2;    // binding of inhibitory Ca2+
    const double a5 = 20;     // binding of activation Ca2+

    const double b1 = 52.;    // removal of IP3 if inhibitory Ca2+ is not bound
    const double b2 = 0.2098; // removal of inhibitory Ca2+ if IP3 is bound
    const double b3 = 377.2;  // removal of IP3 if inhibitory Ca2+ is bound
    const double b4 = 0.0289; // removal of inhibitory Ca2+ if IP3 is not bound
    const double b5 = 1.6468; // removal of activation Ca2+

    double ip3 = std::stof(argv[2]);

    const double dt = 0.00001;
    double t = 0;
    int out_count = 0;
    double xi1;
    double xi2;

    const std::vector<bool> channel_open = {true, true, false};
    std::vector<bool> channel_sub_state = {false, false, false};
    std::vector<std::vector<bool>> channel_state (3, channel_sub_state);
    std::vector<std::vector<std::vector<bool>>> channel_states(N_total, channel_state);

    // ------------------------ Parameters for output file -------------------------------------------------------------
    std::string path = "../out/";
    char parameters[100];
    snprintf(parameters, sizeof(parameters), "_N%d_ip%.2f", N_total, ip3);

    std::string out_file;
    out_file = path + "dyk_model_ca_waves" + parameters;
    std::ofstream file;
    file.open(out_file);
    if (!file.is_open()) {
        std::cout << "Could not open file at: " << out_file << std::endl;
        struct passwd *pw = getpwuid(getuid());
        const char *homedir = pw->pw_dir;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }

    file << "# index puff_time puff_amplitude spike_triggered" << std::endl;
    //------------------------------------------------------------------------------------------------------------------

    double ca = 0.0;
    double j1;
    double j2;
    double j3;

    double p_ip3; /// probability to change the state of the IP3 binding site (occupied/unoccupied)
    double p_caa; // probability to change the state of the activation Ca2+ binding site (occupied/unoccupied)
    double p_cai; // probability to change the state of the inhibition Ca2+ binding site (occupied/unoccupied)
    double p; // probability to change the channel state

    while(t < 600){
        float N_open_channel = 0;
        for(auto &sub_unit: channel_states){
            // For every channel (with state x) in a cluster
            int N_open_subunits = 0;
            for(auto &x : sub_unit){
                // For every subunit (with state x) in a channel

                bool ip3_is_bound = x[0];
                bool caa_is_bound = x[1];
                bool cai_is_bound = x[2];

                if(!ip3_is_bound){
                    // Binding probability of Ip3
                    p_ip3 = a1*ip3;
                }
                if(!caa_is_bound){
                    // Binding probability of activation Ca2+
                    p_caa = a5*ca;
                }
                if(!cai_is_bound){
                    // Binding probability of inhibition Ca2+
                    p_cai = a2*ca;
                }

                // Unbinding probabilities
                if(ip3_is_bound){
                    if(!cai_is_bound){
                        // Unbinding probability of Ip3 if inhibition Ca2+ is unbound
                        p_ip3 = b1;
                    }
                    else{
                        // Unbinding probability of Ip3 if inhibition Ca2+ is bound
                        p_ip3 = b3;
                    }
                }
                if(caa_is_bound) {
                    // Unbinding probability of activation Ca2+
                    p_caa = b5;
                }
                if(cai_is_bound) {
                    if (!ip3_is_bound) {
                        // Unbinding probability of inhibition Ca2+ if Ip3 is unbound
                        p_cai = b4;
                    } else {
                        // Unbinding probability of inhibition Ca2+ if Ip3 is bound
                        p_cai = b2;
                    }
                }

                p = p_ip3 + p_caa + p_cai;

                xi1 = uniform_dist(generator);
                if (xi1 < p*dt){
                    xi2 = uniform_dist(generator);
                    if(xi2 < p_ip3/p){
                        x[0] = !x[0];
                    }
                    else if(xi2 < (p_caa+p_ip3)/p){
                        x[1] = !x[1];
                    }
                    else{
                        x[2] = !x[2];
                    }
                }
                if(x == channel_open){
                    N_open_subunits += 1;
                }
            }
            if(N_open_subunits >= 3){
                N_open_channel += 1;
            }
        }

        j1 = v1*N_open_channel/N_total*((c0 - ca*(1.+c1)));
        j2 = v2*(c0 - ca*(1.+c1));
        j3 = v3*pow(ca, 2)/(pow(k3, 2) + pow(ca, 2));

        ca += (j1 + j2 - j3)*dt;
        t += dt;
        out_count += 1;
        if(out_count==1000){
            out_count = 0;
            file << t << " " << ca << " " << j1 << " " << j2 << " " << j3 << " " << N_open_channel/N_total <<"\n";
        }
    }

    return 0;
}
