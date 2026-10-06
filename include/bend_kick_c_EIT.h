#ifndef BEND_KICK_UPDATER_H
#define BEND_KICK_UPDATER_H
#include <vector>
#include <complex>

//void bend_kick(double dt, double Bz, Ensemble ensemble, Force forces, int num_steps);

extern double _epsilon0;
extern double _k;
extern double rabi_eit;
extern double det1_eit;
extern double det2_eit;
extern double gamma1_eit;
extern double gamma2_eit;
extern std::vector<double> k1_eit;
extern std::vector<double> k2_eit;
extern double k1_eit_mag;
extern double k2_eit_mag;
extern double beam_w_eit;
extern std::complex<double> II;

void ca_bend_kick_update_vector(double time);

void run_trap_potential(double kz);

void run_doppler_only();

void run_trap_coulomb(double kz);

void update_dens_mat_half_step();

void eit_force_midpoint();

void update_dens_mat();

void update_dens_mat_OLD(int step);

void eit_force();

void trap_force(double kz);

void coulomb_force_fmm();

void damping_force();

void coulomb_force();

void rotating_wall_force();

static void coulomb_force_one_pair(int i, const double* r0, const double* r1, double kij);

static double distance(const double *r/*,double delta*/);

void run_no_forces();

void write_to_outfile();

void write_to_outfile_state();

void write_beam_params();

void write_EIT_params();

void reset_forces();

void doppler_force();

void Radiation_Pressure(int beam, int ion, double dt);

double gaussian_intensity(int beam, int ion);

double detuning(int beam, int ion);

double compute_nbar(double dt, double inten, double det);

double scattering_rate(double inten, double det);

void add_radiation_pressure(int beam, int ion, struct CARandCtx* ctx, double nbar);

static void add_radiation_pressure_one(int beam, int ion, struct CARandCtx* ctx, double hbar_k_nrm, double nbar);

static void add_radiation_pressure_small_n(int beam, int ion,
        struct CARandCtx* ctx,
        double hbar_k_nrm,
        int n);


#endif
