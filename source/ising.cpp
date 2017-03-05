#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>

#define UP 0
#define RIGHT 1
#define LEFT 2
#define DOWN 3

#define L 16
#define SIZE L*L

#define DATA 9

#define MAG 0
#define MAG2 1
#define MAG4 2
#define MAGERR 3
#define SUSERR 4
#define ENE 5
#define ENE2 6
#define ENE4 7
#define ENERR 8
#define CHERR 9

using namespace std;


void initialize(bool spins[SIZE], mt19937& gen, uniform_int_distribution<int>& brandom);
void get_neighbors(int neighs[SIZE][4]);


void do_step(bool spins[SIZE],  int neighs[SIZE][4], double tstar, int N, double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos, double m[DATA]);
void do_step_wolff(bool spins[SIZE],  int neighs[SIZE][4], double tstar, int N, double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos, double m[DATA]);
double magnetization(bool spins[SIZE]);

void flip_spin(bool spins[SIZE], int neighs[SIZE][4], double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos);
void add_to_cluster(bool spins[SIZE], int neighs[SIZE][4], int pos, double& energy, double p, mt19937& gen, uniform_real_distribution<double>& ran_u);
double get_energy(bool spins[SIZE], int neighs[SIZE][4]);

void write(bool spins[SIZE]);

void w_output(ofstream& file, double tstar, int N, double m[DATA]);


int main(void)
{
    bool spins[SIZE]; //Stores the spins

	int i,j; //Counters

	int N; //Number of averages done into the system

	double h[5]; //Values of the exp(-J/kT)

	double energy; //Value of the energy of the system

    double m[DATA]; //This array contains several moments of the magnetization and energy

	int neighs[SIZE][4]; //To store nearest neighbours

    mt19937 gen(958431198); //Mersenne Twister RNG
	uniform_int_distribution<int> brandom(0, 1); //Get any random integer
	uniform_int_distribution<int> ran_pos(0, SIZE-1); //Get any random integer
	uniform_real_distribution<double> ran_u(0.0, 1.0); //Our uniform variable generator

	double tstar; //Control parameter
	double deltat, deltat_crit; //Change in control parameter by iteration
	double tmax, tmin, tcrit_up, tcrit_down; //Max and min temperature, and interval where we apply Wolff

    tmax = 5.0;
    tmin = 0.1;
    tcrit_up = 2.4;
    tcrit_down = 2.2;

    deltat = 0.1;
    deltat_crit = 0.01;

	ofstream output; //Output of the stream

    initialize(spins, gen, brandom); //Init randomly
	get_neighbors(neighs); //Get neighbour table
	energy = get_energy(spins, neighs); //Compute initial energy

    N = 1e5;

	output.open("test_alt162.txt");
	for (tstar = tmax; tstar > tcrit_up; tstar -= deltat)
	{
	    do_step(spins, neighs, tstar, N, h, energy, gen, ran_u, ran_pos, m);
	    w_output(output, tstar, N, m);
	    cout << tstar << endl;
	}
    for (tstar = tcrit_up - deltat_crit; tstar > tcrit_down; tstar -= deltat_crit)
	{
	    do_step_wolff(spins, neighs, tstar, N, energy, gen, ran_u, ran_pos, m);
	    w_output(output, tstar, N, m);
	    cout << tstar << endl;
        //do_step(spins, neighs, tstar, N, h, energy, gen, ran_u, ran_pos, m);
        //output << 1.0/tstar << " " << m[0] << " "<< m[1] - m[0] * m[0] << " " << 1.0 - m[2]/(3.0 * m[1] * m[1]) << endl;
	}
    for (tstar = tcrit_down - deltat; tstar >= tmin; tstar -= deltat)
	{
	    do_step(spins, neighs, tstar, N, h, energy, gen, ran_u, ran_pos, m);
	    w_output(output, tstar, N, m);
	    cout << tstar << endl;
        //output << tstar << 1.0/tstar << " " << m[0] << " " << m[1] - m[0] * m[0] << " " << 1.0 - m[2]/(3.0 * m[1] * m[1]) << " " << m[4] - m[3] * m[3] << endl;
	}
	output.close();

	return 0;

}

void write(bool spins[SIZE])
{
    int i;
    ofstream output;
    output.open("check.txt");
    for (i=0; i < SIZE; i++)
    {
        if (i % SIZE == 0 and i != 0) output << endl;
        output << spins[i] << " ";
    }
    return;
}

//Initialices the grid in which we are going to do the Ising, using random values
void initialize(bool spins[SIZE], mt19937& gen, uniform_int_distribution<int>& brandom)
{
	int i,j;

	//Init spins with a random distribution
	for (i=0; i < SIZE; i++)
	{
        spins[i] = brandom(gen); //Generate numbers
	}

	return;
}

//Fills the neigbour table
void get_neighbors(int neighs[SIZE][4])
{
	int i,j;
	int u,d,r,l;

	for (i=0; i < L; i++)
	{
		for (j=0; j < L; j++)
		{
		    //Get the (x,y) with periodic boundaries
			u = j+1 == L ? 0 : j+1;
			d = j-1 == -1 ? L-1 : j-1;
			r = i+1 == L ? 0 : i+1;
			l = i-1 == -1 ? L-1 : i-1;

            //(x,y) to index notation and store in table
			neighs[i+j*L][UP] = i+u*L;
			neighs[i+j*L][DOWN] = i+d*L;
			neighs[i+j*L][RIGHT] = r+j*L;
			neighs[i+j*L][LEFT] = l+j*L;
		}
	}

	return;
}

//Do all the things needed for a certain temperature
void do_step(bool spins[SIZE],  int neighs[SIZE][4], double tstar, int N, double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos, double m[DATA])
{

	int i,j; //Counters
	double sum;//To compute the sum of spins
    double energysum;
    double chi, heat;
	double old_sum, old_chi, old_heat, old_energy;

    for (i=0; i < DATA; i++) m[i] = 0.0; //Init the values

	//Compute the factors of exp(-dH/kT)
	for (i=-4; i <= 4; i += 2)
	{
		h[(i+4)/2] =  min(1.0, exp(- 2.0 * i / tstar));
	}

	//Thermalize the state
	for (j=0; j < 1100; j++)
	{
		flip_spin(spins, neighs, h,  energy, gen, ran_u, ran_pos);
	}

	///----- TODO: optimize the number of steps for thermalization/measures
	old_sum = 0.0;
	old_chi = 0.0;
	old_heat = 0.0;
	old_energy = 0.0;
	for (i=0; i < N; i++)
	{
		//Make changes and then average
		for (j=0; j < 1100; j++)
		{
			flip_spin(spins, neighs, h,  energy, gen, ran_u, ran_pos);
		}

        //Compute quantities at time j
        sum = abs(magnetization(spins));
        chi = sum * sum;
        heat = energy * energy;

        //Add all the quantities
		m[MAG] += sum; //Magnetization
		m[MAG2] += chi; //For the susceptibility
		m[MAG4] += chi * chi; //For the Binder cumulant and also variance of susceptibility
		m[ENE] += energy; //Energy
		m[ENE2] += heat; //For specific heat
		m[ENE4] += heat * heat; //For the variance of specific heat
		//This are used for errors,
		m[MAGERR] += old_sum * sum; //in magnetization
		m[SUSERR] += old_chi * chi; //in susceptibility
		m[ENERR] += old_energy * energy; //in energy
		m[CHERR] += old_heat * heat; //in specific heat

		//Get the value for the next iteration
		old_sum = sum;
		old_energy = energy;
		old_chi = chi;
		old_heat = heat;
	}

    //Finish the average
    for (i=0; i < DATA; i++)  m[i] /= (1.0 * N);

	return;
}


//Do all the things needed for a certain temperature
void do_step_wolff(bool spins[SIZE],  int neighs[SIZE][4], double tstar, int N, double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos, double m[DATA])
{

	int i,j; //Counters
	double sum; //To compute the sum of spins
    double chi, heat; //To compute magnetic susceptibility
	//Note: we use directly variable energy for the energy
    //To remember last results
	double old_sum, old_chi, old_heat, old_energy;


	double pa = 1.0 - exp(- 2.0 / tstar); //TODO change

    for (i=0; i < DATA; i++) m[i] = 0.0; //Init the values

    //Make changes and then average
    for (j=0; j < 15; j++) add_to_cluster(spins, neighs, ran_pos(gen), energy, pa, gen, ran_u);

	///----- TODO: optimize the number of steps for thermalization/measures
	old_sum = 0.0;
	old_chi = 0.0;
	old_heat = 0.0;
	old_energy = 0.0;


	for (i=0; i < N; i++)
	{
		//Make changes and then average
        for (j=0; j < 12; j++) add_to_cluster(spins, neighs, ran_pos(gen), energy, pa, gen, ran_u);

        //Compute quantities at time j
        sum = abs(magnetization(spins));
        chi = sum * sum;
        heat = energy * energy;

        //Add all the quantities
		m[MAG] += sum; //Magnetization
		m[MAG2] += chi; //For the susceptibility
		m[MAG4] += chi * chi; //For the Binder cumulant and also variance of susceptibility
		m[ENE] += energy; //Energy
		m[ENE2] += heat; //For specific heat
		m[ENE4] += heat * heat; //For the variance of specific heat
		//This are used for errors,
		m[MAGERR] += old_sum * sum; //in magnetization
        m[ENERR] += old_energy * energy; //in energy
		m[SUSERR] += old_chi * chi; //in susceptibility
		m[CHERR] += old_heat * heat; //in specific heat

		//Get the value for the next iteration
		old_sum = sum;
        old_energy = energy;
		old_chi = chi;
		old_heat = heat;
	}

    //Finish the average
    for (i=0; i < DATA; i++)  m[i] /= (1.0 * N);

	return;
}

//Flip a spin via Metropolis
void flip_spin(bool spins[SIZE], int neighs[SIZE][4], double h[5], double& energy, mt19937& gen, uniform_real_distribution<double>& ran_u, uniform_int_distribution<int>& ran_pos)
{
    int index = ran_pos(gen); //Get a random position to flip
    //Compute the sum of neighbours
    int sum_neigh = spins[neighs[index][UP]] + spins[neighs[index][DOWN]] + spins[neighs[index][RIGHT]] + spins[neighs[index][LEFT]];
    //Use this to get the energy change (depending on the value of my spin)
    int change = spins[index] ? 2.0 * (sum_neigh) - 4.0 : 4.0 - 2.0 * (sum_neigh);

    //Apply Metropolis skim
    if (ran_u(gen) < h[(change+4)/2])
    {
        spins[index] = !spins[index];
        energy += (2.0*change)/(1.0*SIZE);
        //cout << change << "  " << (2.0*change)/(1.0*SIZE) << "  " << energy << endl;
    }

    return;
}

//Compute the magnetization
double magnetization(bool spins[SIZE])
{
    int i;

    double sum = 0.0;
    //Sum all the values of the spins
    for (i=0; i < SIZE; i++)
    {
        sum += spins[i];
    }
    //And then return them
    return (2.0*sum - SIZE)/(SIZE);
}

//Computes the energy of the system
double get_energy(bool spins[SIZE], int neighs[SIZE][4])
{
    int i; //Counters
    int sum_neigh;

    int energy = 0; //Sum

    //For every spin,
    for (i=0; i < SIZE; i++)
    {
        //Get sum of the neighbours
        sum_neigh = spins[neighs[i][UP]] + spins[neighs[i][DOWN]] + spins[neighs[i][RIGHT]] + spins[neighs[i][LEFT]];
        //And compute the energy change
        energy += spins[i] ? 2.0 * (sum_neigh) - 4.0 : 4.0 - 2.0 * (sum_neigh);
    }

    return 2.0*energy/(1.0*SIZE); //Return the energy
}

//Executes Wolff algorithm in a recursive way
void add_to_cluster(bool spins[SIZE], int neighs[SIZE][4], int pos, double& energy, double p, mt19937& gen, uniform_real_distribution<double>& ran_u)
{
    int i; //Counter
    int n; //Neighbour position


    //Compute sum of neighbours and change in energy
    int sum_neigh = spins[neighs[pos][UP]] + spins[neighs[pos][DOWN]] + spins[neighs[pos][RIGHT]] + spins[neighs[pos][LEFT]];
    int delta_energy = spins[pos] ? 2.0 * (sum_neigh) - 4.0 : 4.0 - 2.0 * (sum_neigh);

    //Then modify the energy
    energy += (2.0*delta_energy)/(1.0*SIZE);


    spins[pos] = !spins[pos]; //Flip the spin

    //For every neighbour,
    for (i=0; i < 4; i++)
    {
        n = neighs[pos][i]; //Get its position
        if (spins[n] != spins[pos]) //Check if it the same (remember we flipped)
        {
            //Add it to the cluster with certain probability.
            if (ran_u(gen) < p)
            {
                add_to_cluster(spins, neighs, n, energy, p, gen, ran_u);
            }
        }
    }

    return;
}

void w_output(ofstream& file, double tstar, int N, double m[DATA])
{
    file << tstar << " " << 1.0/tstar << " "; //Write T and B

    //We here take in account residual errors, which, for low T, makes the quantities chi, ch, etc.
    //to diverge. This must be substracted. That's why we use an abs for correlation time and also
    //a check to avoid zero value of variances.

    //Then write the quantities and the corresponding errors to a file. The four things are equal,
    //but each one referred to a different quantity.

    double chi = m[MAG2] - m[MAG] * m[MAG]; //Magnetic susceptibility (up to T factor)
    double rhom = chi != 0 ? (m[MAGERR] - m[MAG] * m[MAG]) / chi : 0.0; //Rho magnetization, computed if chi != 0
    double taugm = rhom != 1.0 ? rhom / (1.0 - rhom) : 0.0; //Taug magnetization, computed if rhom != 0
    file << m[MAG] << " " << sqrt(chi * abs(2.0 * taugm + 1) / (1.0*N)) << " "; //Write everything

    double fourth = m[MAG4] - m[MAG2] * m[MAG2]; //Susceptibility variance
    double rhos = fourth != 0.0 ? (m[SUSERR] - m[MAG2] * m[MAG2]) / fourth : 0.0; //Rho susceptibility
    double taugs = rhos != 1.0 ? rhos /(1.0 - rhos) : 0.0; //Taug susceptibility
    double error_sq = sqrt(fourth * abs(2.0 * taugs + 1) / (1.0*N));
    file << " " << chi << " " << error_sq << " ";

    double heat = m[ENE2] - m[ENE] * m[ENE]; //Specific heat (up to T^2 factor)
    double rhoe = heat != 0.0 ? (m[ENERR] - m[ENE]*m[ENE]) / heat : 0.0;
    double tauge = rhoe != 1.0 ? rhoe / (1.0 - rhoe) : 0.0;
    file << " " << m[ENE] << " " << sqrt(heat * abs(2.0 * tauge + 1) / (1.0*N)) << " ";

    double fourth_ene = m[ENE4] - m[ENE2] * m[ENE2];
    double rhoc = fourth_ene != 0.0 ? (m[CHERR] - m[ENE2] * m[ENE2]) / fourth_ene : 0.0;
    double taugc = rhoc != 1.0 ? rhoc / (1.0 - rhoc) : 0.0;
    file << " " << heat << " " << sqrt(fourth_ene * abs(2.0 * taugc + 1) / (1.0*N)) << " ";

    //Binder cumulant
    double binder = 1.0 - m[MAG4]/(3.0 * m[MAG2] * m[MAG2]); //Computes 4th cumulant minus one, b-1.
    file << binder << " " << 2.0 * (1.0 - binder) * (error_sq / m[MAG2]) << endl;
    return;
}
