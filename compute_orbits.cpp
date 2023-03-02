// Assessment 3: n-body gravitational solver

// To avoid warnings tell the compiler to use a recent standard of C++:
// g++ -std=c++17 vector3d.cpp body.cpp compute_orbits.cpp -o compute_orbits
// ./compute_orbits sun_earth.csv test_case.csv 1e-4 110
// ./compute_orbits sun_earth.csv test_case.csv 1e-3 1 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include "vector3d.hpp"
#include "body.hpp"

using std::cout, std::endl;

// *** Function declarations ***

// -----------------------------------------------------------
// ADD YOUR FUNCTION DECLARATIONS HERE
void compute_energy_L(std::vector<body> &system);
void update_acc(std::vector<body> &system);
void vel_verlet(std::vector<body> &system, double dt);


// -----------------------------------------------------------
// Read input data from file
void read_init(std::string input_file, std::vector<body> &system);
// Read the components of a 3d vector from a line
void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z);
// Save the data to file
void save_data(std::ofstream& savefile, const std::vector<body> &system, double t);

int main(int argc, char* argv[])
{
  // Checking if number of arguments is equal to 4:
  if (argc != 6) {
    cout << "ERROR: need 4 arguments - compute_orbits <input_file> <output_file> <dt> <T> <Tsave>" << endl;
    return EXIT_FAILURE;
  }
  // Process command line inputs:
  std::string input_file = argv[1];
  std::string output_file = argv[2];
  double dt = atof(argv[3]); // Time step
  int T = atoi(argv[4]); // Total number of time steps
  int Tsave = atoi(argv[5]); // Number of steps between each save

  std::vector<body> system; // Create an empty vector container for bodies
  read_init(input_file, system); // Read bodies from input file into system
  int N = system.size(); // Number of bodies in system

  cout << "--- Orbital motion simulation ---" << endl;
  cout << " number of bodies N: " << N << endl;
  for(int p=0; p < N; p++)
  {
    cout << "- " << system[p].get_name() << endl; // Display names
  }
  cout << "       time step dt: " << dt << endl;
  cout << "  number of steps T: " << T << endl;
  cout << "   save steps Tsave: " << Tsave << endl;

  std::ofstream savefile (output_file); // Open save file
  if (!savefile.is_open()) {
    cout << "Unable to open file: " << output_file << endl; // Exit if save file has not opened
    return EXIT_FAILURE;
  }
  savefile << std::setprecision(16); // Set the precision of the output to that of a double
  // Write a header for the save file
  savefile << dt << "," << T << "," << Tsave << "," << N << endl;
  for(int p=0; p < (N-1); p++)
  {
    savefile << system[p].get_name() << ",";
  }
  savefile << system[N-1].get_name() << endl;

  // -----------------------------------------------------------
  // ADD YOUR CODE HERE
  // Initialise variables
  double total_time;
  int Tsave_counter = 0;

  // Loop steps until it reaches total number for time steps 'T'.
  for (int t=1; t<=T; t++)
  {
    // Runs Vel_verlet function.
    vel_verlet(system, dt);
    
    // Adds 1 to Tsave step counter.
    Tsave_counter++;

    // Checks if step counter is equal to Tsave.
    if (Tsave_counter == Tsave)
    {
      // Calculates the total time at this timestep.
      total_time = t*dt;

      // Computes energies and momentum and saves data to file.
      compute_energy_L(system);
      save_data(savefile, system, total_time);
      
      // Reset Tsave counter to 0
      Tsave_counter = 0;   
    }
  }
  //compute_energy_L(system);
  //update_acc(system);
  //save_data(savefile, system, 0); // Example of saving data (for initial state here) 

  // -----------------------------------------------------------
  
  savefile.close();
  return EXIT_SUCCESS; 
}

// *** Function implementations ***


// -----------------------------------------------------------
// ADD YOUR FUNCTION IMPLEMENTATIONS HERE
void vel_verlet(std::vector<body> &system, double dt)
{
  int N = system.size();
  int i;
  vec pos;
  vec velocity;
  vec new_velocity;
  vec acceleration; 
  vec new_acceleration;
  vec add_acceleration;
  vec new_pos;
  update_acc(system);
  std::vector<body> new_system = system;

  // Loop through planets in the system updating velocity and acceleration.
  for (i = 0; i < N; i++)
    {
      //Initialise planet variables
      
      pos = system[i].get_pos();
      velocity = system[i].get_vel();
      acceleration = system[i].get_acc();
      
      // Calculate components of the equation.
      velocity = velocity *=(dt);
      acceleration = acceleration*=(0.5 * dt * dt);

      // This could cause an error
      // Calculates the new position at the new time.
      new_pos = pos +=(velocity)+=(acceleration);
      
      // Sets new positions of planet
      new_system[i].set_pos(new_pos);

      // Calculates the new acceleration of planet
      update_acc(new_system);

      // Calculates new velocity of planet using old and new variables
      new_acceleration = new_system[i].get_acc();
      add_acceleration = new_acceleration +=(acceleration);
      new_velocity = add_acceleration*=(dt*0.5);
      new_velocity = new_velocity +=(velocity);

      // Saves new velocity.
      new_system[i].set_vel(new_velocity);
    }
  
  // Sets all the new values back to the origonal object 'system'
  system = new_system;
}

void update_acc(std::vector<body> &system)
{
  int N = system.size();
  int i;
  int j;
  vec acceleration;
  vec answer;
  vec pos1;
  vec pos2;
  double mass1;
  vec distance;
  double length;
  double length_cubed;

  for (i = 0; i < N; i++)
    {
      // Initialise variables for planet 1.
      vec acceleration(0,0,0);
      pos1 = system[i].get_pos();
      mass1 = system[i].get_mass();

      for (j=0; j < N; j++)
      {
        if (i!=j)
        {
          // Initialise variables for planet 2.
          pos2 = system[j].get_pos();

          // Calculates the vector between them.
          distance = pos1 -(pos2);
          
          // Calculates the length of the vector.
          // It will always be positive as its been squared.
          length = distance.length();
          length_cubed = pow(length, 3);
          // This could cause an error
          answer = distance/=(length_cubed)*=(mass1);
          //answer = answer*=(mass1);

          // Sums the different accelerations on the planet.
          acceleration +(answer);
        }
      }
    }
    // Sets acceleration variable of planet.
    system[i].set_acc(acceleration);
}

void compute_energy_L(std::vector<body> &system)
{
    int N = system.size();
    int i;
    int j;
    int k;
    vec pos;
    vec pos1; 
    vec pos2;
    double mass;
    double mass1; 
    double mass2; 
    vec velocity;
    vec velocity1; 
    vec velocity2; 
    vec distance;
    double length;
    double gpe = 0.;
    double ke;
    vec L;

    for (i = 0; i < N; i++)
    {
      // Resets the value of gpe for each planet
      gpe = 0.0;
      // Initialise variables for planet 1
      pos1 = system[i].get_pos();
      mass1 = system[i].get_mass();
      velocity1 = system[i].get_vel();

      for (j=0; j < N; j++)
      {
        if (i!=j)
        {
          // Initialise variables for planet 2.
          pos2 = system[j].get_pos();
          mass2 = system[j].get_mass();
          velocity2 = system[j].get_vel();

          // Calculates the vector between them.
          distance = pos1 -(pos2);
          
          // Calculates the length of the vector.
          // It will always be positive as its been squared.
          length = distance.length();
          
          // Calculates the gpe.
          
          gpe += -(mass1*mass2)/length;
          //std::cout << "gpe" << gpe << std::endl;
          
          // Sets the gpe value to the system 1.
          system[i].set_gpe(gpe);
          
        }
      }
    }

    for (k = 0; k < N; k++) 
    {
      // Calculates the kinetic energy 'ke'
      velocity = system[k].get_vel();
      mass = system[k].get_mass();
      ke = velocity.dot(velocity);
      ke *= 0.5 * mass;
     
      // Sets the variable for kinetic energy
      system[k].set_ke(ke);
      
      // Calculates the angular momentum vector 'L'
      pos = system[k].get_pos();
      L = pos*(mass);
      L = L.cross(velocity);

      // Sets the angular momentum variable for this planet.
      system[k].set_L(L);
    }
    
    
    
    //return EXIT_SUCCESS;
}

// -----------------------------------------------------------
void read_init(std::string input_file, std::vector<body> &system)
{
  std::string line; // Declare a string to store each line
  std::string name; // String to store body name
  double m, x, y, z, vx, vy, vz; // Doubles to store vector components
  int line_cnt = 0; // Line counter

  // Declare and initialise an input file stream object
  std::ifstream data_file(input_file); 

  while (getline(data_file, line)) // Read the file line by line
  {
    line_cnt++;
    std::stringstream data_line(line); // Create a string stream from the line
    switch (line_cnt)
    {
      case 1:
        name = line;
        break;
      case 2:
        m = std::stod(line); // Convert string line into double
        break;            
      case 3:
        read_vector3d(data_line, x, y, z); // Read the 3 components of the vector on the line
        break;
      case 4:
        read_vector3d(data_line, vx, vy, vz); // Read the 3 components of the vector on the line
        break;                   
      }
      if (line_cnt==4) // Data for one body has been extracted
      {
        line_cnt = 0; // Reset line counter
        body b(name,m,vec(x,y,z),vec(vx,vy,vz)); // Package data into body
        system.push_back(b); // Add body to system
      }
    }
    // Close the file
    data_file.close();
}

void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z)
{
  std::string value; // Declare a string to store values in a line
  int val_cnt = 0; // Value counter along the line
  
  while (getline(data_line, value, ','))
  {
    val_cnt++;
    switch (val_cnt)
    {
      case 1:
        x = std::stod(value); 
        break;
      case 2:
        y = std::stod(value); 
        break;
      case 3:
        z = std::stod(value);
        break;
    }
  } 
}

void save_data(std::ofstream& savefile, const std::vector<body> &system, double t)
{
    // Function for saving the simulation data to file.

    vec L; // Total angular momentum
    double E = 0.0; // Total energy
    for(int p = 0; p < system.size(); p++)
    { 
      E += system[p].get_ke() + 0.5*system[p].get_gpe();
      L += system[p].get_L();
    }
    double Lmag = L.length(); // Magnitude of total angular momentum

    // Write a header for this time-step with time, total energy and total mag of L:
    savefile << t << "," << E << "," << Lmag << endl;

    // Loop over the bodies:
    for(int p = 0; p < system.size(); p++)
    { 
      // Output position and velocity for each body:
      savefile << system[p].get_pos() << "," << system[p].get_vel() << endl;
    }
}