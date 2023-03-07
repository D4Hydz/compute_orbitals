// Assessment 3: n-body gravitational solver

// To avoid warnings tell the compiler to use a recent standard of C++:
// g++ -std=c++17 vector3d.cpp body.cpp compute_orbits.cpp -o compute_orbits

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


// ######################################
// ####### Function declarations ########
// ######################################
// -----------------------------------------------------------
// ADD YOUR FUNCTION DECLARATIONS HERE

// Compute the gravitational potential energy, the kinetic energy and the angular momentum. 
void compute_energy_L(std::vector<body> &system);
// Calculate the accelerations of planets in the system at a timestep.
void update_acc(std::vector<body> &system);
// Calculate the new positions and velocites of planets in the system at a given timestep.
void vel_verlet(std::vector<body> &system, double dt);
// -----------------------------------------------------------
// Read input data from file
void read_init(std::string input_file, std::vector<body> &system);
// Read the components of a 3d vector from a line
void read_vector3d(std::stringstream& data_line, double& x, double& y, double& z);
// Save the data to file
void save_data(std::ofstream& savefile, const std::vector<body> &system, double t);


// ######################################
// ################ MAIN ################
// ######################################
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
  double total_time = 0.;
  int Tsave_counter = 0;
  
  // Calculate and save data for timestep t = 0. 
  compute_energy_L(system);
  save_data(savefile, system, total_time);
    
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
      total_time = t * dt;

      // Computes energies and momentum and saves data to file.
      compute_energy_L(system);
      save_data(savefile, system, total_time);
      
      // Reset Tsave counter to 0.
      Tsave_counter = 0;   
    }
  }
  // -----------------------------------------------------------
  savefile.close();
  return EXIT_SUCCESS; 
}


// ######################################
// ###### Function implementations ######
// ######################################
// -----------------------------------------------------------
// ADD YOUR FUNCTION IMPLEMENTATIONS HERE
void vel_verlet(std::vector<body> &system, double dt)
{
  // Initialise variables.
  int i;
  vec pos;
  vec velocity;
  vec acceleration;
  vec new_acceleration;
  vec new_position;
  vec new_velocity;
  int N = system.size();
  
  // Creates a copy of the vector class system.
  std::vector<body> new_system = system;

  // Updates the acceleration of the system.
  update_acc(system);

  // Loop through planets in the system updating position.
  for (i = 0; i < N; i++)
  {
    // Initialise variables for this planet.
    pos = system[i].get_pos();
    velocity = system[i].get_vel();
    acceleration = system[i].get_acc();
    
    // Calculates and sets new position to the new system.
    new_position = pos + (velocity * (dt)) + (acceleration * (dt*dt*0.5));
    new_system[i].set_pos(new_position);
  }
  
  // Calculates new accelerations based on all the new positions.
  update_acc(new_system);

  // Loop to update the new velocities of planets.
  for (i=0; i<N; i++)
  {
    // Initialise variables for planets in this loop.
    new_acceleration = new_system[i].get_acc();
    acceleration = system[i].get_acc();
    velocity = system[i].get_vel();

    //Calculates the new velocities based on the new acceleration and old velocity.
    new_velocity = new_acceleration + (acceleration);
    new_velocity *= (dt*0.5);
    new_velocity += (velocity);
    
    // Sets new velocity to new system.
    new_system[i].set_vel(new_velocity);

  }
  // Returns the updated values back to original system.
  system = new_system;
} 

void update_acc(std::vector<body> &system)
{
  // Initalise variables.
  int p;
  int j;
  vec answer;
  vec posp;
  vec posj;
  vec distance;
  double massj;
  double length;
  double length_cubed;
  int N = system.size();

  // Loop through first planet.
  for (p=0; p<N; p++)
    {
      // Initialise variables for first planet.
      vec acceleration(0.,0.,0.);
      posp = system[p].get_pos();

      // Loop through second planet.
      for (j=0; j<N; j++)
      {
          // Planets must be different!
          if (p != j)
          {
          // Initialise variables for second planet.
          posj = system[j].get_pos();
          massj = system[j].get_mass();

          // Calculates the vector between them.
          distance = posj - (posp);
          
          // Calculates the length of the vector.
          // It will always be positive as its been squared.
          length = distance.length();
          length_cubed = pow(length, 3);
              
          // Futher calculation for acceleration.
          answer = distance / (length_cubed);
          answer *= (massj);

          // Sums the different accelerations on the planet.
          acceleration += (answer);
          }
      }
      // Sets acceleration variable of planet.
      system[p].set_acc(acceleration);
    }
}

void compute_energy_L(std::vector<body> &system)
{
    // Initialise variables.
    int i;
    int j;
    int k;
    vec pos;
    vec pos1; 
    vec pos2;
    vec velocity; 
    vec distance;
    vec L;
    double mass;
    double mass1; 
    double mass2; 
    double length;
    double ke;
    double gpe;
    int N = system.size();

    // Loop Planets to calculate gpe.
    for (i=0; i<N; i++)
    {
      // Resets the value of gpe for each planet
      gpe = 0.0;
      
      // Initialise variables for first planet.
      pos1 = system[i].get_pos();
      mass1 = system[i].get_mass();
      
      // Loop other planets.
      for (j=0; j<N; j++)
      {
        // Skip if it is the same planet.
        if (i != j)
        {
          // Initialise variables for second planet.
          pos2 = system[j].get_pos();
          mass2 = system[j].get_mass();

          // Calculates the vector between them.
          distance = pos1 - (pos2);
          
          // Calculates the length of the vector.
          // It will always be positive as its been squared.
          length = distance.length();
          
          // Sums the gpe based on every planet interaction.
          gpe += (mass1 * mass2) / length;   
        }
      }
      // Negates the value.
      gpe *= -1;
      
      // Sets the final gpe value to the planet.
      system[i].set_gpe(gpe);
    }

    // Loop planets to calculate kinetic energy and angular momentum.
    for (i=0; i<N; i++) 
    {
      // Calculates the kinetic energy 'ke'
      velocity = system[i].get_vel();
      mass = system[i].get_mass();
      ke = velocity.dot(velocity);
      ke *= (0.5 * mass);
     
      // Sets the variable for kinetic energy
      system[i].set_ke(ke);
      
      // Calculates the angular momentum vector 'L'
      pos = system[i].get_pos();
      L = pos.cross(velocity);
      L *= (mass);

      // Sets the angular momentum variable for this planet.
      system[i].set_L(L);
    }
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