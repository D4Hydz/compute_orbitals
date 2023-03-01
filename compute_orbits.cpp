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

#include "vector3d.hpp"
#include "body.hpp"

using std::cout, std::endl;

// *** Function declarations ***

// -----------------------------------------------------------
// ADD YOUR FUNCTION DECLARATIONS HERE
double compute_energy_L(std::vector<body> &system);


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
  compute_energy_L(system);
  


  save_data(savefile, system, 0); // Example of saving data (for initial state here) 



  // -----------------------------------------------------------
  
  savefile.close();
  return EXIT_SUCCESS; 
}

// *** Function implementations ***


// -----------------------------------------------------------
// ADD YOUR FUNCTION IMPLEMENTATIONS HERE
double compute_energy_L(std::vector<body> &system)
{
    // Creates variables for the positions of the 2 bodies.
    
    vec pos1 = system[0].get_pos();
    vec pos2 = system[1].get_pos();
    double mass1 = system[0].get_mass();
    double mass2 = system[1].get_mass();
    vec velocity1 = system[0].get_vel();
    vec velocity2 = system[1].get_vel();
    
    // Calculates the vector between them.
    vec distance = pos1 -(pos2);
    std::cout << "Vector distance " << distance << std::endl;
    
    // Calculates the length of the vector.
    // I need to make sure it is positive though...
    double length = distance.length();
    std::cout << "Length " << length << std::endl;
    
    // Calculates the gpe.
    double gpe = -(mass1*mass2)/length;
    //std::cout << "gpe" << gpe << std::endl;
    
    // Sets the gpe value to the system 1.
    system[1].set_gpe(gpe);
    system[0].set_gpe(gpe);
    
    // Output and check it is a valid value.
    double a = system[1].get_gpe();
    //std::cout << "gpe " << a << std::endl;
    
    // Calculates the kinetic energy
    double ke2 = velocity2.dot(velocity2);
    ke2 *= 0.5 * mass2;
    
    double ke1 = velocity1.dot(velocity1);
    ke1 *= 0.5 * mass1;
    
    // Sets the variable for kinetic energy
    system[1].set_ke(ke2);
    system[0].set_ke(ke1);
    
    // Calculates the angular momentum vector
    vec L1 = pos1*(mass1);
    vec L2 = pos2*(mass2);
    L1 = L1.cross(velocity1);
    L2 = L2.cross(velocity2);
    
    system[0].set_L(L1);
    system[1].set_L(L2);
    
    
    
    
    return EXIT_SUCCESS;
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