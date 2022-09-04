#pragma once
#include <unordered_map>
#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <windows.h>
#include <fileapi.h>
#include <random>
#include <chrono>

#include <cuda_runtime.h>

#include "cublas_v2.h"
#include "cusparse_v2.h"

#include "cublas_wrapper.h"
#include "misc.h"
#include "preconditioner.h"
#include "solver.h"

#include "tools/globals.h"
#include "Filament.h"


class Plasma_Ball
{
private:

	/*
	*
	*/
	bool _debug;

	/*
	*
	*/
	bool _grad_descent;

	/*
	* threshold for convergence
	*/
	double _convergence;

	/*
	* [0, 1] range of random numbers for DBM
	*/
	double _randomness;
	
	/*
	* [1, 32] how jagged and branching - pre novel
	*/
	double _eta;

	/*
	* Whether new points are found along whole filament or only at tip
	*/
	bool _front;

	/*
	* Serial execution is filaments independently
	*/
	bool _serial;

	/*
	* Whether the potentials get updated during filament growth
	*/
	bool _updating;

	/*
	* [0, 1] amount of branching that occurs
	*/
	double _branching;

	/*
	* [0, 1] range for determining field strength
	*/
	double _field_strength;

	/*
	* Number of filaments and branches
	*/
	unsigned _number_of_filaments;

	/*
	* Number of roots placed
	*/
	unsigned _number_of_roots;
	
	/*
	* True while a filament is still alive
	*/
	bool _alive;

	/*
	* Applied voltage
	*/
	double _Vapp;

	/*
	*
	*/
	Point _centre;

	/*
	* max stepsize of filament
	*/
	int _stepsize;

	/*
	* Actual voltage
	*/
	double _Vact;

	/*
	* How many cells in each dimension
	*/
	int _resolution;

	/*
	* CIURCULAR or SPERICAL
	*/
	Geometry _dim;

	/*
	* radius of ball in m
	*/
	double _real_size;

	/*
	* Size of inner electrode in m
	*/
	double _real_electr_size;

	/*
	* Size of ball cells
	*/
	double _radius;

	/*
	* _radius + _buffer
	*/
	double _grid_size;

	/*
	* Size of electrode in cells
	*/
	double _electrode_radius;

	/*
	* buffer at edges
	*/
	double _buffer;

	/*
	* _real_size / _resolution
	*/
	double _dxyz;

	/*
	*  real inner span of ball: _real_size - _real_electr_size
	*/
	double _min_fil_length;

	/*
	* _dxyz * _dxyz
	*/
	double _dxyz2;

	/*
	* time step size
	*/
	double _dt;

	/*
	*
	*/
	double _epsilon0 = 0.000000000008854187;

	/*
	* Map of regions of Materials
	*/
	GRID _Map;

	/*
	* _Map with plasma regions
	*/
	GRID _Mask;

	/*
	* how many steps have been executed
	*/
	int _step = 0;

	/*
	* cumulative vertical shift
	*/
	double _shift = 0;

	/*
	* amount of charge left by plasma relative to Vact
	*/
	double _residual;

	/*
	*
	*/
	int _fil_step;

	/*
	* Map containing electrical potentials as calculated by Poisson
	*/
	GRID _potentials;

	/*
	* Map containing positions and magnitude of charges placed in ball (i.e. plasma)
	*/
	GRID _charge_map;

	/*
	* Copy of charge map edited for root placement
	*/
	GRID _root_charge_map;

	GRID _gradients;

	/*
	* Tool to visualise directions of E-Field
	*/
	GRID _gradient_directions;

	GRID _Gx;

	GRID _Gy;
	
	GRID _Gz;

	/*
	*
	*/
	std::vector<GRID> _leaders;

	//GRID _Sparse_sol;

	/*
	* all maps of potentials
	*/
	std::vector<GRID> _potentials_maps;

	/*
	* Keep track of named maps to be shown
	*/
	std::unordered_map<std::string, int> show_counters_;

	Filament _root_points;

	std::vector<Filament> _filaments;

	std::vector<Filament> _last_step_filaments;

	std::vector<Filament> _branches;
	/*
	* list of all filaments and whether they are active filaments still growing
	*/
	int _active_filaments;

	int _electrode_count;
	int _plasmas; //plasma count

	/*
	*
	*/
	std::string _rootdir;
	
	/*
	*
	*/
	std::string _dir;

	/*
	*
	*/
	std::string _logname = "log.txt";

	/*
	*
	*/
	std::ofstream _logfile;

	/*
	*
	*/
	std::string _name;

	/*
	*
	*/
	float *_cuda_potentials;
	
	float *_cuda_mask;

	float *_cuda_charge_map;

	void Poisson(int its = 1);

	void Poisson_Neuman(int its = 1);

	void create_dir();

	void init_map();

	void reset_maps();

	void create_det_root_points();

	/*
	* to set init roots as you wish (evenly spaced)
	*/
	void select_init_roots(int indx);

	void select_root();

	void perform_DBM();

	void add_charge(double val = 0);

	bool find_living_filament();

	void clean_maps();

public:



	Plasma_Ball(Geometry geo = CIRCULAR, double realsize = 0.15, double resolution = 64, double voltage = -10000);

	void init_DBM(unsigned fils = 0, int stepsize = 1, double branching = 0.0, int fps = 1, double residual = 0.01);

	void get_gradient();

	void get_gradient_directions();
	
	void step(bool debug = false);

	void finish_step();

	void perform_grad_descent();

	void set_voltage(double v);
	
	void place_charges();

	void vert_drift();

	void set_updating(bool b);

	void set_grad_descent(bool b);

	GRID get_potentials();

	GRID get_map();

	GRID get_mask();

	void show(GRID &Map, std::string name = "");

	void name_dir(std::string name);
};
