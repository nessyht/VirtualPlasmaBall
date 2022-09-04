#pragma once
#include <random>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "tools/globals.h"
# define M_PI           3.14159265358979323846  /* pi */
class Filament
{

	std::vector<Point> _electrode_root_points; //ONLY FOR THE CENTRAL ELECTRODE

public:

	bool _is_branch = false;
	bool _is_parent = false;
	bool _found_branch = false;
	double _branch_pre = 100000; // previous branching length
	
	double _min_length; // distance from inner to outer
	double _resolution; // squares per m
	
	Point _root;
	Point _shifted_root;
	Point _new_point;
	Point _next_new_point;
	Point _new_branch_point;
	bool _front = false; // is there a new point
	std::vector<Point> _points;
	std::vector<Point> _ordered_points; //if there is an order to the points 
	std::vector<Point> _nodes;
	std::vector<Point> _new_point_nodes;
	std::vector<Point> _charges;
	std::vector<double> _phis;
	std::vector<double> _new_point_phis;
	std::vector<double> _conductivities;

	std::vector<VEC> _descent_points;


	bool _alive;
	bool _valid; // whether filament made it to end or should be deleted
	double _charge;
	VEC _ball_centre;

	Filament();
	Filament(Point root, VEC ball_centre, double min_length, double resolution, bool is_branch = false);

	void extend();

	void descent_steps(const GRID &Gx, const GRID &Gy, const GRID &Gz);

	/*
	Iterate through all the nodes around the new point
	Any valid (not a point of the Filament or already a node and inside the map) points are added to the nodes
	*/
	void find_nodes(const GRID &mask);

	void get_front_nodes(const GRID &map);

	/*
	* Return 8 all points in contact with passed point (window size 3)
	*/
	std::vector<Point> get_nodes(Point point, int window_size, int dim);

	void fill_phis(const GRID & map, bool punish_diagonals = false);

	
	void grad_find_next(const GRID & Gx, const GRID & Gy, const GRID & Gz, const GRID &map, double branching = 0.0);
	void grad_find_stepped_next(const GRID & Gx, const GRID & Gy, const GRID & Gz, const GRID & map, const int stepsize, double branching = 0);
	void grad_find_next(const GRID &grad_dir, const GRID &map, double branching = 0.0);

	void init_find_next(const GRID &grad_dir, const GRID &map);

	void init_find_next(const GRID & Gx, const GRID & Gy, const GRID & Gz, const GRID &map);

	void grad_find_next_root(GRID grad, GRID grad_dir);

	/*
	Using the phis calculate next point to be selected
	*/
	void find_next(double dim = 1.0, bool front = false, double rand = 1.0, double branching = 0.0);

	Point get_next_electrode_root_point(int indx);	
	void set_next_electrode_root_point(std::vector<Point> selected);
};

