#include "Filament.h"

static std::vector<double> generate_data(size_t size, double rand = 0)
{
	using value_type = double;
	// Use static in order to instantiate the random engine
	// and the distribution once only.
	// Had some thread-safety issues at 8 CPU threads.
	static std::uniform_real_distribution<value_type> distribution(0.0, 1.0);
	static std::default_random_engine generator;

	std::vector<value_type> data(size + 1); // extra is placeholder for later use
	std::generate(data.begin(), data.end(), []() { return distribution(generator); });
	return data;
}

Filament::Filament() : _alive(true)
{
}

Filament::Filament(Point root, VEC ball_centre, double min_length, double resolution, bool is_branch) : 
	_root(root), _ball_centre(ball_centre), _is_branch(is_branch), _min_length(min_length), _resolution(resolution)
{
	_alive = true;
	_front = true;
	_points.push_back(root);
	_new_point = root;
	_descent_points.push_back(VEC(root));
}

void Filament::extend()
{
	if (_alive)
	{
		_points.push_back(_new_point);
		_new_point = _next_new_point;
		//_points.push_back(_descent_points[_descent_points.size() -1]); // CHANGED HERE
	}
}

// Old stepping
void Filament::descent_steps(const GRID &Gx, const GRID &Gy, const GRID &Gz)
{

	VEC last_point = _descent_points[_descent_points.size() - 1];
	double x = Gx.interpolate(last_point);
	double y = Gy.interpolate(last_point);
	double z = (Gx.deps() == 1) ? 0 : Gz(last_point);
	VEC dir(x, y, z);
	dir = dir.normalise();
	VEC next_point = last_point - dir;
	_new_point = next_point;
	if (std::find(_descent_points.begin(), _descent_points.end(), next_point) != _descent_points.end()
		|| _descent_points.size() > 150)
		this->_alive = false;
	_descent_points.push_back(next_point);
}

// Find ALL nodes
void Filament::find_nodes(const GRID &map) // To Do: Parellise
{
	int dim = (map.deps() == 1) ? 2 : 3;

	_new_point_nodes = std::vector<Point>();
	auto new_point_nodes = (_front) ? get_nodes(_new_point, 3, dim) : std::vector<Point>();
	for (auto point : new_point_nodes)
		if (map(point) != ELECTRODE && map(point) != PLASMA) // && _new_point.l2(point) == 1.) // second condition of no diagonals
			_new_point_nodes.push_back(point);

	for (Point point : _points)
	{
		std::vector<Point> nodes = get_nodes(point, 3, dim);
		for (auto node : nodes)
		{
			if ((map.deps() == 1 && node.z() == 0) // in the plane for 2D
				|| map.deps() > 1)					// or 3D
				if (map(node) != -1 && !(node == point)) // inside the map, not the current point
				{
					if (!(std::find(_points.begin(), _points.end(), node) != _points.end()))  // a point
					{
						if (!(std::find(_nodes.begin(), _nodes.end(), node) != _nodes.end())) // already a node
						{
							if (map(node) == AIR)
							{
								_new_point = point;
								_alive = false;
								_nodes.clear();
								return;
							}
							if (map(node) != ELECTRODE && map(node) != BORDER && map(node) != GLASS)
							{
								_nodes.push_back(node);
								//check if node is on wavefront
							}
						}
					}
				}
		}
	}
	
}

// Get viable front nodes 
void Filament::get_front_nodes(const GRID &map)
{
	int dim = (map.deps() == 1) ? 2 : 3;

	_new_point_nodes = std::vector<Point>();
	auto new_point_nodes = (_front) ? get_nodes(_new_point, 3, dim) : std::vector<Point>();
	for (auto point : new_point_nodes)
		if (map(point) != ELECTRODE && map(point) != PLASMA) // && _new_point.l2(point) == 1.) // second condition of no diagonals
		{
			_new_point_nodes.push_back(point);
			if (map(point) == AIR)
				_alive = false;
		}
}

// Find viable nodes for point
std::vector<Point> Filament::get_nodes(Point point, int window_size, int dim)
{
	std::vector<Point> nodes;
	if(dim == 3)
		for (int z = point.z() - int(window_size / 2); z < point.z() + int(window_size / 2) + 1; ++z) //-1 in y to +1 in z
		{
			for (int y = point.y() - int(window_size / 2); y < point.y() + int(window_size / 2) + 1; ++y) //-1 in y to +1 in y
			{
				for (int x = point.x() - int(window_size / 2); x < point.x() + int(window_size / 2) + 1; ++x) // -1 in x to +1 in x
				{
					if (x != y && y != z)
						nodes.push_back(Point(x, y, z));
				}
			}
		}
	else
		for (int y = point.y() - int(window_size / 2); y < point.y() + int(window_size / 2) + 1; ++y) //-1 in y to +1 in y
		{
			for (int x = point.x() - int(window_size / 2); x < point.x() + int(window_size / 2) + 1; ++x) // -1 in x to +1 in x
			{
				if (!(x == point.x() && y == point.y()))
					nodes.push_back(Point(x, y, 0));
			}
		}
		
	return nodes;
}

void Filament::fill_phis(const GRID &map, bool punish_diagonals)
{
	_phis.clear();
	_new_point_phis.clear();
	int dim = (map.deps() == 1) ? 2 : 3;
	for (auto node : _nodes)
	{
		/*
		std::vector<Point> punishment = get_nodes(node, 3, dim);
		double pun_amount = 0;
		for (auto p : punishment)
			if (map(p) == ELECTRODE)
				++pun_amount;
		*/
		_phis.push_back(map(node));
	}

	for (auto node : _new_point_nodes)
	{
		double dist = node.l2(_new_point);
		auto modifier = dist > 1.;
		//if(!(dist > 1.)) // to take adjacents only
		_new_point_phis.push_back(map(node));
	}
}

static bool abs_compare(int a, int b)
{
	return (std::abs(a) < std::abs(b));
}

// To see grad directions
Point get_direction(double val, Point center)
{
	if (-7. * M_PI / 8. < val && val <= -5. * M_PI / 8.)
		return center.ne();
	if (-5. * M_PI / 8. < val && val <= -3. * M_PI / 8.)
		return center.n();
	if (-3. * M_PI / 8. < val && val <= -1. * M_PI / 8.)
		return center.nw();
	if (-1. * M_PI / 8. < val && val <= 1. * M_PI / 8.)
		return center.w();
	if (1. * M_PI / 8. < val && val <= 3. * M_PI / 8.)
		return center.sw();
	if (3. * M_PI / 8. < val && val <= 5. * M_PI / 8.)
		return center.s();
	if (5. * M_PI / 8. < val && val <= 7. * M_PI / 8.)
		return center.se();
	if (7. * M_PI / 8. < val || val <= -7. * M_PI / 8.)
		return center.e();
}

//naive grad stepping
void Filament::grad_find_next(const GRID &grad_dir, const GRID &map, double branching)
{

	get_front_nodes(map);
	Point dir = get_direction(grad_dir(_new_point), _new_point);
	if (std::find(_new_point_nodes.begin(), _new_point_nodes.end(), dir) != _new_point_nodes.end())
	{
		_next_new_point = dir;
		double length = _root.l2(_new_point) * _resolution;
		if (branching > 0.0 && !_is_branch && !_is_parent && length >= _min_length*0.5) // only 1 branch for now
		{
			double current = grad_dir(_new_point);
			double lower = grad_dir(_new_point) - branching;
			double upper = grad_dir(_new_point) + branching;
			Point bcurrent = get_direction(current, _new_point);
			Point blower = get_direction(lower, _new_point);
			Point bupper = get_direction(upper, _new_point);
			if (blower != dir || bupper != dir)
			{
				_found_branch = true;
				_next_new_point = get_direction(grad_dir(_new_point) - 2 * M_PI / 8., _new_point);
				_new_branch_point = get_direction(grad_dir(_new_point) + 2 * M_PI / 8., _new_point);
			}
		}
	}
	else
	{
		if (std::find(_nodes.begin(), _nodes.end(), dir) != _nodes.end() || dir == _new_point)
			_valid = false;
		_alive = false;
	}

}

//naive grad steppingxy
void Filament::grad_find_next(const GRID & Gx, const GRID & Gy, const GRID & Gz, const GRID &map, double branching)
{

	VEC dir(Gx(_new_point), Gy(_new_point), (map.deps() > 1) ? Gz(_new_point) : 0);

	Point dirpoint = get_direction(std::atan2f(dir.y(), dir.x()), _new_point);

	dir = dir.normalise(true);
	_next_new_point = dirpoint; // _new_point - Point(dir.x() * 1.5, dir.y() * 1.5, dir.z() * 1.5);

	auto nodes = get_nodes(_next_new_point, 3, 2); // new filament joining feature
	//for (const auto & n : nodes)
	//	if (map(n) == PLASMA)
	//	{
	//		_next_new_point = n;
	//		_alive = false;
	//	}
	if (_new_point == _next_new_point || map(_next_new_point) != GAS) 
	{
		_alive = false;
		std::cout << "filament killed early" << std::endl;
	}

}

// Grad based stepping
void Filament::grad_find_stepped_next(const GRID & Gx, const GRID & Gy, const GRID & Gz, const GRID &map, const int stepsize, double branching)
{
	bool branch = false;
	double b_check = ((_min_length / _resolution) / 58.);
	double branch_chance = (_new_point.l2(_root) > _branch_pre) ? _new_point.l2(_root) - _branch_pre : 1.0;
	if (_new_point.l2(_root) > (_min_length / _resolution) * 0.6
		&& !_is_parent && !_is_branch
		&& generate_data(1)[0] * branch_chance > 1 - (branching / ((_min_length / _resolution) / 64.))
		) // random at first
	{
		_is_parent = true;
		_found_branch = true;
		Point last_point = _points.back();
		Point direction = _new_point - last_point; // only in 2D
		if (direction.l2() == 1) // (0, 1), (1, 0), (-1, 0), (0, -1)
		{
			_new_point = last_point + direction;
			_new_branch_point = _new_point;
			_new_point = _new_point + Point(direction.y() * -1, direction.x() * -1, 0);
			_new_branch_point = _new_branch_point + Point(direction.y(), direction.x(), 0);
		}
		else
		{
			_new_point = Point(last_point.x() + direction.x(), last_point.y(), 0);
			_new_branch_point = Point(last_point.x(), last_point.y() + direction.y(), 0);
		}
		_descent_points.push_back(_new_point);

	}

	VEC dir(Gx(_new_point), Gy(_new_point), (map.deps() > 1) ? Gz(_new_point) : 0);

	dir = dir.normalise(false);

	VEC pos = _descent_points.back();

	double step = _new_point.l2(_ball_centre);

	double steplength = step + stepsize;
	
	// Calc steps
	while(step < steplength)
	{

		pos = pos - dir;

		Point disc_pos(pos.x() + 0.5 - (pos.x()<0), pos.y() + 0.5 - (pos.y() < 0), pos.z() + 0.5 - (pos.z() < 0));
		
		if (_new_point != disc_pos) // made a step on discrete level
		{

			_next_new_point = disc_pos;
			auto nodes = get_nodes(_next_new_point, 3, 2); // new filament joining feature
			for (const auto & n : nodes)
				if (map(n) == AIR)
				{
					_next_new_point = n;
					extend();
					_new_point = n;
					_alive = false;
					return;
				}

			extend();
			step = _next_new_point.l2(_ball_centre);
		}
	}
	_descent_points.push_back(pos);
}

// old stepping
void Filament::init_find_next(const GRID &grad_dir, const GRID &map)
{
	Point dir = get_direction(grad_dir(_new_point), _new_point);
	get_front_nodes(map);
	if (std::find(_new_point_nodes.begin(), _new_point_nodes.end(), dir) != _new_point_nodes.end())
	{
		_next_new_point = dir;

	}
	else
		_alive = false;
}

// find first steps
void Filament::init_find_next(const GRID &Gx, const GRID &Gy, const GRID &Gz, const GRID &map)
{
	VEC dir(Gx(_new_point), Gy(_new_point), (map.deps() > 1) ? Gz(_new_point) : 0);

	dir = dir.normalise(true);

	Point dirpoint = get_direction(std::atan2f(dir.y(), dir.x()), _new_point);
	_new_point = dirpoint; // _new_point - Point(dir.x() * 1.5, dir.y() * 1.5, dir.z() * 1.5);

	//prep second point as first time

	VEC dir2(Gx(_new_point), Gy(_new_point), (map.deps() > 1) ? Gz(_new_point) : 0);

	dirpoint = get_direction(std::atan2f(dir2.y(), dir2.x()), _new_point);
	dir2 = dir2.normalise(true);
	_next_new_point = dirpoint; // _new_point - Point(dir2.x()* 1.5, dir2.y() * 1.5, dir2.z() * 1.5);

	auto nodes = get_nodes(_next_new_point, 3, 2); // new filament joining feature

	for (const auto & n : nodes)
		if (map(n) == PLASMA)
		{
			_next_new_point = n;
			_alive = false;
		}
	if (_new_point == _next_new_point)
	{
		_alive = false;
		std::cout << "filament killed early" << std::endl;
	}
}

// root attribution
void Filament::grad_find_next_root(GRID grad, GRID grad_dir)
{
	double max = 0;
	for (auto node : _nodes)
	{
		if (grad(node) >= max)
		{
			max = grad(node);
			_new_point = node;
		}
	}
}

// old find next
void Filament::find_next(double eta, bool front, double rand, double branching) // no longer DBM based
{
	//srand(time(NULL));

	std::vector<double> phis = _phis;

	std::vector<int> maxes;

	int max_idx = 0;
	int old_max_idx = 0;
	double max =-10000;
	double old_max = -10000;
	for (int i = 0; i < phis.size(); i++)
	{
		double val = phis[i];
		if (val == max)
		{
			maxes.push_back(i);
			
		}
		else if	(val > max) 
		{ 
			maxes.clear();
			maxes.push_back(i);
			max = val;
		}
		// phis.push_back(val); // normalise between 0 <= phi < 1
	}

	max_idx = maxes[0];

	if (maxes.size() > 1)
		max_idx = maxes[(maxes.size()/2)];

	_new_point = _nodes[max_idx];

	_phis.clear();
	_new_point_phis.clear();
}

Point Filament::get_next_electrode_root_point(int indx)
{
	return _electrode_root_points[indx];
}

void Filament::set_next_electrode_root_point(std::vector<Point> selected)
{
	_electrode_root_points = selected;

}
