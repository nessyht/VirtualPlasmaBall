#include "Plasma_Ball.h"

Plasma_Ball::Plasma_Ball(Geometry geo, double realsize, double resolution, double voltage) :
_dim(geo), _real_size(realsize), _resolution(resolution), _radius(resolution/2.), _Vapp(voltage)
{
	_rootdir = std::filesystem::current_path();
	_buffer = (int)(_radius/2);
	_grid_size = _resolution + _buffer;

	_dxyz = _real_size / _radius;
	_dxyz2 = _dxyz * _dxyz;
	//_dt = 1. / 5.;

	_real_electr_size = 0.02;
	_min_fil_length = _real_size - _real_electr_size;

	_electrode_radius = _real_electr_size * 1 / _dxyz;

}

void Plasma_Ball::init_map()
{
	switch (_dim)
	{
	case CIRCULAR:
		_Map = GRID(_grid_size, _grid_size, 1, BORDER);
		_Map.circleBlock(ELECTRODE, 0, _electrode_radius - 1);
		_Map.circleBlock(ELECTRODE, _electrode_radius - 1, _electrode_radius);
		_Map.circleBlock(GAS, _electrode_radius, _radius);
		_Map.circleBlock(AIR, _radius, _radius + _buffer/2. - 2);
		_electrode_count = _Map.count_instance(ELECTRODE);
		_Map.block(5, 50, 0, 5, 5, 1, CONDUCTOR);
		_Gx = _Map;
		_Gy = _Map;
		break;
	case SPERICAL:
		_Map = GRID(_grid_size, _grid_size, _grid_size, AIR);
		_Map.sphereBlock(ELECTRODE, 0, _electrode_radius);
		_Map.sphereBlock(GAS, _electrode_radius, _radius);
		_Map.sphereBlock(AIR, _radius, _radius + _buffer / 2. - 2);

		_Gx = _Map;
		_Gy = _Map;
		_Gz = _Map;
		break;
	}

}
// BASE
void Plasma_Ball::reset_maps()
{
	_Mask = GRID(_Map);
	_potentials = GRID(_Map, true);
	_charge_map = GRID(_Map, true);
	_potentials.maskBlock(_Mask, ELECTRODE, _Vapp);
	_charge_map = _potentials;

}

void Plasma_Ball::create_det_root_points()
{

	double r = _electrode_radius + 2;
	double rd = (_dim == SPERICAL) ?  r : 0;

	double w = _grid_size / 2.;
	double h = _grid_size / 2.;
	double d = (_dim == SPERICAL) ? _grid_size / 2. : 0.0;

	int zterm = (d + rd > 1.0) ? d + rd : 1.0;

	_centre = Point(w, h, d);

	for (int z = d - rd; z < zterm; ++z)
		for (int y = h - r; y < h + r; ++y)
			for (int x = w - r; x < w + r; ++x)
			{
				Point point(x, y, z);
				double pr = point.l2(_centre);
				if (pr >= _electrode_radius - 1 && pr < _electrode_radius)
					_root_points._points.push_back(Point(x, y, z));
			}

	_root_points.find_nodes(_Mask);

	if (_dim == CIRCULAR)
	{
		int offset = 2;
		Point midcentre = Point(_centre.x(), _root_points._nodes[0].y(), 0); // centre of circle x, y from first node (bottom left)
		std::vector<Point> points = _root_points._nodes;
		std::vector<Point> ordered_points;
		
		ordered_points.push_back(midcentre);

		points.erase(std::find(points.begin(), points.end(), midcentre));

		//Point leftcentre = Point(midcentre.x() - 1, midcentre.y(), 0);
		//Point rightcentre = Point(midcentre.x() + 1, midcentre.y(), 0);


		//ordered_points.push_back(leftcentre);
		//ordered_points.push_back(rightcentre);

		//points.erase(std::find(points.begin(), points.end(), leftcentre));
		//points.erase(std::find(points.begin(), points.end(), rightcentre)); // prep first 3 as base has 2 points, rest 1

		while (points.size() > 0)
		{
			Point candidate;
			int ind = 0;
			double candl2 = 200;

			for (int i = 0; i < points.size(); ++i)
			{

				double l2 = ordered_points[ordered_points.size() - 1].l2(points[i]); // centremid has no more neighbours start at [1]
				if (l2 == 1.0) // sqrt(2) or 1
				{
					ind = i;
					candl2 = l2;
					break;
				}
			}

			ordered_points.push_back(points[ind]);
			points.erase(points.begin() + ind);
		}
		std::vector<Point> selected;
	
		for (int i = 2; i < ordered_points.size(); ++i)
		{
			if (i % ((ordered_points.size() / _number_of_roots) - 1) == 0)
			{
				selected.push_back(ordered_points[i]);
			}
		}
		std::random_device _rd;
		std::mt19937 g(_rd());

		std::shuffle(selected.begin(), selected.end(), g);
		_root_points._ordered_points = ordered_points;
		_root_points.set_next_electrode_root_point(selected);


		//////// testing
		//Point top = _root_points._ordered_points[_root_points._ordered_points.size() / 2];

		//double max_dist = _root_points._ordered_points[0].l2(top);



		//for (auto p : _root_points._ordered_points)
		//{

		//	double root_distance = p.l2(top);
		//	double inv_root_dist = (max_dist - root_distance);
		//	static std::exponential_distribution<double>d(1);
		//	double exp_res = d(g);
		//	std::cout << "Root_dist: " << root_distance << ", exp_res: " << exp_res << ", res: " << (root_distance < (exp_res - 0.5)) << std::endl;

		//}
		//		//////
	}
	else // Spherical
	{

		std::vector<Point> points = _root_points._nodes;
		std::vector<Point> ordered_points;

		int offset = 2;
		ordered_points.push_back(points[offset]);
		points.erase(points.begin() + offset);
		while (points.size() > 0)
		{
			Point candidate;
			int ind = 0;
			double candl2 = 200;

			for (int i = 0; i < points.size(); ++i)
			{

				double l2 = ordered_points[ordered_points.size() - 1].l2(points[i]);
				if (l2 == 1.0) // sqrt(2) or 1
				{
					ind = i;
					candl2 = l2;
					break;
				}
			}

			ordered_points.push_back(points[ind]);
			points.erase(points.begin() + ind);
		}
		std::vector<Point> selected;

		for (int i = 0; i < ordered_points.size(); ++i)
		{
			if (i % (ordered_points.size() / _number_of_roots) == 0)
			{
				selected.push_back(ordered_points[i]);
			}
		}
		std::random_device _rd;
		std::mt19937 g(_rd());

		std::shuffle(selected.begin(), selected.end(), g);
		_root_points._ordered_points = ordered_points;
		_root_points.set_next_electrode_root_point(selected);
	}

}

void Plasma_Ball::select_init_roots(int indx)
{
	Filament new_fil(_root_points.get_next_electrode_root_point(indx), _centre, _min_fil_length, _dxyz);

	//new_fil.descent_steps(_Gx, _Gy, _Gz);


	new_fil.grad_find_stepped_next(_Gx, _Gy, _Gz, _Mask, _stepsize);

	_filaments.push_back(new_fil);



	++_active_filaments;
	++_number_of_filaments;
	++_plasmas;
	//add_charge();

}

void Plasma_Ball::select_root()
{
	if (_number_of_filaments > 0)
	{
		_root_points.fill_phis(_root_charge_map);
		_root_points.find_next(0, false, 0);
		std::vector<Point> rooterase = _root_points.get_nodes(_root_points._new_point, 3, _dim);
		rooterase.push_back(_root_points._new_point);
		_root_charge_map.block(rooterase, 0);
		for (auto p : _root_points._ordered_points)
		{
			if (_Mask(p) != PLASMA)
				_root_points._phis.push_back(sqrt(_Gx(p)*_Gx(p) + _Gy(p)*_Gy(p)));
			else
				auto pp = p;
		}
		--_number_of_filaments;
	}
	else
	{
		_root_points.fill_phis(_potentials);
		_root_points.find_next(0, false, 0);
	}
	//show(_potentials, "root_pots");

	Filament new_fil(_root_points._new_point, _centre, _min_fil_length, _dxyz);

	show(_root_charge_map, "root_charge_map");

	for (const auto &fil : _last_step_filaments)
	{
		if (new_fil._root == fil._shifted_root)
		{
			new_fil._branch_pre = fil._new_branch_point.l2(fil._root);
		}
	}
	new_fil.grad_find_stepped_next(_Gx, _Gy, _Gz, _Mask, 3);

	_filaments.push_back(new_fil);
	++_active_filaments;
	//++_number_of_filaments;
	++_plasmas;

}



void Plasma_Ball::set_voltage(double v)
{
	_Vapp = v;
}

void Plasma_Ball::place_charges()
{
	_last_step_filaments = _filaments;

	std::random_device _rd;
	std::mt19937 g(_rd());
	static std::uniform_real_distribution<float> distribution(0.0, 1.0);

	_charge_map = GRID(_Map, true);

	GRID temp(_Mask, true);
	// TODO: ROOT FIXES HERE
	//temp.maskBlock(_Mask, PLASMA, (-_Vact * 0.01)/ _dxyz2);

	Point top = _root_points._ordered_points[_root_points._ordered_points.size() / 2];

	bool deleted = false;

	for (const auto &fil : _filaments)
	{
		double root_distance = fil._root.l2(top);

		float val = distribution(g);

		float expval = std::abs(std::log(val) / root_distance);

		if (root_distance > _electrode_radius / 2. || expval < 0.5 || deleted) //TODO: dying at top killed here
		{
			if(fil._root.l2(_centre)  < _electrode_radius + 2)
				++_number_of_filaments;
			for (const auto &p : fil._points)
			{

				temp(p) = (-_Vact * _residual * _dxyz2);
			}
		}
		else
		{
			std::cout << "root deleted" << std::endl;
			//if(_stepsize == 5)
				//_debug = true;
			//deleted = true;
			//Point root(_root_points._ordered_points[0]);
			//int left = (fil._root.x() - _centre.x() < 0) ? -1 : 1;
			//for (int i = 0; i < _radius - _electrode_radius - 5; ++i)
			//{
			//	if (i == 0)
			//		std::cout << "deleted filament seeded on centre.x() " << left << std::endl;
			//	temp(root.x() + left, root.y() - i, root.z()) = (-_Vact * _residual* _dxyz2);
			//}
		}
	}
	_charge_map = temp;
	vert_drift();

	show(_charge_map, "pre_smoothed");
	bool _smooth_charge = true;
	if (_smooth_charge)
	{
		if (_dim == CIRCULAR)
			for (int i = 0; i < 1; ++i)
			{
				GRID pre_charge(_charge_map);
				_charge_map.maskBlock(_Mask, GAS, [pre_charge, this](float p, int i, int j, int k)
				{
					return (1. / 256.)*((pre_charge(i - 2, j - 2, 0) + pre_charge(i + 2, j + 2, 0) + pre_charge(i - 2, j + 2, 0) + pre_charge(i + 2, j - 2, 0)) +
						4.*(pre_charge(i - 1, j - 2, 0) + pre_charge(i + 1, j + 2, 0) + pre_charge(i - 1, j + 2, 0) + pre_charge(i + 1, j - 2, 0) +
							pre_charge(i - 2, j - 1, 0) + pre_charge(i + 2, j + 1, 0) + pre_charge(i - 2, j + 1, 0) + pre_charge(i + 2, j - 1, 0)) +
						6.*(pre_charge(i, j - 2, 0) + pre_charge(i, j + 2, 0) + pre_charge(i - 2, j, 0) + pre_charge(i + 2, j, 0)) +
						16.*(pre_charge(i - 1, j - 1, 0) + pre_charge(i + 1, j + 1, 0) + pre_charge(i - 1, j + 1, 0) + pre_charge(i + 1, j - 1, 0)) +
						24.*(pre_charge(i - 1, j, 0) + pre_charge(i + 1, j, 0) + pre_charge(i, j + 1, 0) + pre_charge(i, j - 1, 0)) +
						36.*(pre_charge(i, j, 0)));
				});
				_charge_map.maskBlock(_Mask, PLASMA, [pre_charge, this](float p, int i, int j, int k)
				{
					return (1. / 256.)*((pre_charge(i - 2, j - 2, 0) + pre_charge(i + 2, j + 2, 0) + pre_charge(i - 2, j + 2, 0) + pre_charge(i + 2, j - 2, 0)) +
						4.*(pre_charge(i - 1, j - 2, 0) + pre_charge(i + 1, j + 2, 0) + pre_charge(i - 1, j + 2, 0) + pre_charge(i + 1, j - 2, 0) +
							pre_charge(i - 2, j - 1, 0) + pre_charge(i + 2, j + 1, 0) + pre_charge(i - 2, j + 1, 0) + pre_charge(i + 2, j - 1, 0)) +
						6.*(pre_charge(i, j - 2, 0) + pre_charge(i, j + 2, 0) + pre_charge(i - 2, j, 0) + pre_charge(i + 2, j, 0)) +
						16.*(pre_charge(i - 1, j - 1, 0) + pre_charge(i + 1, j + 1, 0) + pre_charge(i - 1, j + 1, 0) + pre_charge(i + 1, j - 1, 0)) +
						24.*(pre_charge(i - 1, j, 0) + pre_charge(i + 1, j, 0) + pre_charge(i, j + 1, 0) + pre_charge(i, j - 1, 0)) +
						36.*(pre_charge(i, j, 0)));
				});
			}
		else
			std::cout << "No 3D smoothing in place\n";
	}
	_root_charge_map = _charge_map;
	
	show(_charge_map, "_smooth_charge_map");
	/* // add past maps
	for (int i = _potentials_maps.size() - 1; i >= 0; --i)
	{
		GRID temp(_Mask, true);
		temp.maskBlock(_Mask, GAS, 0.0);
		temp.maskBlock(_Mask, PLASMA, _Vact / _number_of_roots);
		temp.block(1, 1, 0, temp.cols() - 1, temp.cols() - 1, 1, [threshold](float val, int i, int j, int k) {return (std::abs(val) > threshold) ? val : 0; });
		_charge_map -= temp * multiplier;
		multiplier *= 0.5;
	}
	*/
}

void Plasma_Ball::vert_drift()
{
	bool incremental_shift = false;
	double increment = 0.1;
	double counter = 0;

	show(_charge_map, "pre_drift");
	GRID temp(_charge_map, true);
	double velocity = 1; //TODO: use velocity to kill drift here
	_shift += _dt * velocity; 
	
	double shift = _shift;
	if (incremental_shift) {
		std::cout << "softshifting" << std::endl;
		temp.maskBlock(_Mask, PLASMA, [&counter, increment, this](float pot, int i, int j, int k)
		{
			Point south = Point(i, j, k).s();
			if (_Mask(south) == GAS || _Mask(south) == PLASMA)
			{
				counter += increment;
				if (counter == 1.)
				{
					counter = 0.0;
					return (float)_charge_map(south);
				}
				else
					return (float)_charge_map(i, j, k);
			}
		});
		temp.maskBlock(_Mask, GAS, [&counter, increment, this](float pot, int i, int j, int k)
		{
			Point south = Point(i, j, k).s();
			if (_Mask(south) == GAS || _Mask(south) == PLASMA)
			{
				counter += increment;
				if (counter == 1.)
				{
					counter = 0.0;
					return (float)_charge_map(south);
				}
				else
					return (float)_charge_map(i, j, k);
			}
		});
		show(temp, "mid_drift");
		std::cout << "softshiftingsuccess" << std::endl;
	}
	if (shift >= 1.0)
	{
		std::cout << "shifting" << std::endl;
		while (shift >= 1.0)
		{

			temp.maskBlock(_Mask, PLASMA, [this](float pot, int i, int j, int k)
			{
				Point south = Point(i, j, k).s();
				if (_Mask(south) == GAS || _Mask(south) == PLASMA)
					return (float)(_charge_map(south));
			});
			temp.maskBlock(_Mask, GAS, [this](float pot, int i, int j, int k)
			{
				Point south = Point(i, j, k).s();
				if (_Mask(south) == GAS || _Mask(south) == PLASMA)
					return (float)(_charge_map(south));
			});
			show(temp, "mid_drift");
			temp.block(_root_points._ordered_points, [this](float pot, Point p)
			{
				int size = _root_points._ordered_points.size();
				Point centremid(_centre.x(), _root_points._ordered_points[size / 2].y(), 0);

				for (int a = 0; a < size; ++a)
				{
					if (_root_points._ordered_points[a] == p)
					{
						if (_root_points._ordered_points[a] == centremid)
							return 0.f;
						if (a == size - 1) // bot centre left
						{
							return (float)_charge_map(_root_points._ordered_points[0] * 0.45);
						}
						if (a == 0)
						{
							return (float)_charge_map(_root_points._ordered_points[0] * 0.1);
						}
						if (a == 1)
						{
							return (float)_charge_map(_root_points._ordered_points[0] * 0.45);
						}
						if (a < (size / 2.))
						{
							return (float)_charge_map(_root_points._ordered_points[a - 1])*10;
						}
						else
						{
							return (float)_charge_map(_root_points._ordered_points[a + 1])*10;

						}
					}
				}
			});
			int size = _root_points._ordered_points.size();
			Point centremid(_centre.x(), _root_points._ordered_points[size / 2].y(), 0);// at the top
			for (auto &fil : _last_step_filaments)
			{
				auto root = fil._root;

				for (int a = 0; a < size; ++a)
				{
					if (_root_points._ordered_points[a] == root)
					{
						if (_root_points._ordered_points[a] == centremid)
						{
							fil._shifted_root = root;
							break;
						}
						if (a == 0) // special case
						{
							fil._shifted_root = _root_points._ordered_points[a + 1]; // TODO: fix side
							break;

						}
						if (a < (size / 2.))
						{
							fil._shifted_root = _root_points._ordered_points[a + 1];
							break;

						}
						else
						{
							fil._shifted_root = _root_points._ordered_points[a - 1];
							break;

						}
					}
				}
			}

			shift--; //step completed
			_shift--;
		}
	}
	else
	{
		for (auto &fil : _last_step_filaments)
		{
			fil._shifted_root = fil._root;
		}
		temp = _charge_map;
	}
	//else if (_shift > 0.5 && _shift < 1.0)
	//{
	//	temp.maskBlock(_Mask, PLASMA, [this](float pot, int i, int j, int k)
	//	{
	//		Point south = Point(i, j, k).s();
	//		if (_Mask(south) == GAS || _Mask(south) == PLASMA)
	//			return (float)(_charge_map(south) * 0.5 + _charge_map(i, j, k) * 0.5);
	//	});
	//	temp.maskBlock(_Mask, GAS, [this](float pot, int i, int j, int k)
	//	{
	//		Point south = Point(i, j, k).s();
	//		if (_Mask(south) == GAS || _Mask(south) == PLASMA)
	//			return (float)(_charge_map(south) * 0.5 + _charge_map(i, j, k) * 0.5);
	//	});
	//}


	if (false && _dim == CIRCULAR)
	{
		int min_y = 10000;
		int max_y = 0;

		int mid_bot = 0;
		int mid_top = 0;
		for (auto p : _root_points._ordered_points)
		{
			std::vector<Point> botrow;
			std::vector<Point> toprow;
			if (p.y() <= min_y)
			{
				if (p.y() == min_y)
					botrow.push_back(p);
				if (p.y() < min_y)
				{
					min_y = p.y();
					botrow.clear();
					botrow.push_back(p);
				}
			}

		}
	}
	_charge_map = temp;
	_root_charge_map = temp;
	show(_charge_map, "post_drift");
}

void Plasma_Ball::set_updating(bool b)
{
	_updating = b;
}

void Plasma_Ball::set_grad_descent(bool b)
{
	_grad_descent = b;
}

void Plasma_Ball::init_DBM(unsigned fils, int stepsize, double branching, int fps, double residual)
{
	_randomness = 0;
	_number_of_filaments = 0;
	_number_of_roots = fils;
	_eta = 2.0;
	_serial = false;
	_updating = true;
	_front = true;
	_dt = 1. / (float)fps;
	_residual = residual;
	_branching = branching;
	_stepsize = stepsize;

	init_map();

	reset_maps();

	create_det_root_points();

	create_dir();
}


void Plasma_Ball::get_gradient()
{
	_gradients = GRID(_potentials, true);
	_gradient_directions = _gradients;

	GRID Gx = _gradients;
	GRID Gy = _gradients;
	
	Gy.maskBlock(_Mask, GAS, [this](float current, int i, int j, int k)
	{
		return  _potentials(i - 1, j - 1, k) + 2 * _potentials(i, j - 1, k) + _potentials(i + 1, j - 1, k)
			- _potentials(i - 1, j + 1, k) - 2 * _potentials(i, j + 1, k) - _potentials(i - 1, j + 1, k);
	}
	);
	Gx.maskBlock(_Mask, GAS, [this](float current, int i, int j, int k)
	{
		return  _potentials(i - 1, j - 1, k) + 2 * _potentials(i - 1, j, k) + _potentials(i - 1, j + 1, k)
			- _potentials(i + 1, j - 1, k) - 2 * _potentials(i + 1, j, k) - _potentials(i + 1, j + 1, k);
	}
	);
	_gradients = Gx.power(2) * Gy.power(2);
	_gradients.sqrt();
	_gradient_directions.maskBlock(_Mask, GAS, [Gx, Gy](float c, int i, int j, int k) {return std::atan2f(Gy(i, j, k), Gx(i, j, k)); });
	if (_debug)
	{
		GRID directions = _gradient_directions;

		directions.maskBlock(_Mask, GAS, [](float val, int i, int j, int k)
		{

			if (-7. * M_PI / 8. < val && val <= -5. * M_PI / 8.)
				return 8;
			if (-5. * M_PI / 8. < val && val <= -3. * M_PI / 8.)
				return 1;
			if (-3. * M_PI / 8. < val && val <= -1. * M_PI / 8.)
				return 2;
			if (-1. * M_PI / 8. < val && val <= 1. * M_PI / 8.)
				return 3;
			if (1. * M_PI / 8. < val && val <= 3. * M_PI / 8.)
				return 4;
			if (3. * M_PI / 8. < val && val <= 5. * M_PI / 8.)
				return 5;
			if (5. * M_PI / 8. < val && val <= 7. * M_PI / 8.)
				return 6;
			if (7. * M_PI / 8. < val || val <= -7. * M_PI / 8.)
				return 7;
		});
		show(directions, "directions");
	}
	

	get_gradient_directions();

	_Gx = Gx;
	_Gy = Gy;
	_Gx.maskBlock(_Mask, GAS, [this](float current, int i, int j, int k)
	{
		return (_potentials(i - 1, j, k) - _potentials(i + 1, j, k));
	});
	_Gy.maskBlock(_Mask, GAS, [this](float current, int i, int j, int k)
	{
		return (_potentials(i, j - 1, k) - _potentials(i, j + 1, k));
	});
	
	show(Gx, "Gx");
	show(_Gx, "Gx2");
	show(Gy, "Gy");
	show(_Gy, "Gy2");
	show(_potentials, "pots");

	if (_dim == SPERICAL)
	{
		_Gz = _gradients;
		_Gz.maskBlock(_Mask, GAS, [this](float current, int i, int j, int k)
		{
			return (_potentials(i, j, k - 1) - _potentials(i, j, k + 1));
		});
	}



	if (_dim == SPERICAL)
	{
		Gx.maskBlock(_Mask, GAS, [this](float current, int i, int j, int k)
		{
			return (_potentials(i, j, k - 1) - _potentials(i, j, k + 1));
		});
		
	}
}

void Plasma_Ball::get_gradient_directions()
{
	_gradients = _Gx.power(2) * _Gy.power(2);
	_gradients.sqrt();
	_gradient_directions = _gradients;
	_gradient_directions.maskBlock(_Mask, GAS, [this](float c, int i, int j, int k) {return std::atan2f(_Gy(i, j, k), _Gx(i, j, k)); });
	if (_debug)
	{
		GRID directions = _gradient_directions;

		directions.maskBlock(_Mask, GAS, [](float val, int i, int j, int k)
		{

			if (-7. * M_PI / 8. < val && val <= -5. * M_PI / 8.)
				return 8;
			if (-5. * M_PI / 8. < val && val <= -3. * M_PI / 8.)
				return 1;
			if (-3. * M_PI / 8. < val && val <= -1. * M_PI / 8.)
				return 2;
			if (-1. * M_PI / 8. < val && val <= 1. * M_PI / 8.)
				return 3;
			if (1. * M_PI / 8. < val && val <= 3. * M_PI / 8.)
				return 4;
			if (3. * M_PI / 8. < val && val <= 5. * M_PI / 8.)
				return 5;
			if (5. * M_PI / 8. < val && val <= 7. * M_PI / 8.)
				return 6;
			if (7. * M_PI / 8. < val || val <= -7. * M_PI / 8.)
				return 7;
		});
		auto gx = _Gx * _Gx;
		auto gy = _Gy * _Gy;
		GRID e_field = gx + gy;
		e_field = e_field.sqrt() * -1;
		for (int i = 0; i < _grid_size; ++i)
		{
			_logfile << "V: " << _potentials(i, _centre.y(), 0) << std::endl;
			_logfile << "E: " << e_field(i, _centre.y(), 0) << std::endl;
		}
		show(e_field, "e_field");
		show(directions, "directions");
	}
}

void Plasma_Ball::clean_maps()
{
	// TODO: adjust for charge on outside of ball
	_Mask = _Map;
	if (_potentials_maps.size() == 0)
	{
		show(_Map, "_Map_init_");
		show(_Mask, "_Mask_init_");
		//_potentials / _dxyz2;
		//_potentials.block(0, 0, 0, _grid_size, _grid_size, 1, [this](float current, int i, int j , int k)
		//{
		//	if (_Mask(i, j, k) == GAS || _Mask(i, j, k) == AIR)
		//	{
		//		return _Vact * (_radius - _center.l2(Point(i, j, k))) / _radius;
		//	}
		//	else if (_Mask(i, j, k) == ELECTRODE)
		//		return _Vact / _dxyz2;
		//	else
		//		return 0.0;

		//});
		_potentials = _potentials / _dxyz2;
		_potentials.maskBlock(_Mask, GAS, 0);
		_potentials.maskBlock(_Mask, AIR, 0);
		
		_potentials_maps.push_back(_potentials);
		//GRID _potentials_indexed = _potentials_maps[0];
		//int index = 0;
		//_potentials_indexed.block(0, 0, 0, _grid_size, _grid_size, 1, [&index](float current, int i, int j, int k)
		//{
		//	if (current == GAS || current == AIR)
		//	{
		//		int val = index;
		//		++index;
		//		return val;
		//	}
		//	else
		//		return -1;
		//});

		show(_potentials_maps[0], "_boundary_conditions_init_");
		Poisson(_radius * 2);

		//_potentials_maps.push_back(_potentials);
	}
	else
	{
		_potentials = _potentials_maps[0];
		show(_potentials, "_potentials_init_");

		//_potentials = _potentials/_dxyz2;
		//show(_potentials, "_potentials_init_div_dxyz2");

		_potentials = _potentials + _charge_map;
		show(_potentials, "_potentials_init_w_charge_map");

		_potentials_maps.push_back(_potentials);
	}
	show(_potentials, "_potentials_init_");

	_filaments.clear();
	_active_filaments = 0;
}



void Plasma_Ball::step(bool debug)
{	
	_step++;

	_logfile << "Step: " << _step << std::endl;
	_logfile << "Degug: " << std::to_string(debug) << std::endl;

	auto begin = std::chrono::high_resolution_clock::now();

	std::cout << "Performing step: " << _step << std::endl;
	_debug = debug;
	_active_filaments = 0;
	_plasmas = 0;
	_Vact = _Vapp;
	_number_of_filaments = 0;
	_fil_step = 0;



	if (_step > 1)
	{
		place_charges();
	}

   	clean_maps(); // previous filaments cleared here

	show(_charge_map, "_charges_");

	// execute filament by filament
	if(_serial)
		for (int f = 0; f < _number_of_roots; ++f)
		{
			_potentials = _potentials_maps[0];
			_potentials.maskBlock(_Mask, PLASMA, _Vact);
			_alive = true;
			//Poisson(2 *_radius);
			show(_potentials, "root_pots");
			//_root_points.fill_phis(_root_charge_map);
			if (_step > 1)
				select_root();
			else
				select_init_roots(f);
			show(_Mask, "mask_with_roots");
			if (_grad_descent)
			{
				perform_grad_descent();
			}
			else
				perform_DBM();


		}
	//execute filaments in parallel
	else
	{
		Poisson();
		for (int f = 0; f < _number_of_roots; ++f)
		{
			_alive = true;

			if (_step == 1) // first step deterministic roots
			{		
				select_init_roots(f);
				add_charge();
			}
			else
			{			
				Poisson();
				select_root();
			}
		}
		_number_of_filaments = _number_of_roots;
		show(_potentials, "pre_filaments_");
		_fil_step += 1;
		perform_grad_descent();

	}

	finish_step();
	auto end = std::chrono::high_resolution_clock::now();
	auto time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

	_logfile << "Time taken to executute step: " << _step << " : " << std::to_string(time/(1e9)) << "s" << std::endl;

}

void Plasma_Ball::finish_step()
{
	add_charge();
	show(_Mask, "Leaders_");

	_leaders.push_back(_Mask);
	if (_dim == CIRCULAR)
	{
		GRID vol_leader(_grid_size, _grid_size, 5, 0.0);
		GRID gauss(_potentials_maps[0], true);
		gauss.maskBlock(_Mask, PLASMA, _Vact);
		show(gauss, "Leaders_");

		for (int i = 0; i < 1; ++i)
		{
	
			//add_charge(_Vact);
			GRID gauss_temp(gauss);
			gauss.maskBlock(_Map, GAS, [gauss_temp, this](float p, int i, int j, int k)
			{
				return (1. / 256.)*((gauss_temp(i - 2, j - 2, 0) + gauss_temp(i + 2, j + 2, 0) + gauss_temp(i - 2, j + 2, 0) + gauss_temp(i + 2, j - 2, 0)) +
					4.*(gauss_temp(i - 1, j - 2, 0) + gauss_temp(i + 1, j + 2, 0) + gauss_temp(i - 1, j + 2, 0) + gauss_temp(i + 1, j - 2, 0) +
						gauss_temp(i - 2, j - 1, 0) + gauss_temp(i + 2, j + 1, 0) + gauss_temp(i - 2, j + 1, 0) + gauss_temp(i + 2, j - 1, 0)) +
					6.*(gauss_temp(i, j - 2, 0) + gauss_temp(i, j + 2, 0) + gauss_temp(i - 2, j, 0) + gauss_temp(i + 2, j, 0)) +
					16.*(gauss_temp(i - 1, j - 1, 0) + gauss_temp(i + 1, j + 1, 0) + gauss_temp(i - 1, j + 1, 0) + gauss_temp(i + 1, j - 1, 0)) +
					24.*(gauss_temp(i - 1, j, 0) + gauss_temp(i + 1, j, 0) + gauss_temp(i, j + 1, 0) + gauss_temp(i, j - 1, 0)) +
					36.*(gauss_temp(i, j, 0)));
			});
			gauss_temp = gauss;
		}
		//gauss.maskBlock(_Mask, ELECTRODE, _Vapp);
		show(gauss, "Leaders_");
		vol_leader.fill2Dto3D(gauss);
		show(vol_leader, "Leaders_vol_");

	}
}

void Plasma_Ball::perform_grad_descent()
{

	if (false && !_updating)
	{
		Poisson(2 * _radius);
		if (_grad_descent)
		{
			show(_gradients, "gradients");
			show(_gradient_directions, "grad_dir");
		}
	}
	_alive = true;

	while (_alive)
	{
		std::cout << "Filament step: " << _fil_step << ", of time step " << _step << std::endl;
		if (_plasmas > _electrode_count)
			_debug = false;
		int active = 0;
		_active_filaments = _number_of_filaments;

		for (int i = 0; i < _number_of_filaments; ++i)
		{
			Filament &fil = _filaments[i];
			bool step_occured = false;
			if (fil._alive)
			{
				++active;

				if (_stepsize == 0)
				{
					fil.extend();
					++_plasmas;

					fil.grad_find_next(_Gx, _Gy, _Gz, _Mask, _branching);
					if (fil._alive == false)
						break;
				}
				else if (fil._new_point.l2(_centre) - _electrode_radius < _fil_step + 2 * _stepsize) 
					// only if fil is not overgrown
				{
					step_occured = true;
					fil.grad_find_stepped_next(_Gx, _Gy, _Gz, _Mask, _stepsize, _branching);
					if (fil._found_branch)
					{

						std::cout << "Filament " << i << " has spawned a branch\n";
						fil._found_branch = false;
						Filament branch_fil(Filament(fil._new_branch_point, _centre, _min_fil_length, _dxyz, true));
						branch_fil.grad_find_stepped_next(_Gx, _Gy, _Gz, _Mask, _stepsize, _branching);
						_filaments.push_back(branch_fil);

						_number_of_filaments++;
						_active_filaments++;
					}
				}
				else
					std::cout << "fil " << i << " paused in step " << _fil_step << std::endl;
				//fil.descent_steps(_Gx, _Gy, _Gz);
				if (_Mask(_filaments[i]._new_point) == ELECTRODE || _Mask(_filaments[i]._new_point) == AIR)
				{
					std::cout << "fil " << i << " terminated in step " << _fil_step <<
						" because new_point was " << Material((int)_Mask(_filaments[i]._new_point)) << std::endl;

					fil._alive = false;
				}
				if (fil._alive == false && fil._new_point.l2(_centre) < _radius - _electrode_radius - 3)
					std::cout << "fil " << i << " died in step " << _fil_step 
					<< "at length " <<  fil._new_point.l2(fil._root) << std::endl;
			}
			if (_updating && step_occured) // only if step actually occured
			{
				Poisson();
				show(_potentials, "pot_steps");
			}
			if (!fil._alive)
				--_active_filaments;
		}
		_fil_step += _stepsize; // total steps each filament has taken
		if (!(_active_filaments > 0) || active == 0)
			_alive = false;
	}
	/*for (auto fil : _filaments)
	{
		for (auto v : fil._descent_points)
		{
			Point pointver(v.x(), v.y(), v.z());
			_Mask(pointver) == PLASMA;
			_potentials(pointver) == _Vact;
		}
		show(_Mask, "filamement");
	}*/
}

void Plasma_Ball::perform_DBM()
{
	if (!_updating) // if poisson not updating per step, converge once before
	{
		Poisson(2 * _radius);
		if (_grad_descent)
		{
			show(_gradients, "gradients");
			show(_gradient_directions, "grad_dir");
		}
	}
	int counter = 0;
	while (_alive)
	{
		if (_updating && !_grad_descent) // only if stepwise poisson solve
		{
			Poisson(_radius * 0.25);
		}
			/*
		GRID next_map = _potentials + _charge_map;

		show(next_map, "_next_map");
		// */
		int active = 0;
		for (auto &fil : _filaments)
		{
			if (fil._alive)
			{

				fil.find_nodes(_Mask);
				fil.fill_phis(_potentials);
				if (!fil._alive)
				{
					--_active_filaments;
					continue;
				}
				///
				// fil.fill_phis(next_map);
				///

				if (_grad_descent)
					fil.grad_find_next(_gradient_directions, _Mask);
				else
					fil.find_next(_eta, _front, _randomness, _branching);
				fil.extend();
				std::cout << fil._new_point << std::endl;
				if (_branching != 0.0) // only if branching is active
					if (fil._new_branch_point >= Point(0, 0, 0) && !fil._is_branch)
					{
						_filaments.push_back(Filament(fil._new_branch_point, _min_fil_length, _dxyz, true));
						++_active_filaments;
						add_charge();
					}
				_potentials(fil._new_point) = _Vact;
				_Mask(fil._new_point) = PLASMA;
			}
			if (fil._alive)
				++active;
		}
		if (!(active > 0))
			_alive = false;
		show(_potentials, "pot_steps");
		if (_updating && _grad_descent) // only if stepwise poisson solve
		{
			Poisson(_radius * 0.25);
			if (_grad_descent)
			{
				show(_gradients, "gradients");
				show(_gradient_directions, "grad_dir");
			}
		}
		++counter;
	}
	
}

void Plasma_Ball::add_charge(double val)
{
	//double v = (val != 0) ? val : _Vact * 10;// (_grid_size * _Vact * _electrode_count) / _plasmas;
	for (auto &Fil : _filaments)
	{
		//if (Fil._alive)
		//{
			for (const auto &p : Fil._points)
			{
				_potentials(p) = _Vact / _dxyz2; // *p.l2(Point(_center));
				_Mask(p) = PLASMA;
			}
		//}
	}

}

bool Plasma_Ball::find_living_filament()
{
	for (auto &Fil : _filaments)
		if (Fil._alive)
			return true;
	return false;
}

GRID Plasma_Ball::get_map()
{
	return _Map;
}

GRID Plasma_Ball::get_mask()
{
	return _Mask;
}

GRID Plasma_Ball::get_potentials()
{
	return _potentials;
}

extern "C" void perform_conjugate_gradient_solve(const int m, float *res);

void Plasma_Ball::Poisson(int its)
{
	bool _cuda = true;
	if (_cuda) //GPU FEM
	{
		if (_dim == CIRCULAR)
		{
			if (_potentials_maps.size() > 0)
			{
				_potentials = _potentials_maps[_step-1];
				add_charge();
				//test
				//_Mask.block(5, 50, 0, 10, 55, 1, CONDUCTOR); 
				//_potentials.maskBlock(_Mask, CONDUCTOR, 100000/_dxyz2);
				//show(_Mask, "external_map");
			}
			GRID Mask(_Mask);
			if (_fil_step - _stepsize - 1 > 0) // block off done filaments to reduce interference
			{
				int electrode_layer = _fil_step - _stepsize * 0.5;
				Mask.circleBlock(FAKEELECTRODE, _electrode_radius, _electrode_radius + electrode_layer);
				//double sum = _electrode_count * _Vact; //inside band
				//double count = _electrode_count;
				//_potentials.maskBlock(Mask, FAKEELECTRODE, [&sum, &count](float current, int i, int j, int k)
				//{++count; sum += current; return current; });
				float vact = (_Vact / _dxyz2);
				//vact = (_Vact * (_electrode_radius/(_electrode_radius + electrode_layer)))/_dxyz2;
				_potentials.maskBlock(Mask, FAKEELECTRODE, vact); // keep potential accurate? sum/count?
				show(_potentials, "pre_pots_masked"); // = _potentials; // / (float)_dxyz;
				show(Mask, "fakeelectrode_Mask"); // = _potentials; // / (float)_dxyz;

			}
			//_potentials = GRID(_Mask.cols(), _Mask.rows(), 1, -1.602e-7 / 8.85418782);
			show(_potentials, "pre_pots"); // = _potentials; // / (float)_dxyz;

			const int m = _grid_size;

			float *res = (float*)malloc(m * m * sizeof(float));


			//perform_conjugate_gradient_solve(m, res);
			cublasHandle_t cublas_handle;
			cusparseHandle_t cusparse_handle;
			cublasCreate(&cublas_handle);
			cusparseCreate(&cusparse_handle);

			size_t free_mem = 0;
			size_t total_mem = 0;
			int cudaStat1 = cudaMemGetInfo(&free_mem, &total_mem);
			assert(cudaSuccess == cudaStat1);
			//std::cout << "Free memory: " << free_mem << " bytes, Total memory: " << total_mem << " bytes\n";


			//const int    m = 1028;
			const float		alpha = 0.01;
			int				max_iter = 1000; // 5 * _stepsize; // m * 0.6;
			float			tolerance = 0.0000000001;
			auto copystart = std::chrono::high_resolution_clock::now();

			float *b, *mask, *x, *xprime, *yprime;
			cudaMalloc((void **)&b, (m*m) * sizeof(float));
			cudaMalloc((void **)&mask, (m*m) * sizeof(float));
			cudaMalloc((void **)&x, (m*m) * sizeof(float));
			cudaMalloc((void **)&xprime, (m*m) * sizeof(float));
			cudaMalloc((void **)&yprime, (m*m) * sizeof(float));
			//device_memset<float>(b, 1.0, m*m);

			cudaMemcpy(mask, Mask.data(), (m*m) * sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(b, _potentials.data(), (m*m) * sizeof(float), cudaMemcpyHostToDevice);

			int num_iter = 0;
			const float dxyz = _dxyz;

			auto kernelstart = std::chrono::high_resolution_clock::now();

			num_iter = solve_with_conjugate_gradient<float>(cublas_handle, cusparse_handle,
				m, alpha, dxyz, b, mask, x, xprime, yprime, max_iter, tolerance, NULL);
			//std::cout << num_iter << " iterations" << std::endl;

			auto kernelend = std::chrono::high_resolution_clock::now();

			auto kerneltime = std::chrono::duration_cast<std::chrono::nanoseconds>(kernelend - kernelstart).count();

			_logfile << "Kernel time: " << std::to_string(kerneltime / (1e9)) << "s" << std::endl;


			cudaMemcpy(res, x, (m*m) * sizeof(float), cudaMemcpyDeviceToHost);
			_potentials.GRID(res);
			cudaMemcpy(res, xprime, (m*m) * sizeof(float), cudaMemcpyDeviceToHost);
			_Gx.GRID(res);
			cudaMemcpy(res, yprime, (m*m) * sizeof(float), cudaMemcpyDeviceToHost);
			_Gy.GRID(res);

			show(_potentials, "pre_norm_pots");
			//_potentials.mask_clamp(_Mask, GAS, _Vact, 0);
			//_potentials.normalise(_Vact);
			//if (_debug)
			//{
			//	if (_potentials_maps.size() == 1)
			//		_potentials_maps.push_back(_potentials);
			//	_potentials += _potentials_maps[1];
			//}
			//show(_potentials, "pots");
			show(_Gx, "gx");
			show(_Gy, "gy");

			cudaFree(mask);
			cudaFree(b);
			cudaFree(x);
			cudaFree(xprime);
			cudaFree(yprime);
			cublasDestroy(cublas_handle);
			cusparseDestroy(cusparse_handle);

			free(res);

			auto copyend = std::chrono::high_resolution_clock::now();

			auto copytime = std::chrono::duration_cast<std::chrono::nanoseconds>(copyend - copystart).count();

			_logfile << "Copy time: " << std::to_string((copytime-kerneltime) / (1e9)) << "s" << std::endl;

			//get_gradient();
			if(_debug)
				get_gradient_directions();
		} 
		else // SPHERICAL - N/A
		{
			const int m = _grid_size;
			const int mmm = m * m * m;
			float *res = (float*)malloc(m * m * m * sizeof(float));

			/*
					for (int i = 0; i < _potentials.size(); ++i)
						mask[i] = _Mask(i);
			*/

			//perform_conjugate_gradient_solve(m, res);
			cublasHandle_t cublas_handle;
			cusparseHandle_t cusparse_handle;
			cublasCreate(&cublas_handle);
			cusparseCreate(&cusparse_handle);

			//const int		m = 1028;
			const float		alpha = _dxyz2;
			int				max_iter = 100;
			float			tolerance = 0.0000001;

			float *b, *mask, *x, *xprime, *yprime, *zprime;
			cudaMalloc((void **)&b, (mmm) * sizeof(float));
			cudaMalloc((void **)&mask, (mmm) * sizeof(float));
			cudaMalloc((void **)&x, (mmm) * sizeof(float));
			cudaMalloc((void **)&xprime, (mmm) * sizeof(float));
			cudaMalloc((void **)&yprime, (mmm) * sizeof(float));
			cudaMalloc((void **)&zprime, (mmm) * sizeof(float));
			//device_memset<float>(b, 1.0, m*m);

			cudaMemcpy(mask, _Mask.data(), (mmm) * sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(b, _potentials.data(), (mmm) * sizeof(float), cudaMemcpyHostToDevice);

			int num_iter = 0;
			num_iter = solve_with_conjugate_gradient3D<float>(cublas_handle, cusparse_handle,
				m, alpha, b, mask, x, xprime, yprime, zprime, max_iter, tolerance, NULL);
			// std::cout << num_iter << " iterations" << std::endl;

			cudaMemcpy(res, x, (mmm) * sizeof(float), cudaMemcpyDeviceToHost);
			_potentials.GRID(res);
			cudaMemcpy(res, xprime, (mmm) * sizeof(float), cudaMemcpyDeviceToHost);
			_Gx.GRID(res);
			cudaMemcpy(res, yprime, (mmm) * sizeof(float), cudaMemcpyDeviceToHost);
			_Gy.GRID(res);
			cudaMemcpy(res, zprime, (mmm) * sizeof(float), cudaMemcpyDeviceToHost);
			_Gz.GRID(res);

			show(_potentials, "pots");
			show(_Gx, "gx");
			show(_Gy, "gy");
			show(_Gy, "gz");

		}

	}
	else // CPU FD
	{
		_convergence = 1;
		double convold = 0;
		for (int i = 0; i < its; ++i)
		{
			GRID pre_pots(_potentials);
			if (_dim == CIRCULAR)
			{
				_potentials.maskBlock(_Mask, GAS, [pre_pots, this](float p, int i, int j, int k)
				{
					double res = (pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) / 4. - _charge_map(i, j, k) - pre_pots(i, j, k);
					return pre_pots(i, j, k) + res;
					//return (_dxyz2*(pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) - _charge_map(i, j, k) * _dxyz2 * _dxyz2) / (4 * (_dxyz2));
					//return pre_pots(i, j, k) + ((2 * (pre_pots(i - 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i + 1, j, k) + pre_pots(i, j + 1, k))
					//	+ (pre_pots(i - 1, j - 1, k) + pre_pots(i + 1, j - 1, k) + pre_pots(i + 1, j + 1, k) + pre_pots(i - 1, j + 1, k))) / 12. - pre_pots(i, j, k));
				});
				_potentials.maskBlock(_Mask, AIR, [pre_pots, this](float p, int i, int j, int k)
				{
					double res = (pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) / 4. - pre_pots(i, j, k);
					return pre_pots(i, j, k) + res;
					//return (_dxyz2*(pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) - _charge_map(i, j, k) * _dxyz2 * _dxyz2) / (4 * (_dxyz2));
					//return pre_pots(i, j, k) + ((2 * (pre_pots(i - 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i + 1, j, k) + pre_pots(i, j + 1, k))
					//	+ (pre_pots(i - 1, j - 1, k) + pre_pots(i + 1, j - 1, k) + pre_pots(i + 1, j + 1, k) + pre_pots(i - 1, j + 1, k))) / 12. - pre_pots(i, j, k));
				});
			}
			else
			{
				show(_potentials, "pre_pots");
				_potentials.maskBlock(_Mask, GAS, [pre_pots, this](float p, int i, int j, int k)
				{
					double res = (pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k) + pre_pots(i, j, k - 1) + pre_pots(i, j, k + 1)) / 6. - _charge_map(i, j, k) - pre_pots(i, j, k);
					return pre_pots(i, j, k) + res;
					//return (_dxyz2*(pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) - _charge_map(i, j, k) * _dxyz2 * _dxyz2) / (4 * (_dxyz2));
					//return pre_pots(i, j, k) + ((2 * (pre_pots(i - 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i + 1, j, k) + pre_pots(i, j + 1, k))
					//	+ (pre_pots(i - 1, j - 1, k) + pre_pots(i + 1, j - 1, k) + pre_pots(i + 1, j + 1, k) + pre_pots(i - 1, j + 1, k))) / 12. - pre_pots(i, j, k));
				});
				show(_potentials, "pre_pots");
				/*
				_potentials.maskBlock(_Mask, AIR, [pre_pots, this](float p, int i, int j, int k)
				{
					double res = (pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k) + pre_pots(i, j, k - 1) + pre_pots(i, j, k + 1)) / 6. - pre_pots(i, j, k);
					return pre_pots(i, j, k) + res;
					//return (_dxyz2*(pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) - _charge_map(i, j, k) * _dxyz2 * _dxyz2) / (4 * (_dxyz2));
					//return pre_pots(i, j, k) + ((2 * (pre_pots(i - 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i + 1, j, k) + pre_pots(i, j + 1, k))
					//	+ (pre_pots(i - 1, j - 1, k) + pre_pots(i + 1, j - 1, k) + pre_pots(i + 1, j + 1, k) + pre_pots(i - 1, j + 1, k))) / 12. - pre_pots(i, j, k));
				});
				*/
			}
			GRID temp = (_potentials.power(2) - pre_pots.power(2));
			temp.abs();
			temp.sqrt();
			double conv = temp.sum();
			std::cout << i << ": " << conv << std::endl;
			if (conv < _convergence || i == 9999)
			{
				std::cout << i << ": " << conv << ", " << convold << std::endl;
				temp = temp;
				break;
			}
			convold = conv;


		if (_grad_descent)
			get_gradient();
		}
	}
}

// FD neumann test
void Plasma_Ball::Poisson_Neuman(int its)
{
	_convergence = 1;
	double convold = 0;
	for (int i = 0; i < its; ++i)
	{
		GRID pre_pots(_potentials);
		_potentials.maskBlock(_Mask, GAS, [pre_pots, this](float p, int i, int j, int k)
		{
			double nw = pre_pots(i - 1, j - 1, k);
			double n = pre_pots(i - 1, j, k);
			double ne = pre_pots(i - 1, j + 1, k);
			double w = pre_pots(i, j - 1, k);
			double c = pre_pots(i, j, k);
			double e = pre_pots(i, j + 1, k);
			double sw = pre_pots(i + 1, j - 1, k);
			double s = pre_pots(i + 1, j, k);
			double se = pre_pots(i + 1, j + 1, k);

			// double res = 2 * (n + w + e + s) + (nw + ne + sw + se) / 12. - c;
			double res = (n + w + e + s) / 4. - c;

			//double res = (pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k)  + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k))/4. - pre_pots(i, j, k);
			//return pre_pots(i, j, k) + res;
			//return (_dxyz2*(pre_pots(i - 1, j, k) + pre_pots(i + 1, j, k) + pre_pots(i, j - 1, k) + pre_pots(i, j + 1, k)) - _charge_map(i, j, k) * _dxyz2 * _dxyz2) / (4 * (_dxyz2));
			return c + res;
		});
		GRID temp = (_potentials.power(2) - pre_pots.power(2));
		temp.abs();
		temp.sqrt();
		double conv = temp.sum(); 
		std::cout << i << ": " << conv << std::endl;
		if (conv < _convergence || i == 9999)
		{
			std::cout << i << ": " << conv << ", " << convold << std::endl;
			temp = temp;
			break;
		}
		convold = conv;

	}

}

void Plasma_Ball::create_dir()
{
	
	std::string dir = "\\_r-" + std::to_string((int)_radius)
		+ "_fils-" + std::to_string(_number_of_roots)
		+ "_fps-" + std::to_string(1./_dt)
		+ "_stepsize-" + std::to_string(_stepsize)
		+ "_residual-" + std::to_string(_residual)
		+ "_branching-" + std::to_string(_branching);


	_dir = _rootdir + dir + "_" + _name;
	CreateDirectory(const_cast<char *>(_dir.c_str()), NULL);
	_logname = _dir + "\\" + _logname;

	_logfile.open(_logname);
	_logfile << "Dimensions: " << std::to_string(_grid_size) << " x " << std::to_string(_grid_size) << std::endl;
	_logfile << "Radius in pixels: " << std::to_string((int)_radius) << std::endl;
	_logfile << "Number of seeded filaments: " << std::to_string(_number_of_roots) << std::endl;
	_logfile << "Stepsize: " << std::to_string(_stepsize) << std::endl;
	_logfile << "Residual multiplier: " << std::to_string(_residual) << std::endl;
	_logfile << "Branching multiplier: " << std::to_string(_branching) << std::endl;
}

void Plasma_Ball::show(GRID & Map, std::string name)
{
	if (!_debug) // only save vtk in debug
		if (name != "Leaders_")
		{	
			if (name != "_smooth_charge_map")
				if (name != "Leaders_vol_")
					return;

		}
	show_counters_[name]++;
	std::string it;

	std::string filename = name + "_r-" + std::to_string(_radius)
		+ "_eta-" + std::to_string(_eta)
		+ "_convergence-" + std::to_string(_convergence)
		+ "_fils-" + std::to_string(_number_of_roots)
		+ "_updating-" + std::to_string(_updating)
		+ "_serial-" + std::to_string(_serial)
		+ "step_size-" + std::to_string(_stepsize);
	std::vector<char> writearray;

	it = std::to_string(show_counters_[name]);
	it.insert(it.begin(), 4 - it.size(), '0');
	//std::cout << "Writing " + name + " it " + it + " as vtk" << std::endl;
	for (int n = 0; n < 1; n++)
	{
		//std::string dirname = (name == "Leaders_") ? "vtk_images/DBM/" : "vtk_images/DBM/debug/";
		std::string fname = _dir + "\\" + name + "_it-" + it + ".vtk";

		int POINTS = Map.cols()  * Map.rows() * Map.deps();
		int CELLS = Map.cols() - 1 * Map.rows() - 1 * (1u > Map.deps() - 1) ? 1u : Map.deps() - 1;
		std::ofstream file(fname, std::ios::out | std::ios::trunc);
		if (file)
		{
			// /* Structured grid
			file << "# vtk DataFile Version 2.0\n";
			file << "Plasma Ball, radius = " << _radius << ", eta = " << _eta << "\n";
			file << "ASCII" << "\n";
			file << "DATASET STRUCTURED_POINTS\n";
			file << "DIMENSIONS " << Map.cols() << " " << Map.rows() << " " << Map.deps() << "\n";
			//file << "POINTS " << Map.cols() * Map.rows() << "\n";
			file << "ASPECT_RATIO 1 1 1" << "\n";
			file << "ORIGIN 0 0 0" << "\n";
			file << "POINT_DATA " << POINTS << "\n";
			file << "SCALARS volume_scalars float 1\n";
			file << "LOOKUP_TABLE default";
			// */

			int rowindent = 0;
			for (unsigned k = 0; k < Map.deps(); ++k)
			{
				for (unsigned i = 0; i < Map.rows(); ++i)
				{
					for (unsigned j = 0; j < Map.cols(); ++j)
					{
						if (rowindent % 8 == 0) // 1 points per row (1 cell)
							file << "\n";
						file << Map(j, i, k) << " ";
						rowindent++;
					}
				}
			}
			file.close();
		}
		else {
			file.clear();
			file.close();
			std::cout << "Error writing to vtk" << std::endl;
		}
	}
}

void Plasma_Ball::name_dir(std::string name)
{
	_name = name;
}

/*
switch (_dim)
{
case CIRCULAR:

case SPERICAL:

default:
	break;
}
*/