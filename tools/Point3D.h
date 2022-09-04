#pragma once
#include <cmath>
#include <ostream>



template <typename T>
class Point3D
{
public:
	T _x;
	T _y;
	T _z;

	Point3D(T X, T Y, T Z) : _x(X), _y(Y), _z(Z)
	{}

	Point3D(T X, T Y) : _x(X), _y(Y), _z(0)
	{}

	Point3D(T X) : _x(X), _y(0), _z(0)
	{}

	template <typename TT>
	Point3D(const Point3D<TT> &copy) : _x(copy.x()), _y(copy.y()), _z(copy.z())
	{}

	Point3D()
	{}

	T x() const
	{
		return _x;
	}

	T y() const
	{
		return _y;
	}

	T z() const
	{
		return _z;
	}

	Point3D nw()
	{
		return Point3D<int>(_x - 1, _y + 1, _z);
	}

	Point3D n()
	{
		return Point3D<int>(_x, _y + 1, _z);
	}

	Point3D ne()
	{
		return Point3D<int>(_x + 1, _y + 1, _z);
	}

	Point3D w()
	{
		return Point3D<int>(_x - 1, _y, _z);
	}

	Point3D e()
	{
		return Point3D<int>(_x + 1, _y, _z);
	}

	Point3D sw()
	{
		return Point3D<int>(_x - 1, _y - 1, _z);
	}

	Point3D s()
	{
		return Point3D<int>(_x, _y - 1, _z);
	}

	Point3D se()
	{
		return Point3D<int>(_x + 1, _y - 1, _z);
	}

	const bool operator==(const Point3D &other)
	{
		return this->_x == other.x() && this->_y == other.y() && this->_z == other.z();
	}

	const bool operator!=(const Point3D &other)
	{
		return !operator==(other);
	}

	const bool operator>=(const Point3D &other)
	{
		return this->_x >= other.x() && this->_y >= other.y() && this->_z >= other.z();
	}

	Point3D operator+(const Point3D &other)
	{
		return Point3D(_x + other.x(), _y + other.y(), _z + other.z());
	}

	Point3D operator-(const Point3D &other)
	{
		return Point3D(_x - other.x(), _y - other.y(), _z - other.z());
	}
	
	Point3D operator*(const double mul)
	{
		return Point3D(_x * mul, _y * mul, _z * mul);
	}

	friend std::ostream& operator<<(std::ostream& os, const Point3D& dt)
	{
		os << "(" << dt.x() << ", " << dt.y() << ", " << dt.z() << ")";
		return os;
	}

	double l2() const
	{
		return sqrt(_x*_x + _y * _y + _z * _z);
	}

	double l2(const Point3D p) const
	{
		return sqrt((_x - p.x())*(_x - p.x()) + (_y - p.y())* (_y - p.y()) + (_z - p.z()) * (_z - p.z()));
	}
	
	Point3D normalise(bool max_1 = false)
	{
		if (max_1)
		{
			// scale to 1 step in a direction
			T xt = std::abs(_x);
			T yt = std::abs(_y);
			T zt = std::abs(_z);
			T max = (xt > yt) ? xt : yt;
			max = (max > zt) ? max : zt;
			return Point3D(_x / max, _y / max, _z / max);
		}
		else
		{
			double length = this->l2();

			return Point3D(_x / length, _y / length, _z / length);
		}

	}


	bool is_neighbour(const Point3D &other);
};

typedef Point3D<double> VEC;

typedef Point3D<int> Point;