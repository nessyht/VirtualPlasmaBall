#include "Point3D.h"




int Point3D::x() const
{
	return _x;
}

int Point3D::y() const
{
	return _y;
}

int Point3D::z() const
{
	return _z;
}

Point3D Point3D::nw()
{
	return Point3D(_x - 1, _y + 1, _z);
}
Point3D Point3D::n()
{
	return Point3D(_x, _y + 1, _z);
}
Point3D Point3D::ne()
{
	return Point3D(_x + 1, _y + 1, _z);
}
Point3D Point3D::w()
{
	return Point3D(_x - 1, _y, _z);
}
Point3D Point3D::e()
{
	return Point3D(_x + 1, _y, _z);
}
Point3D Point3D::sw()
{
	return Point3D(_x - 1, _y - 1, _z);
}
Point3D Point3D::s()
{
	return Point3D(_x, _y - 1, _z);
}
Point3D Point3D::se()
{
	return Point3D(_x + 1, _y - 1, _z);
}

const bool Point3D::operator==(const Point3D & other)
{
	return this->_x == other.x() && this->_y == other.y() && this->_z == other.z();
}

const bool Point3D::operator!=(const Point3D & other)
{
	return !operator==(other);
}

const bool Point3D::operator>=(const Point3D & other)
{
	return this->_x >= other.x() && this->_y >= other.y() && this->_z >= other.z();
}

double Point3D::l2() const 
{
	return sqrt(_x*_x + _y*_y + _z*_z);
}

double Point3D::l2(const Point3D p) const 
{
	return sqrt((_x-p.x())*(_x - p.x()) + (_y - p.y())* (_y - p.y()) + (_z - p.z()) * (_z - p.z()));
}

std::ostream & operator<<(std::ostream & os, const Point3D & p)
{
	os << "(" << p.x() << ", " << p.y() << ", " << p.z() << ")";
	return os;
}
