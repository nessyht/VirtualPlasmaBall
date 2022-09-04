#pragma once
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cmath>
#include <cassert>
#include "Point3D.h"

// Use (void) to silent unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

template <typename T>
class GridXd
{
private:
	int _x;
	int _y;
	int _z;
	T _init;
	std::vector<T> _GRID;

public:

	GridXd(int y, int x, int z, T fill);
	GridXd(const GridXd<T> &other);
	GridXd(const GridXd<T> &other, bool reset);
	GridXd();

	T& operator()(const unsigned &, const unsigned &, const unsigned &);
	T& operator()(Point point);
	T& operator()(int index); // returns value at this index in GRID warning

	const T& operator()(const unsigned &, const unsigned &, const unsigned &) const;
	const T& operator()(const Point &point) const;
	GridXd operator+(const GridXd<T> &other);
	GridXd operator-(const GridXd<T> &other);
	GridXd operator*(const GridXd<T> &other);
	GridXd operator/(const GridXd<T> &other);

	const bool operator==(const GridXd<T> &other);
	const void operator=(const GridXd<T> &other);
	void operator+=(const GridXd<T> &other);
	void operator-=(const GridXd<T> &other);

	unsigned rows() const;
	unsigned cols() const;
	unsigned deps() const;
	T init() const;
	const std::vector<T> mat() const;

	int pos(int col, int row, int dep);
	int pos(Point point);
	T sum() const;
	void pow(const T &p);
	GridXd power(const T &p );
	T maxCoeff(); 
	T minCoeff(); 
	void abs();
	GridXd sqrt();
	T* data() noexcept;
	unsigned size() const;
	unsigned length() const; // size() for wierd edgcases 
	void normalise(T val);
	void mask_clamp(const GridXd<T>& mask, const T value, T lower, T upper);
	int count_instance(T fill) const;
	void fill2Dto3D(const GridXd<T> &other);

	void GRID(T *other)
	{
		std::vector<T> temp;
		for (int i = 0; i < _GRID.size(); ++i)
		{
			temp.push_back(*other);
			//if (*other != 0)
			//	std::cout << *other << std::endl;
			++other;
			
		}
		_GRID = temp;

	}

	GridXd operator+(const T input);
	GridXd operator-(const T input);
	GridXd operator/(const T input);
	GridXd operator*(const T input);

	T interpolate(double x, double y, double z) const; 
	T interpolate(Point3D<double> p) const;
	T interpolate(Point3D<float> p) const; 
	T interpolate1D(double a, double b, double x) const; 
	T interpolate2D(double a, double b, double c, double d, double x, double y) const; 
	T interpolate3D(double a, double b, double c, double d,
				double e, double f, double g, double h, 
				double x, double y, double z) const; 

	GridXd multiply(const T input) const;


	GridXd operator/=(const T input);
	
	friend std::ostream& operator<<(std::ostream& os, const GridXd<T>& lhs);


	void block(int x, int y, int z, int cols, int rows, int dep, T fill);
	void block(int x, int y, int z, int cols, int rows, int dep, const GridXd<T>& fill);
	void block(int x, int y, int z, int cols, int rows, int dep, std::function<T(T)> const& lambda);
	void block(int x, int y, int z, int cols, int rows, int dep, std::function<T(T, int, int, int)> const& lambda);
	void block(std::vector<Point> points, std::function<T(T, Point)> const& lambda);
	void block(std::vector<Point> points, T fill);
	GridXd block(int x, int y, int z, int cols, int rows, int dep);
	void blockweightedneighbours(std::vector<VEC> points); //TODO: FILL IN REVERSE INTERPOLATION STYLE

	/*
	* 3D version as no z value is given
	*/
	void sphereBlock(const T fill, const float radius_inner, const float radius_outer = -1);
	/*
	* 2D version circle at z = depth;
	*/
	void circleBlock(const T fill, const float radius_inner,  const float radius_outer = -1, const unsigned depth = 0);
	/*
	* 2D circle Block with accurate distribution inside cell
	*/
	void circleBlockprecise(const T fill, const float radius_inner, const float radius_outer = -1, const unsigned depth = 0);
	/*
	* Where mask = value, let this = fill
	*/
	void maskBlock(const GridXd<T>& mask, const T value, const float fill);
	/*
	* Where mask = value, let this = fill by lambda
	*/
	void maskBlock(const GridXd<T>& mask, const T value, std::function<T(T, int, int, int)> const& lambda);
	   
	/*
	* Where mask = value from list, let this = fill by lambda
	*/
	//void maskBlock(const GridXd<T>& mask, const T * values, std::function<T(T, int, int, int)> const & lambda);

#ifdef cl
/*/**
	 * @brief Convert volume data to UCHAR format and generate OpenCL image textue memory object.
	 * @param volumeData The raw volume data.
	 */
	template<class T>
	cl::Image3D volDataToCLmem()
	{
		// reinterpret raw data (char) to input format
		auto s = reinterpret_cast<const T *>(_GRID.data());
		auto e = reinterpret_cast<const T *>(_GRID.data() + _GRID.size());
		// convert input vector to the desired output precision
		std::vector<unsigned char> convertedData(s, e);
		try
		{
			cl::ImageFormat format;
			format.image_channel_order = CL_R;
			format.image_channel_data_type = CL_UNORM_INT8;

			cl::Image3D image(cl::Image3D(_contextCL,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				format,
				_x,
				_y,
				_z,
				0, 0,
				const_cast<unsigned char*>(convertedData.data())));
		}
		catch (cl::Error err)
		{
			throw std::runtime_error("ERROR: " + std::string(err.what()) + "("
				+ getCLErrorString(err.err()) + ")");
		}
	}*/
#endif // cl


};

template<typename T>
 GridXd<T>::GridXd(int y, int x, int z, T fill) : _y(y), _x(x), _z(z), _init(fill)
{
	_GRID = std::vector<T>(y*x*z, fill);
}

 template<typename T>
 GridXd<T>::GridXd(const GridXd<T> &other) : _x(other.cols()), _y(other.rows()), _z(other.deps()), _GRID(other.mat()), _init(other.init())
 {
 }

template<typename T>
GridXd<T>::GridXd(const GridXd<T> &other, bool reset) : _x(other.cols()), _y(other.rows()), _z(other.deps()), _GRID(other.mat()), _init(other.init())
{
	_GRID = (reset) ? std::vector<T>(_y*_x*_z, 0.0) : other.mat();
}

template<typename T>
GridXd<T>::GridXd()
{
}

template<typename T>
T& GridXd<T>::operator()(const unsigned &col, const unsigned &row, const unsigned &dep) {
	return _GRID[_x*_y*dep + _x * row + col];
}

template<typename T>
T& GridXd<T>::operator()(Point point)
{
	assert(point.x() >= 0);
	assert(point.y() >= 0);
	assert(point.z() >= 0);
	return _GRID[_x*_y*point.z() + _x * point.y() + point.x()];
}

template<typename T>
inline T & GridXd<T>::operator()(int index)
{
	return _GRID[index];
}

template<typename T>
const T& GridXd<T>::operator()(const unsigned &col, const unsigned &row, const unsigned &dep) const {
	return _GRID[_x*_y*dep + _x * row + col];
}

template<typename T>
const T & GridXd<T>::operator()(const Point &point) const
{
	return _GRID[_x*_y*point.z() + _x * point.y() + point.x()];
}

template<typename T>
 GridXd<T> GridXd<T>::operator+(const GridXd<T> &other)
{
	assert(this->rows() == other.rows());
	assert(this->cols() == other.cols());
	assert(this->deps() == other.deps());
	GridXd<T> res(other);
	for (int k = 0; k < _z; ++k)
		for (int j = 0; j < _y; ++j)
			for (int i = 0; i < _x; ++i)
				res(i, j, k) += this->operator()(i, j, k);
	return res;
}

template<typename T>
GridXd<T> GridXd<T>::operator-(const GridXd<T> &other)
{
	assert(this->rows() == other.rows());
	assert(this->cols() == other.cols());
	assert(this->deps() == other.deps());
	GridXd<T> res(other);
	for (int k = 0; k < _z; ++k)
		for (int j = 0; j < _y; ++j)
			for (int i = 0; i < _x; ++i)
				res(i, j, k) -= this->operator()(i, j, k);
	return res;
}

template<typename T>
GridXd<T> GridXd<T>::operator*(const GridXd<T>& other)
{
	assert(this->rows() == other.rows());
	assert(this->cols() == other.cols());
	assert(this->deps() == other.deps());
	GridXd<T> res(other);
	for (int k = 0; k < _z; ++k)
		for (int j = 0; j < _y; ++j)
			for (int i = 0; i < _x; ++i)
				res(i, j, k) *= this->operator()(i, j, k);
	return res;
}

template<typename T>
GridXd<T> GridXd<T>::operator/(const GridXd<T>& other)
{
	assert(this->rows() == other.rows());
	assert(this->cols() == other.cols());
	assert(this->deps() == other.deps());
	GridXd<T> res(other);
	for (int k = 0; k < _z; ++k)
		for (int j = 0; j < _y; ++j)
			for (int i = 0; i < _x; ++i)
				res(i, j, k) /= this->operator()(i, j, k);
	return res;
}


template<typename T>
const bool GridXd<T>::operator==(const GridXd<T> & other)
{
	assert(this->rows() == other.rows());
	assert(this->cols() == other.cols());
	assert(this->deps() == other.deps());
	return _GRID == other.mat();
}

template<typename T>
const void GridXd<T>::operator=(const GridXd<T> & other)
{
	_x = other.cols();
	_y = other.rows();
	_z = other.deps();
	_init = other.init();
	_GRID = other.mat();
}

template<typename T>
void GridXd<T>::operator+=(const GridXd<T>& other)
{
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				this->operator()(i, j, k) = this->operator()(i, j, k) + other(i, j, k);
}

template<typename T>
void GridXd<T>::operator-=(const GridXd<T>& other)
{
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				this->operator()(i, j, k) = this->operator()(i, j, k) - other(i, j, k);
}

template<typename T>
unsigned GridXd<T>::rows() const
{
	return _y;
}

template<typename T>
unsigned GridXd<T>::cols() const
{
	return this->_x;
}

template<typename T>
unsigned GridXd<T>::deps() const
{
	return _z;
}

template<typename T>
T GridXd<T>::init() const
{
	return _init;
}

template<typename T>
const std::vector<T> GridXd<T>::mat() const
{
	return _GRID;
}

template<typename T>
inline int GridXd<T>::pos(int col, int row, int dep)
{
	return _x * _y*dep + _x * row + col;
}

template<typename T>
inline int GridXd<T>::pos(Point point)
{
	return _x * _y*point.z() + _x * point.y() + point.x();
}

template<typename T>
T GridXd<T>::sum() const
{
	T sum = 0.0;
	for (T el : _GRID)
		sum += el;
	return sum;
}

template<typename T>
void GridXd<T>::pow(const T & p)
{
	std::for_each(_GRID.begin(), _GRID.end(), [p](T &el) {std::pow(el, p); });
}

template<typename T>
GridXd<T> GridXd<T>::power(const T &p)
{
	GridXd<T> res(*this);
	res.pow(p);
	return res;
}

template<typename T>
 T GridXd<T>::maxCoeff()
{
	return *std::max_element(_GRID.begin(), _GRID.end());
}

template<typename T>
 T GridXd<T>::minCoeff()
{
	return *std::min_element(_GRID.begin(), _GRID.end());
}

 template<typename T>
 inline void GridXd<T>::abs()
 {
	 for (auto &el : _GRID)
		 el = std::abs(el);
	 //std::for_each(_GRID.begin(), _GRID.end(), [](T &el) {std::abs(el); });
 }

template<typename T>
GridXd<T> GridXd<T>::sqrt()
{
	GridXd<T> res(*this);
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				res(i, j, k) = std::sqrt(this->operator()(i, j, k));
	return res;
	//std::for_each(_GRID.begin(), _GRID.end(), [](T &el) {std::sqrtf(el); });
}

template<typename T>
T * GridXd<T>::data() noexcept
{
	return _GRID.data();
}

template<typename T>
unsigned GridXd<T>::size() const
{
	return _GRID.size();
}

template<typename T>
unsigned GridXd<T>::length() const
{
	return _GRID.size();
}

template<typename T>
inline void GridXd<T>::normalise(T val)
{
	val = std::abs(val);
	T max = this->maxCoeff();
	T min = this->minCoeff();
	max = (std::abs(max) > std::abs(min)) ? std::abs(max) : std::abs(min);
	std::cout << max << std::endl; //debug
	T fact = val / max;
	for (T &i : _GRID)
		i *= fact;
}

template<typename T>
inline void GridXd<T>::mask_clamp(const GridXd<T> &mask, const T value, T lower, T upper)
{
	lower = (lower < upper) ? lower : upper; // idiot proof
	upper = (lower > upper) ? lower : upper;

	T max = this->maxCoeff();
	T min = this->minCoeff();
	if (min >= lower && max <= upper) // nothing to be done
		return;

	for (int i = 0; i < this->size(); ++i)
		if (mask(i) == value)
		{
			if (_GRID[i] < lower)
				_GRID[i] = lower;
			else if (_GRID[i] > upper)
				_GRID[i] = upper;
		}
}

template<typename T>
inline int GridXd<T>::count_instance(T fill) const
{
	int count = 0;
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				if (this->operator()(i, j, k) == fill)
					++count;
	return count;
}

template<typename T>
inline void GridXd<T>::fill2Dto3D(const GridXd<T>& other)
{
	int cols = this->cols();
	int rows = this->rows();
	int deps = this->deps();
	int plane_size = other.length();

	for (int k = 0; k < deps; ++k)
		for (int i = 0; i < plane_size; ++i)
			_GRID[k *plane_size + i] = other(i);
}

template<typename T>
GridXd<T> GridXd<T>::operator+(const T input)
{
	GridXd<T> res(*this);
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				res(i, j, k) = this->operator()(i, j, k) + input;
	return res;
}

template<typename T>
 GridXd<T> GridXd<T>::operator-(const T input)
{
	GridXd<T> res(*this);
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				res(i, j, k) = this->operator()(i, j, k) - input;
	return res;
}

template<typename T>
 GridXd<T> GridXd<T>::operator/(const T input)
{
	GridXd<T> res(*this);
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				res(i, j, k) = this->operator()(i, j, k) / input;
	return res;
}

 template<typename T>
 GridXd<T> GridXd<T>::operator*(const T input)
 {
	 GridXd<T> res(this->cols(), this->rows(), this->deps(), 0.0);
	 auto first = res(0, 0, 0);
	 for (int k = 0; k < res.deps(); ++k)
		 for (int j = 0; j < res.rows(); ++j)
			 for (int i = 0; i < res.cols(); ++i)
				 res(i, j, k) = this->operator()(i, j, k) * input;
	 return res;
}

 template<typename T>
 inline T GridXd<T>::interpolate(double x, double y, double z) const
 {
	 int row0 = (int)x;
	 int col0 = (int)y;
	 int dep0 = (int)z;

	 int dim = 0;
	 double inx = std::abs(row0 - x);
	 double iny = std::abs(col0 - y);
	 double inz = std::abs(dep0 - z);

	 std::vector<Point> indexes;
	 int row1 = row0 + 1;
	 int col1 = col0 + 1;
	 int dep1 = dep0 + 1;
	 
	 if (inx > 0)
	 {
		 ++dim;
		 indexes.push_back(Point(row0, col0, dep0));
		 indexes.push_back(Point(row0 + 1, col0, dep0));
	 }
	 if (iny > 0)
	 {
		 ++dim;

		 indexes.push_back(Point(row0, col0 + 1, dep0));
		 indexes.push_back(Point(row0 + 1, col0 + 1, dep0));

	 }
	 if (inz > 0)
	 {
		 ++dim;
		 indexes.push_back(Point(row0, col0, dep0 + 1));
		 indexes.push_back(Point(row0 + 1, col0, dep0 + 1));
		 indexes.push_back(Point(row0, col0 + 1, dep0 + 1));
		 indexes.push_back(Point(row0 + 1, col0 + 1, dep0 + 1));
	 }
	 std::vector<T> values;
	 for (auto p : indexes)
		 values.push_back(this->operator()(p));
	 
	 if (dim == 0)
		 return this->operator()(x, y, z);
	 else if (dim == 1)
		 return interpolate1D(values[0], values[1], x);
	 else if (dim == 2)
		 return interpolate2D(values[0], values[1], values[2], values[3], x, y);
	 else if (dim == 3)
		 return interpolate3D(	values[0], values[1], values[2], values[3],
								values[4], values[5], values[6], values[7],
								x, y, z
							 );
 }

 template<typename T>
 inline T GridXd<T>::interpolate(Point3D<double> p) const
 {
	 return interpolate(p.x(), p.y(), p.z());
 }

 template<typename T>
 inline T GridXd<T>::interpolate1D(double a, double b, double x) const
 {
	 return a * (1 - x) + b * (x);
 }

 template<typename T>
 inline T GridXd<T>::interpolate2D(double a, double b, double c, double d, double x, double y) const
 {
	 return interpolate1D(interpolate1D(a, b, x), interpolate1D(c, d, x), y);
 }

 template<typename T>
 inline T GridXd<T>::interpolate3D(double a, double b, double c, double d,
									 double e, double f, double g, double h,
									 double x, double y, double z) const
 {
	 T m = interpolate2D(a, c, e, g, x, y);
	 T n = interpolate2D(b, d, f, h, x, y);
	 return interpolate1D(m, n, z);
 }



 template<typename T>
 GridXd<T> GridXd<T>::multiply(const T input) const
 {
	 GridXd<T> res(*this);
	 for (int k = 0; k < res.deps(); ++k)
		 for (int j = 0; j < res.rows(); ++j)
			 for (int i = 0; i < res.cols(); ++i)
				 res(i, j, k) = res(i, j, k) * input;
	 return res;
 }

template<typename T>
GridXd<T> GridXd<T>::operator/=(const T input)
{
	for (int k = 0; k < deps(); ++k)
		for (int j = 0; j < rows(); ++j)
			for (int i = 0; i < cols(); ++i)
				this->operator()(i, j, k) = this->operator()(i, j, k) / input;
	return *this;
}

template<typename T>
std::ostream & operator<<(std::ostream & os, const GridXd<T> & mat)
{
	for (int k = 0; k < mat.deps(); ++k)
	{
		for (int j = 0; j < mat.rows(); ++j)
		{
			for (int i = 0; i < mat.cols(); ++i)
				os << mat(i, j, k) << " ";
			os << std::endl;
		}
		os << std::endl;
	}
	return os;
}

template<typename T>
void GridXd<T>::block(int x, int y, int z, int cols, int rows, int deps, T fill)
{
	for (int k = z; k < deps; ++k)
		for (int j = y; j < rows; ++j)
			for (int i = x; i < cols; ++i)
				this->operator()(i, j, k) = fill;
}

template<typename T>
void GridXd<T>::block(int x, int y, int z, int cols, int rows, int deps, const GridXd<T>& fill)
{
	assert(cols - x == fill.cols());
	assert(rows - y == fill.rows());
	assert(deps - z == fill.deps());
	for (int k = z; k < z + deps; ++k)
		for (int j = y; j < y + rows; ++j)
			for (int i = x; i < x + cols; ++i)
				this->operator()(i, j, k) = fill(i - x, j - y, k - z);
}

template<typename T>
void GridXd<T>::block(int x, int y, int z, int cols, int rows, int deps, std::function<T(T)> const& lambda)
{
	for (int k = z; k < z + deps; ++k)
		for (int j = y; j < y + rows; ++j)
			for (int i = x; i < x + cols; ++i)
				this->operator()(i, j, k) = lambda(this->operator()(i, j, k));
}

template<typename T>
void GridXd<T>::block(int x, int y, int z, int cols, int rows, int deps, std::function<T(T, int, int, int)> const& lambda)
{
	for (int k = z; k < z + deps; ++k)
		for (int j = y; j < y + rows; ++j)
			for (int i = x; i < x + cols; ++i)
				this->operator()(i, j, k) = lambda(this->operator()(i, j, k), i, j, k);
}

template<typename T>
void GridXd<T>::block(std::vector<Point> points, std::function<T(T, Point)> const& lambda)
{
	for (const auto &point : points)
		this->operator()(point) = lambda(this->operator()(point), point);
}

template<typename T>
void GridXd<T>::block(std::vector<Point> points, T fill)
{
	for (auto point : points)
		this->operator()(point) = fill;
}

template<typename T>
 GridXd<T> GridXd<T>::block(int x, int y, int z, int cols, int rows, int deps)
{
	GridXd<T> res(cols-x, rows-y, deps-z, 0);
	for (int k = z; k < deps; ++k)
		for (int j = y; j < rows; ++j)
			for (int i = x; i < cols; ++i)
				res(i, j, k) = this->operator()(i, j, k);
	return res;
}

 template<typename T>
 inline void GridXd<T>::blockweightedneighbours(std::vector<VEC> points)
 {
	 for (auto p : points)
	 {
		 return;
		 //TODO: FILL IN REVERSE INTERPOLATION STYLE
	 }

 }

template<typename T>
 void GridXd<T>::sphereBlock(const T fill, const float radius_inner, const float radius_outer)
{
	float real_radius_outer = (radius_outer == -1) ? radius_inner + 1 : radius_outer;

	if (_z > 1)
	{
		const float cx = _x / 2.;
		const float cy = _y / 2.;
		const float cz = _z / 2.;
		this->block(cx - real_radius_outer - 1, cy - real_radius_outer - 1, cz - real_radius_outer - 1,
					cx + real_radius_outer + 1, cy + real_radius_outer + 1, cz - real_radius_outer + 1,
			[radius_inner, real_radius_outer, fill, cx, cy, cz](T current, int i, int j, int k) {
			float R = (i - cx)*(i - cx) + (j - cy)*(j - cy) + (k - cz)*(k - cz);
			if ( R >= radius_inner*radius_inner && R < (real_radius_outer) * (real_radius_outer))
				return fill;
			else
				return current;
		});
	}
}

template<typename T>
 void GridXd<T>::circleBlock(const T fill, const float radius_inner, const float radius_outer, const unsigned depth)
{
	float real_radius_outer = (radius_outer == -1) ? radius_inner + 1 : radius_outer;

	const float cx = _x / 2.;
	const float cy = _y / 2.;
	Point center(cx, cy, 0);
	this->block(cx - real_radius_outer - 1, cy - real_radius_outer - 1, depth, 
				cx + real_radius_outer + 1, cy + real_radius_outer + 1, depth + 1,
		[radius_inner, real_radius_outer, fill, center](T current, int i, int j, int k) {
		double R = Point(i, j, k).l2(center);
		//float R = (i - cx)*(i - cx) + (j - cy)*(j - cy);
		if (R >= radius_inner && R < real_radius_outer)
			return fill;
		else
			return current;
	});
}

 template<typename T>
 inline void GridXd<T>::circleBlockprecise(const T fill, const float radius_inner, const float radius_outer, const unsigned depth)
 {
	 float real_radius_outer = (radius_outer == -1) ? radius_inner + 1 : radius_outer;

	 const float cx = _x / 2.;
	 const float cy = _y / 2.;
	 Point center(cx, cy, 0);
	 this->block(0, 0, depth, _x, _y, depth + 1,
		 [cx, cy, radius_inner, real_radius_outer, fill, center](T current, int i, int j, int k) {
		 double R = Point(i, j, k).l2(center);

		 if (R >= radius_inner && R < real_radius_outer)
		 {
			 double accumulator = 0;
			 double rpw = std::sqrt((i - 0.5 - cx)*(i - 0.5 - cx) + (j - cy)*(j - cy));
			 double rpn = std::sqrt((i + 0.5 - cx)*(i + 0.5 - cx) + (j - cy)*(j - cy));
			 double rpe = std::sqrt((i - cx)*(i - cx) + (j - 0.5 - cy)*(j - 0.5 - cy));
			 double rps = std::sqrt((i - cx)*(i - cx) + (j + 0.5 - cy)*(j + 0.5 - cy));
			 if (rpw >= radius_inner && rpw < real_radius_outer)
				 accumulator += 0.25;
			 if (rpn >= radius_inner && rpn < real_radius_outer)
				 accumulator += 0.25;
			 if (rpe >= radius_inner && rpe < real_radius_outer)
				 accumulator += 0.25;
			 if (rps >= radius_inner && rps < real_radius_outer)
				 accumulator += 0.25;
			return (float)(accumulator * fill);
		 }
		 else
			 return current;
	 });
 }

template<typename T>
 void GridXd<T>::maskBlock(const GridXd<T>& mask, const T value, const float fill)
{
	assert(mask.cols() == this->cols());
	assert(mask.rows() == this->rows());
	assert(mask.deps() == this->deps());

	for (int k = 0; k < mask.deps(); ++k)
		for (int j = 0; j < mask.rows(); ++j)
			for (int i = 0; i < mask.cols(); ++i)
				if (mask(i, j, k) == value)
					this->operator()(i, j, k) = fill;

	/*const auto& maskmat = mask.mat();
	for (int i = 0; i < mask.mat().size(); ++i)
		if (maskmat[i] == value)
			_GRID[i] = fill;
*/
 }

 template<typename T>
 void GridXd<T>::maskBlock(const GridXd<T>& mask, const T value, std::function<T(T, int, int, int)> const& lambda)
 {
	 // std::cout << mask.cols() <<  "  " << this->cols();
	 assert(mask.cols() == this->cols());
	 assert(mask.rows() == this->rows());
	 assert(mask.deps() == this->deps());

	 for (int k = 0; k < mask.deps(); ++k)
		 for (int j = 0; j < mask.rows(); ++j)
			 for (int i = 0; i < mask.cols(); ++i)
				if(mask(i, j, k) == value)
					this->operator()(i, j, k) = lambda(this->operator()(i, j, k), i, j, k);	
 }

 typedef GridXd<float> GRID;