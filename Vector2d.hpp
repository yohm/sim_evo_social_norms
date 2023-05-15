#include <vector>


template <typename T>
class Vector2d {
public:
  Vector2d() : _data(), n_rows(0), n_cols(0) {};
  Vector2d(size_t n_rows, size_t n_cols, T init) : _data(n_rows*n_cols, init), n_rows(n_rows), n_cols(n_cols) {};
  T& operator()(size_t i, size_t j) { return _data[i*n_cols+j]; }
  const T& operator()(size_t i, size_t j) const { return _data[i*n_cols+j]; }
  size_t Rows() const { return n_rows; }
  size_t Cols() const { return n_cols; }
  size_t size() const { return _data.size(); }
  T* data() { return _data.data(); }
  std::vector<T> _data;
  size_t n_rows;
  size_t n_cols;
};

