#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<cstdlib>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>
#include<cmath>
#include<queue>
#include<random>
#include<iomanip>
using namespace std;


class MyException :public exception {
private:
	string _info;
public:
	MyException(string _info_str) {
		_info = _info_str;
	}
	const char* what() const throw() {
		return _info.c_str();
	}
};

class Equation {
public:
	int arg_num;
	int op;
	vector<double> w;
	double b;

	Equation(string raw_equ_str) {
		vector<double> tmp = str_to_vector(raw_equ_str);
		arg_num = tmp.size() - 2;

		w = tmp;
		w.erase(w.end() - 1);
		w.erase(w.end() - 1);

		op = (int)*(tmp.end() - 1);
		b = *(tmp.end() - 2);
	}

	Equation(int _arg_num, int _op, vector<double> _w, double _b) {
		arg_num = _arg_num;
		op = _op;
		w = _w;
		b = _b;
	}

	Equation(int _arg_num, int _op, int idx, double _b) {
		arg_num = _arg_num;
		op = _op;
		b = _b;
		if (idx >= arg_num) {
			throw MyException("Equation::Equation(int, int, int, double) idx is out of bound arg_num");
		}
		//构造系数
		vector<double> tmp(arg_num);
		for (int i = 0; i < arg_num; ++i) {
			tmp[i] = 0;
		}
		tmp[idx] = 1;
		w = tmp;
	}

	vector<double> gen_form_equ() {
		vector<double> tmp = w;
		tmp.push_back(b);
		return tmp;
	}

	static vector<double> str_to_vector(string str) {
		vector<string> tmp = split(str, ' ');
		vector<double> res;
		for (vector<string>::iterator it = tmp.begin(); it != tmp.end(); ++it) {
			res.push_back(stod(*it));
		}
		return res;
	}

	static vector<string> split(string str, char symbol) {
		int m = str.size();
		vector<string> res;
		string buf;
		for (int i = 0; i < m; ++i) {
			if (str[i] != symbol) {
				buf += str[i];
			}
			else {
				res.push_back(buf);
				buf = "";
			}
		}
		if (buf.size()) {
			res.push_back(buf);
		}
		return res;
	}
};

class _base_matrix {
protected:
	vector<vector<double>> _matrix;
	static const int INIT_SIZE = 5;
	static const int OUT_WIDTH = 7;
public:
	_base_matrix() {
		_matrix.clear();
		for (int i = 0; i < INIT_SIZE; ++i) {
			_matrix.push_back(vector<double>(INIT_SIZE));
		}
	}

	_base_matrix(int _row, int _col) {
		if (_row < 0 ||_col < 0) {
			//throw
		}
		_matrix.clear();
		for (int i = 0; i < _row; ++i) {
			_matrix.push_back(vector<double>(_col));
		}
	}

	_base_matrix(const _base_matrix& __matrix) {
		_matrix = __matrix._matrix;
	}

	_base_matrix(vector<vector<double>> __matrix) {
		_matrix = __matrix;
	}

	//设置
	void set(int row, int col, double elem) {
		int row_bound = shape().first;
		int col_bound = shape().second;
		if (row < 0 || col < 0 || row >= row_bound || col >= col_bound) {
			stringstream ss;
			ss << "_base_matrix::set (" << row << "," << col << ") is out of bound (" << row_bound << "," << col_bound << ").";
			throw MyException(ss.str());
		}

		_matrix[row][col] = elem;
	}
	//用他自己设置
	void set(vector<vector<double>> __matrix) {
		_matrix = __matrix;
	}

	//得到矩阵的形状
	pair<int, int> shape() {
		int row = _matrix.size();
		if (_matrix.size() == 0) {
			//throw
		}
		int col = _matrix[0].size();
		pair<int, int> _shape(row, col);
		return _shape;
	}

	//删除矩阵的某一行
	static _base_matrix _del_single_row(_base_matrix __matrix,int row_id) {
		if (row_id < 0 || row_id >= __matrix.shape().first) {
			throw MyException("_base_matrix::_del_single_row row_id is out of bound.");
		}
		_base_matrix tmp = __matrix;
		tmp._matrix.erase(tmp._matrix.begin() + row_id);
		return tmp;
	}

	
	void _del_single_row(int row_id) {
		if (row_id < 0 || row_id >= shape().first) {
			throw MyException("_base_matrix::_del_single_row row_id is out of bound.");
		}
		_matrix.erase(_matrix.begin() + row_id);
	}

	//删除矩阵的某几行
	static _base_matrix _del_row(_base_matrix __matrix, vector<int> rows_id) {
		int m = rows_id.size();
		_base_matrix tmp = __matrix;
		for (int i = 0; i < m; ++i) {
			tmp._del_single_row(rows_id[i]);
		}
		return tmp;
	}

	
	void _del_row(vector<int> rows_id) {
		int m = rows_id.size();
		for (int i = 0; i < m; ++i) {
			_del_single_row(rows_id[i]);
		}
	}

	//删除矩阵的某一列
	static _base_matrix del_single_col(_base_matrix __matrix,int col_id) {
		if (col_id > __matrix.shape().second) {
			throw MyException("_base_matrix::del_single_col col_id is out of bound.");
		}
		_base_matrix tmp = __matrix;
		int row = tmp.shape().first;
		for (int i = 0; i < row; i++) {
			tmp._matrix[i].erase(tmp._matrix[i].begin() + col_id);
		}
		return tmp;
	}

	
	void _del_single_col(int col_id) {
		if (col_id > shape().second) {
			throw MyException("_base_matrix::del_single_col col_id is out of bound.");
		}

		int row = _matrix.size();
		for (int i = 0; i < row; i++) {
			_matrix[i].erase(_matrix[i].begin() + col_id);
		}
	}

	//删除矩阵的某几列
	static _base_matrix del_cols(_base_matrix __matrix,vector<int> cols_id) {
		int max_col = __matrix._matrix[0].size();
		vector<int>::iterator max_cols_iter = max_element(cols_id.begin(), cols_id.end());
		if (*max_cols_iter >= max_col) {
			throw MyException("_base_matrix<_Tp>::del_cols cols_id is out of bounds");
		}
		_base_matrix tmp = __matrix;
		sort(cols_id.begin(), cols_id.end());
		int m = cols_id.size();
		for (int i = 0; i < m; ++i) {
			tmp._del_single_col(cols_id[i]);
			tmp.vector_elem_add_num<int>(cols_id, -1);
		}
		return tmp;
	}

	void _del_cols(vector<int> cols_id) {
		int max_col = _matrix[0].size();
		vector<int>::iterator max_cols_iter = max_element(cols_id.begin(), cols_id.end());
		if (*max_cols_iter >= max_col) {
			throw MyException("_base_matrix<_Tp>::del_cols cols_id is out of bounds");
		}
		sort(cols_id.begin(), cols_id.end());
		int m = cols_id.size();
		for (int i = 0; i < m; ++i) {
			_del_single_col(cols_id[i]);
			vector_elem_add_num<int>(cols_id, -1);
		}
	}

	template<typename T>
	void vector_elem_add_num(vector<T> &v, T elem) {
		int m = v.size();
		for (int i = 0; i < m; ++i) {
			v[i] += elem;
		}
	}
	
	//在矩阵行最后添加一个矩阵
	static _base_matrix append_single_row(_base_matrix __matrix,_base_matrix _add_matrix) {
		if (__matrix.shape().second != _add_matrix.shape().second) {
			throw MyException("_base_matrix<_Tp>::append_single_row matrix col num is unfit");
		}
		_base_matrix tmp = __matrix;
		int m = _add_matrix.shape().first;
		for (int i = 0; i < m; ++i) {
			tmp._matrix.push_back(_add_matrix._matrix[i]);
		}
		return tmp;
	}

	void append_single_row(_base_matrix _add_matrix) {
		if (shape().second != _add_matrix.shape().second) {
			throw MyException("_base_matrix<_Tp>::append_single_row matrix col num is unfit");
		}
		
		int m = _add_matrix.shape().first;
		for (int i = 0; i < m; ++i) {
			_matrix.push_back(_add_matrix._matrix[i]);
		}
	}

	//在矩阵行最后添加一行相同的数字
	static _base_matrix append_single_row_number(_base_matrix __matrix,int _number) {
		int col_num = __matrix.shape().second;
		_base_matrix _tmp(1, col_num);
		for (int i = 0; i < col_num; ++i) {
			_tmp._matrix[0][i] = _number;
		}
		_base_matrix tmp = __matrix;
		tmp.append_single_row(_tmp);
		return tmp;
	}

	//在矩阵列最后添加一个矩阵
	static _base_matrix append_single_col(_base_matrix __matrix,_base_matrix _add_matrix) {

	}

	//在矩阵列最后添加一矩阵
	void _append_single_col(_base_matrix __matrix) {
		if (__matrix.shape().second != 1) {
			throw MyException("_base_matrix<_Tp>::_append_single_col the __matrix shape:col num is invalid.");
		}
		if (__matrix.shape().first != _matrix.size()) {
			throw MyException("_base_matrix<_Tp>::_append_single_col the __matrix shape:row num is invalid.");
		}

		int m = _matrix.size();
		for (int i = 0; i < m; ++i) {
			_matrix[i].push_back(__matrix._matrix[i][0]);
		}
	}

	//在矩阵列最后添加一个相同的数字
	static _base_matrix append_single_col_number(_base_matrix __matrix, int _number) {

	}

	void append_single_col_number(int _number) {

	}

	//添加松弛变量
	void _append_loose_pos(int arg_id) {
		_append_single_col_single_num(arg_id, 1);
	}

	void _append_loose_neg(int arg_id) {
		_append_single_col_single_num(arg_id, -1);
	}

	//添加人工变量
	void _append_arti(int arg_id) {
		_append_single_col_single_num(arg_id, 1);
	}

	void _append_single_col_single_num(int arg_id, double suffix) {
		//根据第几行，加入松弛变量
		int m = _matrix.size();
		if (arg_id < 0 || arg_id >= m) {
			stringstream ss;
			ss << "_base_matrix::append_loose the arg_id ";
			ss << arg_id << "is out of bounds " << m - 1;
			throw MyException(ss.str());
		}
		_base_matrix tmp(m, 1);
		for (int i = 0; i < m; ++i) {
			tmp._matrix[i][0] = 0;
		}
		tmp._matrix[arg_id][0] = suffix;
		_append_single_col(tmp);
	}

	static _base_matrix insert_single_row(_base_matrix __matrix,int row_id, _base_matrix _insert_matrix) {

	}

	void insert_single_row(int row_id, _base_matrix _insert_matrix) {

	}

	static _base_matrix insert_single_col(_base_matrix __matrix,int col_id, _base_matrix _insert_matrix) {

	}

	void insert_single_col(int col_id, _base_matrix _insert_matrix) {

	}

	//取元素
	double get(int row, int col) {
		if (row < 0 ||  row >= _matrix.size()) {
			string tmp = "_base_matrix<_Tp>::get_elem  row range" + to_string(row) + "out of bounds(0-" + to_string(_matrix.size()) + ").";
			throw MyException(tmp);
		}
		if (col < 0 || col >= _matrix[0].size()) {
			string tmp = "_base_matrix<_Tp>::get_elem  col range" + to_string(col) + "out of bounds(0-" + to_string(_matrix[0].size()) + ").";
			throw MyException(tmp);
		}

		return _matrix[row][col];
	}

	//高斯消元
	void gauss_elimination(pair<int, int> choose_id) {
		int id_row = choose_id.first;
		int id_col = choose_id.second;
		int row = _matrix.size();
		int col = _matrix[0].size();
		//vector<vector<_Tp>> _tmp_matrix = _matrix;
		//先将这一行归一
		double axis_elem = _matrix[id_row][id_col];
		for (int i = 0; i < col; ++i) {
			_matrix[id_row][i] /= axis_elem;
		}

		//消元
		for (int i = 0; i < id_row; ++i) {
			double tmp = _matrix[i][id_col];
			for (int j = 0; j < col; j++) {
				_matrix[i][j] = _matrix[i][j] - tmp * _matrix[id_row][j];
			}
		}

		for (int i = id_row + 1; i < row; ++i) {
			double tmp = _matrix[i][id_col];
			for (int j = 0; j < col; j++) {
				_matrix[i][j] = _matrix[i][j] - tmp * _matrix[id_row][j];
			}
		}
	}

	//复制
	void _copy(vector<vector<double>> __matrix) {
		_matrix = __matrix;
	}

	//格式化输出
	void print_form() {
		int row = shape().first;
		int col = shape().second;
		cout << "shape (" << row << "," << col << ")" << endl;
		cout << endl;
		
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				cout.width(4);
				cout << _matrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	//有时间再重载一下运算符

	//全部加/减一个数
	static _base_matrix add_number(_base_matrix __matrix, double _number) {
		int row = __matrix.shape().first;
		int col = __matrix.shape().second;
		_base_matrix _tmp_base_matrix = __matrix;
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; ++j) {
				_tmp_base_matrix._matrix[i][j] += _number;
			}
		}
		return _tmp_base_matrix;
	}

	void add_number(double _number) {
		int row = shape().first;
		int col = shape().second;

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; ++j) {
				_matrix[i][j] += _number;
			}
		}
	}

	//乘法
	static _base_matrix mul(_base_matrix a, _base_matrix b) {
		if (a.shape().second != b.shape().first) {
			stringstream ss;
			ss << "_base_matrix::mul (" << a.shape().first << "," << a.shape().second << ") dosen't fit ("
				<< b.shape().first << "," << b.shape().second << ")";
			throw MyException(ss.str());
		}
		int ar = a.shape().first;
		int ac = a.shape().second;
		int br = b.shape().first;
		int bc = b.shape().second;
		_base_matrix _tmp_matrix(ar, bc);

		for (int i = 0; i < ar; ++i) {
			for (int k = 0; k < bc; ++k) {
				double tmp = 0;
				for (int j = 0; j < ac; ++j) {
					tmp += a._matrix[i][j] * b._matrix[j][i];
				}
				_tmp_matrix._matrix[i][k] = tmp;
			}
		}

		return _tmp_matrix;
	}

	void mul(_base_matrix b) {
		if (shape().second != b.shape().first) {
			stringstream ss;
			ss << "_base_matrix::mul (" << shape().first << "," << shape().second << ") dosen't fit ("
				<< b.shape().first << "," << b.shape().second << ")";
			throw MyException(ss.str());
		}
		int ar = shape().first;
		int ac = shape().second;
		int br = b.shape().first;
		int bc = b.shape().second;
		_base_matrix _tmp_matrix(ar, bc);

		for (int i = 0; i < ar; ++i) {
			for (int k = 0; k < bc; ++k) {
				double tmp = 0;
				for (int j = 0; j < ac; ++j) {
					tmp += _matrix[i][j] * b._matrix[j][i];
				}
				_tmp_matrix._matrix[i][k] = tmp;
			}
		}

		_matrix = _tmp_matrix._matrix;
	}

	//加法
	static _base_matrix add(_base_matrix a, _base_matrix b) {
		if (a.shape().second != b.shape().second && a.shape().first != b.shape().first) {
			stringstream ss;
			ss << "_base_matrix::add (" << a.shape().first << "," << a.shape().second << ") dosen't fit ("
				<< b.shape().first << "," << b.shape().second << ")";
			throw MyException(ss.str());
		}
		int ar = a.shape().first;
		int ac = a.shape().second;
		_base_matrix _tmp_matrix(ar, ac);
		for (int i = 0; i < ar; ++i) {
			for (int j = 0; j < ac; ++j) {
				_tmp_matrix._matrix[i][j] = a._matrix[i][j] + b._matrix[i][j];
			}
		}

		return _tmp_matrix;
	}

	void add(_base_matrix b) {
		if (shape().second != b.shape().second && shape().first != b.shape().first) {
			stringstream ss;
			ss << "_base_matrix::add (" << shape().first << "," << shape().second << ") dosen't fit ("
				<< b.shape().first << "," << b.shape().second << ")";
			throw MyException(ss.str());
		}
		int ar = shape().first;
		int ac = shape().second;
		for (int i = 0; i < ar; ++i) {
			for (int j = 0; j < ac; ++j) {
				_matrix[i][j] += b._matrix[i][j];
			}
		}
	}

	//减法
	static _base_matrix sub(_base_matrix a, _base_matrix b) {
		if (a.shape().second != b.shape().second && a.shape().first != b.shape().first) {
			stringstream ss;
			ss << "_base_matrix::sub (" << a.shape().first << "," << a.shape().second << ") dosen't fit ("
				<< b.shape().first << "," << b.shape().second << ")";
			throw MyException(ss.str());
		}
		int ar = a.shape().first;
		int ac = a.shape().second;
		_base_matrix _tmp_matrix(ar, ac);
		for (int i = 0; i < ar; ++i) {
			for (int j = 0; j < ac; ++j) {
				_tmp_matrix._matrix[i][j] = a._matrix[i][j] - b._matrix[i][j];
			}
		}

		return _tmp_matrix;
	}

	void sub(_base_matrix b) {
		if (shape().second != b.shape().second && shape().first != b.shape().first) {
			stringstream ss;
			ss << "_base_matrix::sub (" << shape().first << "," << shape().second << ") dosen't fit ("
				<< b.shape().first << "," << b.shape().second << ")";
			throw MyException(ss.str());
		}
		int ar = shape().first;
		int ac = shape().second;
		for (int i = 0; i < ar; ++i) {
			for (int j = 0; j < ac; ++j) {
				_matrix[i][j] -= b._matrix[i][j];
			}
		}
	}

	//取反
	static _base_matrix neg(_base_matrix a) {
		_base_matrix _tmp_matrix = a;
		int row = _tmp_matrix.shape().first;
		int col = _tmp_matrix.shape().second;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				_tmp_matrix._matrix[i][j] = -_tmp_matrix._matrix[i][j];
			}
		}
		return _tmp_matrix;
	}

	void neg() {
		int row = shape().first;
		int col = shape().second;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				_matrix[i][j] = -_matrix[i][j];
			}
		}
	}

	bool operator==(_base_matrix& b) {
		if (shape() != b.shape()) {
			return false;
		}

		int row = shape().first;
		int col = shape().second;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				if (fabs(_matrix[i][j] - b._matrix[i][j]) > 1e-6) {
					return false;
				}
			}
		}
		return true;
	}

	_base_matrix operator+(_base_matrix& b) {
		return add(*this, b);
	}

	_base_matrix operator-(_base_matrix& b) {
		return sub(*this, b);
	}

	_base_matrix operator*(_base_matrix& b) {
		return mul(*this, b);
	}

	_base_matrix operator+=(_base_matrix& b) {
		this->add(b);
		return *this;
	}

	_base_matrix operator-=(_base_matrix& b) {
		this->sub(b);
		return *this;
	}

	_base_matrix operator*=(_base_matrix& b) {
		this->mul(b);
		return *this;
	}

	bool operator!=(_base_matrix& b) {
		return !(*this == b);
	}

	ostream& operator<<(ostream& cout) {
		int row = shape().first;
		int col = shape().second;
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				cout.width(OUT_WIDTH);
				cout << fixed << setprecision(3) << _matrix[i][j];
			}
			cout << endl;
		}
		return cout;
	}

};

class Simplex_Matrix :public _base_matrix {
public:
	vector<int> nonbase;		//非基变量下标
	vector<int> base;			//基变量下标
	vector<int> arti;			//人工变量下标
	vector<int> loose;			//松弛变量下标
	int arg_num;				//总的变量数

	Simplex_Matrix(){
		
	}

	Simplex_Matrix(const Simplex_Matrix& __sim_matrix):_base_matrix(__sim_matrix._matrix){

	}

	Simplex_Matrix(vector<Equation> equs):_base_matrix(equs.size(), equs[0].arg_num + 1) {
		//通过Eqution初始化Simplex_Matrix
		int m = equs.size();
		//在这里加入检测是否需要添加人工变量
		if (m <= 0) {
			throw MyException("Simplex_Matrix::Simplex_Matrix equations is null.");
		}
		arg_num = equs[0].arg_num;
		//先将初始的矩阵写入
		for (int i = 0; i < m; ++i) {
			_matrix[i] = equs[i].w;
		}

		//在根据符号添加松弛变量
		for (int i = 0; i < m; ++i) {
			if (equs[i].op == 1) {
				//小于等于号，+松弛变量
				_append_loose_pos(i);
				loose.push_back(arg_num++);

			}
			else if (equs[i].op == 2) {
				//大于等于号，-松弛变量
				_append_loose_neg(i);
				loose.push_back(arg_num++);
			}
		}
		//cout << "添加松弛变量后" << endl;
		//print_form();
		//松弛变量添加完毕
		//现在检查矩阵中是否有可行的单位矩阵，如果没有，则添加人工变量
		int こんにちは = 0;
		check_normal_matrix();

		//cout << "添加人工变量后" << endl;
		//print_form();
		//计算出非基变量
		nonbase = get_nonbase();

		//添加完毕，最后在加上b变为增广矩阵即可
		for (int i = 0; i < m; ++i) {
			_matrix[i].push_back(equs[i].b);
		}
	}

	void check_normal_matrix() {
		int row_size = _matrix.size();
		int col_size = _matrix[0].size();

		bool *exist = new bool[row_size];
		for (int i = 0; i < row_size; ++i) {
			exist[i] = false;
		}
		for (int i = 0; i < row_size; ++i) {
			for (int j = 0; j < col_size; ++j) {
				if (_matrix[i][j] == 1) {
					int cnt = 0;
					for (int k = 0; k < row_size; ++k) {
						if (_matrix[k][j] == 0) cnt++;
					}
					if (cnt == row_size-1 && !exist[i]) {
						//刚好这一列只有一个1,其余全为0
						base.push_back(j);
						exist[i] = true;
					}
				}
			}
		}
		sort(base.begin(), base.end());
		bool need_arti = (base.size() != row_size);
		if (need_arti) {
			//需要添加人工变量
			vector<int> tmp_to_append_raw_idx = get_to_append_raw_idx(exist);
			for (vector<int>::iterator it = tmp_to_append_raw_idx.begin(); it != tmp_to_append_raw_idx.end(); it++) {
				_append_arti(*it);//在最后一列加上一列第*it行元素为1的列向量
				base.push_back(arg_num);
				arti.push_back(arg_num);
				arg_num++;
			}
		}
		delete[]exist;
	}

	vector<int> get_to_append_raw_idx(bool* exist) {
		vector<int> v;
		int raw_size = shape().first;
		for (int i = 0; i < raw_size; ++i) {
			if (!exist[i]) {
				v.push_back(i);
			}
		}
		return v;
	}

	vector<int> get_nonbase() {
		return get_vector_not_in_vector(base, 0, _matrix[0].size());
	}

	vector<int> get_vector_not_in_vector(vector<int> v, int begin, int end) {
		vector<int> tmp;
		int m = v.size();
		//左闭右开
		for (int i = begin; i < end; ++i) {
			if (!elem_is_in_vector(v, i)) {
				tmp.push_back(i);//v中没有该元素，则添加该元素
			}
		}
		return tmp;
	}

	vector<vector<double>> get_matrix() {
		return _matrix;
	}

	template<typename T>
	bool elem_is_in_vector(vector<T> v, T elem) {
		int m = v.size();
		for (int i = 0; i < m;++i) {
			if (v[i] == elem) {
				return true;
			}
		}
		return false;
	}

	void print() {
		cout << "Base:" << endl;
		print_vector(base, 3);
		cout << "Nonbase:" << endl;
		print_vector(nonbase,3);
		cout << "Loose:" << endl;
		print_vector(loose,3);
		cout << "Arti:" << endl;
		print_vector(arti,3);
		print_form();
	}

	template<typename T>
	static void print_vector(vector<T> v, int width, char split = '\0',bool with_huanhang = true) {
		int m = v.size();
		const int presicions = 2;
		for (int i = 0; i < m; ++i) {
			cout.width(width);
			cout <<fixed << setprecision(presicions) << v[i];
			cout << split;
		}
		if(with_huanhang)
		cout << endl;
	}
};




enum Res_Type{Only_optima,Alternative_optima,Unbounded,Infeasible,Other};

class Result {
public:
	Res_Type res_type;					//结果类型
	Simplex_Matrix Wxb;	//约束条件系数矩阵
	vector<double> Wz;			//目标函数系数
	vector<int> nonbase;		//非基变量下标
	vector<int> base;			//基变量下标
	vector<int> arti;			//人工变量下标
	vector<int> loose;			//松弛变量下标
	vector<int> init;			//初始时就有的变量的下标
	vector<double> sigma;		//检验数
	vector<double> theta;		//增量
	vector<double> X;			//当前解
	double Z;

	Result() {
		res_type = Other;
	}
	Result(Res_Type _type) {
		res_type = _type;
	}
	Result(Res_Type _res_type, Simplex_Matrix _Wxb, vector<double> _Wz,
		vector<int> _nonbase, vector<int> _base, vector<int> _arti, vector<int> _loose, vector<int> _init,
		vector<double> _sigma, vector<double> _theta) {
		res_type = _res_type;
		Wxb = _Wxb;
		Wz = _Wz;
		nonbase = _nonbase;
		base = _base;
		arti = _arti;
		loose = _loose;
		init = _init;
		sigma = _sigma;
		theta = _theta;
		Z = 0;

		//解析出当前的解
		int m = Wz.size();
		X = _Wz;
		for (int i = 0; i < m; ++i) {
			X[i] = 0;
		}
		int lcol = Wxb.shape().second - 1;
		for (int i = 0; i < Wxb.shape().first; ++i) {
			X[base[i]] = Wxb.get(i,lcol);
		}

		if (res_type == Only_optima || res_type == Alternative_optima) {
			int X_size = X.size();
			for (int i = 0; i < X_size; ++i) {
				Z += X[i] * Wz[i];
			}
		}
	}
};

class Simplex_Model
{
public:
	Simplex_Matrix Wxb;	//约束条件系数矩阵
	vector<double> Wz;			//目标函数系数
	vector<int> nonbase;		//非基变量下标
	vector<int> base;			//基变量下标
	vector<int> arti;			//人工变量下标
	vector<int> loose;			//松弛变量下标
	vector<int> init;			//初始时有的变量的下标
	vector<double> sigma;		//检验数
	vector<double> theta;		//增量
	Result _res;				//储存结果
	vector<Equation> equtions;	//储存初始化该模型时用的方程组
	int depth;
public:
	//打印
	void fomulate_printout(vector<double> __Wz) {
		const int presicions = 2;
		int arg_num = __Wz.size();
		cout << "-------------------------------------------------------------------------"<<endl;
		cout.width(15);
		cout << "C     ";
		cout << "| ";
		Simplex_Matrix::print_vector(__Wz,6);
		cout << "-------------------------------------------------------------------------" << endl;
		cout << " Cb | Xb |  B  |";
		for (int i = 0; i < arg_num; i++) {
			cout.width(5);
			cout << "X" << i;
		}
		cout << endl;
		cout << "-------------------------------------------------------------------------" << endl;
		int row_size = Wxb.shape().first;
		for (int i = 0; i < row_size; ++i) {
			cout.width(4);
			cout <<fixed<<setprecision(presicions) << __Wz[base[i]];
			cout << "| X" << base[i] << " |";
			cout.width(5);
			cout << fixed << setprecision(presicions) << Wxb.get(i, arg_num);
			cout << "|";
			for (int j = 0; j < arg_num; j++) {
				cout.width(6);
				cout.right;
				cout <<fixed << setprecision(presicions) << Wxb.get(i, j);
			}
			cout << endl;
		}
		cout << "-------------------------------------------------------------------------" << endl;
		cout.width(15);
		cout << "Sigma";
		cout << "|";
		Simplex_Matrix::print_vector(sigma,6);
		cout << "-------------------------------------------------------------------------" << endl;
	}



	//解决问题
	void solve() {
		Result res;
		if (need_two_phase()) {
			res = two_phase(Wz,Wxb);	//两阶段法
		}
		else {
			res = iter_cal(Wz,Wxb);		//单纯形法迭代
		}
		_res = res;
		print_res(res);
	}

	void print_res(Result res) {
		Res_Type res_type = res.res_type;
		switch (res_type)
		{
		case Only_optima:
			res_only(res);	//第0种情况，有唯一最优解
			break;

		case Alternative_optima:
			res_alter(res);	//第1种情况，有无穷多最优解
			break;

		case Unbounded:
			res_unbound(res);	//第2种情况，无界解
			break;

		case Infeasible:
			res_infeasi(res);	//第3种情况，无可行解
			break;

		default:
			res_null(res);
			break;	//否则继续
		}
	}

	void print_x_z_res(Result res) {
		cout << "X* = [";
		Simplex_Matrix::print_vector(res.X, 5, ',', false);
		cout << "]" << endl;
		cout << "Z* = " << res.Z << endl;
	}

	void res_only(Result res) {
		cout << "唯一最优解" << endl;
		print_x_z_res(res);
	}

	void res_alter(Result res) {
		cout << "无穷多最优解" << endl;
		print_x_z_res(res);
	}

	void res_unbound(Result res) {
		cout << "无界解" << endl;
		cout << "X* is unbound" << endl << "Z* is unbound." << endl;
	}

	void res_infeasi(Result res) {
		cout << "无可行解" << endl;
	}

	void res_null(Result res) {
		cout << "Null" << endl;
	}

	bool need_two_phase() {
		return (arti.size() > 0);
	}

	Result two_phase(vector<double>_Wz, Simplex_Matrix _Wbx) {
		vector<double> two_phase_Wz = init_two_phase_Wz();
		cout << "需要用两阶段法\n第一阶段" << endl;
		Result res1 = iter_cal(two_phase_Wz, _Wbx);
		if (detect_res(res1)) {
			//不满足,即人工变量不为0
			//res_3();	//无可行解
			return Result(Infeasible, _Wbx, _Wz, nonbase, base, arti, loose, init, vector<double>(), vector<double>());
		}	

		res1.Wz = _Wz;
		update_base_res(res1);	//更新基可行解
		cout << "第二阶段" << endl;
		Wz = res1.Wz;
		Wxb = res1.Wxb;
		Result res2 = iter_cal(res1.Wz, res1.Wxb);
		return res2;
	}
	
	Result iter_cal(vector<double>_Wz, Simplex_Matrix _Wbx) {
		if (_Wz.size() <= 0 || _Wbx.shape().first <= 0) {
			string tmp = "Simplex_Model::iter_cal _Wz.size or _Wbx.shape is non-positive.";
			throw MyException(tmp);
		}
		int last_out_arg_id = -1;
		int cnt = 0;
		while (1) {
			sigma = cal_sigma(_Wz,_Wbx);	//计算检验数
			Res_Type res_type = detect_cases(_Wbx,sigma);
			fomulate_printout(_Wz);
			
			switch (res_type)
			{
			case Only_optima:
				//res_0();	//第0种情况，有唯一最优解
				return Result(Only_optima,_Wbx,_Wz,nonbase,base,arti,loose,init,sigma,theta);
				break;

			case Alternative_optima:
				//res_1();	//第1种情况，有无穷多最优解
				return Result(Alternative_optima, _Wbx, _Wz, nonbase, base, arti, loose, init,sigma, theta);
				break;

			case Unbounded:
				//res_2();	//第2种情况，无界解
				return Result(Unbounded, _Wbx, _Wz, nonbase, base, arti, loose, init,sigma, theta);
				break;

			case Infeasible:
				//res_3();	//第3种情况，无可行解
				return Result(Infeasible, _Wbx, _Wz, nonbase, base, arti, loose, init,sigma, theta);
				break;

			default:
				break;	//否则继续
			}
			theta = cal_theta(_Wz, _Wbx, sigma);
			int max_theta_idx = find_vector_max_idx(theta);
			if (theta[max_theta_idx] < 0) {
				//所有的theta都小于0，问题为无界解
				return Result(Unbounded, _Wbx, _Wz, nonbase, base, arti, loose, init, sigma, theta);
			}
			pair<int, int> change_id = choose_change_id(sigma,theta);
			if (change_id.second == last_out_arg_id) {
				//这一次换进去的变量是上一次被换出的变量
				if (get_vector_pos_num(theta) > 1) {
					//如果还有次大的，选择次大的
					theta = cal_theta(_Wz, _Wbx, sigma, true);
					change_id = choose_change_id(sigma, theta, true);
				}
				else {
					return Result(Alternative_optima, _Wbx, _Wz, nonbase, base, arti, loose, init, sigma, theta);
				}
			}
			last_out_arg_id = base[change_id.first];//上一次被换出的id
			update_base(change_id);	//更新基变量
			_Wbx.gauss_elimination(change_id);
			Wxb = _Wbx;
			cnt++;
			if (cnt > 100) {
				return Result(Infeasible, _Wbx, _Wz, nonbase, base, arti, loose, init, sigma, theta);
			}
		}
	}

	int get_vector_pos_num(vector<double> v) {
		int cnt = 0;
		for (vector<double>::iterator it = v.begin(); it != v.end(); it++) {
			if (*it > 0) {
				cnt++;
			}
		}
		return cnt;
	}

	int find_vector_second_max_idx(vector<double> v) {
		int m = v.size();
		double max = v[0];
		int max_id = 0;
		int second_max_id = 0;
		double second_max = 0;
		for (int i = 1; i < m; ++i) {
			if (v[i] > max) {
				max = v[i];
				max_id = i;
			}
		}
		for (int i = 1; i < m; ++i) {
			if (v[i] > second_max && i != max_id) {
				second_max = v[i];
				second_max_id = i;
			}
		}
		return second_max_id;
	}

	vector<double> cal_theta(vector<double>&_Wz, Simplex_Matrix& _Wbx, vector<double>& sigma, bool second_max = false) {
		int sigma_max_idx;
		if (second_max) {
			sigma_max_idx = find_vector_second_max_idx(sigma);
		}
		else {
			sigma_max_idx = find_vector_max_idx(sigma);
		}
		int row = _Wbx.shape().first;
		int col = _Wz.size();
		vector<double> tmp(row);
		for (int i = 0; i < row; ++i) {
			if (fabs(_Wbx.get(i, sigma_max_idx)) <= 1e-6) {
				tmp[i] = -1;
			}
			else {
				tmp[i] = _Wbx.get(i, col) / _Wbx.get(i, sigma_max_idx);
			}
		}

		return tmp;
	}

	vector<double> cal_sigma(vector<double> _Wz, Simplex_Matrix _Wbx) {
		int col = _Wz.size();
		int row = _Wbx.shape().first;
		vector<double> tmp = _Wz;
		for (int i = 0; i < col; ++i) {
			for (int j = 0; j < row; ++j) {
				tmp[i] -= _Wz[base[j]] * _Wbx.get(j,i);
			}
		}
		return tmp;
	}

	pair<int, int> choose_change_id(vector<double> sigma, vector<double> theta, bool second_max = false) {
		int max_sigma_idx;
		if (second_max) {
			max_sigma_idx = find_vector_second_max_idx(sigma);
		}
		else {
			max_sigma_idx = find_vector_max_idx(sigma);
		}
		vector<int> theta_pos_idx;
		int theta_size = theta.size();
		for (int i = 0; i < theta_size;++i) {
			if (theta[i] >= 0) {
				theta_pos_idx.push_back(i);
			}
		}
		if (!theta_pos_idx.size()) {
			return pair<int, int>(-1, -1);
		}
		double min = theta[theta_pos_idx[0]];
		int theta_pos_min_idx = 0;
		for (int i = 1; i < theta_pos_idx.size(); ++i) {
			if (theta[theta_pos_idx[i]] < min) {
				min = theta[theta_pos_idx[i]];
				theta_pos_min_idx = i;
			}
		}

		return pair<int, int>(theta_pos_idx[theta_pos_min_idx], max_sigma_idx);//行列
	}

	void update_base(pair<int, int> change_id) {
		//换出基，换入基
		int tmp = base[change_id.first];
		base[change_id.first] = change_id.second;
		for (vector<int>::iterator it = nonbase.begin(); it != nonbase.end(); it++) {
			if (*it == change_id.second) {
				nonbase.erase(it);
				break;
			}
		}
		nonbase.push_back(tmp);
	}

	Res_Type detect_cases(Simplex_Matrix _Wbx, vector<double> theta) {
		if (detect_only_optima_case(_Wbx, theta)) {
			return Only_optima;
		}
		else if (detect_alternative_optima_case(_Wbx, theta)){
			return Alternative_optima;
		}
		else if(detect_unbounded_case(_Wbx, theta)){
			return Unbounded;
		}
		else if (detect_infeasible_case(_Wbx, theta)) {
			return Infeasible;
		}
		return Other;
	}

	bool detect_only_optima_case(Simplex_Matrix _Wbx, vector<double> theta) {
		int m = theta.size();
		
		//确保检验数都小于等于0
		for (int i = 0; i < m; ++i) {
			if (theta[i] > 0) {
				return false;
			}
		}

		//确保非基变量的检验数不等于0
		int nbm = nonbase.size();
		for (int i = 0; i < nbm; ++i) {
			if (theta[nonbase[i]] == 0) {
				return false;
			}
		}

		return true;
	}

	bool detect_alternative_optima_case(Simplex_Matrix _Wbx, vector<double> theta) {
		int m = theta.size();
		
		//确保检验数都小于等于0
		for (int i = 0; i < m; ++i) {
			if (theta[i] > 0) {
				return false;
			}
		}

		//确保非基变量的检验数有等于0的
		int bm = base.size();
		for (int i = 0; i < bm; ++i) {
			if (theta[nonbase[i]] == 0) {
				return true;
			}
		}

		return false;
	}

	bool detect_unbounded_case(Simplex_Matrix _Wbx, vector<double> theta) {
		int max_idx = find_vector_max_idx<double>(theta);
		int row = _Wbx.shape().first;
		if (theta[max_idx] <= 0) {
			return false;
		}
		
		//最大的检验数大于0，作为换入变量，检查该变量的所有系数是否为负
		for (int i = 0; i < row; i++) {
			if (_Wbx.get(i, max_idx) > 0) {
				return false;
			}
		}

		//系数均为非正
		return true;
	}

	bool detect_infeasible_case(Simplex_Matrix _Wbx, vector<double> theta) {
		
		return false;
	}

	//找到vector中最大的元素，返回其下标.
	template<typename M_TP>
	int find_vector_max_idx(vector<M_TP> v) {
		if (v.size() <= 0) {
			throw MyException("Simplex_Model::find_vector_max_idx v.size() <= 0 wrong.");
		}
		int max_idx = 0;
		int max = v[0];
		int m = v.size();
		for (int i = 1; i < m; ++i) {
			if (v[i] > max) {
				max = v[i];
				max_idx = i;
			}
		}
		return max_idx;
	}
	vector<double> init_two_phase_Wz() {
		vector<double> tmp(Wz.size());
		for (int i = 0; i < tmp.size(); ++i) {
			tmp[i] = 0;
		}
		for (int i = 0; i < arti.size(); ++i) {
			tmp[arti[i]] = -1.0;
		}
		return tmp;
	}

	bool detect_res(Result res) {
		for (int i = 0; i < res.arti.size(); ++i) {
			if (res.X[arti[i]] != 0) {
				return true;
			}
		}
		return false;
	}

	bool vector_one_row_is_all_this_num(vector<double> v,double num) {
		int m = v.size();
		for (int i = 0; i < m; ++i) {
			if (fabs(v[i] - num) > 1e-6) {
				return false;
			}
		}
		return true;
	}

	bool check_mxb_arti_row_is_zero(Simplex_Matrix Mxb, int row) {
		if (row < 0 || row >= Mxb.shape().first) {
			throw MyException("Simplex_Model::check_mxb_arti_row_is_zero row is out of range!");
		}
		return vector_one_row_is_all_this_num(Mxb.get_matrix()[row], 0);
	}

	void update_base_res(Result& res) {
		vector<int> arti_in_base_id;//里面的是arti在base里的id
		//人工变量是否在基变量中
		for (int i = 0; i < res.arti.size();++i) {
			for (int j = 0; j < res.base.size();++j) {
				if (arti[i] == base[j]) {
					arti_in_base_id.push_back(j);
				}
			}
		}
		if (arti_in_base_id.size()) {
			//有人工变量是基变量，且当前取值为0
			for (int i = 0; i < arti_in_base_id.size(); ++i) {
				if (check_mxb_arti_row_is_zero(res.Wxb, arti_in_base_id[i])) {
					//该行全为0
					res.Wxb._del_single_row(arti_in_base_id[i]);
				}
				else {
					//该行不全为0
					//强制让一个非基，非人工变量入基
					update_base(pair<int, int>(arti_in_base_id[i], res.nonbase[rand() % res.nonbase.size()]));
				}
			}
		}
		
		//计算得出非人工变量的下标
		vector<int> nonarti;
		int m = res.Wz.size();
		for (int i = 0; i < m; ++i) {
			if (elem_in_vector<int>(res.arti, i)) {

			}
			else{
				nonarti.push_back(i);
			}
		}
		vector<double> new_Wz;
		for (int i = 0; i < nonarti.size(); ++i) {
			new_Wz.push_back(res.Wz[nonarti[i]]);
		}
		//得到没有人工变量的Wz
		res.Wz = new_Wz;

		//得到没有人工变量的base
		int base_size = base.size();
		int nonarti_size = nonarti.size();
		vector<int> new_base;
		for (int i = 0; i < base_size; ++i) {
			for (int j = 0; j < nonarti_size; ++j) {
				if (base[i] == nonarti[j]) {
					new_base.push_back(base[i]);
					break;
				}
			}
		}

		res.base = new_base;
		base = new_base;
		
		//Simplex_Matrix<double> tmp_matrix = res.Wxb;
		res.Wxb._del_cols(arti);

		//修改Nonbase
		vector<int> tmp_nonbase;
		for (int i = 0; i < nonbase.size(); ++i) {
			if (elem_in_vector(nonarti, nonbase[i])) {
				tmp_nonbase.push_back(nonbase[i]);
			}
		}
		nonbase = tmp_nonbase;

		//修改人工变量
		arti.clear();
	}
	template<typename T>
	bool elem_in_vector(vector<T>& v, T elem){
		int m = v.size();
		for (int i = 0; i < m; ++i) {
			if (elem == v[i]) {
				return true;
			}
		}
		return false;
	}
};

class Simplex_Reader {
public:
	string file_path;

	Simplex_Reader(string _str) {
		file_path = _str;
	}

	void readin_info(Simplex_Model& _simplex_model) {
		vector<Equation> equs;
		vector<string> res;
		ifstream infile(file_path.c_str());
		if (!infile.is_open()) {
			throw MyException("_base_reader::readlines file open failed.");
		}
		int optimise;
		infile >> optimise;
		string Wz_str;
		getline(infile, Wz_str);//去掉换行符
		getline(infile, Wz_str);
		vector<double> tmp_wz_v = Equation::str_to_vector(Wz_str);
		if (optimise == 0) {
			for (int i = 0; i < tmp_wz_v.size(); ++i) {
				tmp_wz_v[i] = -tmp_wz_v[i];
			}
		}
		_simplex_model.Wz = tmp_wz_v;
		int arg_num, equ_num;
		infile >> arg_num >> equ_num;
		//更新最初的变量下标
		for (int i = 0; i < arg_num; i++) {
			_simplex_model.init.push_back(i);
		}
		string buf;
		getline(infile, buf);//去掉换行符
		for (int i = 0; i < equ_num; ++i) {
			string equ_str;
			getline(infile, equ_str);
			Equation tmp_equ(equ_str);
			equs.push_back(tmp_equ);
		}
		Simplex_Matrix tmp_simplex_matrix(equs);
		_simplex_model.Wxb = tmp_simplex_matrix;
		_simplex_model.arti = tmp_simplex_matrix.arti;
		_simplex_model.base = tmp_simplex_matrix.base;
		_simplex_model.nonbase = tmp_simplex_matrix.nonbase;
		_simplex_model.loose = tmp_simplex_matrix.loose;
		_simplex_model.equtions = equs;
		int total = _simplex_model.loose.size() + _simplex_model.arti.size();
		for (int i = 0; i < total; ++i) {
			_simplex_model.Wz.push_back(0);
		}
	}

	static void readin_branch_info(Simplex_Model& _simplex_model, vector<Equation> equations, vector<double> _Wz) {
		_simplex_model.Wz = _Wz;
		Simplex_Matrix tmp_simplex_matrix(equations);
		_simplex_model.Wxb = tmp_simplex_matrix;
		_simplex_model.arti = tmp_simplex_matrix.arti;
		_simplex_model.base = tmp_simplex_matrix.base;
		_simplex_model.nonbase = tmp_simplex_matrix.nonbase;
		_simplex_model.loose = tmp_simplex_matrix.loose;
		_simplex_model.equtions = equations;

		//更新最初的变量下标
		for (int i = 0; i < equations[0].arg_num; i++) {
			_simplex_model.init.push_back(i);
		}

		int total = tmp_simplex_matrix.arg_num - _Wz.size();
		for (int i = 0; i < total; ++i) {
			_simplex_model.Wz.push_back(0);
		}
	}
};

class Integer_Model {
private:
	queue<Simplex_Model> branch_queue;
	vector<Simplex_Model> res_storage[100];
	string file_path;
	Result _integer_optimise_res;
	static const int MAX_SEARCH_DEPTH = 5;//最多搜索层数
public:
	Integer_Model() {
		
	}

	Integer_Model(string _str) {
		file_path = _str;
		
	}

	void print_res() {
		Simplex_Model sm;
		sm.print_res(_integer_optimise_res);
	}

	void clear_queue() {
		while (!branch_queue.empty()) {
			branch_queue.pop();
		}
	}

	void solve() {
		clear_queue();//先将队列清空
		Simplex_Model _simplex_model1;
		_simplex_model1.depth = 0;
		Simplex_Reader _simplex_reader(file_path);
		_simplex_reader.readin_info(_simplex_model1);
		_simplex_model1.solve();
		Result _first_res = _simplex_model1._res;
		if (_first_res.res_type == Unbounded 
			|| _first_res.res_type == Infeasible) {
			_integer_optimise_res = _first_res;
			return;
		}
		if (res_is_integer(_first_res)) {
			_integer_optimise_res = _first_res;
			return;
		}
		double low_bound = -1e308;
		double up_bound = _first_res.Z;
		Result _optimise_res;
		create_branch(_simplex_model1);

		while (!branch_queue.empty()) {
			Simplex_Model _tmp_model = branch_queue.front();
			branch_queue.pop();
			cout << "queue num : " << branch_queue.size() << endl;
			if (_tmp_model.depth >= MAX_SEARCH_DEPTH) {
				continue;
			}
			_tmp_model.solve();
			Result _tmp_res = _tmp_model._res;
			
			if (_tmp_res.res_type != Only_optima && _tmp_res.res_type != Alternative_optima) {
				//无可行解
				continue;
			}

			if (_tmp_res.Z > up_bound) {
				continue;
			}

			//结果是唯一解
			if (_tmp_res.Z <= low_bound) {
				continue;
			}

			if (res_is_integer(_tmp_res) && _tmp_res.Z > low_bound) {
				low_bound = _tmp_res.Z;
				_optimise_res = _tmp_res;//更新最优解
				continue;
			}

			create_branch(_tmp_model);
		}
		_integer_optimise_res = _optimise_res;
	}

	void create_branch(Simplex_Model _model) {
		
		//找到待分割的变量noninteger_idx的上和下取整
		int noninteger_idx = find_random_non_integer_num_idx_init_res(_model._res);
		int upper = (int)_model._res.X[noninteger_idx] + 1;
		int lower = (int)_model._res.X[noninteger_idx];

		//开始分别添加条件
		Simplex_Model _upper_model;
		_upper_model.depth = _model.depth + 1;
		vector<Equation> _upper_equs = _model.equtions;
		_upper_equs.push_back(Equation(_model.init.size(), 2, noninteger_idx, upper));
		Simplex_Reader::readin_branch_info(_upper_model, _upper_equs, _model.Wz);

		Simplex_Model _lower_model;
		_upper_model.depth = _model.depth + 1;
		vector<Equation> _lower_equs = _model.equtions;
		_lower_equs.push_back(Equation(_model.init.size(), 1, noninteger_idx, lower));
		Simplex_Reader::readin_branch_info(_lower_model, _lower_equs, _model.Wz);

		branch_queue.push(_upper_model);
		branch_queue.push(_lower_model);
	}

	int find_random_non_integer_num_idx_init_res(Result res) {
		int m = res.init.size();
		vector<int> non_integer_idx;
		for (int i = 0; i < m; ++i) {
			if (!double_elem_is_interger(res.X[res.init[i]])) {
				non_integer_idx.push_back(res.init[i]);
			}
		}

		//从不是整数的下标里面随机选一个
		int rand_idx = non_integer_idx[0];
		return rand_idx;
	}

	int find_first_non_integer_num_idx_init_res(Result res) {
		
		int m = res.init.size();
		for (int i = 0; i < m;++i) {
			if (!double_elem_is_interger(res.X[res.init[i]])) {
				return res.init[i];
			}
		}
		throw MyException("Integer_Model::find_first_non_integer_num_idx_init_res no non-integer number in the vector.");
	}

	bool res_is_integer(Result res) {
		return vector_elem_in_vector_is_interger(res.X, res.init);
	}

	bool vector_elem_in_vector_is_interger(vector<double>v,vector<int> list) {
		int mv = v.size();
		int ml = list.size();
		for (int i = 0; i < ml; ++i) {
			if (list[i] < mv) {
				//没有越界才进行判断，否则略过
				if (!double_elem_is_interger(v[list[i]])) {
					return false;
				}
			}
		}
		return true;
	}

	//左闭右开
	bool range_elem_in_vector_is_interger(vector<double> v, int begin_idx, int end_idx) {
		int m = v.size();
		int i = begin_idx < 0 ? 0 : begin_idx;
		int _end_idx = end_idx <= m ? end_idx : m;
		for (; i < _end_idx; ++i) {
			if (!double_elem_is_interger(v[i])) {
				return false;
			}
		}
		return true;
	}

	vector<int> get_non_interger_elem_idx_vector(vector<double> v) {
		vector<int> tmp;
		int m = v.size();
		for (int i = 0; i < m; ++i) {
			if (!double_elem_is_interger(v[i])) {
				tmp.push_back(i);
			}
		}
		return tmp;
	}

	bool all_elem_in_vector_is_interger(vector<double> v) {
		int m = v.size();
		for (int i = 0; i < m; ++i) {
			if (!double_elem_is_interger(v[i])) {
				return false;
			}
		}
		return true;
	}

	bool double_elem_is_interger(double elem) {
		int tmp = (int)elem;
		if (elem - tmp <= 1e-6) {
			return true;
		}
		return false;
	}


};
