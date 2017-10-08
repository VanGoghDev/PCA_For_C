#include <iostream> 
#include <cmath> 

using namespace std;

// класс метода главных компонент 
class PCA {
private:
	int row; // количество строк 
	int column; // количество столбцов 
	int ii, jj;  // индексы максимального элемента в матрице
	double f, cosinus, sinus;
	double amax, min;
	double **arr; // динамический массив исходной выборки 
	double *mean_values; // динамический массив средних значений факторов 
	double **centered; // динамический массив центрированных значений 
	double *expectation; // математическое ожидание 
	double *dispersion; // дисперсия 
	double **covmat; // ковариационная матрица 
	double **a;
	double **a1;
	double **u;
	double **ut, **u1, **u2;
	double EPSILON;  // погрешность измерений
	double *eigvalue;  // собственные значения
public:
	void init(int i, int j) { // конструктор класса 
		EPSILON = 0.001;
		min = 10;
		row = i; // инициализируем параметры матрицы 
		column = j;

		// создаем двумерный массив исходных данных 
		arr = new double *[row]; // количество строк в массиве 
		for (int count = 0; count < row; count++)
			arr[count] = new double[column]; // и количество столбцов 

											 // создаем двумерный массив центрированных значений (так же как и исходных данных) 
		centered = new double*[row];
		for (int count = 0; count < row; count++)
			centered[count] = new double[column];

		// создаем двумерный массив представляющий собой ковариационную матрицу 
		covmat = new double*[column];
		for (int count = 0; count < column; count++)
			covmat[count] = new double[column];

		// создаем двумерный массив a
		a = new double*[column];
		for (int count = 0; count < column; count++)
			a[count] = new double[column];

		// создаем двумерный массив a1
		a1 = new double*[column];
		for (int count = 0; count < column; count++)
			a1[count] = new double[column];

		// создаем двумерный массив u
		u = new double*[column];
		for (int count = 0; count < column; count++)
			u[count] = new double[column];

		// создаем двумерный массив u1
		u1 = new double*[column];
		for (int count = 0; count < column; count++)
			u1[count] = new double[column];

		// создаем двумерный массив u2
		u2 = new double*[column];
		for (int count = 0; count < column; count++)
			u2[count] = new double[column];

		// создаем двумерный массив ut
		ut = new double*[column];
		for (int count = 0; count < column; count++)
			ut[count] = new double[column];

		// создаем массив средних значений 
		for (int count = 0; count < row; count++)
			mean_values = new double[column];

		// создаем массив мат. ожиданий 
		for (int count = 0; count < row; count++)
			expectation = new double[column];

		// создаем массив дисперсий 
		for (int count = 0; count < row; count++)
			dispersion = new double[column];

		// создаем массив собственных значений
		for (int count = 0; count < row; count++)
			eigvalue = new double[column];

		//временный фрагмент кода, отвечающий за ввод данных в матрицу-------— 
		for (int i = 0; i < row; i++)
			for (int j = 0; j < column; j++)
				cin >> arr[i][j];
		//-------------------------------------------------------------------— 

		//временный фрагмент кода, отвечающий за вывод данных —-------------— 
		cout << "Matrix:" << endl;
		for (int i = 0; i < row; i++)
		{
			cout << endl;
			for (int j = 0; j < column; j++)
				cout << arr[i][j] << " ";
		}
		cout << endl;
		//-------------------------------------------------------------------— 
	}
	double sqr(double a) {
		return a * a;
	}
	void mean_func() {//метод, считающий средние значения для каждого фактора 
		double sum;
		for (int i = 0; i < column; i++) {
			sum = 0;
			for (int j = 0; j < row; j++)
				sum += arr[j][i];
			mean_values[i] = sum / row;
		}
		//временный код---------------------------------------------------— 
		for (int i = 0; i < column; i++) {
			cout << endl;
			cout << "mean_value [" << i << "]" << "= " << mean_values[i];
		}
		cout << endl;
		//----------------------------------------------------------------— 
	}
	void centered_func() {

		for (int i = 0; i < row; i++) // центрирование значений 
			for (int j = 0; j < column; j++)
				centered[i][j] = arr[i][j];
		for (int i = 0; i < row; i++)
			for (int j = 0; j < column; j++)
				centered[i][j] -= mean_values[j];

		for (int i = 0; i < column; i++) { // мат ожидание 
			expectation[i] = 0;
			for (int j = 0; j < row; j++)
				expectation[i] += centered[j][i];
			expectation[i] = expectation[i] / row;

		}
		//временный код---------------------------------------------------— 
		cout << endl << "centered mass:";
		for (int i = 0; i < row; i++) {
			cout << endl;
			for (int j = 0; j < column; j++)
				cout << centered[i][j] << " ";
		}
		cout << endl << "expectation:";
		for (int i = 0; i < column; i++) {
			cout << expectation[i] << " ";
		}
		cout << endl;
		//----------------------------------------------------------------— 
	}
	void dispersion_func() {
		for (int i = 0; i < column; i++) {
			dispersion[i] = 0;
			for (int j = 0; j < row; j++)
				dispersion[i] += pow((centered[j][i] - expectation[i]), 2);
			dispersion[i] = dispersion[i] / (row - 1);
		}
		//временный код 
		cout << endl;
		cout << "dispersion:" << endl;
		for (int i = 0; i <
			column; i++) {
			cout << dispersion[i] << " ";
		}
		cout << endl;
		//-----------— 
	}
	void covmat_func() {
		for (int i = 0; i < column; i++)
			for (int j = 0; j < column; j++) {
				covmat[i][j] = 0;
				for (int n = 0; n < row; n++)
					covmat[i][j] += centered[n][i] * centered[n][j];
				covmat[i][j] = covmat[i][j] / (row - 1);
			}
		//временный код 
		cout << endl;
		cout << "covmat:" << endl;
		for (int i = 0; i < column; i++) {
			cout << endl;
			for (int j = 0; j < column; j++)
				cout << covmat[i][j] << " ";
		}
		//-----------— 
	}
	// процедура, которая ищет максимальный элемент в матрице a
	// и считает f, cosiuns, sinus
	/*void max_element() {
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++)
	{
	if (j > i) {
	if (fabs(this->amax) < fabs(this->a[i][j])) {
	this->amax = this->a[i][j];
	this->ii = i;
	this->jj = j;
	}
	else continue;
	}
	else continue;
	}
	this->f = 0.5 * atan((2 * this->amax) / (this->covmat[ii][ii] - this->covmat[jj][jj]));
	this->cosinus = cos(this->f);
	this->sinus = sin(this->f);
	}
	void u_init() {
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++)
	{
	this->u[i][j] = 0;
	}
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++) {
	if (i == j) { this->u[i][j] = 1; }
	if (i == ii) { this->u[i][i] = cosinus; }
	if (i == jj) { this->u[i][jj] = cosinus; }
	if ((j == jj) && (i == ii)) {
	this->u[i][j] = -sinus;
	this->u[j][i] = sinus;
	}
	}
	}
	void ut_constructor() {
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++) {
	this->ut[j][i] = this->u[i][j];
	}
	}
	void min_value() {
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++) {
	if (j > i) {
	if ((this->u[i][j] != 0) && (this->u[i][j] != 1) && (fabs(this->min) > fabs(this->u[i][j]))) {
	this->min = this->u[i][j];
	}
	else continue;
	}
	}
	}
	void ut_a_u() {
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++) {
	this->a1[i][j] = 0;
	for (int k = 0; k < column; k++) {
	this->a1[i][j] = this->a1[i][j] + this->ut[i][k] * this->a[k][j];
	}
	}
	for (int i = 0; i < column; i++)
	for (int j = 0; j < column; j++) {
	this->a[i][j] = 0;
	for (int k = 0; k < column; k++) {
	this->a[i][j] = this->a[i][j] + this->a1[i][k] * this->u[k][j];
	}
	}
	}*/

public:
	void eigenvec() {
		// приравниваем матрицу a к матрице covmat
		for (int i = 0; i < column; i++)
			for (int j = 0; j < column; j++) {
				this->a[i][j] = this->covmat[i][j];
			}
		int count = 0;
		do {
			amax = 0.1 * EPSILON;
			count += 1;
			//max_element();
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++)
				{
					if (j > i) {
						if (fabs(this->amax) < fabs(this->a[i][j])) {
							this->amax = this->a[i][j];
							this->ii = i;
							this->jj = j;
						}
						else continue;
					}
					else continue;
				}
			//cout << endl << "amax:" << this->amax;
			this->f = 0.5 * atan((2 * this->amax) / (this->a[ii][ii] - this->a[jj][jj]));
			this->cosinus = cos(this->f);
			this->sinus = sin(this->f);
			//u_init();
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++)
				{
					this->u[i][j] = 0;
				}
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++) {
					if (i == j) { this->u[i][j] = 1; }
					if (i == ii) { this->u[i][i] = cosinus; }
					if (i == jj) { this->u[i][jj] = cosinus; }
					if ((j == jj) && (i == ii)) {
						this->u[i][j] = -sinus;
						this->u[j][i] = sinus;
					}
				}
			//временный код 
			/*cout << endl;
			cout << "u:" << endl;
			for (int i = 0; i < column; i++) {
				cout << endl;
				for (int j = 0; j < column; j++)
					cout << this->u[i][j] << " ";
			} */
			//-----------— 
			//min_value();
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++) {
					if (j > i) {
						if ((this->u[i][j] != 0) && (this->u[i][j] != 1) && (fabs(this->min) > fabs(this->u[i][j]))) {
							this->min = this->u[i][j];
						}
					}
					else continue;
				}
			//ut_constructor();
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++) {
					this->ut[j][i] = this->u[i][j];
				}
			//ut_a_u();
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++) {
					this->a1[i][j] = 0;
					for (int k = 0; k < column; k++) {
						this->a1[i][j] = this->a1[i][j] + this->ut[i][k] * this->a[k][j];
					}
				}
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++) {
					this->a[i][j] = 0;
					for (int k = 0; k < column; k++) {
						this->a[i][j] = this->a[i][j] + this->a1[i][k] * this->u[k][j];
					}
				}
			// записываем собственные значения
			for (int i = 0; i < column; i++)
				for (int j = 0; j < column; j++) {
					if (i == j){
						eigvalue[i] = a[i][j];
					}
				}
			//--------------------------------
			if (count == 1) {
				for (int i = 0; i < column; i++)
					for (int j = 0; j < column; j++) {
						u1[i][j] = u[i][j];
					}
			}
			else {
				for (int i = 0; i < column; i++)
					for (int j = 0; j < column; j++) {
						u2[i][j] = 0;
						for (int k = 0; k < column; k++) {
							u2[i][j] = u2[i][j] + u1[i][k] * u[k][j];
						}
					}
				for (int i = 0; i < column; i++)
					for (int j = 0; j < column; j++) {
						u1[i][j] = u2[i][j];  // u2 - массив собственных векторов (не сортированный)
					}
			}
		} while (fabs(min) > EPSILON);

		//временный фрагмент кода
		cout << endl << "Eigenvectors: ";
		for (int i = 0; i < column; i++)
		{
			cout << endl;
			for (int j = 0; j < column; j++)
			{
				cout << u2[i][j] << " ";
			}
		}

		cout << endl << "EigenValues: " << endl;
		for (int i = 0; i < column; i++) {
			cout << eigvalue[i] << " ";
		}
		//-----------------------
		bubbleSort(eigvalue, column);
		cout << endl << "EigenValues: " << endl;
		for (int i = 0; i < column; i++) {
			cout << eigvalue[i] << " ";
		}
		// ПРАВИЛО СЛОМАННОЙ ТРОСТИ НАЧИНАЕТСЯ ЗДЕСЬ=============
		double trС = 0;  // след ковариационной матрицы 
		for (int i = 0; i < column; i++)
			for (int j = 0; j < column; j++) {
				if (i == j) {
					trС += this->covmat[i][j];
				}
				else continue;
			}
		int n = 0;
		int factors = 0;  // количество факторов
		for (int i = 0; i < column; i++) {
			double l = 0;
			double koef = 0;
			koef = this->eigvalue[i] / trС;
			n = column - i;
			for (int number = 0; number < n; number++) {
				double val = number;
				l += 1/(val + 1);
			}
			if (koef > l) {
				factors += 1;
			}
		}
		//=====================================================
	}
	void bubbleSort(double* arrayPtr, int length_array)  // сортировка пузырьком
	{
		double temp = 0; // временная переменная для хранения элемента массива
		bool exit = false; // болевая переменная для выхода из цикла, если массив отсортирован

		while (!exit) // пока массив не отсортирован
		{
			exit = true;
			for (int int_counter = 0; int_counter < (length_array - 1); int_counter++) // внутренний цикл
																					   //сортировка пузырьком по возрастанию - знак >
																					   //сортировка пузырьком по убыванию - знак <
				if (arrayPtr[int_counter] < arrayPtr[int_counter + 1]) // сравниваем два соседних элемента
				{
					// выполняем перестановку элементов массива
					temp = arrayPtr[int_counter];
					arrayPtr[int_counter] = arrayPtr[int_counter + 1];
					arrayPtr[int_counter + 1] = temp;
					exit = false; // на очередной итерации была произведена перестановка элементов
				}
		}
	}
};
void main()
{
	PCA pca_obj;
	PCA* p_pca_obj = &pca_obj; // указатель на объект 
							   //создаем объект (сначала строки, поток столбцы) 
	pca_obj.init(5, 5);
	pca_obj.mean_func();
	pca_obj.centered_func();
	pca_obj.dispersion_func();
	pca_obj.covmat_func();
	//double b = pca_obj.max_element();
	//cout << "max element:";
	//cout << b;
	pca_obj.eigenvec();
	int a;
	cin >> a;
}
